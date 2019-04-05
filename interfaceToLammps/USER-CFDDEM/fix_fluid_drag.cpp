/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_fluid_drag.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixFluidDrag::FixFluidDrag(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR, "Illegal fix fdrag command");

  // int myrank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // printf("++++=A> fluid_drag created! %5d\n", myrank);

  ffluiddrag = NULL;
  DuDt = NULL;
  foamCpuId = NULL;
  vOld = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  force_reneighbor = 0; // when adding particles, set this value to 1.
  create_attribute = 1;

  if (narg == 3) {
    carrier_rho = 0;
  }
  else if (narg == 4) {
    carrier_rho = atoi(arg[3]);
  }
}

/* ---------------------------------------------------------------------- */

 FixFluidDrag::~FixFluidDrag()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);

//  delete [] foamCpuId;
//  delete [] vOld;
//  delete [] DuDt;
  memory->destroy(foamCpuId);
  memory->destroy(vOld);
  memory->destroy(DuDt);
  memory->destroy(ffluiddrag);
}

/* ---------------------------------------------------------------------- */

int FixFluidDrag::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFluidDrag::init()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      ffluiddrag[i][0] = 0.;
      ffluiddrag[i][1] = 0.;
      ffluiddrag[i][2] = 0.;

      DuDt[i][0] = 0.;
      DuDt[i][1] = 0.;
      DuDt[i][2] = 0.;

      foamCpuId[i] = 0;

      vOld[i][0] = 0.;
      vOld[i][1] = 0.;
      vOld[i][2] = 0.;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixFluidDrag::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFluidDrag::post_force(int vflag)
{
  // apply drag force to atoms in group
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double timeStep = update->dt;
  double rho = 0.;
  double accX = 0.;
  double accY = 0.;
  double accZ = 0.;

  double **v = atom->v;
  // get the mass of particles
  double *rmass = atom->rmass;
  double *r = atom->radius;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

     rho = 3.0*rmass[i]/(4.0*3.14159265358917323846*r[i]*r[i]*r[i]);
     accX = ((v[i][0] - vOld[i][0])/timeStep),
     accY = ((v[i][1] - vOld[i][1])/timeStep),
     accZ = ((v[i][2] - vOld[i][2])/timeStep),

     f[i][0] += ffluiddrag[i][0] +
                carrier_rho/rho*0.5*rmass[i]*(DuDt[i][0] - accX);
     f[i][1] += ffluiddrag[i][1] +
                carrier_rho/rho*0.5*rmass[i]*(DuDt[i][1] - accY);
     f[i][2] += ffluiddrag[i][2] +
                carrier_rho/rho*0.5*rmass[i]*(DuDt[i][2] - accZ);

     vOld[i][0] = v[i][0];
     vOld[i][1] = v[i][1];
     vOld[i][2] = v[i][2];

       //printf("%i, f = %e  %e  %e \n", i, ffluiddrag[i][0], ffluiddrag[i][1],ffluiddrag[i][2]);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixFluidDrag::memory_usage()
{
  double bytes = 3*atom->nmax * sizeof(double);
  bytes += atom->nmax * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixFluidDrag::grow_arrays(int nmax)
{
  memory->grow(ffluiddrag,nmax,3,"fluid_drag:ffluiddrag");
  memory->grow(DuDt,nmax,3,"fluid_drag:DuDt");
  memory->grow(foamCpuId,nmax,"fluid_drag:foamCpuId");
  memory->grow(vOld,nmax,3,"fluid_drag:vOld");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixFluidDrag::copy_arrays(int i, int j, int delflag)
{
  ffluiddrag[j][0] = ffluiddrag[i][0];
  ffluiddrag[j][1] = ffluiddrag[i][1];
  ffluiddrag[j][2] = ffluiddrag[i][2];
  DuDt[j][0] = DuDt[i][0];
  DuDt[j][1] = DuDt[i][1];
  DuDt[j][2] = DuDt[i][2];
  foamCpuId[j] = foamCpuId[i];
  vOld[j][0] = vOld[i][0];
  vOld[j][1] = vOld[i][1];
  vOld[j][2] = vOld[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixFluidDrag::pack_exchange(int i, double *buf)
{
  buf[0] = ffluiddrag[i][0];
  buf[1] = ffluiddrag[i][1];
  buf[2] = ffluiddrag[i][2];
  buf[3] = DuDt[i][0];
  buf[4] = DuDt[i][1];
  buf[5] = DuDt[i][2];
  buf[6] = foamCpuId[i];
  buf[7] = vOld[i][0];
  buf[8] = vOld[i][1];
  buf[9] = vOld[i][2];
  return 10;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixFluidDrag::unpack_exchange(int nlocal, double *buf)
{
  ffluiddrag[nlocal][0] = buf[0];
  ffluiddrag[nlocal][1] = buf[1];
  ffluiddrag[nlocal][2] = buf[2];
  DuDt[nlocal][0] = buf[3];
  DuDt[nlocal][1] = buf[4];
  DuDt[nlocal][2] = buf[5];
  foamCpuId[nlocal] = buf[6];
  vOld[nlocal][0] = buf[7];
  vOld[nlocal][1] = buf[8];
  vOld[nlocal][2] = buf[9];
  return 10;
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixFluidDrag::set_arrays(int i)
{
	foamCpuId[i] = 0;

	ffluiddrag[i][0] = 0.;
	ffluiddrag[i][1] = 0.;
	ffluiddrag[i][2] = 0.;

    DuDt[i][0] = 0.;
    DuDt[i][1] = 0.;
    DuDt[i][2] = 0.;

    vOld[i][0] = 0.;
    vOld[i][1] = 0.;
    vOld[i][2] = 0.;
}
