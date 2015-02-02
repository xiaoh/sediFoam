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
  if (narg != 3) error->all(FLERR, "Illegal fix fdrag command");

  // int myrank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // printf("++++=A> fluid_drag created! %5d\n", myrank);

  ffluiddrag = NULL;
  foamCpuId = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  force_reneighbor = 1;
}

/* ---------------------------------------------------------------------- */

 FixFluidDrag::~FixFluidDrag()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);

  delete [] foamCpuId;
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

      foamCpuId[i] = 0;
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

  // int myrank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // printf("++++=A> myrank = %5d, nlocal = %5d, nghost = %5d,
  //    ndrag = %5d\n", myrank, nlocal, nghost, ndrag);
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
       f[i][0] += ffluiddrag[i][0];
       f[i][1] += ffluiddrag[i][1];
       f[i][2] += ffluiddrag[i][2];
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
  memory->grow(foamCpuId,nmax,"fluid_drag:foamCpuId");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixFluidDrag::copy_arrays(int i, int j, int delflag)
{
  ffluiddrag[j][0] = ffluiddrag[i][0];
  ffluiddrag[j][1] = ffluiddrag[i][1];
  ffluiddrag[j][2] = ffluiddrag[i][2];
  foamCpuId[j] = foamCpuId[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixFluidDrag::pack_exchange(int i, double *buf)
{
  buf[0] = ffluiddrag[i][0];
  buf[1] = ffluiddrag[i][1];
  buf[2] = ffluiddrag[i][2];
  buf[3] = foamCpuId[i];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixFluidDrag::unpack_exchange(int nlocal, double *buf)
{
  ffluiddrag[nlocal][0] = buf[0];
  ffluiddrag[nlocal][1] = buf[1];
  ffluiddrag[nlocal][2] = buf[2];
  foamCpuId[nlocal] = buf[3];
  return 4;
}
