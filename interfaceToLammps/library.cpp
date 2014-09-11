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

// C or Fortran style library interface to LAMMPS
// new LAMMPS-specific functions can be added

#include "mpi.h"
#include "library.h"
#include "lammps.h"
#include "update.h"
#include "input.h"
#include "atom.h"
#include "fix_fluid_drag.h" // added JS
#include "modify.h"
#include "string.h"
#include "stdio.h"
#include "math.h"
#include "math_const.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   create an instance of LAMMPS and return pointer to it
   pass in command-line args and MPI communicator to run on
------------------------------------------------------------------------- */

void lammps_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{

  LAMMPS *lammps = new LAMMPS(argc,argv,communicator);
  *ptr = (void *) lammps;
}

/* ----------------------------------------------------------------------
   destruct an instance of LAMMPS
------------------------------------------------------------------------- */

void lammps_close(void *ptr)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  delete lammps;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void lammps_file(void *ptr, char *str)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  lammps->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *lammps_command(void *ptr, char *str)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  return lammps->input->one(str);
}

// This is for use with gdb debugging only;
void lammps_sync(void *ptr)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  int idb=0; while(idb);
  MPI_Barrier(lammps->world);
}

/* ----------------------------------------------------------------------
   add LAMMPS-specific library functions
   all must receive LAMMPS pointer as argument
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

int lammps_get_global_n(void *ptr)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lammps->atom->natoms);
  return natoms;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
// Originaly designed to provide coordinate and velocity information to Foam
// Not extened to provide particle diameter, density, and type as well.
// The name will, however, not be updated.
void lammps_get_initial_info(void* ptr, double* coords, double* velos, double* diam,
                     double* rho_, int* tag_, int* lmpCpuId_, int* type_)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  LAMMPS *lammps = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lammps->atom->natoms);

  double *copyx = new double[3*natoms];
  double *copyv = new double[3*natoms];
  double *copyd = new double[natoms];
  double *copyRho = new double[natoms];
  int *copytype = new int[natoms];
  int *copytag = new int[natoms];
  int *copylmpCpuId = new int[natoms];

  for (int i = 0; i < 3*natoms; i++) {
    copyx[i] = 0.0;
    copyv[i] = 0.0;
  }

  for (int i = 0; i < natoms; i++) {
    copyd[i] = 0.0;
    copyRho[i] = 0.0;
    copytype[i] = 0;
    copytag[i] = 0;
    copylmpCpuId[i] = 0;
  }

  //Info provided to Foam Cloud.
  double **x = lammps->atom->x;
  double **v = lammps->atom->v;
  double *r = lammps->atom->radius;
  // double *rho = lammps->atom->density;
  // commented by rui
  // lammps cannot access density but rmass;
  double *rho = lammps->atom->rmass;
  int *type = lammps->atom->type;

  // Other info used in this function
  int *tag = lammps->atom->tag;
  int nlocal = lammps->atom->nlocal;

  int id,offset;
  for (int i = 0; i < nlocal; i++) {
    id = tag[i];
    offset = 3*(id-1);
    // Get coordinates
    copyx[offset+0] = x[i][0];
    copyx[offset+1] = x[i][1];
    copyx[offset+2] = x[i][2];
    // Copy velocities
    copyv[offset+0] = v[i][0];
    copyv[offset+1] = v[i][1];
    copyv[offset+2] = v[i][2];
    // Copy diameters
    copyd[id-1] = r[i] * 2.0;
    // Copy density
    // commented by rui
    // modify the density
    copyRho[id-1] = 3.0*rho[i]/(4.0*3.14159265358917323846*r[i]*r[i]*r[i]);
    // Copy type
    copytype[id-1] = type[i];
    copytag[id-1] = tag[i];
    copylmpCpuId[id-1] = myrank;
  }

  MPI_Allreduce(copyx, coords, 3*natoms, MPI_DOUBLE, MPI_SUM,
                lammps->world);
  MPI_Allreduce(copyv, velos, 3*natoms, MPI_DOUBLE, MPI_SUM, lammps->world);
  MPI_Allreduce(copyd, diam, natoms, MPI_DOUBLE, MPI_SUM, lammps->world);
  MPI_Allreduce(copyRho, rho_, natoms, MPI_DOUBLE, MPI_SUM, lammps->world);
  MPI_Allreduce(copytag, tag_, natoms, MPI_INT, MPI_SUM, lammps->world);
  MPI_Allreduce(copylmpCpuId, lmpCpuId_, natoms, MPI_INT, MPI_SUM, lammps->world);
  MPI_Allreduce(copytype, type_, natoms, MPI_INT, MPI_SUM, lammps->world);

  delete [] copyx;
  delete [] copyv;
  delete [] copyd;
  delete [] copyRho;
  delete [] copytag;
  delete [] copylmpCpuId;
  delete [] copytype;
}

/* ---------------------------------------------------------------------- */
// Provide particle number.
int lammps_get_local_n(void* ptr)
{
  LAMMPS *lammps = (LAMMPS *) ptr;

  int nlocal = lammps->atom->nlocal;

  // printf("++++=A> nlocal is: %5d\n", nlocal);

  return nlocal;
}


/* ---------------------------------------------------------------------- */
// Provide particle info (incl. coordinate, velocity etc.) to Foam
// Note: No MPI communication occurs! Info is sent to the same processor.
void lammps_get_local_info(void* ptr, double* coords_, double* velos_,
                           int* foamCpuId_, int* lmpCpuId_, int* tag_)
{

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  LAMMPS *lammps = (LAMMPS *) ptr;

  // the pointer to the fix_fluid_drag class
  class FixFluidDrag *drag_ptr = NULL;

  int i;
  for (i = 0; i < (lammps->modify->nfix); i++)
    if (strcmp(lammps->modify->fix[i]->style,"fdrag") == 0) break;

  if (i < lammps->modify->nfix)
    //initialize the pointer
    drag_ptr = (FixFluidDrag *) lammps->modify->fix[i];


  // Coordinate, velocity etc. on Lammps *local* processor
  double **x = lammps->atom->x;
  double **v = lammps->atom->v;
  int *tag = lammps->atom->tag;
  int *type = lammps->atom->type;
  int nlocal = lammps->atom->nlocal;


  double *copyx = new double[3*nlocal];
  double *copyv = new double[3*nlocal];
  int *copyfoamCpuId = new int[nlocal];
  int *copytag = new int[nlocal];

  for (int i = 0; i < 3*nlocal; i++) {
    copyx[i] = 0.0;
    copyv[i] = 0.0;
  }

  for (int i = 0; i < nlocal; i++) {
    copyfoamCpuId[i] = 0;
    copytag[i] = 0;
  }

  for (int i = 0; i < nlocal; i++) {
    coords_[3*i+0] = x[i][0];
    coords_[3*i+1] = x[i][1];
    coords_[3*i+2] = x[i][2];
    velos_[3*i+0] = v[i][0];
    velos_[3*i+1] = v[i][1];
    velos_[3*i+2] = v[i][2];

    foamCpuId_[i] = drag_ptr->foamCpuId[i];

    lmpCpuId_[i] = myrank;
    tag_[i] = tag[i];
  }

  delete [] copyx;
  delete [] copyv;
  delete [] copyfoamCpuId;
  delete [] copytag;
}


/* ---------------------------------------------------------------------- */
// Provide particle info (incl. coordinate, velocity etc.) to Foam
// Note: MPI communication occurs!
void lammps_put_local_info(void* ptr, int nLocalIn, double* fdrag, 
                           int* foamCpuIdIn, int* tagIn)
{

  LAMMPS *lammps = (LAMMPS *) ptr;

  // the pointer to the fix_fluid_drag class
  class FixFluidDrag *drag_ptr = NULL;

  int i;
  for (i = 0; i < (lammps->modify->nfix); i++)
    if (strcmp(lammps->modify->fix[i]->style,"fdrag") == 0) break;

  int *tag = lammps->atom->tag;
  int *type = lammps->atom->type;
  int nlocal = lammps->atom->nlocal;

  if (i < lammps->modify->nfix)
    //initialize the pointer
    drag_ptr = (FixFluidDrag *) lammps->modify->fix[i];

  if (nLocalIn != nlocal)
    {
      printf("Incoming drag not consisent with local particle number.");
    }

  //initialize the tag pair to sort
  std::vector<tagpair> lmptagpair (nlocal);
  std::vector<tagpair> intagpair (nlocal);
  for (int i = 0; i < nlocal; i++) {
      lmptagpair[i].tag = tag[i];
      lmptagpair[i].index = i;
      intagpair[i].tag = tagIn[i];
      intagpair[i].index = i;
  }

  std::sort(lmptagpair.begin(), lmptagpair.end(), by_number());
  std::sort(intagpair.begin(), intagpair.end(), by_number());

  for (int j = 0; j < nlocal; j++) {

    int fromfoamid = intagpair[j].index;
    int tolmpid = lmptagpair[j].index;

    drag_ptr->foamCpuId[tolmpid] = foamCpuIdIn[fromfoamid];

    drag_ptr->ffluiddrag[tolmpid][0] = fdrag[3*fromfoamid+0];
    drag_ptr->ffluiddrag[tolmpid][1] = fdrag[3*fromfoamid+1];
    drag_ptr->ffluiddrag[tolmpid][2] = fdrag[3*fromfoamid+2];
  }
}

/* ---------------------------------------------------------------------- */

// Evolve n steps forward without overhead
void lammps_step(void *ptr, int n)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  // lammps -> Run -> command(5, n, "pre", "no", "post", "no");
  char cmd[50];
  int nc;

  nc = sprintf(cmd, "run %d pre no post no", n);
  if(nc <= 0) {
    printf("lammps_step: Command print failure.");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  lammps ->input->one(cmd);
}

/* ---------------------------------------------------------------------- */

double lammps_get_timestep(void *ptr)
{
   LAMMPS *lammps = (LAMMPS *) ptr;
   return lammps->update->dt;
}

/* ---------------------------------------------------------------------- */

void lammps_set_timestep(void *ptr, double dt_i)
{
  LAMMPS *lammps = (LAMMPS *) ptr;
  lammps->update->dt = dt_i;
}

