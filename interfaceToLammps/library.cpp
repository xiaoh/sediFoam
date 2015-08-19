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
#include "atom_vec.h"
#include "fix_fluid_drag.h" // added JS
#include "modify.h"
#include "string.h"
#include "stdio.h"
#include "math.h"
#include "math_const.h"
#include "random_park.h"
#include "memory.h"
#include "group.h"
#include "domain.h"

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
void lammps_get_initial_np(void* ptr, int* np_)
{
  printf("start getting info..");
  int myrank;
  int mysize;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &mysize);

  LAMMPS *lammps = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lammps->atom->natoms);

  printf("creating array of the number of particle..");
  int *copyNp = new int[mysize];

  printf("initializing the number of particles in each processor");
  for (int i = 0; i < mysize; i++) {
    copyNp[i] = 0;
  }

  int nlocal = lammps->atom->nlocal;

  copyNp[myrank] = nlocal;

  printf("calling MPI_Allreduce");
  MPI_Allreduce(copyNp, np_, mysize, MPI_INT, MPI_SUM, lammps->world);

  delete [] copyNp;
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
// Originaly designed to provide coordinate and velocity information to Foam
// Not extened to provide particle diameter, density, and type as well.
// The name will, however, not be updated.
void lammps_get_initial_info(void* ptr, double* coords, double* velos, double* diam,
                     double* rho_, int* tag_, int* lmpCpuId_, int* type_)
{
  printf("start getting info..");
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  LAMMPS *lammps = (LAMMPS *) ptr;
  int natoms = static_cast<int> (lammps->atom->natoms);
  int nlocal = lammps->atom->nlocal;

  printf("initializing position and velocity");
  for (int i = 0; i < 3*nlocal; i++) {
    coords[i] = 0.0;
    velos[i] = 0.0;
  }

  printf("initializing other values");
  for (int i = 0; i < nlocal; i++) {
    diam[i] = 0.0;
    rho_[i] = 0.0;
    type_[i] = 0;
    tag_[i] = 0;
    lmpCpuId_[i] = 0;
  }

  printf("find the pointers");
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

  int id,offset;
  printf("start the loop!");
  for (int i = 0; i < nlocal; i++) {

    offset = 3*i;
    // Get coordinates
    coords[offset+0] = x[i][0];
    coords[offset+1] = x[i][1];
    coords[offset+2] = x[i][2];
    // Copy velocities
    velos[offset+0] = v[i][0];
    velos[offset+1] = v[i][1];
    velos[offset+2] = v[i][2];
    // Copy diameters
    diam[i] = r[i] * 2.0;
    // Copy density
    // commented by rui
    // modify the density
    rho_[i] = 3.0*rho[i]/(4.0*3.14159265358917323846*r[i]*r[i]*r[i]);
    // Copy type
    type_[i] = type[i];
    tag_[i] = tag[i];
    lmpCpuId_[i] = myrank;
  }
}

/* ---------------------------------------------------------------------- */
// Provide particle number.
int lammps_get_local_n(void* ptr)
{
  LAMMPS *lammps = (LAMMPS *) ptr;

  int nlocal = lammps->atom->nlocal;

  return nlocal;
}


/* ---------------------------------------------------------------------- */
// Provide the local domain of each processor
void lammps_get_local_domain(void* ptr, double* domain_)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  LAMMPS *lammps = (LAMMPS *) ptr;

  double *sublo = lammps->domain->sublo; 
  double *subhi = lammps->domain->subhi; 

  domain_[0] = sublo[0];
  domain_[1] = subhi[0];
  domain_[2] = sublo[1];
  domain_[3] = subhi[1];
  domain_[4] = sublo[2];
  domain_[5] = subhi[2];
  printf("++++=A> domain in processor %5d\n is: x0 %f, x1 %f, y0 %f; y1 %f z0 %f, z1 %f.", 
          myrank, domain_[0], domain_[1], domain_[2], domain_[3], domain_[4], domain_[5]);
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
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    printf("Incoming drag not consistent with local particle number.\n");
    printf("Incoming drag is: %5d, local particle number is: %5d.", nLocalIn, nlocal);
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

  lammps->input->one(cmd);
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

/* ---------------------------------------------------------------------- */

void lammps_create_particle(void* ptr, int npAdd, double* position, double* tag, 
                            double diameter, double rho, int type)
{

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  LAMMPS *lammps = (LAMMPS *) ptr;

  int natom = static_cast<int> (lammps->atom->natoms);
  double **x = lammps->atom->x;
  int ntype = type;
  double *xNew = new double[2];

  tagint max = 0;
  tagint maxtag_all = 0;

  tagint *tagOld = lammps->atom->tag;
  for (int i = 0; i < lammps->atom->nlocal; i++) max = MAX(max,tagOld[i]);
  // MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,lammps->world);

  for (int m = 0; m < npAdd; m++) {
    xNew[0] = position[0+m*3];
    xNew[1] = position[1+m*3];
    xNew[2] = position[2+m*3];

    lammps->atom->avec->create_atom(ntype,xNew);
    int nlocal = lammps->atom->nlocal;

    int n = lammps->atom->nlocal - 1;
    // lammps->atom->tag[n] = maxtag_all + 1;
    lammps->atom->tag[n] = tag[m];
 
    // // Setting PBC?
    // for (int i = 0; i < nlocal; i++) {
    //   double delx = xNew[0] - x[i][0];
    //   double dely = xNew[1] - x[i][1];
    //   double delz = xNew[2] - x[i][2];
    //   lammps->domain->minimum_image(delx,dely,delz);
    // }

    // TODO: may have problem here
    int igroup = lammps->group->find("active");
    int groupbit = lammps->group->bitmask[igroup];

    lammps->atom->mask[n] = 1 | groupbit;
    lammps->atom->image[n] = ((imageint) IMGMAX << IMG2BITS) |
                             ((imageint) IMGMAX << IMGBITS) | IMGMAX;
    lammps->atom->v[n][0] = 0.0;
    lammps->atom->v[n][1] = 0.0;
    lammps->atom->v[n][2] = 0.0;

    double radtmp = 0.5*diameter;
    lammps->atom->radius[n] = radtmp;
    lammps->atom->rmass[n] = 4.0*3.14159265358917323846/3.0 * radtmp*radtmp*radtmp * rho;

    for (int j = 0; j < lammps->modify->nfix; j++)
    {
      if (lammps->modify->fix[j]->create_attribute)
      {
          lammps->modify->fix[j]->set_arrays(n);
      }
    }
  }

  int nGlobal = 0;
  int nlocal = lammps->atom->nlocal;
  MPI_Allreduce(&nlocal,&nGlobal,1,MPI_INT,MPI_SUM,lammps->world);
  lammps->atom->natoms = nGlobal;

  if (lammps->atom->map_style) {
    lammps->atom->nghost = 0;
    lammps->atom->map_init();
    lammps->atom->map_set();
  }

  for (int j = 0; j < lammps->modify->nfix; j++)
  {
    bigint nt = lammps->update->ntimestep;
    lammps->modify->fix[j]->next_reneighbor = nt + 1;
  }

  delete [] xNew;

}

void lammps_delete_particle(void *ptr, int* deleteList, int nDelete)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  LAMMPS *lammps = (LAMMPS *) ptr;
  class RanPark *random = new RanPark(lammps,100);

  int i,j,m,iwhichglobal,iwhichlocal;
  int ndel,ndeltopo[4];
  int *list,*mark;

  // grow list and mark arrays if necessary

  int nmax = 0;
  if (lammps->atom->nlocal > nmax) {
    lammps->memory->destroy(list);
    lammps->memory->destroy(mark);
    nmax = lammps->atom->nmax;
    lammps->memory->create(list,nmax,"evaporate:list");
    lammps->memory->create(mark,nmax,"evaporate:mark");
  }

  // ncount = # of deletable atoms in region that I own
  // nall = # on all procs // nbefore = # on procs before me
  // list[ncount] = list of local indices of atoms I can delete

  double **x = lammps->atom->x;
  int *mask = lammps->atom->mask;
  tagint *tag = lammps->atom->tag;
  int nlocal = lammps->atom->nlocal;

  int ncount = 0;

  for (i = 0; i < nlocal; i++)
  {
    for (j = 0; j < nDelete; j++)
    {
       if (deleteList[j] == tag[i]) // delete the particle when assigned in openfoam
       {
           list[ncount++] = i;
       }
    }
  }

  int nall,nbefore;
  MPI_Allreduce(&ncount,&nall,1,MPI_INT,MPI_SUM,lammps->world);

  MPI_Scan(&ncount,&nbefore,1,MPI_INT,MPI_SUM,lammps->world);
  nbefore -= ncount;

  // ndel = total # of atom deletions, in or out of region
  // ndeltopo[1,2,3,4] = ditto for bonds, angles, dihedrals, impropers
  // mark[] = 1 if deleted

  ndel = 0;

  for (i = 0; i < nlocal; i++) mark[i] = 0;

  double flux = nDelete;

  double totalFlux = 0;
  MPI_Allreduce(&flux,&totalFlux,1,MPI_DOUBLE,MPI_SUM,lammps->world);

  totalFlux -= 0.5;

  while (nall && ndel < totalFlux) {
    iwhichglobal = static_cast<int> (nall*random->uniform());
    if (iwhichglobal < nbefore)
    {
      nbefore--;
    }
    else if (iwhichglobal < nbefore + ncount)
    {
      iwhichlocal = iwhichglobal - nbefore;
      mark[list[iwhichlocal]] = 1;
      list[iwhichlocal] = list[ncount-1];
      ncount--;
    }
    ndel++;
    nall--;
  }

  // atomic deletions
  // choose atoms randomly across all procs and mark them for deletion
  // shrink eligible list as my atoms get marked
  // keep ndel,ncount,nall,nbefore current after each atom deletion

  // delete my marked atoms
  // loop in reverse order to avoid copying marked atoms

  AtomVec *avec = lammps->atom->avec;

  for (i = nlocal-1; i >= 0; i--) {
    if (mark[i]) {
      avec->copy(lammps->atom->nlocal-1,i,1);
      lammps->atom->nlocal--;
    }
  }

  // reset global natoms and bonds, angles, etc
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  lammps->atom->natoms -= ndel;

  if (ndel && lammps->atom->map_style) {
    lammps->atom->nghost = 0;
    lammps->atom->map_init();
    lammps->atom->map_set();
  }

  // // statistics
  for (int j = 0; j < lammps->modify->nfix; j++)
  {
    bigint nt = lammps->update->ntimestep;
    lammps->modify->fix[j]->next_reneighbor = nt + 1;
    // lammps->modify->fix[j]->force_reneighbor = 1;
  }

  // ndeleted += ndel;
  // next_reneighbor = update->ntimestep + nevery;
  // delete [] random;
}
