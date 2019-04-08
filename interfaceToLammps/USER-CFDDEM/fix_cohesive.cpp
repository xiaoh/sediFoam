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

#include "string.h"
#include "stdlib.h"
#include "fix_cohesive.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "mpi.h"
#include "comm.h"
#include "memory.h"


using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 10000;
/* ---------------------------------------------------------------------- */

FixCohe::FixCohe(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix cohesive command");

  ah = atof(arg[3]);
  lam = atof(arg[4]);
  smin = atof(arg[5]);
  smax = atof(arg[6]);
  opt = atoi(arg[7]);

  nmax = 0;
  nvalues = 7;   // Number of output columns : PID FX FY FZ NX NY NZ
  laststep = -1;
  size_local_cols = nvalues;
  local_flag = 1;
  laststep_local = -1; // last time step for compute_local()
  // compute_local_flag = 1; // Calling compute_local flag
}

/* ---------------------------------------------------------------------- */

int FixCohe::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCohe::init()
{ 

 int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  

if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixCohe::init_list(int id, NeighList *ptr)
{
      list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixCohe::compute_local()
{
     if(laststep_local == update->ntimestep)
       return;
        
      int npairs = count_pairs(0);
      
     //nmax = nconts
     if (npairs > nmax)     
       reallocate(npairs);
     size_local_rows = npairs;
     size_local_cols = nvalues;
    if(npairs) 
        calc_pairs();

     laststep_local = update->ntimestep;
        
}



/* ---------------------------------------------------------------------- */

void FixCohe::setup()
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(1,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }

}

/* ---------------------------------------------------------------------- */

void FixCohe::min_setup()
{
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixCohe::post_force(int vflag)
{
  int i,ii=0,j=0,jj=0,k,*numneigh,inum,**firstneigh;
  int jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r, rinv;
  double ccel,ccelx,ccely,ccelz;
  int *jlist,*ilist;
  double del;
  double **f = atom->f;
  double **x = atom->x;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *mask = atom->mask;

  double PInv = 0.25/atan(1.0);
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  if (opt == 0){
   printf(" inum %i", inum);
  
  for (ii = 0; ii < nlocal; ii++){
	  i = ilist[ii];
	  if (!(mask[i] & groupbit)) continue;
	  xtmp = x[i][0];
	  ytmp = x[i][1];
	  ztmp = x[i][2];
	  radi = radius[i];
	  jlist = firstneigh[i];
	  jnum = numneigh[i];

	  for (jj = 0; jj < jnum; jj++) {
		  j = jlist[jj];
		  delx = xtmp - x[j][0];
		  dely = ytmp - x[j][1];
		  delz = ztmp - x[j][2];
		  rsq = delx*delx + dely*dely + delz*delz;
		  radj = radius[j];
		  radsum = radi + radj;

	  if (rsq < (radsum + smax)*(radsum + smax)) {
	    r = sqrt(rsq);
	    del = r - radsum;
	    if (del > lam*PInv)
	      ccel = - ah*radsum*lam*
		(6.4988e-3 - 4.5316e-4*lam/del + 1.1326e-5*lam*lam/del/del)/del/del/del;
	    else if (del > smin)
	      ccel = - ah * (lam + 22.242*del)*radsum*lam/24.0/(lam + 11.121*del)
		/(lam + 11.121*del)/del/del;
	    else 
	      ccel = - ah * (lam + 22.242*smin)*radsum*lam/24.0/(lam + 11.121*smin)
	    /(lam + 11.121*smin)/smin/smin;
	    rinv = 1/r;

	    ccelx = delx*ccel*rinv ;
	    ccely = dely*ccel*rinv ;
	    ccelz = delz*ccel*rinv ;
	    f[i][0] += ccelx;
	    f[i][1] += ccely;
	    f[i][2] += ccelz;

	    if (newton_pair || j < nlocal) {
	      f[j][0] -= ccelx;
	      f[j][1] -= ccely;
	      f[j][2] -= ccelz;
	    }
	  }
	}
  }
}else if (opt ==1){      


   for (ii = 0; ii < nlocal; ii++){
	i = ilist[ii];
	if (!(mask[i] & groupbit)) continue;
	xtmp = x[i][0];
	ytmp = x[i][1];
	ztmp = x[i][2];
	radi = radius[i];
	jlist = firstneigh[i];
	jnum = numneigh[i];
                  

	for (jj = 0; jj < jnum; jj++) {
	  j = jlist[jj];
	  delx = xtmp - x[j][0];
	  dely = ytmp - x[j][1];
	  delz = ztmp - x[j][2];
	  rsq = delx*delx + dely*dely + delz*delz;
	  radj = radius[j];
	  radsum = radi + radj;
	  
	  if (rsq < (radsum + smax)*(radsum + smax)){
	    r = sqrt(rsq);
	    del = r - radsum;
	    if (del > smin)
	      ccel = - ah*pow(radsum,6)/6.0/del/del/(r + radsum)/(r + radsum)/r/r/r;
	    else 
	      ccel = 0;
	     // ccel = - ah*pow(radsum,6)/6.0/smin/smin/(smin+ 2.0*radsum)/(smin + 2.0*radsum)
	//	/(smin + radsum)/(smin + radsum)/(smin + radsum);
	    rinv = 1/r;

	    ccelx = delx*ccel*rinv;
	    ccely = dely*ccel*rinv;
	    ccelz = delz*ccel*rinv;
	    f[i][0] += ccelx;
	    f[i][1] += ccely;
	    f[i][2] += ccelz;

		if (newton_pair || j < nlocal){
	      f[j][0] -= ccelx;
	      f[j][1] -= ccely;
	      f[j][2] -= ccelz;
	    }
	  }
	}
  }
}else error->all(FLERR,"invalid option for cohesive force model");
}


/* ---------------------------------------------------------------------- */

void FixCohe::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCohe::min_post_force(int vflag)
{
  post_force(vflag);
}

void FixCohe::extract_cohe(int *p_opt, double *p_ah, double *p_lam, 
				   double *p_smin,double *p_smax)
{
  *p_opt = opt;
  *p_ah = ah;
  *p_lam = lam;
  *p_smin = smin;
  *p_smax = smax;
}

/*---------------------------------------------------------------------------------------------*/

/* Count the number of pairs of the contact first. */

int FixCohe::count_pairs(int flag)
{ 

  
  int i,j,m,ii,jj,inum,jnum,itype;
  double xtmp,ytmp,ztmp,delx,dely,delz,radi,radsum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq,factor_lj,factor_coul;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;


 // double *special_coul = force->special_coul;
 // double *special_lj = force->special_lj;
// Invoke half neighbor list (will copy or build if necessary)  
  if (flag == 0) neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Loop over neighbours of my atoms 
  // Skip if I or J are not in group 
  // for flag = 0 just count the pair interactions within force cutoff 
  // for flag = 1, calculate requested output fields. 
  m = 0;
  for (ii = 0; ii<inum; ii++)
     {i = ilist[ii];
      if(!(mask[i] & groupbit)) continue;
       
     xtmp = x[i][0];
     ytmp = x[i][1];
     ztmp = x[i][2];
     radi = radius[i];
     itype = type[i];
     jlist = firstneigh[i];
     jnum = numneigh[i];

   for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
     // factor_lj = special_lj[sbmask(j)];
     // factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;
      if (newton_pair == 0 && j >= nlocal) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radius[j];
       if( rsq < (radsum+smax)*(radsum+smax))
          m++;
       }
   }
     m = 2*m ;  // Counting each contact twice for outputting it as a single particle contact

 return m;

}

void FixCohe::calc_pairs()
{  
   int i,j=0,m,n,ii=0,jj=0,k,inum,*numneigh,**firstneigh;  
   int jnum;
   double xtmp,ytmp,ztmp,delx,dely,delz,radi,radj,radsum;
   double facti,factj;
   double rsq,rinv,r,del;
   int *jlist, *ilist;
   double ccel,ccelx,ccely,ccelz;
   double **x = atom->x;
   double **v = atom->v;
   double f[3]; 
   int *tag = atom->tag;
   double **omega = atom->omega;
   double torque[3];
   double *radius = atom->radius;
   double *rmass = atom->rmass;
   int *mask = atom->mask;
   int nlocal = atom->nlocal;
   
  double PInv = 0.25/atan(1.0); 
  inum = list->inum;
  ilist = list->ilist;
  numneigh=list->numneigh;
  firstneigh=list->firstneigh;

  n = 0;
  nmax = size_local_rows;
  
  if (opt == 0)
  {
    
  for ( ii = 0; ii < nlocal; i++)
     {
         i = ilist[ii];
         xtmp = x[i][0];
         ytmp = x[i][1];
         ztmp = x[i][2];
         radi = radius[i];
         jlist = firstneigh[i];
         jnum = numneigh[i];

       	  for (jj = 0; jj < jnum; jj++)
          {
		  j = jlist[jj];
		  delx = xtmp - x[j][0];
		  dely = ytmp - x[j][1];
		  delz = ztmp - x[j][2];
		  rsq = delx*delx + dely*dely + delz*delz;
		  radj = radius[j];
		  radsum = radi + radj;

	         if (rsq < (radsum + smax)*(radsum + smax)) 
                  {
                    f[0] = 0.0 ; f[1] = 0.0 ; f[2] = 0.0;
	            r = sqrt(rsq);
	            del = r - radsum;
	            if (del > lam*PInv)
	            ccel = - ah*radsum*lam*(6.4988e-3 - 4.5316e-4*lam/del + 1.1326e-5*lam*lam/del/del)/del/del/del;
	            else if (del > smin)
	            ccel = - ah * (lam + 22.242*del)*radsum*lam/24.0/(lam + 11.121*del)/(lam + 11.121*del)/del/del;
	            else 
	            ccel = - ah * (lam + 22.242*smin)*radsum*lam/24.0/(lam + 11.121*smin)/(lam + 11.121*smin)/smin/smin;
	    
                    rinv = 1/r;

	            ccelx = delx*ccel*rinv ;
	            ccely = dely*ccel*rinv ;
	            ccelz = delz*ccel*rinv ;
	            f[0] += ccelx;
	            f[1] += ccely;
	            f[2] += ccelz;

               //     facti = (radi + (del*0.5))/2 ;
               //     factj = (radj + (del*0.5))/2 ;
                     
                         array_local[n][0] = tag[i];
                         array_local[n][1] = f[0];
                         array_local[n][2] = f[1];
                         array_local[n][3] = f[2];
                         array_local[n][4] = delx;
                         array_local[n][5] = dely;
                         array_local[n][6] = delz;  
			n++;
                         array_local[n][0] = tag[j];
                         array_local[n][1] = -f[0];
                         array_local[n][2] = -f[1];
                         array_local[n][3] = -f[2];
                         array_local[n][4] = -delx;
                         array_local[n][5] = -dely;
                         array_local[n][6] = -delz;  
                        n++;

	          }
            }
        
     }
 }



else if (opt ==1){

   for (ii = 0; ii < nlocal; ii++){
	i = ilist[ii];
	if (!(mask[i] & groupbit)) continue;
	xtmp = x[i][0];
	ytmp = x[i][1];
	ztmp = x[i][2];
	radi = radius[i];
	jlist = firstneigh[i];
	jnum = numneigh[i];

	for (jj = 0; jj < jnum ; jj++) {
	  j = jlist[jj];
	  delx = xtmp - x[j][0];
	  dely = ytmp - x[j][1];
	  delz = ztmp - x[j][2];
	  rsq = delx*delx + dely*dely + delz*delz;
	  radj = radius[j];
	  radsum = radi + radj;
	  
	  if (rsq < (radsum + smax)*(radsum + smax)){
	    r = sqrt(rsq);
	    del = r - radsum;
            f[0] = 0.0 ; f[1]= 0.0; f[2]= 0.0;
	    if (del > smin)
	      ccel = - ah*pow(radsum,6)/6.0/del/del/(r + radsum)/(r + radsum)
		/r/r/r;
	    else 
	      ccel = 0;
//	      ccel = - ah*pow(radsum,6)/6.0/smin/smin/(smin+ 2.0*radsum)/(smin + 2.0*radsum)
//		/(smin + radsum)/(smin + radsum)/(smin + radsum);
	    rinv = 1/r;

	    ccelx = delx*ccel*rinv;
	    ccely = dely*ccel*rinv;
	    ccelz = delz*ccel*rinv;
	    f[0] += ccelx;
	    f[1] += ccely;
	    f[2] += ccelz;
            
                     
                         array_local[n][0] = tag[i];
                         array_local[n][1] = f[0];
                         array_local[n][2] = f[1];
                         array_local[n][3] = f[2];
                         array_local[n][4] = delx;
                         array_local[n][5] = dely;
                         array_local[n][6] = delz;  
			n++;
                         array_local[n][0] = tag[j];
                         array_local[n][1] = -f[0];
                         array_local[n][2] = -f[1];
                         array_local[n][3] = -f[2];
                         array_local[n][4] = -delx;
                         array_local[n][5] = -dely;
                         array_local[n][6] = -delz;  
                       n++;

	  }
	}
  }
}else error->all(FLERR,"invalid compute option for this fix cohesion");



}




/*-------------------------------------------------------------------------------------------*/

void FixCohe::reallocate(int n)
{  

 // grow vector or array and indices array
  while (nmax < n) nmax += DELTA ;
  
  memory->destroy(array_local);
  memory->create(array_local,nmax,nvalues,"Cohe:array");

  size_local_rows = n;
  size_local_cols = nvalues;

}




