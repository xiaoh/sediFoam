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

#ifdef FIX_CLASS

FixStyle(cohesive,FixCohe)

#else

#ifndef LMP_FIX_COHESIVE_H
#define LMP_FIX_COHESIVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCohe : public Fix {
 public:
  FixCohe(class LAMMPS *, int, char **);
//  ~FixCohe();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup();
  void min_setup();
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
  virtual void min_post_force(int);

  void extract_cohe(int *, double *, double *, double *, double *);

  double smax; //maximum seperatioin for force cutoff
   void compute_local();

 private:
  double ah; //Hammaker constant
  double lam; // London retardation wavelength
  double smin; //minimum separation
  int opt; //option for cohesive force model
  int nlevels_respa;
  int nmax, nvalues;
   int npairs;
  bigint laststep;
  bigint laststep_local;

  class NeighList *list;
  int count_pairs(int flag);  // Count number of cohesive interactions for a particle 
  void reallocate(int n);
  void calc_pairs();
};

}

#endif
#endif


