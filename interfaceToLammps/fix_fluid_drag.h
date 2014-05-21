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

// added by Rui Sun
#ifdef FIX_CLASS

FixStyle(fdrag,FixFluidDrag)

#else

#ifndef LMP_FIX_FLUID_DRAG_H
#define LMP_FIX_FlUID_DRAG_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFluidDrag : public Fix {
 public:
  double **ffluiddrag;
  int *foamCpuId;

  FixFluidDrag(class LAMMPS *, int, char **);
  ~FixFluidDrag();

  int setmask();
  void init();
  void setup(int);
  virtual void post_force(int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 private:

};

}

#endif
#endif
