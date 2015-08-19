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

#ifdef COMPUTE_CLASS

ComputeStyle(gran/local,ComputeGranLocal)

#else

#ifndef LMP_COMPUTE_GRAN_LOCAL_H
#define LMP_COMPUTE_GRAN_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGranLocal : public Compute {
 public:
  ComputeGranLocal(class LAMMPS *, int, char **);
  ~ComputeGranLocal();
  void init();
  void init_list(int, class NeighList *);
  void compute_local();
  double memory_usage();

 private:
  int nvalues,dflag,eflag,fflag;
  int ncount;

  int *pstyle;              // style of each requested output
  int *pindex;              // for pI, index of the output (0 to M-1)
  int singleflag;

  int nmax;
  double *vector;
  double **array;

  class NeighList *list;

  int compute_pairs(int);
  void reallocate(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid keyword in compute gran/local command

Self-explanatory.

E: No gran style is defined for compute gran/local

Self-explanatory.

E: Gran style does not support compute gran/local

The gran style does not have a single() function, so it can
not be invoked by compute gran/local.

E: Gran style does not have extra field requested by compute gran/local

The gran style does not support the pN value requested by the compute
gran/local command.

*/
