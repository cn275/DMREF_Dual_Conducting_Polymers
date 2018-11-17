/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(addforcefield,FixAddForceField)

#else

#ifndef LMP_FIX_ADDFORCEFIELD_H
#define LMP_FIX_ADDFORCEFIELD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAddForceField : public Fix {
 public:
  FixAddForceField(class LAMMPS *, int, char **);
  ~FixAddForceField();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();
  int calc_grid(double, double, double, double, double, double, double, double, double, int, int, int); 


 private:

  enum{LOOKUP,LINEAR,SPLINE,BITMAP};
  int tabstyle,tablength;
  int xgrid,ygrid,zgrid, NX, NY, NZ, tablex, tabley, tablez;
  struct Table {
    int ninput,rflag,fpflag,match,ntablebits, nx, ny, nz, nxflag, nyflag, nzflag;
    int nshiftbits,nmask;
    double rlo,rhi,fplo,fphi,cut;
    double *rfile,*efile,*ffile, *fxfile, *fyfile, *fzfile;
    double *e2file,*f2file;
    double innersq,delta,invdelta,deltasq6;
    double *rsq,*drsq,*e,*de,*f,*df,*e2,*f2;
    double *fx, *fy, *fz;
  };

  int **tabindex;

  int ntables;
  Table *tables;


  void allocate();
  void null_table(Table *tb);
  
  void read_table(Table *, char *, char *);
  void param_extract(Table *, char *);
  void bcast_table(Table *);
  void compute_table(Table *);

  double lx, ly, lz; 
  double xspace, yspace, zspace, xinv, yinv, zinv;


  double xvalue,yvalue,zvalue, scale_factor;
  int varflag,iregion;
  char *xstr,*ystr,*zstr,*estr, *filestr, *tablename;
  char *idregion;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle,estyle;
  double foriginal[4],foriginal_all[4];
  int force_flag;
  int ilevel_respa;



  int maxatom;
  double **sforce;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix addforce does not exist

Self-explanatory.

E: Variable name for fix addforce does not exist

Self-explanatory.

E: Variable for fix addforce is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix addforce

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix addforce

Must define an energy vartiable when applyting a dynamic
force during minimization.

*/
