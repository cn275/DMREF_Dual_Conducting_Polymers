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

#include <mpi.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_addforce_field.h"
#include "atom.h"
#include "atom_masks.h"
#include "accelerator_kokkos.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "universe.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"

#define MAXLINE 1024
#define EPSILONR 1.0e-6

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM, BMP, RSQ, RLINEAR};

/* ---------------------------------------------------------------------- */

FixAddForceField::FixAddForceField(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), estr(NULL), idregion(NULL), sforce(NULL)

{
  if (narg < 6) error->all(FLERR,"Illegal fix addforce command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  scale_factor = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  xstr = ystr = zstr = NULL;

  ntables = 0;
  tables = NULL;


  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    xvalue = force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else {
    yvalue = force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  // optional args

  nevery = 1;
  iregion = -1;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix addforce command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix addforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix addforce command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        estr = new char[n];
        strcpy(estr,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal fix addforce command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix addforce command");
      int n = strlen(&arg[iarg+1][2]) + 1;
      filestr = new char[n];
      strcpy(filestr,&arg[iarg+1][0]);

      n=strlen(&arg[iarg+2][2]) + 1;        
      tablename = new char[n];

       
      strcpy(tablename,&arg[iarg+2][0]);
      iarg += 3;
//      int me;
//      MPI_Comm_rank(world,&me);
//      tables = (Table *)
//      memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");
//      Table *tb = &tables[ntables];
//      null_table(tb);
//      if (me == 0) read_table(tb,arg[7],arg[8]);
    } else if (strcmp(arg[iarg],"scale_factor") == 0) {

//      fprintf(universe->uscreen," scale=: %s\n", arg[iarg+1]);     
//      scale_factor = atoi(arg[iarg+1]);
      scale_factor = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Illegal fix addforce command");
  }



  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,4,"addforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixAddForceField::~FixAddForceField()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] filestr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixAddForceField::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddForceField::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  }
  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0)
      error->all(FLERR,"Variable name for fix addforce does not exist");
    if (input->variable->atomstyle(evar)) estyle = ATOM;
    else error->all(FLERR,"Variable for fix addforce is invalid style");
  } else estyle = NONE;

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix addforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (varflag == CONSTANT && estyle != NONE)
    error->all(FLERR,"Cannot use variable energy with "
               "constant force in fix addforce");
  if ((varflag == EQUAL || varflag == ATOM) &&
      update->whichflag == 2 && estyle == NONE)
    error->all(FLERR,"Must use variable energy with fix addforce");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }



//  fprintf(universe->uscreen," Keyword: %s\n", tablename);
//  fprintf(universe->uscreen," Keyword: %s\n", tablename);
  if (1==1){
    int me;
    MPI_Comm_rank(world,&me);
//    fprintf(universe->uscreen," Keyword: %s\n", tablename);
    tables = (Table *)  memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");
//    fprintf(universe->uscreen," Keyword: %s\n", tablename);
    Table *tb = &tables[ntables];
    null_table(tb);
    if (me == 0) read_table(tb,filestr,tablename);
    bcast_table(tb);



  }


  if (0==1){
  error->all(FLERR,"force end");
  }


}

/* ---------------------------------------------------------------------- */

void FixAddForceField::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForceField::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForceField::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;

  Table *tb;
  tb = &tables[0];

  int me;
  MPI_Comm_rank(world,&me);

//  if (me != 0){
//    fprintf(universe->uscreen," me: %s\n", me);
//  }


//  int tablex,tabley,tablez;


  int tablex=tb->nx;
  int tabley=tb->ny;
  int tablez=tb->nz;

//  tablex=NX;
//  tabley=NY;
//  tablez=NZ;

//  fprintf(universe->uscreen,"tablex/y/z: %d %d %d\n",tablex,tabley,tablez); 

//  double xspace, yspace, zspace, xinv, yinv, zinv;
//  double 
  lx=(domain->xprd), ly=(domain->yprd) ,lz=(domain->zprd);
  xspace=lx/(float)(tablex);
  yspace=ly/(float)(tabley);
  zspace=lz/(float)(tablez);

//  xinv=(float)(tablex)/lx;
//  yinv=(float)(tabley)/ly;
//  zinv=(float)(tablez)/lz;




//  int xgrid,ygrid,zgrid, gridind;

//  fprintf(universe->uscreen,"x/ly/lz: %f %f %f\n",xspace,yspace,zspace);

  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  if (update->ntimestep % nevery) return;

  if (lmp->kokkos)
    atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
                      (unsigned int) F_MASK);

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if ((varflag == ATOM || estyle == ATOM) && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,4,"addforce:sforce");
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    double unwrap[3];
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);

        xgrid=static_cast<int> (fmod((fmod(x[i][0],lx)+lx),lx)/xspace);
        ygrid=static_cast<int> (fmod((fmod(x[i][1],ly)+ly),ly)/yspace);
        zgrid=static_cast<int> (fmod((fmod(x[i][2],lz)+lz),lz)/zspace);
//        int gridind=(static_cast<int> (fmod((fmod(x[i][0],lx)+lx),lx)/xspace))+tablex*(static_cast<int> (fmod((fmod(x[i][1],ly)+ly),ly)/yspace))+(tabley*tablex)*(static_cast<int> (fmod((fmod(x[i][2],lz)+lz),lz)/zspace));
        int gridind=xgrid+ygrid*tablex+zgrid*(tablex*tabley);

//        fprintf(universe->uscreen,"loc: %d %f %f %f, spacings:  %d %d %d %d  yep force vec: %f %f %f\n",i, x[i][0],x[i][1],x[i][2], tablex, tabley, tablez, gridind, tb->fx[gridind], tb->fy[gridind],tb->fz[gridind]);

        foriginal[0] += tb->efile[gridind];
        foriginal[1] += tb->fx[gridind];
        foriginal[2] += tb->fy[gridind];
        foriginal[3] += tb->fz[gridind];
        f[i][0] += tb->fx[gridind];
        f[i][1] += tb->fy[gridind];
        f[i][2] += tb->fz[gridind];
      }

  // variable force, wrap with clear/add
  // potential energy = evar if defined, else 0.0
  // wrap with clear/add

  }
//    fprintf(universe->uscreen,"energy: %f\n",foriginal[0]);

    } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],4,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],4,0);
    if (estyle == ATOM)
      input->variable->compute_atom(evar,igroup,&sforce[0][3],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        if (estyle == ATOM) foriginal[0] += sforce[i][3];
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        if (xstyle == ATOM) f[i][0] += sforce[i][0];
        else if (xstyle) f[i][0] += xvalue;
        if (ystyle == ATOM) f[i][1] += sforce[i][1];
        else if (ystyle) f[i][1] += yvalue;
        if (zstyle == ATOM) f[i][2] += sforce[i][2];
        else if (zstyle) f[i][2] += zvalue;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForceField::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForceField::min_post_force(int vflag)
{
  post_force(vflag);
}



void FixAddForceField::null_table(Table *tb)
{
  tb->rfile = tb->efile = tb->ffile = NULL;
  tb->e2file = tb->f2file = NULL;
  tb->rsq = tb->drsq = tb->e = tb->de = NULL;
  tb->f = tb->df = tb->e2 = tb->f2 = NULL;
  tb->fx = tb->fy = tb-> fz = NULL;
}


/* ---------------------------------------------------------------------- */



void FixAddForceField::read_table(Table *tb, char *file, char *keyword)
{
  char line[MAXLINE];

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }


//  fprintf(universe->uscreen," Keyword: %s\n", keyword);
//  fprintf(universe->uscreen," Keyword size :%d\n",strlen(keyword));


  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL){
      error->one(FLERR,"Did not find keyword in table file");}
    if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
    if (line[0] == '#') continue;                          // comment
    char *word = strtok(line," \t\n\r");
//    fprintf(universe->uscreen," word: %s\n", word);
    if (strcmp(word,keyword) == 0) break;           // matching keyword
    fgets(line,MAXLINE,fp);                         // no match, skip section
    param_extract(tb,line);
    fgets(line,MAXLINE,fp);
    for (int i = 0; i < tb->ninput; i++) fgets(line,MAXLINE,fp);
  }


   
//  fprintf(universe->uscreen," before create block\n");

  fgets(line,MAXLINE,fp);
  param_extract(tb,line);
  memory->create(tb->rfile,tb->ninput,"pair:rfile");
  memory->create(tb->efile,tb->ninput,"pair:efile");
  memory->create(tb->ffile,tb->ninput,"pair:ffile");
  memory->create(tb->fx,tb->ninput,"pair:fx");
  memory->create(tb->fy,tb->ninput,"pair:fy");
  memory->create(tb->fz,tb->ninput,"pair:fz");

  tablength=tb->ninput;
//  fprintf(universe->uscreen," tablength %d\n",tablength);

  
  tb->ntablebits = 0;
  int masklo,maskhi,nmask,nshiftbits;

  int itmp;
  double rfile,rnew;

  int rerror = 0;
  int cerror = 0;

//  fprintf(universe->uscreen,"ninput= %d\n", tb->ninput);

  fgets(line,MAXLINE,fp);
  for (int i = 0; i < tb->ninput; i++) {
    if (NULL == fgets(line,MAXLINE,fp))
      error->one(FLERR,"Premature end of file in pair table");
    if (7 != sscanf(line,"%d %lg %lg %lg %lg %lg %lg ",
                    &itmp,&rfile,&tb->efile[i],&tb->ffile[i] ,&tb->fx[i] ,&tb->fy[i], &tb->fz[i]))  ++cerror;



    tb->efile[i]=tb->efile[i]*scale_factor;
    tb->ffile[i]=tb->ffile[i]*scale_factor;
    tb->fx[i]=tb->fx[i]*scale_factor;
    tb->fy[i]=tb->fy[i]*scale_factor;
    tb->fz[i]=tb->fz[i]*scale_factor;
//    fprintf(universe->uscreen," Keyword: %f\n",tb->efile[i] );
    


    rnew = rfile;
    if (tb->rflag == RLINEAR)
      rnew = tb->rlo + (tb->rhi - tb->rlo)*i/(tb->ninput-1);
    else if (tb->rflag == RSQ) {
      rnew = tb->rlo*tb->rlo +
        (tb->rhi*tb->rhi - tb->rlo*tb->rlo)*i/(tb->ninput-1);
      rnew = sqrt(rnew);
    } 

    if (tb->rflag && fabs(rnew-rfile)/rfile > EPSILONR) rerror++;

    tb->rfile[i] = rnew;
  }
  fclose(fp);


  double r,e,f,fx,fy,fz,rprev,rnext,eprev,enext,fleft,fright;

  int ferror = 0;
  for (int i = 1; i < tb->ninput-1; i++) {
    r = tb->rfile[i];
    rprev = tb->rfile[i-1];
    rnext = tb->rfile[i+1];
    e = tb->efile[i];
    eprev = tb->efile[i-1];
    enext = tb->efile[i+1];
    f = tb->ffile[i];
//    fprintf(universe->uscreen," grid, fx, fy, fz: %d %f %f %f\n", i, tb->fx[i], tb->fy[i], tb->fz[i]); 


//    fx = tb->fxfile[i];
//    fy = tb->fyfile[i];
//    fz = tb->fzfile[i]; 
//    fleft = - (e-eprev) / (r-rprev);
//    fright = - (enext-e) / (rnext-r);
//    if (f < fleft && f < fright) ferror++;
//    if (f > fleft && f > fright) ferror++;
  }

//  for (int i = 0; i < tb->ninput; i++) {
//    fprintf(universe->uscreen," grid, fx, fy, fz: %d %f %f %f\n", i, tb->fx[i], tb->fy[i], tb->fz[i]);
//  }

  if (ferror) {
    char str[128];
    sprintf(str,"%d of %d force values in table are inconsistent with -dE/dr.\n"
            "  Should only be flagged at inflection points",ferror,tb->ninput);
    error->warning(FLERR,str);
  }


  if (rerror) {
    char str[128];
    sprintf(str,"%d of %d distance values in table with relative error\n"
            "  over %g to re-computed values",rerror,tb->ninput,EPSILONR);
    error->warning(FLERR,str);
  }
  if (cerror) {
    char str[128];
    sprintf(str,"%d of %d lines in table were incomplete\n"
            "  or could not be parsed completely",cerror,tb->ninput);
    error->warning(FLERR,str);
  }
}
 




/* ----------------------------------------------------------------------*/

void FixAddForceField::param_extract(Table *tb, char *line)
{ 
  tb->ninput = 0;
  tb->rflag = NONE;
  tb->fpflag = 0;
  
  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->ninput = atoi(word);
    } else if (strcmp(word,"R") == 0 || strcmp(word,"RSQ") == 0 ||
               strcmp(word,"BITMAP") == 0) { 
      if (strcmp(word,"R") == 0) tb->rflag = RLINEAR;
      else if (strcmp(word,"RSQ") == 0) tb->rflag = RSQ;
      else if (strcmp(word,"BITMAP") == 0) tb->rflag = BMP;
      word = strtok(NULL," \t\n\r\f");
      tb->rlo = atof(word);
      word = strtok(NULL," \t\n\r\f");
      tb->rhi = atof(word);
    } else if (strcmp(word,"FPRIME") == 0) {
      tb->fpflag = 1;
      word = strtok(NULL," \t\n\r\f");
      tb->fplo = atof(word);
      word = strtok(NULL," \t\n\r\f");
      tb->fphi = atof(word);

    } else if (strcmp(word,"NX") == 0) {
      tb->nxflag = 1;
      word = strtok(NULL," \t\n\r\f");
//      fprintf(universe->uscreen,"word: %s\n",word); 
      tb->nx = atof(word);
      NX=atof(word);
//      fprintf(universe->uscreen,"nx= %d\n", tb->nx); 

    } else if (strcmp(word,"NY") == 0) {
      tb->nyflag = 1;
      word = strtok(NULL," \t\n\r\f");
      tb->ny = atof(word);
      NY=atof(word);

    } else if (strcmp(word,"NZ") == 0) {
      tb->nzflag = 1;
      word = strtok(NULL," \t\n\r\f");
      tb->nz = atof(word);
      NZ=atof(word);

    } else {
      printf("WORD: %s\n",word);
      error->one(FLERR,"Invalid keyword in pair table parameters");
    }
    word = strtok(NULL," \t\n\r\f");
  }
  
  if (tb->ninput == 0) error->one(FLERR,"Pair table parameters did not set N");
}






/* ----------------------------------------------------------------------*/





void FixAddForceField::compute_table(Table *tb)
{
  int tlm1 = tablength-1;

//  fprintf(universe->uscreen," tablength %d\n",tablength);

  double inner;
  if (tb->rflag) inner = tb->rlo;
  else inner = tb->rfile[0];
  tb->innersq = inner*inner;
  tb->delta = (tb->cut*tb->cut - tb->innersq) / tlm1;
  tb->invdelta = 1.0/tb->delta;


  if (1==1){
//    fprintf(universe->uscreen," Inside lookup loop\n"); 

//  if (tabstyle == LOOKUP) {

    memory->create(tb->e,tlm1,"pair:e");
    memory->create(tb->f,tlm1,"pair:f");
    fprintf(universe->uscreen," After EF\n"); 
    memory->create(tb->fx,tablength,"pair:fx");
    memory->create(tb->fy,tablength,"pair:fy");
    memory->create(tb->fz,tablength,"pair:fz");

//    fprintf(universe->uscreen," after creates\n"); 

    double r,rsq;
    for (int i = 0; i < tlm1; i++) {
      rsq = tb->innersq + (i+0.5)*tb->delta;
      r = sqrt(rsq);
//      tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
//      tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
    }
  }




  if (tabstyle == LINEAR) {
    memory->create(tb->rsq,tablength,"pair:rsq");
    memory->create(tb->e,tablength,"pair:e");
    memory->create(tb->f,tablength,"pair:f");
    memory->create(tb->fx,tablength,"pair:fx");
    memory->create(tb->fy,tablength,"pair:fy");
    memory->create(tb->fz,tablength,"pair:fz");
    memory->create(tb->de,tlm1,"pair:de");
    memory->create(tb->df,tlm1,"pair:df");

    double r,rsq;
    for (int i = 0; i < tablength; i++) {
      rsq = tb->innersq + i*tb->delta;
      r = sqrt(rsq);
      tb->rsq[i] = rsq;
      if (tb->match) {
        tb->e[i] = tb->efile[i];
        tb->f[i] = tb->ffile[i]/r;
      } else {
//        tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
//        tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
      }
    }

    for (int i = 0; i < tlm1; i++) {
      tb->de[i] = tb->e[i+1] - tb->e[i];
      tb->df[i] = tb->f[i+1] - tb->f[i];
    }
  }


} 







/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddForceField::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
//  fprintf(universe->uscreen," SOUP\n");
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAddForceField::compute_vector(int n)
{
  // only sum across procs one time
  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAddForceField::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}



int FixAddForceField::calc_grid(double xi, double yi, double zi, double xspace, double yspace, double zspace, double lx, double ly, double lz, int tablex, int tabley, int tablez)
{
  int res=0;
  return res;
}





void FixAddForceField::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->rfile,tb->ninput,"pair:rfile");
    memory->create(tb->efile,tb->ninput,"pair:efile");
    memory->create(tb->ffile,tb->ninput,"pair:ffile");
    memory->create(tb->fx,tb->ninput,"pair:fx");
    memory->create(tb->fy,tb->ninput,"pair:fy");
    memory->create(tb->fz,tb->ninput,"pair:fz");
    
  }

  else {
    for (int i = 0; i < tb->ninput; i++){
//        fprintf(universe->uscreen,"me=0 bcast: %d %f\n", i, tb->fx[i]);
    }
  }


  MPI_Bcast(tb->rfile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->ffile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->fx,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->fy,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->fz,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(&tb->nx,1,MPI_INT,0,world);
  MPI_Bcast(&tb->ny,1,MPI_INT,0,world);
  MPI_Bcast(&tb->nz,1,MPI_INT,0,world); 




  if (me > 0) {
    for (int i = 0; i < tb->ninput; i++){
//        fprintf(universe->uscreen,"me!=0 bcast: %d %f\n", i, tb->fx[i]);
    }
  }  



  MPI_Bcast(&tb->rflag,1,MPI_INT,0,world);
  if (tb->rflag) {
    MPI_Bcast(&tb->rlo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->rhi,1,MPI_DOUBLE,0,world);
  }
  MPI_Bcast(&tb->fpflag,1,MPI_INT,0,world);
  if (tb->fpflag) {
    MPI_Bcast(&tb->fplo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->fphi,1,MPI_DOUBLE,0,world);
  }
}


