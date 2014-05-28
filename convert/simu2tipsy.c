#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"
#include "../src/libio_serial/io_serial.h"

#define NGRID   100
//#define ASCII

/*============================================================================
 *                              TIPSY DEFINITION
 *============================================================================*/
#define MAXDIM  3
typedef float   Real;

struct dump {
  double time ;
  int nbodies ;
  int ndim ;
  int nsph ;
  int ndark ;
  int nstar ;
} TIPSYheader;

struct gas_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real rho;
  Real temp;
  Real hsmooth;
  Real metals ;
  Real phi ;
} *gas_particles;

struct dark_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real eps;
  Real phi ;
} *dark_particles;

struct star_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real metals ;
  Real tform ;
  Real eps;
  Real phi ;
} *star_particles;


void get_density   (long unsigned npart, float *x, float *y, float *z, float *rdens, int LGRID);
void assign        (long unsigned npart, float *x, float *y, float *z, int LGRID, float dens[LGRID][LGRID][LGRID]);
void interp        (long unsigned npart, float *x, float *y, float *z, int LGRID, float dens[LGRID][LGRID][LGRID], float *rdens);
int  get_ngas      (partptr fst_part, long npart);
int  get_nstar     (partptr fst_part, long npart);

/*============================================================================
 *                               MAIN
 *============================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
  long unsigned int i,j,k,idim,ipart;
  int               BUFSIZE, SWAPBYTES;
  int               no_timestep;
  float             *fx, *fy, *fz, *fvx, *fvy, *fvz, *fw;
  float             *rdens;
  double            dread[7];
  float             fread[7];
  double            m_fac, x_fac, v_fac;
  double            x2tipsy, v2tipsy, m2tipsy;
  char              indata[MAXSTRING], outdata[MAXSTRING];
  FILE              *fpin, *fpout;
  partptr           cur_part;
  
  /* TIPSY stuff */
  int ndim;
  int nbodies;
  int ngas;
  int ndark;
  int nstar;
  int count;
  struct gas_particle  *gp, *lastgp;
  struct dark_particle *dp, *lastdp;
  struct star_particle *sp, *lastsp;
  
  
  if(argc<2) {
    fprintf(stderr,"Use: %s AMIGA_file\n",*argv);
    exit(1);
  }
  
  printf("\n");
  printf("\n====================================================\n");
  printf("\n  read simulation data and convert to TIPSY format\n");
  printf("\n====================================================\n\n");
  
  strcpy(outdata,argv[1]);
  strcat(outdata,".tipsy");
  
  /*====================================================
   * READ INPUT FILE
   *====================================================*/
  input(argv[1]);
  
  
  /*====================================================
   * SET CONVERSION FACTORS
   *====================================================*/
  x_fac = io.header.boxsize;                              // Mpc/h
  v_fac = io.header.boxsize * H0 / io.header.a_current;   // km/sec
  m_fac = io.header.pmass;                                // Msun/h
  fprintf(stderr,"+ conversion factors (AHF internal -> physical):\n");
  fprintf(stderr,"   positions:   %lf\n",x_fac);
  fprintf(stderr,"   velocities:  %lf\n",v_fac);
  fprintf(stderr,"   masses:      %lf\n",m_fac);
  
  x2tipsy = 1.0;        // Mpc/h  -> TIPSY internal units
  v2tipsy = 1.0;        // km/sec -> TIPSY internal units
  m2tipsy = 1.0;        // Msun/h -> TIPSY internal units
  fprintf(stderr,"+ conversion factors (physical -> TIPSY internal):\n");
  fprintf(stderr,"   positions:   %lf\n",x2tipsy);
  fprintf(stderr,"   velocities:  %lf\n",v2tipsy);
  fprintf(stderr,"   masses:      %lf\n",m2tipsy);
  
  // we will only use one conversion factor:
  x_fac *= x2tipsy;
  v_fac *= v2tipsy;
  m_fac *= m2tipsy;
  
  
  /*====================================================
   * MAYBE CALCULATE DENSITY AT PARTICLE POSITIONS
   *====================================================*/
#ifdef RDENS
  /* get density at particle positions */
  fprintf(stderr,"+ getting density at particles positions ... ");
  rdens = (float *) malloc(io.header.no_part*sizeof(float));
  get_density(io.header.no_part, fx, fy, fz, rdens, NGRID);
  fprintf(stderr,"done\n");
#endif
  
  /*====================================================
   * TIPSY HEADER
   *====================================================*/
  fprintf(stderr,"+ copying information to TIPSY header ... ");
  TIPSYheader.nbodies = (int) io.header.no_part;
  TIPSYheader.nsph    = (int) get_ngas(io.fst_part, io.header.no_part);
  TIPSYheader.nstar   = (int) get_nstar(io.fst_part, io.header.no_part);
  TIPSYheader.ndim    = (int) 3;
  TIPSYheader.time    = (double) io.header.a_current;
  
  ndim    = TIPSYheader.ndim;
  nbodies = TIPSYheader.nbodies;
  ngas    = TIPSYheader.nsph;
  nstar   = TIPSYheader.nstar;
  ndark   = TIPSYheader.ndark = nbodies - nstar - ngas;
  fprintf(stderr," (nbodies=%d ndark=%d ngas=%d nstar=%d) done\n",nbodies,ndark,ngas,nstar);
  
  
  /*====================================================
   * TIPSY PARTICLES
   *====================================================*/
  fprintf(stderr,"+ copying particles to TIPSY structures ... ");
  if(ngas != 0) {
    gas_particles = (struct gas_particle *) calloc(ngas,sizeof(*gas_particles));
    if(gas_particles == NULL) {
      printf("<sorry, no memory for gas particles>\n") ;
      return(0) ;
    }
  }
  else {
    gas_particles = NULL;
  }
  
  if(ndark != 0) {
    dark_particles = (struct dark_particle *) calloc(ndark,sizeof(*dark_particles));
    if(dark_particles == NULL) {
      printf("<sorry, no memory for dark particles>\n") ;
      return(0) ;
    }
  }
  else {
    dark_particles = NULL;
  }
  
  if(nstar != 0){
    star_particles = (struct star_particle *) calloc(nstar,sizeof(*star_particles));
    if(star_particles == NULL) {
      printf("<sorry, no memory for star particles>\n") ;
      return(0) ;
    }
  }
  else {
    star_particles = NULL;
  }
  
  /*====================================================
   * CONVERSION OF PARTICLES
   *====================================================*/
  /* transfer to TIPSY arrays */
  dp = dark_particles;
  gp = gas_particles;
  sp = star_particles;
    
  for(ipart=0; ipart<nbodies; ipart++) {
    
    // particle pointer
    cur_part = io.fst_part+ipart;
    
    // gas particle
    if (isgreaterequal(cur_part->u, PGAS)){
      
      // position, velocity & mass
      gp->pos[0] = cur_part->pos[X] * x_fac;
      gp->pos[1] = cur_part->pos[Y] * x_fac;
      gp->pos[2] = cur_part->pos[Z] * x_fac;
      gp->vel[0] = cur_part->mom[X] * v_fac;
      gp->vel[1] = cur_part->mom[Y] * v_fac;
      gp->vel[2] = cur_part->mom[Z] * v_fac;
      gp->mass   = cur_part->weight * m_fac;
      
      // additional stuff
      gp->phi    = 0.0;
      gp->hsmooth= 0.0;
      gp->rho    = 0.0;
      gp->temp   = 0.0;
      gp->metals = 0.0;
      
      // next gas particle
      gp++;
    }
    
    // star particle
    else if (fabs(cur_part->u - PSTAR) < ZERO) {
      
      // position, velocity & mass
      sp->pos[0] = cur_part->pos[X] * x_fac;
      sp->pos[1] = cur_part->pos[Y] * x_fac;
      sp->pos[2] = cur_part->pos[Z] * x_fac;
      sp->vel[0] = cur_part->mom[X] * v_fac;
      sp->vel[1] = cur_part->mom[Y] * v_fac;
      sp->vel[2] = cur_part->mom[Z] * v_fac;
      sp->mass   = cur_part->weight * m_fac;
      
      // additional stuff
      sp->metals = 0.0;
      sp->eps    = 0.0;
      sp->phi    = 0.0;
      sp->tform  = 0.0;
      
      // next star particle
      sp++;
    }

    // dm particle
    else {
      
      // position, velocity & mass
      dp->pos[0] = cur_part->pos[X] * x_fac;
      dp->pos[1] = cur_part->pos[Y] * x_fac;
      dp->pos[2] = cur_part->pos[Z] * x_fac;
      dp->vel[0] = cur_part->mom[X] * v_fac;
      dp->vel[1] = cur_part->mom[Y] * v_fac;
      dp->vel[2] = cur_part->mom[Z] * v_fac;
      dp->mass   = cur_part->weight * m_fac;
      
      // additional stuff
      dp->phi    = 0.0;
      dp->eps    = 0.001;
      
      // next dm particle
      dp++;
    }
    
  }
  fprintf(stderr,"done\n");
  
  
  
  /*====================================================
   * WRITE TIPSY FILE
   *====================================================*/
  if((fpout=fopen(outdata,"w"))==NULL) {
    printf("I cannot open %s\n", outdata);
    exit(1);
  }
  
  fwrite((char *)&TIPSYheader,sizeof(TIPSYheader),1,fpout) ;
  fwrite((char *)gas_particles,sizeof(struct gas_particle),  ngas,fpout) ;
  fwrite((char *)dark_particles,sizeof(struct dark_particle),ndark,fpout) ;
  fwrite((char *)star_particles,sizeof(struct star_particle),nstar,fpout) ;
  
  fclose(fpout);
  
  
  /*====================================================
   * FREE ARRAYS
   *====================================================*/
  free(fx);
  free(fy);
  free(fz);
  free(fvx);
  free(fvy);
  free(fvz);
#ifdef RDENS
  free(rdens);
#endif
  
}

/*=============================================================================
 * get_ngas()
 *============================================================================*/
int get_ngas(partptr fst_part, long npart)
{
  long    ipart;
  int     ngas;
  partptr cur_part;
  
  ngas = 0;
  for(ipart=0; ipart<npart; ipart++) {
    cur_part = fst_part+ipart;
    if (isgreaterequal(cur_part->u, PGAS)) {
      ngas++;
    }
  }
  return(ngas);
}

/*=============================================================================
 * get_nstar()
 *============================================================================*/
int get_nstar(partptr fst_part, long npart)
{
  long    ipart;
  int     nstar;
  partptr cur_part;
  
  nstar = 0;
  for(ipart=0; ipart<npart; ipart++) {
    cur_part = fst_part+ipart;
    if (fabs(cur_part->u - PSTAR) < ZERO) {
      nstar++;
    }
  }
  return(nstar);
}

/*=============================================================================
 * get_density.c:   calculate the density at each particle position
 *                  using a grid with LGRID^3 grid points
 *
 *                  particle positions are assumed to lie in
 *                  the range [0,1] !
 *
 *
 *     INPUT:    npart, x[npart], y[npart], z[npart], LGRID
 *
 *     OUTPUT:   rdens[npart]
 *
 *============================================================================*/
void get_density(long unsigned npart, float *x, float *y, float *z, float *rdens, int LGRID)
{
  float dens[LGRID][LGRID][LGRID], rho_mean;
  int   ix, iy, iz;
  
  rho_mean = (float)npart / (float) (LGRID*LGRID*LGRID);
  
  /* assign particles to the grid */
  assign(npart, x, y, z, LGRID, dens);
  
  /* get density in terms of mean density */
  for(ix=0; ix<LGRID; ix++)
    for(iy=0; iy<LGRID; iy++)
      for(iz=0; iz<LGRID; iz++)
        dens[ix][iy][iz] /= rho_mean;
  
  /* interpolate density contrast back to particle positions */
  interp(npart, x, y, z, LGRID, dens, rdens);
}

/*
 * assign:   assign particles to the grid points using a TSC scheme
 */
void assign(long unsigned npart, float *x, float *y, float *z, int LGRID,
            float dens[LGRID][LGRID][LGRID])
{
  int           ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1;
  long unsigned ipart;
  float         rrx, rry, rrz;
  float         hx, hy, hz;
  float         hx0, hy0, hz0, hxp1, hyp1, hzp1, hxm1, hym1, hzm1;
  
  for(ix=0; ix<LGRID; ix++)
    for(iy=0; iy<LGRID; iy++)
      for(iz=0; iz<LGRID; iz++)
        dens[ix][iy][iz]=0.0;
  
  for(ipart=0; ipart<npart; ipart++)
   {
    /* coordinates in grid units */
    rrx = x[ipart] * (float)(LGRID);
    rry = y[ipart] * (float)(LGRID);
    rrz = z[ipart] * (float)(LGRID);
    
    /* index of nearest grid point */
    ix  = (int)(rrx+0.5);
    iy  = (int)(rry+0.5);
    iz  = (int)(rrz+0.5);
    
    /* distance to nearest grid point */
    hx  = rrx - (float)ix;
    hy  = rry - (float)iy;
    hz  = rrz - (float)iz;
    
    /* keep track of peridoc boundaries */
    ix=(int)fmod(ix,LGRID);
    iy=(int)fmod(iy,LGRID);
    iz=(int)fmod(iz,LGRID);
    
    /* calculate TSC weights */
    hx0=0.75 - hx*hx;
    hxp1=0.5* pow2(0.5 + hx);
    hxm1=0.5* pow2(0.5 - hx);
    hy0=0.75 - hy*hy;
    hyp1=0.5* pow2(0.5 + hy);
    hym1=0.5* pow2(0.5 - hy);
    hz0= 0.75 - hz*hz;
    hzp1=0.5* pow2(0.5 + hz);
    hzm1=0.5* pow2(0.5 - hz);
    
    /* keep track of peridoc boundaries */
    ixp1=(int)fmod(ix+1,LGRID);
    iyp1=(int)fmod(iy+1,LGRID);
    izp1=(int)fmod(iz+1,LGRID);
    ixm1=(int)fmod(ix-1+LGRID,LGRID);
    iym1=(int)fmod(iy-1+LGRID,LGRID);
    izm1=(int)fmod(iz-1+LGRID,LGRID);
    
    /* assign particle according to weights to 27 neighboring nodes */
    dens[ixm1][iym1][izm1]   = dens[ixm1][iym1][izm1]+ hxm1*hym1 *hzm1;
    dens[ix  ][iym1][izm1]   = dens[ix  ][iym1][izm1]+ hx0 *hym1 *hzm1;
    dens[ixp1][iym1][izm1]   = dens[ixp1][iym1][izm1]+ hxp1*hym1 *hzm1;
    dens[ixm1][iy  ][izm1]   = dens[ixm1][iy  ][izm1]+ hxm1*hy0  *hzm1;
    dens[ix  ][iy  ][izm1]   = dens[ix  ][iy  ][izm1]+ hx0 *hy0  *hzm1;
    dens[ixp1][iy  ][izm1]   = dens[ixp1][iy  ][izm1]+ hxp1*hy0  *hzm1;
    dens[ixm1][iyp1][izm1]   = dens[ixm1][iyp1][izm1]+ hxm1*hyp1 *hzm1;
    dens[ix  ][iyp1][izm1]   = dens[ix  ][iyp1][izm1]+ hx0 *hyp1 *hzm1;
    dens[ixp1][iyp1][izm1]   = dens[ixp1][iyp1][izm1]+ hxp1*hyp1 *hzm1;
    dens[ixm1][iym1][iz  ]   = dens[ixm1][iym1][iz  ]+ hxm1*hym1 *hz0;
    dens[ix  ][iym1][iz  ]   = dens[ix  ][iym1][iz  ]+ hx0 *hym1 *hz0;
    dens[ixp1][iym1][iz  ]   = dens[ixp1][iym1][iz  ]+ hxp1*hym1 *hz0;
    dens[ixm1][iy  ][iz  ]   = dens[ixm1][iy  ][iz  ]+ hxm1*hy0  *hz0;
    dens[ix  ][iy  ][iz  ]   = dens[ix  ][iy  ][iz  ]+ hx0 *hy0  *hz0;
    dens[ixp1][iy  ][iz  ]   = dens[ixp1][iy  ][iz  ]+ hxp1*hy0  *hz0;
    dens[ixm1][iyp1][iz  ]   = dens[ixm1][iyp1][iz  ]+ hxm1*hyp1 *hz0;
    dens[ix  ][iyp1][iz  ]   = dens[ix  ][iyp1][iz  ]+ hx0 *hyp1 *hz0;
    dens[ixp1][iyp1][iz  ]   = dens[ixp1][iyp1][iz  ]+ hxp1*hyp1 *hz0;
    dens[ixm1][iym1][izp1]   = dens[ixm1][iym1][izp1]+ hxm1*hym1 *hzp1;
    dens[ix  ][iym1][izp1]   = dens[ix  ][iym1][izp1]+ hx0 *hym1 *hzp1;
    dens[ixp1][iym1][izp1]   = dens[ixp1][iym1][izp1]+ hxp1*hym1 *hzp1;
    dens[ixm1][iy  ][izp1]   = dens[ixm1][iy  ][izp1]+ hxm1*hy0  *hzp1;
    dens[ix  ][iy  ][izp1]   = dens[ix  ][iy  ][izp1]+ hx0 *hy0  *hzp1;
    dens[ixp1][iy  ][izp1]   = dens[ixp1][iy  ][izp1]+ hxp1*hy0  *hzp1;
    dens[ixm1][iyp1][izp1]   = dens[ixm1][iyp1][izp1]+ hxm1*hyp1 *hzp1;
    dens[ix  ][iyp1][izp1]   = dens[ix  ][iyp1][izp1]+ hx0 *hyp1 *hzp1;
    dens[ixp1][iyp1][izp1]   = dens[ixp1][iyp1][izp1]+ hxp1*hyp1 *hzp1;
   }
}

/*
 * interp:  interpolate the density the grid points back to the
 *          particle positions using again a TSC scheme
 */
void interp(long unsigned npart, float *x, float *y, float *z, int LGRID, float dens[LGRID][LGRID][LGRID], float *rdens)
{
  int           ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1;
  long unsigned ipart;
  float         rrx, rry, rrz;
  float         hx, hy, hz;
  float         hx0, hy0, hz0, hxp1, hyp1, hzp1, hxm1, hym1, hzm1;
  float         ac;
  
  for(ipart=0; ipart<npart; ipart++)
   {
    /* coordinates in grid units */
    rrx = x[ipart] * (float)(LGRID);
    rry = y[ipart] * (float)(LGRID);
    rrz = z[ipart] * (float)(LGRID);
    
    /* index of nearest grid point */
    ix  = (int)(rrx+0.5);
    iy  = (int)(rry+0.5);
    iz  = (int)(rrz+0.5);
    
    /* distance to nearest grid point */
    hx  = rrx - (float)ix;
    hy  = rry - (float)iy;
    hz  = rrz - (float)iz;
    
    /* keep track of peridoc boundaries */
    ix=(int)fmod(ix,LGRID);
    iy=(int)fmod(iy,LGRID);
    iz=(int)fmod(iz,LGRID);
    
    /* calculate TSC weights */
    hx0=0.75 - hx*hx;
    hxp1=0.5* pow2(0.5 + hx);
    hxm1=0.5* pow2(0.5 - hx);
    hy0=0.75 - hy*hy;
    hyp1=0.5* pow2(0.5 + hy);
    hym1=0.5* pow2(0.5 - hy);
    hz0= 0.75 - hz*hz;
    hzp1=0.5* pow2(0.5 + hz);
    hzm1=0.5* pow2(0.5 - hz);
    
    /* keep track of peridoc boundaries */
    ixp1=(int)fmod(ix+1,LGRID);
    iyp1=(int)fmod(iy+1,LGRID);
    izp1=(int)fmod(iz+1,LGRID);
    ixm1=(int)fmod(ix-1+LGRID,LGRID);
    iym1=(int)fmod(iy-1+LGRID,LGRID);
    izm1=(int)fmod(iz-1+LGRID,LGRID);
    
    
    ac =   dens[ixm1][iym1][izm1]* hxm1*hym1 *hzm1
    +dens[ix  ][iym1][izm1]* hx0 *hym1 *hzm1
    +dens[ixp1][iym1][izm1]* hxp1*hym1 *hzm1
    +dens[ixm1][iy  ][izm1]* hxm1*hy0  *hzm1
    +dens[ix  ][iy  ][izm1]* hx0 *hy0  *hzm1
    +dens[ixp1][iy  ][izm1]* hxp1*hy0  *hzm1
    +dens[ixm1][iyp1][izm1]* hxm1*hyp1 *hzm1
    +dens[ix  ][iyp1][izm1]* hx0 *hyp1 *hzm1
    +dens[ixp1][iyp1][izm1]* hxp1*hyp1 *hzm1
    +dens[ixm1][iym1][iz  ]* hxm1*hym1 *hz0
    +dens[ix  ][iym1][iz  ]* hx0 *hym1 *hz0
    +dens[ixp1][iym1][iz  ]* hxp1*hym1 *hz0
    +dens[ixm1][iy  ][iz  ]* hxm1*hy0  *hz0
    +dens[ix  ][iy  ][iz  ]* hx0 *hy0  *hz0;
    ac +=  dens[ixp1][iy  ][iz  ]* hxp1*hy0  *hz0
    +dens[ixm1][iyp1][iz  ]* hxm1*hyp1 *hz0
    +dens[ix  ][iyp1][iz  ]* hx0 *hyp1 *hz0
    +dens[ixp1][iyp1][iz  ]* hxp1*hyp1 *hz0
    +dens[ixm1][iym1][izp1]* hxm1*hym1 *hzp1
    +dens[ix  ][iym1][izp1]* hx0 *hym1 *hzp1
    +dens[ixp1][iym1][izp1]* hxp1*hym1 *hzp1
    +dens[ixm1][iy  ][izp1]* hxm1*hy0  *hzp1
    +dens[ix  ][iy  ][izp1]* hx0 *hy0  *hzp1
    +dens[ixp1][iy  ][izp1]* hxp1*hy0  *hzp1
    +dens[ixm1][iyp1][izp1]* hxm1*hyp1 *hzp1
    +dens[ix  ][iyp1][izp1]* hx0 *hyp1 *hzp1
    +dens[ixp1][iyp1][izp1]* hxp1*hyp1 *hzp1;
    rdens[ipart] = ac;
    
   }
}
