#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libio_serial/io_serial.h"
#include "../src/libutility/utility.h"

/*=================================================================================================
 *                                         DEFINES
 *=================================================================================================*/
//#define NPARTDENS // number density rather than mass density
//#define WRITE_PDENSFILE // also write a file containing the density at each particle position

#define IDX(ix,iy,iz)  ((iz)*NGRID*NGRID + (iy)*NGRID + (ix))

/*=================================================================================================
 *                                     COMMON VARIABLES
 *=================================================================================================*/
// a negative NGRID will place the particles into the largest possible cube defined by the positions itself
double BoxSize;
int    MinMaxScale;
double MinPos, MaxPos;
long   NumPart;
float  frac;
int    iseed;

/*=================================================================================================
 *                                       PROTOTYPES
 *=================================================================================================*/
void    get_density(long unsigned npart, double nvpart, float *x, float *y, float *z, float *weight, float *rdens, int NGRID);
void    check_density(int NGRID, float *dens);
void    assign(long unsigned npart, float *x, float *y, float *z, float *weight, int NGRID, float *dens);
void    interp(long unsigned npart, float *x, float *y, float *z, int NGRID, float *dens, float *rdens);
partptr ReadSimuAndCalcPdens(char *infile, int NGRID);
int     qcompareIDs(const void *, const void *);
int     bcompareIDs(const void *, const void *);
void    dump_parts(partptr p, char *out);
void    get_ijk(long unsigned id, int L, int *iL, int *jL, int *kL);

/*=================================================================================================
 *                                        STRUCTURES
 *=================================================================================================*/


/*=================================================================================================
 *                                           MAIN
 *=================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
  char          indata1[MAXSTRING], indata2[MAXSTRING], outdata[MAXSTRING];
  FILE          *fpout;
  int           NGRID, ilen, i, j;
  partptr       Psimu1;
  partptr       Psimu2;
  long unsigned ipart1, ipart2, id1, id2;
  double        dx,dy,dz;
  double        d12, d1L, d2L;
  int           iL, jL, kL, L;
  
  if(argc<6)
   {
    fprintf(stderr,"usage: %s simu_file1 simu_file2 NGRID fraction iseed\n",*argv);
    exit(1);
   }
  
  fprintf(stderr,"===============================================================================\n");
  fprintf(stderr,"  Read 2x simulation binaries and calculate density at each particle position\n");
  fprintf(stderr,"            Cross-Correlate both simulations using particle IDs\n");
  fprintf(stderr,"             (only write a random fraction of the particles)\n");
  fprintf(stderr,"===============================================================================\n");
  
  /*-------------------------------------------------------
   * deal with command line
   *-------------------------------------------------------*/
  NGRID = (int) (atoi(argv[3]));
  if(NGRID < 0)
    MinMaxScale = 1;
  else
    MinMaxScale = 0;
  NGRID = abs(NGRID);
  
  ilen = strlen(argv[1]);
  for(i=ilen-1; i>=0; i--)
   /* chop filename at last '/' */
    if(argv[1][i] == '/') break;
  for(j=i+1; j<ilen; j++)
    indata1[j-(i+1)] = argv[1][j];
  indata1[j-(i+1)] = '\0';
  
  ilen = strlen(argv[2]);
  for(i=ilen-1; i>=0; i--)
   /* chop filename at last '/' */
    if(argv[2][i] == '/') break;
  for(j=i+1; j<ilen; j++)
    indata2[j-(i+1)] = argv[2][j];
  indata2[j-(i+1)] = '\0';
  sprintf(outdata,"CroCo-%d_%s_%s",NGRID,indata1,indata2);
  fprintf(stderr,"\n-> writing result to file %s\n\n",outdata);
  
  if((fpout=fopen(outdata,"w"))==NULL)
   {
    printf("I cannot open %s\n", outdata);
    exit(1);
   }
  
  frac  = (float) atof(argv[4]);
  iseed = (int)   atoi(argv[5]);
  
  /*-------------------------------------------------------
   * get the two relevant particle structure arrays
   *-------------------------------------------------------*/
  NumPart = -1;
  Psimu1  = ReadSimuAndCalcPdens(argv[1], NGRID);
  Psimu2  = ReadSimuAndCalcPdens(argv[2], NGRID);
  
  /*-------------------------------------------------------
   * sorting with respects to particle id
   *-------------------------------------------------------*/
  fprintf(stderr,"  o qsorting Psimu1[] structures according to id ... ");
  qsort((void *)Psimu1, NumPart, sizeof(part), qcompareIDs);
  fprintf(stderr,"done\n");
  fprintf(stderr,"  o qsorting Psimu2[] structures according to id ... ");
  qsort((void *)Psimu2, NumPart, sizeof(part), qcompareIDs);
  fprintf(stderr,"done\n");

//  dump_parts(Psimu1,"Psimu1.dat");
//  dump_parts(Psimu2,"Psimu2.dat");
//  exit(0);

  /*-------------------------------------------------------
   * locate all particles from simu1 in simu2
   *-------------------------------------------------------*/
  L = (int)(pow((double)NumPart,0.333333333333333333333)+0.5);
  
  fprintf(fpout,"# dist12(1) dist1L(2) dist2L(3) dens1(4) dens2(5)\n");
  for(ipart1=0; ipart1<NumPart; ipart1++) {
    
    // id1
    id1 = Psimu1[ipart1].id;
    
    // find position of particle with id1 in array Psimu2[]
    ipart2 = (long unsigned) ((partptr)bsearch(&(id1), Psimu2, NumPart, sizeof(part), bcompareIDs) - Psimu2);
    
    // id2
    id2 = Psimu2[ipart2].id;
    
    // distance between both paticles
    dx = fabs(Psimu1[ipart1].pos[X]-Psimu2[ipart2].pos[X]);
    dy = fabs(Psimu1[ipart1].pos[Y]-Psimu2[ipart2].pos[Y]);
    dz = fabs(Psimu1[ipart1].pos[Z]-Psimu2[ipart2].pos[Z]);
    if(dx > BoxSize/2) dx=BoxSize-dx;
    if(dy > BoxSize/2) dy=BoxSize-dy;
    if(dz > BoxSize/2) dz=BoxSize-dz;
    d12  = sqrt(pow2(dx)+pow2(dy)+pow2(dz));
    
    // distance of ipart1 to its Lagrangian position
    get_ijk(id1, L, &iL, &jL, &kL);
    dx = fabs(Psimu1[ipart1].pos[X]-((double)(iL)/(double)L)*BoxSize);
    dy = fabs(Psimu1[ipart1].pos[Y]-((double)(jL)/(double)L)*BoxSize);
    dz = fabs(Psimu1[ipart1].pos[Z]-((double)(kL)/(double)L)*BoxSize);
    if(dx > BoxSize/2) dx=BoxSize-dx;
    if(dy > BoxSize/2) dy=BoxSize-dy;
    if(dz > BoxSize/2) dz=BoxSize-dz;
    d1L  = sqrt(pow2(dx)+pow2(dy)+pow2(dz));
    
    // distance of ipart2 to its Lagrangian position
    get_ijk(id2, L, &iL, &jL, &kL);
    dx = fabs(Psimu2[ipart2].pos[X]-((double)(iL)/(double)L)*BoxSize);
    dy = fabs(Psimu2[ipart2].pos[Y]-((double)(jL)/(double)L)*BoxSize);
    dz = fabs(Psimu2[ipart2].pos[Z]-((double)(kL)/(double)L)*BoxSize);
    if(dx > BoxSize/2) dx=BoxSize-dx;
    if(dy > BoxSize/2) dy=BoxSize-dy;
    if(dz > BoxSize/2) dz=BoxSize-dz;
    d2L  = sqrt(pow2(dx)+pow2(dy)+pow2(dz));
    
    
    fprintf(fpout,"%g %g %g %g %g\n",d12,d1L,d2L,Psimu1[ipart1].mom[0],Psimu2[ipart2].mom[0]);
    //fprintf(fpout,"%ld %d %d %d %lf %lf %lf\n",id1,iL,jL,kL,((double)(iL)/(double)L)*BoxSize,((double)(jL)/(double)L)*BoxSize,((double)(kL)/(double)L)*BoxSize);
    fflush(fpout);
  }
  
  /*-------------------------------------------------------
   * finish
   *-------------------------------------------------------*/
  fclose(fpout);
  free(Psimu1);
  free(Psimu2);
}

/*=============================================================================
 *
 *============================================================================*/
partptr ReadSimuAndCalcPdens(char *infile, int NGRID)
{
  double min_weight, max_weight, nvpart;
  long   ipart, npart;
  double Xmin, Ymin, Zmin, Xmax, Ymax, Zmax, MinPos, MaxPos;
  float  *x, *y, *z, *weight, *rdens;
#ifdef WRITE_PDENSFILE
  char   indata[MAXSTRING], outdata[MAXSTRING];
  int    ilen,i,j;
  FILE   *fpout;

  ilen = strlen(infile);
  for(i=ilen-1; i>=0; i--)
   /* chop filename at last '/' */
    if(infile[i] == '/') break;
  for(j=i+1; j<ilen; j++)
    indata[j-(i+1)] = infile[j];
  indata[j-(i+1)] = '\0';
  sprintf(outdata,"Pdens-%d_%s",NGRID,indata);
  fprintf(stderr,"-> writing Pdens file %s\n\n",outdata);

  if((fpout=fopen(outdata,"w"))==NULL)
   {
    printf("I cannot open %s\n", outdata);
    exit(1);
   }
#endif /* WRITE_PDENSFILE */
  
  /*-------------------------------------------------------
   * read data
   *-------------------------------------------------------*/
  fprintf(stderr,"o reading %s\n",infile);
  input(infile);
  
  BoxSize    = io.header.boxsize;
  npart      = io.header.no_part;
  nvpart     = io.header.no_vpart;
  min_weight = io.header.min_weight;
  max_weight = io.header.max_weight;
  
  if(NumPart != -1 && NumPart != npart) {
    fprintf(stderr,"The two simulations do not contain the same number of particles: %ld vs. %ld.\nABORTING\n",NumPart,npart);
    exit(0);
   }
  else {
    NumPart = npart;
  }
  
  /*-------------------------------------------------------
   * be verbose
   *-------------------------------------------------------*/
  fprintf(stderr,"Simulation #1:\n");
  fprintf(stderr,"==============\n");
  fprintf(stderr,"Boxsize    = %g\n",BoxSize);
  fprintf(stderr,"no_part    = %ld\n",npart);
  fprintf(stderr,"no_vpart   = %g\n",nvpart);
  fprintf(stderr,"min_weight = %g\n",min_weight);
  fprintf(stderr,"max_weight = %g\n\n",max_weight);
  
  Xmin     = 2*BoxSize;
  Ymin     = 2*BoxSize;
  Zmin     = 2*BoxSize;
  Xmax     = -1;
  Ymax     = -1;
  Zmax     = -1;
  
  /*-------------------------------------------------------
   * copy over to temporary arrays (for legacy reasons)
   *-------------------------------------------------------*/
  x       = (float *)  calloc(npart, sizeof(float));
  y       = (float *)  calloc(npart, sizeof(float));
  z       = (float *)  calloc(npart, sizeof(float));
  weight  = (float *)  calloc(npart, sizeof(float));
  rdens   = (float *)  calloc(npart+1, sizeof(float));
  for(ipart = 0; ipart < npart; ipart++)
   {
    x[ipart]      = fmod((io.fst_part+ipart)->pos[X]+1.0, 1.0);
    y[ipart]      = fmod((io.fst_part+ipart)->pos[Y]+1.0, 1.0);
    z[ipart]      = fmod((io.fst_part+ipart)->pos[Z]+1.0, 1.0);
#ifdef MULTIMASS
    weight[ipart] = (io.fst_part+ipart)->weight;
#else
    weight[ipart] = 1.0;
#endif
    
    if(MinMaxScale)
     {
      if(x[ipart] < Xmin) Xmin=x[ipart];
      if(y[ipart] < Ymin) Ymin=y[ipart];
      if(z[ipart] < Zmin) Zmin=z[ipart];
      if(x[ipart] > Xmax) Xmax=x[ipart];
      if(y[ipart] > Ymax) Ymax=y[ipart];
      if(z[ipart] > Zmax) Zmax=z[ipart];
     }
   }
  
  if(MinMaxScale)
   {
    fprintf(stderr,"  o re-scaling particles into Min-Max cube...");
    MinPos = MIN(Xmin,Ymin);
    MinPos = MIN(MinPos,Zmin);
    MaxPos = MAX(Xmax,Ymax);
    MaxPos = MAX(MaxPos,Zmax);
    for(ipart = 0; ipart < npart; ipart++)
     {
      x[ipart] = fmod(((x[ipart]-MinPos)/(MaxPos-MinPos))+1.0, 1.0);
      y[ipart] = fmod(((y[ipart]-MinPos)/(MaxPos-MinPos))+1.0, 1.0);
      z[ipart] = fmod(((z[ipart]-MinPos)/(MaxPos-MinPos))+1.0, 1.0);
     }
    fprintf(stderr,"done\n");
   }
  
  /*-------------------------------------------------------
   * get density at particle positions
   *-------------------------------------------------------*/
  get_density(npart,nvpart,x,y,z,weight,rdens,NGRID);
  
  /*-------------------------------------------------------
   * scale coordinates to physical coordinates
   *-------------------------------------------------------*/
  for(ipart = 0; ipart < npart; ipart++)
   {
    (io.fst_part+ipart)->pos[X]   = ((io.fst_part+ipart)->pos[X])*BoxSize;
    (io.fst_part+ipart)->pos[Y]   = ((io.fst_part+ipart)->pos[Y])*BoxSize;
    (io.fst_part+ipart)->pos[Z]   = ((io.fst_part+ipart)->pos[Z])*BoxSize;
    
    /* use not needed mom storage for the density at particle position */
    (io.fst_part+ipart)->mom[0]   = rdens[ipart];
   }
  
  /* free memory */
  if(x)      free(x);
  if(y)      free(y);
  if(z)      free(z);
  if(weight) free(weight);
  if(rdens)  free(rdens);
  
#ifdef WRITE_PDENSFILE
  for(ipart = 0; ipart < npart; ipart++) {
    if(ran3(&iseed) <= frac) {
      fprintf(fpout,"%g %g %g  %g  %ld\n",
              (io.fst_part+ipart)->pos[X],
              (io.fst_part+ipart)->pos[Y],
              (io.fst_part+ipart)->pos[Z],
              (io.fst_part+ipart)->mom[0],
              (io.fst_part+ipart)->id);
      fflush(fpout);
    }
  }
  fclose(fpout);
#endif /* WRITE_PDENSFILE */

  /* keep particles though */
  return(io.fst_part);
}


/*=============================================================================
 * get_density.c:   calculate the density at each particle position
 *                  using a grid with NGRID^3 grid points
 *
 *                  particle positions are assumed to lie in
 *                  the range [0,1] !
 *
 *
 *     INPUT:    npart, x[npart], y[npart], z[npart], NGRID
 *
 *     OUTPUT:   rdens[npart]
 *
 *============================================================================*/
void get_density(long unsigned npart, double nvpart, float *x, float *y, float *z, float *weight,
                 float *rdens, int NGRID)
{
  double rho_mean;
  int   ix, iy, iz;
  float *dens;
  
  dens = (float *) calloc(NGRID*NGRID*NGRID, sizeof(float));
  
  rho_mean = nvpart / (double) (NGRID*NGRID*NGRID);
  
  /* assign particles to the grid */
#ifdef NPARTDENS
  fprintf(stderr,"  o assigning particles to grid (number density)...");
#else
  fprintf(stderr,"  o assigning particles to grid (mass density)...");
#endif
  assign(npart, x, y, z, weight, NGRID, dens);
  fprintf(stderr,"done\n");
  
  /* get density in terms of mean density */
  for(ix=0; ix<NGRID; ix++)
    for(iy=0; iy<NGRID; iy++)
      for(iz=0; iz<NGRID; iz++)
        dens[IDX(ix,iy,iz)] /= rho_mean;
  
  /* interpolate density contrast back to particle positions */
  fprintf(stderr,"  o interpolating to particle positions...");
  interp(npart, x, y, z, NGRID, dens, rdens);
  fprintf(stderr,"done\n");
  free(dens);
}

/*==================================================================
 * assign:   assign particles to the grid points using a TSC scheme
 *==================================================================*/
void assign(long unsigned npart, float *x, float *y, float *z, float *weight,
            int NGRID, float *dens)
{
  int           ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1;
  long unsigned ipart;
  float         rrx, rry, rrz, www;
  float         hx, hy, hz;
  float         hx0, hy0, hz0, hxp1, hyp1, hzp1, hxm1, hym1, hzm1;
  
  for(ix=0; ix<NGRID; ix++)
    for(iy=0; iy<NGRID; iy++)
      for(iz=0; iz<NGRID; iz++)
        dens[IDX(ix,iy,iz)]=0.0;
  
  for(ipart=0; ipart<npart; ipart++)
   {
    /* coordinates in grid units */
    rrx = x[ipart] * (float)(NGRID);
    rry = y[ipart] * (float)(NGRID);
    rrz = z[ipart] * (float)(NGRID);
    
    /* particle weight */
#ifdef NPARTDENS
    www = 1.0;
#else
    www = weight[ipart];
#endif
    
    /* index of nearest grid point */
    ix  = (int)(rrx+0.5);
    iy  = (int)(rry+0.5);
    iz  = (int)(rrz+0.5);
    
    /* distance to nearest grid point */
    hx  = rrx - (float)ix;
    hy  = rry - (float)iy;
    hz  = rrz - (float)iz;
    
    /* keep track of peridoc boundaries */
    ix=(int)fmod(ix,NGRID);
    iy=(int)fmod(iy,NGRID);
    iz=(int)fmod(iz,NGRID);
    
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
    ixp1=(int)fmod(ix+1,NGRID);
    iyp1=(int)fmod(iy+1,NGRID);
    izp1=(int)fmod(iz+1,NGRID);
    ixm1=(int)fmod(ix-1+NGRID,NGRID);
    iym1=(int)fmod(iy-1+NGRID,NGRID);
    izm1=(int)fmod(iz-1+NGRID,NGRID);
    
    /* assign particle according to weights to 27 neighboring nodes */
    dens[IDX(ixm1,iym1,izm1)] += hxm1*hym1 *hzm1 * www;
    dens[IDX(ix,  iym1,izm1)] += hx0 *hym1 *hzm1 * www;
    dens[IDX(ixp1,iym1,izm1)] += hxp1*hym1 *hzm1 * www;
    dens[IDX(ixm1,  iy,izm1)] += hxm1*hy0  *hzm1 * www;
    dens[IDX(  ix,  iy,izm1)] += hx0 *hy0  *hzm1 * www;
    dens[IDX(ixp1,  iy,izm1)] += hxp1*hy0  *hzm1 * www;
    dens[IDX(ixm1,iyp1,izm1)] += hxm1*hyp1 *hzm1 * www;
    dens[IDX(  ix,iyp1,izm1)] += hx0 *hyp1 *hzm1 * www;
    dens[IDX(ixp1,iyp1,izm1)] += hxp1*hyp1 *hzm1 * www;
    dens[IDX(ixm1,iym1,  iz)] += hxm1*hym1 *hz0 * www;
    dens[IDX(  ix,iym1,  iz)] += hx0 *hym1 *hz0 * www;
    dens[IDX(ixp1,iym1,  iz)] += hxp1*hym1 *hz0 * www;
    dens[IDX(ixm1,  iy,  iz)] += hxm1*hy0  *hz0 * www;
    dens[IDX(  ix,  iy,  iz)] += hx0 *hy0  *hz0 * www;
    dens[IDX(ixp1,  iy,  iz)] += hxp1*hy0  *hz0 * www;
    dens[IDX(ixm1,iyp1,  iz)] += hxm1*hyp1 *hz0 * www;
    dens[IDX(  ix,iyp1,  iz)] += hx0 *hyp1 *hz0 * www;
    dens[IDX(ixp1,iyp1,  iz)] += hxp1*hyp1 *hz0 * www;
    dens[IDX(ixm1,iym1,izp1)] += hxm1*hym1 *hzp1 * www;
    dens[IDX(  ix,iym1,izp1)] += hx0 *hym1 *hzp1 * www;
    dens[IDX(ixp1,iym1,izp1)] += hxp1*hym1 *hzp1 * www;
    dens[IDX(ixm1,  iy,izp1)] += hxm1*hy0  *hzp1 * www;
    dens[IDX(  ix,  iy,izp1)] += hx0 *hy0  *hzp1 * www;
    dens[IDX(ixp1,  iy,izp1)] += hxp1*hy0  *hzp1 * www;
    dens[IDX(ixm1,iyp1,izp1)] += hxm1*hyp1 *hzp1 * www;
    dens[IDX(  ix,iyp1,izp1)] += hx0 *hyp1 *hzp1 * www;
    dens[IDX(ixp1,iyp1,izp1)] += hxp1*hyp1 *hzp1 * www;
   }
}

/*==============================================================
 * interp:  interpolate the density the grid points back to the
 *          particle positions using again a TSC scheme
 *==============================================================*/
void interp(long unsigned npart, float *x, float *y, float *z,
            int NGRID, float *dens, float *rdens)
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
    rrx = x[ipart] * (float)(NGRID);
    rry = y[ipart] * (float)(NGRID);
    rrz = z[ipart] * (float)(NGRID);
    
    /* index of nearest grid point */
    ix  = (int)(rrx+0.5);
    iy  = (int)(rry+0.5);
    iz  = (int)(rrz+0.5);
    
    /* distance to nearest grid point */
    hx  = rrx - (float)ix;
    hy  = rry - (float)iy;
    hz  = rrz - (float)iz;
    
    /* keep track of peridoc boundaries */
    ix=(int)fmod(ix,NGRID);
    iy=(int)fmod(iy,NGRID);
    iz=(int)fmod(iz,NGRID);
    
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
    ixp1=(int)fmod(ix+1,NGRID);
    iyp1=(int)fmod(iy+1,NGRID);
    izp1=(int)fmod(iz+1,NGRID);
    ixm1=(int)fmod(ix-1+NGRID,NGRID);
    iym1=(int)fmod(iy-1+NGRID,NGRID);
    izm1=(int)fmod(iz-1+NGRID,NGRID);
    
    
    ac =dens[IDX(ixm1,iym1,izm1)]* hxm1*hym1 *hzm1
    +dens[IDX(ix  ,iym1,izm1)]* hx0 *hym1 *hzm1
    +dens[IDX(ixp1,iym1,izm1)]* hxp1*hym1 *hzm1
    +dens[IDX(ixm1,iy  ,izm1)]* hxm1*hy0  *hzm1
    +dens[IDX(ix  ,iy  ,izm1)]* hx0 *hy0  *hzm1
    +dens[IDX(ixp1,iy  ,izm1)]* hxp1*hy0  *hzm1
    +dens[IDX(ixm1,iyp1,izm1)]* hxm1*hyp1 *hzm1
    +dens[IDX(ix  ,iyp1,izm1)]* hx0 *hyp1 *hzm1
    +dens[IDX(ixp1,iyp1,izm1)]* hxp1*hyp1 *hzm1
    +dens[IDX(ixm1,iym1,iz  )]* hxm1*hym1 *hz0
    +dens[IDX(ix  ,iym1,iz  )]* hx0 *hym1 *hz0
    +dens[IDX(ixp1,iym1,iz  )]* hxp1*hym1 *hz0
    +dens[IDX(ixm1,iy  ,iz  )]* hxm1*hy0  *hz0
    +dens[IDX(ix  ,iy  ,iz  )]* hx0 *hy0  *hz0;
    ac+=dens[IDX(ixp1,iy  ,iz  )]* hxp1*hy0  *hz0
    +dens[IDX(ixm1,iyp1,iz  )]* hxm1*hyp1 *hz0
    +dens[IDX(ix  ,iyp1,iz  )]* hx0 *hyp1 *hz0
    +dens[IDX(ixp1,iyp1,iz  )]* hxp1*hyp1 *hz0
    +dens[IDX(ixm1,iym1,izp1)]* hxm1*hym1 *hzp1
    +dens[IDX(ix  ,iym1,izp1)]* hx0 *hym1 *hzp1
    +dens[IDX(ixp1,iym1,izp1)]* hxp1*hym1 *hzp1
    +dens[IDX(ixm1,iy  ,izp1)]* hxm1*hy0  *hzp1
    +dens[IDX(ix  ,iy  ,izp1)]* hx0 *hy0  *hzp1
    +dens[IDX(ixp1,iy  ,izp1)]* hxp1*hy0  *hzp1
    +dens[IDX(ixm1,iyp1,izp1)]* hxm1*hyp1 *hzp1
    +dens[IDX(ix  ,iyp1,izp1)]* hx0 *hyp1 *hzp1
    +dens[IDX(ixp1,iyp1,izp1)]* hxp1*hyp1 *hzp1;
    rdens[ipart] = (float) ac;
    
   }
}


void check_density(int NGRID, float *dens)
{
  int ix, iy, iz;
  
  for(iz=0; iz<NGRID; iz++)
    for(iy=0; iy<NGRID; iy++)
      for(ix=0; ix<NGRID; ix++)
       {
        if(dens[IDX(ix,iy,iz)] < 1.0)
          fprintf(stderr,"underdens cell: %12d %12d %12d = %16.8g\n",ix,iy,iz,dens[IDX(ix,iy,iz)]);
       }
}


/*==============================================================================
 *  compare particle ids (used with qsort and bsearch)
 *==============================================================================*/
int qcompareIDs(const void *p1, const void *p2)
{
	long unsigned id1, id2;
  
	id1 = ((partptr)p1)->id;
	id2 = ((partptr)p2)->id;
  
	return id1 < id2 ? -1 : (id1 > id2 ? 1 : 0);
}

int bcompareIDs(const void *id, const void *p)
{
	const long *i    = (const long *)id;
  const long ipart = ((partptr)p)->id;
  
  return *i < ipart ? -1 : (*i > ipart ? 1 : 0);
}

/*==============================================================================
 *  dump p[] array to file for debugging purposes
 *==============================================================================*/
void dump_parts(partptr p, char *out)
{
  FILE *fpout;
  long i;
  
  fpout = fopen(out,"w");
  for(i=0; i<NumPart; i++) {
    fprintf(fpout,"%g %g %g %ld\n",p[i].pos[X],p[i].pos[Y],p[i].pos[Z],p[i].id);
  }
  
  fclose(fpout);
}

/*==============================================================================
 *  find Lagrangian position of particle with id
 *  -> this routine heavily assumes that NumPart is a power of 2!
 *==============================================================================*/
void get_ijk(long unsigned id, int L, int *iL, int *jL, int *kL)
{
  int i,j,k;
  long unsigned id_tmp;
  
  // get i
  i = (int) ((double)(id-1)/(double)(L*L));
  
  // get j
  id_tmp = (id-1) - i*L*L;
  j = (int) ((double)id_tmp/(double)L);
  
  // get k
  k = id_tmp - j*L;
  
  *iL = i;
  *jL = j;
  *kL = k;
}



