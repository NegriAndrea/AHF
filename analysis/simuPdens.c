#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libio_serial/io_serial.h"
#include "../src/libutility/utility.h"

// DEFINE FEATURES
//=================
//#define STEREO2              /* write output file compatible with stereo2 */
//#define NPARTDENS            /* number density rather than mass density */
//#define AHF_HALOS            /* work on simulation data or halo catalogues */
//#define MMFOCUS              /* only dump high-resolution particles */
//#define WRITE_GRID_DENSITY   /* also write a file that contains the density on the grid */

//#define RFOCUS
#define Xc 23.32017853
#define Yc 33.24143512
#define Zc 27.43802337
#define Rc 0.3

#define IDX(ix,iy,iz)  ((iz)*NGRID*NGRID + (iy)*NGRID + (ix))

// a negative NGRID will place the particles into the largest possible cube defined by the positions itself
int    MinMaxScale;
double MinPos, MaxPos;
char   indata[MAXSTRING];

void get_density(long unsigned npart, double nvpart, float *x, float *y, float *z, float *weight,
                 double *rdens, int NGRID, float *dens);

void check_density(int NGRID, float *dens);

void assign(long unsigned npart, float *x, float *y, float *z, float *weight,
            int NGRID, float *dens);

void interp(long unsigned npart, float *x, float *y, float *z,  
            int NGRID, float *dens, double *rdens);

long read_halos(char *infile, double BoxSize);
int  check_rfocus(float x, float y, float z, float BoxSize);

int main(argc,argv)
int argc;
char **argv;
{
  char          outdata[MAXSTRING], cdummy[20];
  FILE          *fpin, *fpout;
  long unsigned npart, ipart;
  double        nvpart;
  float         *x, *y, *z, *weight, *dens, min_weight, max_weight;
  float         fw[7];
  double        dw[7];
  double        *rdens, logdens;
  long unsigned *index;
  int           NGRID, ilen, i, j, idummy;
  float         BoxSize, dummy;
  float         MaxDens, MinDens, CurDens, r_color, g_color, b_color;
  int           BUFSIZE;
  float         frac;
  int           iseed;
  int           SWAPBYTES;
  long          nhalos;
  double        Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
  
#ifdef AHF_HALOS
  if(argc<4)
   {
    fprintf(stderr,"usage: %s *.AHF_halos NGRID BoxSize\n",*argv);
    exit(1);
   }
  
  fprintf(stderr,"=====================================================================\n");
  fprintf(stderr,"  Read AHF catalogues and calculate density at each halo position\n");
  fprintf(stderr,"           (only write a random fraction of the halos)\n");
  fprintf(stderr,"=====================================================================\n");
#else
  if(argc<5)
   {
    fprintf(stderr,"usage: %s simu_file NGRID fraction iseed\n",*argv);
    exit(1);
   }
  
  fprintf(stderr,"==========================================================================\n");
  fprintf(stderr,"  Read simulation binary and calculate density at each particle position\n");
  fprintf(stderr,"             (only write a random fraction of the particles)\n");
  fprintf(stderr,"==========================================================================\n");
#endif
  
  NGRID = (int) (atoi(argv[2]));
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
    indata[j-(i+1)] = argv[1][j];
  indata[j-(i+1)] = '\0';
#ifdef STEREO2
  sprintf(outdata,"Pdens-%d_%s.geom",NGRID,indata);
#else
  sprintf(outdata,"Pdens-%d_%s",NGRID,indata);
#endif
  fprintf(stderr,"\n-> writing result to file %s\n\n",outdata);
  
  if((fpout=fopen(outdata,"w"))==NULL)
   {
    printf("I cannot open %s\n", outdata);
    exit(1);
   }
  
#ifdef AHF_HALOS
  frac       = 1.;
  iseed      = 0;
  BoxSize    = (double) atof(argv[3]);
  
  nhalos = read_halos(argv[1], BoxSize);
  
  npart      = nhalos;
  nvpart     = io.header.no_vpart;
  min_weight = io.header.min_weight;
  max_weight = io.header.max_weight;
#else
  
  frac  = (float) atof(argv[3]);
  iseed = (int)   atoi(argv[4]);
  
  /*-------------------------------------------------------
   * read data
   *-------------------------------------------------------*/
  input(argv[1]);
  
  BoxSize    = io.header.boxsize;
  npart      = io.header.no_part;
  nvpart     = io.header.no_vpart;
  min_weight = io.header.min_weight;
  max_weight = io.header.max_weight;
#endif
  
  fprintf(stderr,"Boxsize    = %g\n",BoxSize);
  fprintf(stderr,"no_part    = %ld\n",npart);
  fprintf(stderr,"no_vpart   = %g\n",nvpart);
  fprintf(stderr,"min_weight = %g\n",min_weight);
  fprintf(stderr,"max_weight = %g\n\n",max_weight);
  x       = (float *)  calloc(npart, sizeof(float));
  y       = (float *)  calloc(npart, sizeof(float));
  z       = (float *)  calloc(npart, sizeof(float));
  weight  = (float *)  calloc(npart, sizeof(float));
  dens    = (float *)  calloc(NGRID*NGRID*NGRID, sizeof(float));
  rdens   = (double *) calloc(npart+1, sizeof(double));
  index   = (long unsigned *) calloc(npart+1, sizeof(long unsigned));
  Xmin    = 2*BoxSize;
  Ymin    = 2*BoxSize;
  Zmin    = 2*BoxSize;
  Xmax    = -1;
  Ymax    = -1;
  Zmax    = -1;
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
  
  /* get density at particle positions */
  get_density(npart,nvpart,x,y,z,weight,rdens,NGRID,dens);
  
  /* look for underdense regions... */
  //check_density(NGRID, dens);
  
  /* sort rdens according to density value */
  fprintf(stderr,"  o sorting particles...");
  indexx(npart,rdens-1,index-1); /* "-1" because indexx returns index in range [1,npart] */
  fprintf(stderr,"done\n");
  
  /* scale coordinates to physical coordinates */
  for(ipart = 0; ipart < npart; ipart++)
   {
    *(x+ipart)   = ((io.fst_part+ipart)->pos[X])*BoxSize;
    *(y+ipart)   = ((io.fst_part+ipart)->pos[Y])*BoxSize;
    *(z+ipart)   = ((io.fst_part+ipart)->pos[Z])*BoxSize;
    index[ipart] = index[ipart] - 1;    /* indexx returns index in range [1,npart] */
   }
  
  /* write DATA to outfile */
  MinDens = *(rdens+index[0]);
  MaxDens = *(rdens+index[npart-1]);
  fprintf(stderr,"  o writing output (minDens=%g maxDens=%g)...",MinDens,MaxDens);
  for(ipart = 0; ipart < npart; ipart++)
   {
    CurDens = (*(rdens+index[ipart]));
    
#ifdef STEREO2
    
    logdens = log10(CurDens-MinDens+1)/log10(MaxDens-MinDens+1.);
    
    r_color =    pow(logdens,0.9);
    g_color =    pow(logdens,1.5);
    b_color = 1.-pow(logdens,0.2);
    
#ifdef MMFOCUS
    if(fabs(*(weight+index[ipart])-min_weight) < ZERO) 
#endif
#ifdef RFOCUS
    if(check_rfocus(*(x+index[ipart]),*(y+index[ipart]),*(z+index[ipart]),BoxSize))
#endif
     {
      if(ran3(&iseed) < frac)
        fprintf(fpout,"%s %g %g %g %g %g %g   %ld\n",
                "p ",
                *(x+index[ipart]),
                *(y+index[ipart]),
                *(z+index[ipart]),
                r_color,g_color,b_color,
                index[ipart]);
     }
#else /* STEREO2 */
	  
    /* only dump high-res particles */
#ifdef MMFOCUS
    if(fabs(*(weight+index[ipart])-min_weight) < ZERO) 
#endif
#ifdef RFOCUS
    if(check_rfocus(*(x+index[ipart]),*(y+index[ipart]),*(z+index[ipart]),BoxSize))
#endif
     {
      if(ran3(&iseed) <= frac)
        fprintf(fpout,"%g %g %g %g   %ld\n",
                *(x+index[ipart]),
                *(y+index[ipart]),
                *(z+index[ipart]),
                CurDens,
                index[ipart]);
     }
#endif
    fflush(fpout);
   }
  fprintf(stderr,"done\n");
  
  /* close files */
  fclose(fpout);
  
  /* free memory */
  if(x)      free(x);
  if(y)      free(y);
  if(z)      free(z);
  if(weight) free(weight);
  if(dens)   free(dens);
  if(rdens)  free(rdens);
  if(index)  free(index);
}

/*=============================================================================
* get_density.c:   calculate the density at each particle position
*                  using a grid with NGRID^3 grid points
*
*                  particle positions are assumed to lie in
*                  the range [0,1] !
*
* 
*     INPUT:    npart, x[npart], y[npart], z[npart], L1DIM
*
*     OUTPUT:   rdens[npart]
*
*============================================================================*/
void get_density(long unsigned npart, double nvpart, float *x, float *y, float *z, float *weight,
                 double *rdens, int NGRID, float *dens)
{
   double rho_mean;
   int   ix, iy, iz;
#ifdef WRITE_GRID_DENSITY
  char outfile[MAXSTRING];
  FILE *fpout;
#endif
   
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
   
#ifdef WRITE_GRID_DENSITY
  sprintf(outfile,"Ldens-%d_%s",NGRID,indata);
  fpout = fopen(outfile,"w");
  for(ix=0; ix<NGRID; ix++)
    for(iy=0; iy<NGRID; iy++)
      for(iz=0; iz<NGRID; iz++)
        fprintf(fpout,"%16.8f %16.8f %16.8f %16.8f\n",
                ((double)ix+0.5)/(double)NGRID * io.header.boxsize,
                ((double)iy+0.5)/(double)NGRID * io.header.boxsize,
                ((double)iz+0.5)/(double)NGRID * io.header.boxsize,
                dens[IDX(ix,iy,iz)]);
  fclose(fpout);
#endif
  
  
   /* interpolate density contrast back to particle positions */ 
   fprintf(stderr,"  o interpolating to particle positions...");
   interp(npart, x, y, z, NGRID, dens, rdens);
   fprintf(stderr,"done\n");
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
            int NGRID, float *dens, double *rdens)
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
      rdens[ipart] = (double) ac;
      
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

/*======================================================================
 * read an AHF halos catalogue
 *======================================================================*/
long read_halos(char *infile, double BoxSize)
{
  long nhalos, ihalo, idummy;
  FILE *fpin;
  char  instream[MAXSTRING];
  long npart;
  double nvpart, xpos, ypos, zpos;
  
  if( (fpin=fopen(infile,"r")) == NULL )
   {
    fprintf(stderr,"Could not open %s\nAborting\n",infile);
    exit(0);
   }
  
  nhalos = 0;
  fgets(instream,MAXSTRING,fpin);
  while(!feof(fpin))
   {
    nhalos++;
    fgets(instream,MAXSTRING,fpin);
   }
  nhalos--;
  rewind(fpin);
  fgets(instream,MAXSTRING,fpin);

#ifdef MULTIMASS
  io.header.min_weight =  1E40;
  io.header.max_weight = -1E40;
  io.header.no_vpart   = 0.;
#else
  io.header.min_weight = 1.;
  io.header.max_weight = 1.;
  io.header.no_vpart   = (double)nhalos;
#endif
  io.fst_part = (partptr) calloc(nhalos,sizeof(part));
  for(ihalo=0; ihalo<nhalos; ihalo++)
   {
    fgets(instream,MAXSTRING,fpin);
    sscanf(instream,"%ld %ld %ld %lf %ld %lf %lf %lf", &idummy, &idumy, &idummy, &nvpart, &npart, &xpos, &ypos, &zpos);
    
    (io.fst_part+ihalo)->pos[X] = xpos/BoxSize;
    (io.fst_part+ihalo)->pos[Y] = ypos/BoxSize;
    (io.fst_part+ihalo)->pos[Z] = zpos/BoxSize;
#ifdef MULTIMASS
    (io.fst_part+ihalo)->weight = nvpart;
    if(nvpart > io.header.max_weight) io.header.max_weight = nvpart;
    if(nvpart < io.header.min_weight) io.header.min_weight = nvpart;
    io.header.no_vpart += nvpart;
#endif
   }
  
  fprintf(stderr,"  o successfully read %ld halos\n\n",nhalos);

  return(nhalos);
}

/*=======================================================================
 * check for particle position to lie in region centred about Xc, Yc, Zc
 *=======================================================================*/
int check_rfocus(float x, float y, float z, float BoxSize)
{
  double dx, dy, dz;
  
  dx = fabs(x-Xc);
  dy = fabs(y-Yc);
  dz = fabs(z-Zc);
  
  if(dx > BoxSize/2.) dx=dx-BoxSize;
  if(dy > BoxSize/2.) dy=dy-BoxSize;
  if(dz > BoxSize/2.) dz=dz-BoxSize;
  
  if((pow2(dx)+pow2(dy)+pow2(dz)) < pow2(Rc))
    return(1);
  else
    return(0);
}

