#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libio_serial/io_serial.h"
#include "../src/libutility/utility.h"

/*==============================
 * switch on/off various things
 *==============================*/
#define USED -10


/*===================
 * COMMON structures
 *===================*/
info_io   io;

int           NBINS; /* number of BINS               for correlation function */
long unsigned NSUB;  /* size of subsample to be used for correlation function */



/*===================
 * functions
 *===================*/
void statistic(double XA, double XE, int ngrid, FILE *fpout);


/*===================
 * main()
 *===================*/
int main(argc,argv)
int argc;
char **argv;
{
  long unsigned i,j,k,ipart;
  int           no_timestep, ngrid;
  double        RA,RE;
  double        x_fac, v_fac;
  char          indata[MAXSTRING], outdata[MAXSTRING], header[HEADERSTRING];
  FILE         *fpout;
  
  printf("=============================================\n");
  printf("  Read simulation binary and calculate Xi(r)\n");
  printf("=============================================\n");
  printf("Please give name of input  file:                        ");
  scanf("%s", indata);
  printf("%s\n",indata);
  printf("Please give NGRID_DOM for linked-list:                  ");
  scanf("%d", &ngrid);
  printf("%d\n",ngrid);
  printf("Please give number of bins:                             ");
  scanf("%d", &NBINS);
  printf("%d\n",NBINS);
  printf("Please give number of particles to be used (-1: all):   ");
  scanf("%d", &NSUB);
  printf("%d\n",NSUB);
  printf("Please give r_min [Mpc/h]:                              ");
  scanf("%lf", &RA);
  printf("%lf\n",RA);
  printf("Please give r_max [Mpc/h:                               ");
  scanf("%lf", &RE);
  printf("%lf\n",RE);
  printf("\nPlease give name for output file:                     ");
  scanf("%s", outdata);
  printf("%s\n",outdata);
  printf("\n");

  if((fpout=fopen(outdata,"w"))==NULL)
    {
     printf("I cannot open output file %s\n", outdata);
     exit(1);
    }
  
  /*=================
   * read input file
   *=================*/
  input(indata);
  
  RA /= io.header.boxsize;
  RE /= io.header.boxsize;
  
  if(NSUB == -1) NSUB = io.header.no_part;
  
#ifdef HEADER 
  fprintf(fpout,"#%s\n",io.header.header);
  fprintf(fpout,"#no_part     = %ld\n",io.header.no_part);
  fprintf(fpout,"#no_vpart    = %g\n",io.header.no_vpart);
  fprintf(fpout,"#boxsize     = %g\n",io.header.boxsize);
  fprintf(fpout,"#Omega0      = %g\n",io.header.omega0);
  fprintf(fpout,"#Lambda0     = %g\n",io.header.lambda0);
  fprintf(fpout,"#z_initial   = %g (a=%g)\n",
          1.0/io.header.a_initial-1.0,io.header.a_initial);
  fprintf(fpout,"#z_current   = %g (a=%g)\n",
          1.0/io.header.a_current-1.0,io.header.a_current);
  fprintf(fpout,"#no_timestep = %d\n",io.header.no_timestep);
  fprintf(fpout,"#NGRID       = %d\n",ngrid);
  fprintf(fpout,"#NBINS       = %d\n",NBINS);
  fprintf(fpout,"#NSUB        = %d\n",NSUB);
#endif


  /*=======================
   * perform your analysis
   *=======================*/
  statistic(RA,RE,ngrid,fpout);

  fclose(fpout);
}

/*-----------------------------------------------------------------------------
 * randomly choose a subsample of all particles
 *-----------------------------------------------------------------------------*/
int subsample(long unsigned nobj, int *ISEED, long unsigned *isub)
{
   double         random;
   long unsigned i;
   
   for(i=0; i<NSUB; i++)
     {
      random  = ran3(ISEED);
      isub[i] = (long unsigned) (random * ((double)nobj-1.) + 1.);
     }
}

/*-----------------------------------------------------------------------------
 * build initial linked lists...afterwards this list is only updated
 *-----------------------------------------------------------------------------*/
void linklist(long unsigned npart, double *x, double *y, double *z, 
              int ngrid, long unsigned *ihc, long unsigned *nhc, long unsigned *ll) 
{
   
   long unsigned ipart, ngrid2, ngrid3, ibox, inst, ins, insp;
   int           ix, iy, iz;
   double        rfac;
   double         Zpos;
   
   fprintf(stderr," o generating linked list ... ");
   
   ngrid2 = pow2(ngrid);
   ngrid3 = pow3(ngrid);
   
   /* set all values to zero */
   for(ibox=0; ibox<ngrid3; ibox++)
     {
      nhc[ibox] = 0;
      ihc[ibox] = 0;
     }
   for(ipart=0; ipart<npart; ipart++)
     {
      ll[ipart] = 0;
     }
   
   rfac = (double)ngrid;
   
   /* loop over all particles starting with fst_part (first particle) */
   for(ipart = 0; ipart < npart; ipart++)
     {
      Zpos = z[ipart];
      
      /* particle position in linklist grid units */
      ix   = (int)(x[ipart]*rfac);
      iy   = (int)(y[ipart]*rfac);
      iz   = (int)(z[ipart]*rfac);
      
      /* mapping to one-dimensional array index */
      ibox = ix + ngrid*iy + ngrid2*iz;
      
      /* add particle to that cell */
      nhc[ibox] += 1;
      
      inst = ihc[ibox];
      ins  = inst;
      insp = inst;
      
      /* re-order particles in ascending z-coordinate */
      while(ins > 0)
        {
         if(Zpos < z[ins])
           {
            insp = ins;
            ins  = ll[ins];
           }
         else
           {
            break;
           }
        }
      
      /* store particle index */
      ll[ipart] = ins;
      
      if(ins == inst)
         ihc[ibox] = ipart;
      else
         ll[insp]  = ipart;
     }
   
   fprintf(stderr,"done\n");
}



/*-----------------------------------------------------------------------------
* actually calculates the correlation function
*-----------------------------------------------------------------------------*/
int ana(double *x, double *y, double *z, long unsigned npart,
        double RA,double RE, double *CR, int *IC, double *WP, double *WC, double RNORM,
        int ngrid, FILE *fpout)
{
   double        lRA,lRE,ldr,dr,dr2,RE2,dx,dy,dz,dRmin,dRmax,dRmin2,dRmax2,rfac;
   double        RNC;
   long unsigned NR, IR, *nhc, *ihc, *ll, *ipp, cellpart, nn;
   long unsigned ix,iy,iz,imin,imax,jmin,jmax,kmin,kmax,ibox,ngrid2;
   long unsigned i,j,k,iix,iiy,iiz,ibin,ipart,npair,ipair,ipair12,ipair21;
   int           no_use;
   double       *xtmp;
   
   /* temporary storage */
   xtmp  = (double *) calloc(npart, sizeof(double));
   
   /* organize particles in linked-list */
   ihc = (long unsigned *) calloc(pow3(ngrid), sizeof(long unsigned));
   nhc = (long unsigned *) calloc(pow3(ngrid), sizeof(long unsigned));
   ll  = (long unsigned *) calloc(npart,       sizeof(long unsigned));
   linklist(npart,x,y,z,ngrid,ihc,nhc,ll);
   
   /* initialize values */
   lRA = log10(RA);
   lRE = log10(RE);
   RE2 = pow2(RE);
   
   for(NR=0; NR<NBINS; NR++) 
      IC[NR] = 0;
   
   ngrid2 = pow2(ngrid);
   rfac   = (double)ngrid;
   
   fprintf(stderr,"\n o monitoring progress...\n");
   
   /* loop over all BINS */
   for(ibin=0; ibin<NBINS-1; ibin++)
     {
      /* radii covered by current BIN */
      dRmin  = (CR[ibin]);
      dRmax  = (CR[ibin+1]);
      dRmin2 = pow2(CR[ibin]);
      dRmax2 = pow2(CR[ibin+1]);
      
      /* temporarily safe coordinates */
      for(ipart=0; ipart<npart; ipart++)
         xtmp[ipart] = x[ipart];
      
      /* loop over all particles */
      for(ipart=0; ipart<npart; ipart++)
        {
         /* get range for linklist */
         imin  = (int)((x[ipart]-dRmax)*rfac);
         jmin  = (int)((y[ipart]-dRmax)*rfac);
         kmin  = (int)((z[ipart]-dRmax)*rfac);
         imax  = (int)((x[ipart]+dRmax)*rfac)+1;
         jmax  = (int)((y[ipart]+dRmax)*rfac)+1;
         kmax  = (int)((z[ipart]+dRmax)*rfac)+1;
         
         /* loop over all those linklist cells */
         for(i=imin; i<=imax; i++)
            for(j=jmin; j<=jmax; j++)
               for(k=kmin; k<=kmax; k++)
                 {
                  /* take care of periodic boundaries */
                  ix = (int) fmod(i+ngrid,ngrid);
                  iy = (int) fmod(j+ngrid,ngrid);
                  iz = (int) fmod(k+ngrid,ngrid);
                  
                  /* transfer to 1d array index */
                  ibox = ix + ngrid*iy + ngrid2*iz;
                  
                  /* are there particles within that linklist cell */
                  if(nhc[ibox] > 0)
                    {
                     /* index to first particle in that cell */
                     cellpart = ihc[ibox];
                     
                     if(x[cellpart] > USED)
                       {
                        /* distance to actual particle */
                        dx = fabs(x[cellpart] - x[ipart]);
                        dy = fabs(y[cellpart] - y[ipart]);
                        dz = fabs(z[cellpart] - z[ipart]);
                        
                        /* take care of periodic boundaries */
                        if(dx > 0.5) dx = 1 - dx;
                        if(dy > 0.5) dy = 1 - dy;
                        if(dz > 0.5) dz = 1 - dz;
                        
                        dr2 = pow2(dx) + pow2(dy) + pow2(dz);
                        
                        /* particle falls into BIN */
                        if(dRmin2 < dr2 && dr2 < dRmax2)
                           IC[ibin] += 1;
                       }
                     
                     /* remaining particles in that linklist cell */
                     for(nn=1; nn<nhc[ibox]; nn++)
                       {
                        /* run through linklist... */
                        cellpart = ll[cellpart];
                        
                        if(x[cellpart] > USED)
                          {
                           /* distance to actual particle */
                           dx = fabs(x[cellpart] - x[ipart]);
                           dy = fabs(y[cellpart] - y[ipart]);
                           dz = fabs(z[cellpart] - z[ipart]);
                           
                           /* take care of periodic boundaries */
                           if(dx > 0.5) dx = 1 - dx;
                           if(dy > 0.5) dy = 1 - dy;
                           if(dz > 0.5) dz = 1 - dz;
                           
                           dr2 = pow2(dx) + pow2(dy) + pow2(dz);
                           
                           /* particle falls into BIN */
                           if(dRmin2 < dr2 && dr2 < dRmax2)
                              IC[ibin] += 1;
                          }
                       }
                    }
                 }
         x[ipart] = USED;
        }
      
      /* get coordinates back from temporarily storage */
      for(ipart=0; ipart<npart; ipart++)
         x[ipart] = xtmp[ipart];
      
      RNC = (double) IC[ibin];
      WC[ibin] = RNC/(WP[ibin]*RNORM)-1.0;
      fprintf(stderr,"    r = %12.8g , xi = %12.4g  (%12.8g)\n",
              CR[ibin]*io.header.boxsize, WC[ibin], CR[NBINS-1]*io.header.boxsize);
      
      fprintf(fpout,"%g %g %d %g\n",CR[ibin]*io.header.boxsize,WC[ibin],IC[ibin],WP[ibin]);
      fflush(fpout);
     }
}
/*-----------------------------------------------------------------------------
* get proper normalisation of correlation function
*-----------------------------------------------------------------------------*/
int norm(double RA, double RE, double *WP)
{
   double lRA, lRE, VBOX, FAK, VA, CR, CV;
   int    i;
   double boxsize;
   
   boxsize = 1.0; // internal units!
   
   /* initialize values */
   lRA  = log10(RA);
   lRE  = log10(RE);
   VBOX = pow3(boxsize);
   FAK  = 4.0*PI/3.0;
   VA   = FAK*pow3(RA);
   
   
   for(i=0; i<NBINS; i++)
     {
      CR    = (lRA + (i+1) * (lRE-lRA)/(double)NBINS);
      CR    = pow(10.,CR);
      CV    = FAK*pow3(CR);
      WP[i] = (CV-VA)/VBOX;
      VA    = CV;
      RA    = CR;
     }
}


/*===========================================================================
* do you own stuff with simulation data given now in proper units...
*===========================================================================*/
void statistic(double XA, double XE, int ngrid, FILE *fpout)
{
   double         RA, RE, lRA, lRE, RNORM, RNC;
   double         *CR,*WP,*WC;
   int            *IC;
   long unsigned  nobj, nused, i, *isub;
   int           ISEED = -1001;
   double          *xsub, *ysub, *zsub;
   partptr        cur_part;
   
   
   /* initialize some values */
   RA   = (double)XA;
   RE   = (double)XE;
   nobj = io.header.no_part;
   
   lRA = log10(RA);
   lRE = log10(RE);
   
   CR = (double*) calloc(NBINS, sizeof(double));
   WP = (double*) calloc(NBINS, sizeof(double));
   WC = (double*) calloc(NBINS, sizeof(double));
   IC = (int*)    calloc(NBINS, sizeof(int));
   
   
   for(i=0; i<NBINS; i++)
      CR[i] = pow(10.0,(lRA + (i) * (lRE-lRA)/(double)NBINS));

   RNORM = (double)NSUB*(double)(NSUB-1)/2.0;
   
   
   /* create subsample if needed */
   if(NSUB < nobj)
     {
      isub  = calloc(NSUB, sizeof(long unsigned));
      nused = NSUB;
      subsample(nobj, &ISEED, isub);
     }
   else
     {
      isub  = calloc(nobj, sizeof(long unsigned));
      nused = nobj;
      for(i=0; i<nused; i++) *(isub+i)=i;
     }
   
   /* create new arrays to hold subsample */
   xsub = (double *) calloc(nused, sizeof(double));
   ysub = (double *) calloc(nused, sizeof(double));
   zsub = (double *) calloc(nused, sizeof(double));
   
   /* move particles to new arrays */
#ifdef WITH_OPENMP
#pragma omp parallel private(i, cur_part) shared(io, isub, xsub, ysub, zsub, nused)
#pragma omp for schedule(static)
#endif
   for(i=0; i<nused; i++)
     {
      cur_part = io.fst_part + isub[i];
      xsub[i]  = cur_part->pos[X];
      ysub[i]  = cur_part->pos[Y];
      zsub[i]  = cur_part->pos[Z];
     }
   
   /* calculate normalisation for correlation function */
   norm(RA, RE, WP);
   
   /* perform correlation function analysis */
   ana(xsub,ysub,zsub,nused,RA,RE,CR,IC,WP,WC,RNORM,ngrid,fpout);
   
}


