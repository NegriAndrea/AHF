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
//#define PK_GAS_ONLY  // use only the gas density
#define PK_DM_ONLY   // use only the dark matter particles


/*===================
 * COMMON structures
 *===================*/
info_io     io;
info_global global;
gridls     *cur_grid;

long unsigned NGRID, NYQUIST;

/*===================
 * functions
 *===================*/
void assign   (long unsigned no_part, flouble *x, flouble *y, flouble *z, flouble *mass, flouble *dens);
void add_gas  (flouble *dens);
void get_power(long unsigned no_part, double no_vpart, flouble *x, flouble *y, flouble *z, flouble *mass, flouble *dens,
               double *dpower, long *jpower);
void shuffle_subcubes(long unsigned no_prt, flouble *x, flouble *y, flouble *z, flouble *mass, int level);


/*===================
 * main()
 *===================*/
int main(argc,argv)
int argc;
char **argv;
{   
   partptr       cur_part;
   long unsigned ipart, no_part, FFTarray_length, *idx, BUFSIZE, ik, istart, iend;
   int           ilen, i, j, k, no_timestep, level, MAXLEVEL;
   char          indata[MAXSTRING], outdata[MAXSTRING], header[HEADERSTRING], tmpstring[MAXSTRING];
   FILE          *fpin, *fpout;
   flouble       *dens;
   flouble       *x,  *y,  *z, *mass;
   double        xmin,xmax,ymin,ymax,zmin,zmax;
   double        SubBox, Praw, Pcorr, pik, rkNy, RhoMean, RhoMeanInv;
   
   double        *dpower, FFTnorm;
   double        *rk, *Pk; 
   long          *jpower;
   
   char           dummyline[MAXSTRING], cdummy;
   
   if(argc<4)
     {
      fprintf(stderr,"usage: %s simu_file NGRID MAXLEVEL\n",*argv);
      exit(1);
     }
   
   strcpy(indata, argv[1]);
   NGRID    = (int) atoi(argv[2]);
   MAXLEVEL = (int) atoi(argv[3]);
   NYQUIST  = NGRID/2;  /* Nyquist frequency in grid units */
  
   /*==================
    * open output file
    *==================*/
   /* prepare filename for output file */
   ilen = strlen(argv[1]);
   for(i=ilen-1; i>=0; i--)
     {
      /* chop filename at last '/' */
      if(argv[1][i] == '/')
         break;
     }
   
   for(j=i+1; j<ilen; j++)
      tmpstring[j-(i+1)] = argv[1][j];
   
   tmpstring[j-(i+1)] = '\0';
   sprintf(outdata,"Pk-%05ld-%03d_%s",NGRID,MAXLEVEL,tmpstring);
   
   if((fpout=fopen(outdata,"w"))==NULL)
     {
      printf("I cannot open %s\n", outdata);
      exit(1);
     }
   
   
   /*=================
    * read input file
    *=================*/
   input(argv[1]);
   
   rkNy       = PI*(double)NGRID*(double)MAXLEVEL/io.header.boxsize;       /* Nyquist frequency in physical units */
   RhoMean    = (double)io.header.no_part/pow3(io.header.boxsize);         /* number density of objects           */
   RhoMeanInv = 1./RhoMean;

#ifdef HEADER
   fprintf(fpout,"#%s\n",io.header.header);
   fprintf(fpout,"#no_part     = %ld\n",io.header.no_part);
   fprintf(fpout,"#no_vpart    = %g\n",io.header.no_vpart);
   fprintf(fpout,"#boxsize     = %g\n",io.header.boxsize);
   fprintf(fpout,"#Omega0      = %g\n",io.header.omega0);
   fprintf(fpout,"#Lambda0     = %g\n",io.header.lambda0);
   fprintf(fpout,"#z_initial   = %g (a=%g)\n",1.0/io.header.a_initial-1.0,io.header.a_initial);
   fprintf(fpout,"#z_current   = %g (a=%g)\n",1.0/io.header.a_current-1.0,io.header.a_current);
   fprintf(fpout,"#no_timestep = %d\n",io.header.no_timestep);
   fprintf(fpout,"#NGRID       = %ld\n",NGRID);
   fprintf(fpout,"#MAXLEVEL    = %d\n",MAXLEVEL);
   fprintf(fpout,"#k[h/Mpc](1)    P[(Mpc/h)^3](2)    Correction[(Mpc/h)^3](3)    Nmodes(4)    Pmodes(5)\n");
   fflush(fpout);
#endif
   
   
  /*====================================================================================
   * ALLOCATE ARRAYS
   *====================================================================================*/
   FFTarray_length = (long unsigned)(2*NGRID*NGRID*NGRID);
   FFTnorm         = (double)(NGRID*NGRID*NGRID);
   no_part         = io.header.no_part;
   
   if((dens   = (flouble*) calloc(FFTarray_length,   sizeof(flouble))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR dens (asked for %g GB)\nABORTING\n",
             (double)FFTarray_length*sizeof(flouble)/1024./1024./1024.);
     exit(0);
    }
   if((x      = (flouble*) calloc(no_part,           sizeof(flouble))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR x\nABORTING\n");
     exit(0);
    }
   if((y      = (flouble*) calloc(no_part,           sizeof(flouble))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR y\nABORTING\n");
     exit(0);
    }
   if((z      = (flouble*) calloc(no_part,           sizeof(flouble))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR z\nABORTING\n");
     exit(0);
    }
   if((mass   = (flouble*) calloc(no_part,           sizeof(flouble))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR mass\nABORTING\n");
     exit(0);
    }
   if((dpower = (double*)  calloc(NGRID+1,           sizeof(double))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR dpower\nABORTING\n");
     exit(0);
    }
   if((jpower = (long*)    calloc(NGRID+1,           sizeof(int))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR jpower\nABORTING\n");
     exit(0);
    }
   if((rk     = (double*)  calloc(MAXLEVEL*(NGRID+1),sizeof(double))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR rk\nABORTING\n");
     exit(0);
    }
   if((Pk     = (double*)  calloc(MAXLEVEL*(NGRID+1),sizeof(double))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR Pk\nABORTING\n");
     exit(0);
    }
   if((idx    = (long unsigned*)calloc(MAXLEVEL*(NGRID+1),sizeof(long unsigned))) == NULL)
    {
     fprintf(stderr,"NOT ENOUGH MEMORY FOR idx\nABORTING\n");
     exit(0);
   }
    
  
  /*====================================================================================
   * MAIN LOOP OVER ALL THE SUB-VOLUMES
   *====================================================================================*/
   ik = 0;
   for(level=1; level <= MAXLEVEL; level++)
     {
      /* sidelength of SubBox */
      SubBox = io.header.boxsize/pow(2.0,(double)level-1.0);
      
      fprintf(stderr,"\nworking on level %d (%g) => %ld subcubes\n",level,SubBox,pow3((long)pow(2.0,(double)level-1.0)));
      
      /* map all particles onto subcube extension */
      shuffle_subcubes(no_part, x, y, z, mass, level);
      
      /* get power spectrum for current subcube map */
      get_power(no_part, io.header.no_vpart, x, y, z, mass, dens, dpower, jpower);
      
      if(level == MAXLEVEL)
        {
         istart = 1;
         iend   = NYQUIST;
        }
      else if(level == 1 && level != MAXLEVEL)
        {
         istart = 1;
         iend   = NYQUIST/4;
        }
      else if(level == MAXLEVEL)
        {
         istart = 4;
         iend   = NYQUIST;
        }
      else
        {
         istart = 4;
         iend   = NYQUIST/4;
        }
      
      /* save power spectrum in array */
      for(i=istart; i<=iend; i++)
        {
         rk[ik] = ( (double)i*TWOPI/SubBox + (double)(i+1)*TWOPI/SubBox ) / 2.;
         Pk[ik] = dpower[i]/(double)jpower[i]/FFTnorm * pow3(io.header.boxsize)/FFTnorm;
         ik++;
        }
      
     }
   ik--;
   
   
  /*====================================================================================
   * WRITING OUTPUT
   *====================================================================================*/
   /* sort according to rk array */
   indexx(ik,rk-1,idx-1); /* "-1" because indexx returns index in range [1,npart] */
   
   fprintf(stderr," o writing result to file %s\n",outdata);

   if(MAXLEVEL == 1)
     {
      for(i=0; i<=ik; i++)
       {
        /* correction term ala Nuza et al. (2012, Eq.10 in 1202.6057), but for TSC assignment */
        pik   = (PI*rk[i]/(2.*rkNy));
        Pcorr = RhoMeanInv * (1. - 2./3.*pow2(sin(pik)) + 2./15.*pow4(sin(pik)));

        fprintf(fpout,"%g %g %g %ld %g\n",
                 rk[i],
                 Pk[i],
                 Pcorr,
                 jpower[i],
                 dpower[i]);
       }
     }
   else
     {
      for(i=0; i<4; i++)
       {
        fprintf(fpout,"%g %g %g %ld %g\n",
                rk[i],
                Pk[i],
                RhoMeanInv*(1. - 2./3.*pow2(sin(PI*rk[i]/(2.*rkNy))) + 2./15.*pow4(sin(PI*rk[i]/(2.*rkNy)))),
                jpower[i],
                dpower[i]);
       }
      
      for(i=5; i<ik; i+=2)
        {
         /* smooth out P(k) a little bit... */
         Praw  = (Pk[idx[i-1]-1]+Pk[idx[i]-1]+Pk[idx[i+1]-1])/3.;
         
         /* correction term ala Nuza et al. (2012), but for TSC assignment */
         pik   = (PI*rk[idx[i]-1]/(2.*rkNy));
         Pcorr = RhoMeanInv * (1. - 2./3.*pow2(sin(pik)) + 2./15.*pow4(sin(pik)));
         
         fprintf(fpout,"%g %g %g %ld %g\n",
                 rk[idx[i]-1],
                 Praw,
                 Pcorr,
                 (jpower[idx[i-1]-1]+jpower[idx[i]-1]+jpower[idx[i+1]-1])/3,
                 (dpower[idx[i-1]-1]+dpower[idx[i]-1]+dpower[idx[i+1]-1])/3);
        }
     }
   fclose(fpout);
   
   free(io.fst_part);
   free(dens);
   free(x);
   free(y);
   free(z);
   free(mass);
   free(dpower);
   free(jpower);
   free(rk);
   free(Pk);
   free(idx);
}

/*=======================================================================
* shuffle_subcubes:  map all particles from various subcubes onto [0,1]
*=======================================================================*/
void shuffle_subcubes(long unsigned no_part, flouble *x, flouble *y, flouble *z, flouble *mass, int level)
{
   partptr       cur_part;
   long unsigned no_1Dsubcubes, modx, mody, modz, ipart;
   double        dx, dy, dz;
   double        xx, yy, zz;
   
   no_1Dsubcubes = (long unsigned) pow((double)2.0, (double)(level-1));
   
   dx = 1.0/(double)no_1Dsubcubes;
   dy = 1.0/(double)no_1Dsubcubes;
   dz = 1.0/(double)no_1Dsubcubes;
   
   fprintf(stderr,"shuffling sub-cubes...");
  
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(dx, dy, dz) private(ipart, cur_part, xx, yy, zz, modx, mody, modz) shared(io, x, y, z)
#pragma omp for schedule(static)
#endif
   for(ipart=0; ipart<no_part; ipart++)
     {
      cur_part    = io.fst_part + ipart;

#ifdef MULTIMASS
      mass[ipart] = cur_part->weight;
#else
      mass[ipart] = 1.0;
#endif

      xx          = cur_part->pos[X];
      yy          = cur_part->pos[Y];
      zz          = cur_part->pos[Z];
      modx        = (int)(xx/dx);
      mody        = (int)(yy/dy);
      modz        = (int)(zz/dz);
      
      x[ipart]    = (xx-modx*dx)/dx;
      y[ipart]    = (yy-mody*dy)/dy;
      z[ipart]    = (zz-modz*dz)/dz;
     }

   fprintf(stderr,"done\n");
}

/*====================================================================
* get_power:   assign particles to grid
*              calculate density contrast
*              forward FFT
*              average delta^2 in k-space
*====================================================================*/
void get_power(long unsigned no_part, double no_vpart, flouble *x, flouble *y, flouble *z, flouble *mass, flouble *dens,
               double *dpower, long *jpower)
{
   long unsigned nn[NDIM], npart;
   long          i, j, k;
   long          ix, iy, iz;
   double        RhoMean;
   int           *kmod, ks;
   double        kx, ky, kz, fcorrect;
   
   fprintf(stderr,"getting power:\n");

   /* this will help extracting information from data^ */
   kmod = (int*) calloc(NGRID,sizeof(long));
   for(i=0; i<NGRID; i++)
     {
      if(i <= NYQUIST)
         kmod[i] = i;
      else
         kmod[i] = -(NGRID-i);
     }
   
   for(i=0; i<NYQUIST; i++)
     {
      dpower[i] = (double)0.0;
      jpower[i] = 0;
     }
   
   /* FFT array dimensions */
   nn[0]    = (long unsigned) NGRID;
   nn[1]    = (long unsigned) NGRID;
   nn[2]    = (long unsigned) NGRID;
   
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(NGRID) private(iz, iy, ix) shared(dens)
#pragma omp for schedule(static)
#endif
   /* zero density array */
   for(iz=0; iz<NGRID; iz++)
      for(iy=0; iy<NGRID; iy++)
         for(ix=0; ix<NGRID; ix++)
           {
            dens[Re(ix,iy,iz,(long)NGRID)]=0.0;
            dens[Im(ix,iy,iz,(long)NGRID)]=0.0;
           }
            
            
   /* assign particles to grid */
   fprintf(stderr," o assigning particles to %ld grid...",NGRID);
   assign(no_part, x, y, z, mass, dens);
   fprintf(stderr,"done\n");
   
   /* get density contrast */
   fprintf(stderr," o transferring to density contrast...");
   RhoMean = no_vpart / (double)pow3(NGRID);
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(RhoMean, NGRID) private(k, j, i) shared(dens)
#pragma omp for schedule(static)
#endif
   for(k=0; k<NGRID; k++)
      for(j=0; j<NGRID; j++)
         for(i=0; i<NGRID; i++)
           {
            dens[Re(i,j,k,NGRID)] = dens[Re(i,j,k,NGRID)] / RhoMean - 1.;
            dens[Im(i,j,k,NGRID)] = 0.0;
           }
   fprintf(stderr,"done\n");
   
   /* perform FFT of density contrast */
   fprintf(stderr," o doing FFT...");
   fourn(dens-1, nn-1, NDIM, 1);
   fprintf(stderr,"done\n");
   
   /* get |delta|^2 (square of a complex number is real !) */
   fprintf(stderr," o calculating delta^2...");
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(NGRID) private(k, j, i) shared(dens)
#pragma omp for schedule(static)
#endif
   for(k=0; k<NGRID; k++)
      for(j=0; j<NGRID; j++)
         for(i=0; i<NGRID; i++)
           {
            dens[Re(i,j,k,NGRID)] = pow2(dens[Re(i,j,k,NGRID)])+pow2(dens[Im(i,j,k,NGRID)]);
            dens[Im(i,j,k,NGRID)] = 0.0;
           }
   fprintf(stderr,"done\n");
   
   /* get P(k) */
   fprintf(stderr," o calculating P(k)...");
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(kmod, NGRID) private(k, j, i, ks, kx, ky, kz, fcorrect) shared(dpower, jpower)
#pragma omp for schedule(static)
#endif
   for(k=0; k<NGRID; k++)
      for(j=0; j<NGRID; j++)
         for(i=0; i<NGRID; i++)
           {
            /* wave number */
            ks = (int) floor( sqrt( (float)pow2(kmod[i]) +
                                    (float)pow2(kmod[j]) +
                                    (float)pow2(kmod[k]) ));
            
            fcorrect = 1.0;
#ifdef TSC_CORRECT
            kx = PI*(double)kmod[i]/(double)NGRID;
            ky = PI*(double)kmod[j]/(double)NGRID;
            kz = PI*(double)kmod[k]/(double)NGRID;
            
            if(kx > 0.) fcorrect *= sin(kx)/kx;
            if(ky > 0.) fcorrect *= sin(ky)/ky;
            if(kz > 0.) fcorrect *= sin(kz)/kz;
            
            fcorrect = 1./pow3(fcorrect);
#endif
#ifdef WITH_OPENMP
#pragma omp critical
#endif
            if(ks >= 1 && ks <= NYQUIST)
              {
               dpower[ks] += fcorrect*(double)dens[Re(i,j,k,NGRID)];               
               jpower[ks] += 1;
              }
           }  
   fprintf(stderr,"done\n");
}

/*====================================================================
 * assign:   assign particles to the grid points using a TSC scheme
 *====================================================================*/
void assign(long unsigned no_part, flouble *x, flouble *y, flouble *z, flouble *mass, flouble *dens)
{
   long          ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1;
   long unsigned ipart;
   double        rrx, rry, rrz;
   double        hx, hy, hz;
   double        hx0, hy0, hz0, hxp1, hyp1, hzp1, hxm1, hym1, hzm1;
   double        pw, f_b;
   
    /* baryon fraction */
   f_b = 0.0;
   
#ifndef PK_GAS_ONLY   
   for(ipart=0; ipart<no_part; ipart++)
     {
      /* particle mass */
      pw  = (1.-f_b) * (double) mass[ipart];
      
      /* coordinates in subcube units [0,NGRID] */
      rrx = (x[ipart]) * (double)(NGRID);
      rry = (y[ipart]) * (double)(NGRID);
      rrz = (z[ipart]) * (double)(NGRID);
      
      /* index of nearest grid point [0,NGRID] */
      ix  = (int)(rrx+0.5);
      iy  = (int)(rry+0.5);
      iz  = (int)(rrz+0.5);
      
      /* distance to nearest grid point */
      hx  = rrx - (double)ix;
      hy  = rry - (double)iy;
      hz  = rrz - (double)iz;
      
      /* keep track of peridoc boundaries -> [0,NGRID-1] ; NGRID=0  */
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
      dens[Re(ixm1,iym1,izm1,NGRID)] += hxm1*hym1 *hzm1 * pw;
      dens[Re(ix,  iym1,izm1,NGRID)] += hx0 *hym1 *hzm1 * pw;
      dens[Re(ixp1,iym1,izm1,NGRID)] += hxp1*hym1 *hzm1 * pw;
      dens[Re(ixm1,  iy,izm1,NGRID)] += hxm1*hy0  *hzm1 * pw;
      dens[Re(  ix,  iy,izm1,NGRID)] += hx0 *hy0  *hzm1 * pw;
      dens[Re(ixp1,  iy,izm1,NGRID)] += hxp1*hy0  *hzm1 * pw;
      dens[Re(ixm1,iyp1,izm1,NGRID)] += hxm1*hyp1 *hzm1 * pw;
      dens[Re(  ix,iyp1,izm1,NGRID)] += hx0 *hyp1 *hzm1 * pw;
      dens[Re(ixp1,iyp1,izm1,NGRID)] += hxp1*hyp1 *hzm1 * pw;
      dens[Re(ixm1,iym1,  iz,NGRID)] += hxm1*hym1 *hz0 * pw;
      dens[Re(  ix,iym1,  iz,NGRID)] += hx0 *hym1 *hz0 * pw;
      dens[Re(ixp1,iym1,  iz,NGRID)] += hxp1*hym1 *hz0 * pw;
      dens[Re(ixm1,  iy,  iz,NGRID)] += hxm1*hy0  *hz0 * pw;
      dens[Re(  ix,  iy,  iz,NGRID)] += hx0 *hy0  *hz0 * pw;
      dens[Re(ixp1,  iy,  iz,NGRID)] += hxp1*hy0  *hz0 * pw;
      dens[Re(ixm1,iyp1,  iz,NGRID)] += hxm1*hyp1 *hz0 * pw;
      dens[Re(  ix,iyp1,  iz,NGRID)] += hx0 *hyp1 *hz0 * pw;
      dens[Re(ixp1,iyp1,  iz,NGRID)] += hxp1*hyp1 *hz0 * pw;
      dens[Re(ixm1,iym1,izp1,NGRID)] += hxm1*hym1 *hzp1 * pw;
      dens[Re(  ix,iym1,izp1,NGRID)] += hx0 *hym1 *hzp1 * pw;
      dens[Re(ixp1,iym1,izp1,NGRID)] += hxp1*hym1 *hzp1 * pw;
      dens[Re(ixm1,  iy,izp1,NGRID)] += hxm1*hy0  *hzp1 * pw;
      dens[Re(  ix,  iy,izp1,NGRID)] += hx0 *hy0  *hzp1 * pw;
      dens[Re(ixp1,  iy,izp1,NGRID)] += hxp1*hy0  *hzp1 * pw;
      dens[Re(ixm1,iyp1,izp1,NGRID)] += hxm1*hyp1 *hzp1 * pw;
      dens[Re(  ix,iyp1,izp1,NGRID)] += hx0 *hyp1 *hzp1 * pw;
      dens[Re(ixp1,iyp1,izp1,NGRID)] += hxp1*hyp1 *hzp1 * pw;
     }       
#endif /* PK_GAS_ONLY */
}

/*====================================================================
 * add_gas:   add the gas density as stored in cur_node->u[Udens]
 *====================================================================*/
void add_gas(flouble *dens)
{
   long          ix, iy, iz;
   pqptr         cur_pquad;
   cqptr         cur_cquad;
   nqptr         cur_nquad;
   nptr          cur_node;
   
#ifndef PK_DM_ONLY
   /* we explicitly assume that cur_grid is a regular grid! */
   cur_pquad = cur_grid->pquad;
   
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(cur_pquad) private(cur_cquad, cur_nquad, cur_node, ix, iy, iz) shared(dens)
#pragma omp for schedule(static)
#endif
   for(iz=0; iz<cur_pquad->length; iz++)
     {
      cur_cquad = cur_pquad->loc + iz;
      iy        = cur_cquad->y;
      for(cur_nquad = cur_cquad->loc;  
          cur_nquad < cur_cquad->loc + cur_cquad->length; 
          cur_nquad++, iy++) 
        { 
         ix = cur_nquad->x;
         for(cur_node = cur_nquad->loc; 
             cur_node < cur_nquad->loc + cur_nquad->length; 
             cur_node++, ix++)
           {
            dens[Re(ix,iy,iz,NGRID)] += cur_node->u[Udens];
           }
        }
     }
   
#endif /* PK_DM_ONLY */
}
