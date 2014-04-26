#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
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
#define VLOG   // use logarithmic binning for velocity distribution


/*===================
 * COMMON structures
 *===================*/
info_io   io;

int       NBINS;     /* number of BINS for velocity distribution  */
double    x_fac, v_fac;



/*===================
 * functions
 *===================*/
void statistic (double vmin, double vmax, FILE *fpout);



/*===================
 * main()
 *===================*/
int main(argc,argv)
int argc;
char **argv;
{
   long unsigned i,j,k,np;
   int           ii,jj,ilen, no_timestep;
   double        weight;
   double        vmean[3], vm, vm3D, vmin, vmax, vtemp, sigV, sigV3D, BUFSIZE;
   double        Lmean[3], Lm, Lm3D;
   char          indata[MAXSTRING], outdata[MAXSTRING], tmpstring[MAXSTRING], header[HEADERSTRING];
   FILE          *fpout;
   partptr       cur_part;
   double        x, y, z, vx, vy, vz;
   
   if(argc<3)
     {
      fprintf(stderr,"usage: %s simu_file NBINS\n",*argv);
      exit(1);
     }
   
   printf("==================================================================\n");
   printf("  Read AMIGA binary and calculate velocity distribution function\n");
   printf("==================================================================\n");

   strcpy(indata, argv[1]);
   NBINS = (int) atoi(argv[2]);
   
   /*==================
    * open output file
    *==================*/
   /* prepare output filename */
   ilen = strlen(argv[1]);
   for(ii=ilen-1; ii>=0; ii--)
     {  
      /* chop filename at last '/' */
      if(argv[1][ii] == '/') 
         break;
     }
   for(jj=ii+1; jj<ilen; jj++)
      tmpstring[jj-(ii+1)] = argv[1][jj];
   
   tmpstring[jj-(ii+1)] = '\0';
   sprintf(outdata,"Vstat-%d_%s",NBINS,tmpstring);
   fprintf(stderr,"\n-> writing result to file %s\n\n",outdata);
         
   if((fpout=fopen(outdata,"w"))==NULL)
     {
      printf("I cannot open %s\n", outdata);
      exit(1);
     }

   
   /*=================
    * read input file
    *=================*/
   input(argv[1]);
   
#ifdef HEADER 
  fprintf(fpout,"#%s\n",io.header.header);
  fprintf(fpout,"#no_part     = %ld\n",io.header.no_part);
  fprintf(fpout,"#no_vpart    = %g\n",io.header.no_vpart);
  fprintf(fpout,"#no_vpart    = %g\n",io.header.no_vpart);
  fprintf(fpout,"#boxsize     = %g\n",io.header.boxsize);
  fprintf(fpout,"#Omega0      = %g\n",io.header.omega0);
  fprintf(fpout,"#Lambda0     = %g\n",io.header.lambda0);
  fprintf(fpout,"#z_initial   = %g (a=%g)\n",
          1.0/io.header.a_initial-1.0,io.header.a_initial);
  fprintf(fpout,"#z_current   = %g (a=%g)\n",
          1.0/io.header.a_current-1.0,io.header.a_current);
  fprintf(fpout,"#no_timestep = %d\n",io.header.no_timestep);
  fprintf(fpout,"#NBINS       = %d\n",NBINS);
#endif

  x_fac  = io.header.boxsize;
  v_fac  = io.header.boxsize/io.header.t_unit / io.header.a_current;

  
  
  
  /*================================================================================
   * calculate means
   *================================================================================*/
  vmean[0] =  0.0;
  vmean[1] =  0.0;
  vmean[2] =  0.0;
  Lmean[0] =  0.0;
  Lmean[1] =  0.0;
  Lmean[2] =  0.0;
  vm       =  0.0;
  Lm       =  0.0;
  vm3D     =  0.0;
  vmin     =  1000000.;
  vmax     = -1000000.;
  for(i=0;i<io.header.no_part;i++)
    {
     cur_part = io.fst_part + i;
     
     
     /* set mass to 1 in case we are dealing with single mass simulation */
#ifdef MULTIMASS
     weight = cur_part->weight;
#else
     weight = 1.0;
#endif
     
     /* convert coordinates and velocities */
     x  = cur_part->pos[X] * x_fac;
     y  = cur_part->pos[Y] * x_fac;
     z  = cur_part->pos[Z] * x_fac;
     vx = cur_part->mom[X] * v_fac;
     vy = cur_part->mom[Y] * v_fac;
     vz = cur_part->mom[Z] * v_fac;
     
     /* velocity statistics */
     vmean[0] += weight * vx;
     vmean[1] += weight * vy;
     vmean[2] += weight * vz;
     vm       += weight * sqrt( pow2(vx) + pow2(vy) + pow2(vz) );
     
     /* angular momentum statistics */
     Lmean[0] += weight * (vy*z - vz*y);
     Lmean[1] += weight * (vz*x - vx*z);
     Lmean[2] += weight * (vx*y - vy*x);
     Lm       += weight * sqrt( pow2(vy*z - vz*y) + pow2(vz*x - vx*z) + pow2(vx*y - vy*x) );
     
     /* this will not give credible results for OpenMP */
     vtemp = sqrt( pow2(vx) + pow2(vy) + pow2(vz) );
     if(vtemp > vmax) vmax = vtemp;
     if(vtemp < vmin) vmin = vtemp;
    }
  
  vmean[0] /= (double)io.header.no_vpart;
  vmean[1] /= (double)io.header.no_vpart;
  vmean[2] /= (double)io.header.no_vpart;
  vm       /= (double)io.header.no_vpart;
  Lmean[0] /= (double)io.header.no_vpart;
  Lmean[1] /= (double)io.header.no_vpart;
  Lmean[2] /= (double)io.header.no_vpart;
  Lm       /= (double)io.header.no_vpart;
  vm3D      = sqrt(pow2(vmean[0])+pow2(vmean[1])+pow2(vmean[2]));
  Lm3D      = sqrt(pow2(Lmean[0])+pow2(Lmean[1])+pow2(Lmean[2]));

  /*================================================================================
   * calculate dispersions
   *================================================================================*/
  sigV   = 0.0;
  sigV3D = 0.0;
  for(i=0;i<io.header.no_part;i++)
    {
     cur_part = io.fst_part + i;
     vx       = cur_part->mom[X] * v_fac;
     vy       = cur_part->mom[Y] * v_fac;
     vz       = cur_part->mom[Z] * v_fac;
     
     sigV3D += pow2(vx-vmean[0]) + pow2(vy-vmean[1]) + pow2(vz-vmean[2]);
     vtemp   = sqrt( pow2(vx) + pow2(vy) + pow2(vz) );
     sigV   += pow2(vtemp-vm);
    }
  sigV   = sqrt(sigV  /(double)io.header.no_part);
  sigV3D = sqrt(sigV3D/(double)io.header.no_part);

  printf("\nvelocity statistic:\n");
  printf("min. v        = %12.3g km/sec\n",vmin);
  printf("max. v        = %12.3g km/sec\n",vmax);
  printf(" sum |v|/N    = %12.3g km/sec  +- %12.3g km/sec (1-sigma)\n",vm,sigV);
  printf("|sum v| /N    = %12.3g km/sec  +- %12.3g km/sec (1-sigma) (vx=%g, vy=%g, vz=%g)\n",
	 vm3D,sigV3D,vmean[0],vmean[1],vmean[2]);
  printf("|sum v|       = %12.3g km/sec\n",vm3D*io.header.no_vpart);
  printf("angular momentum statistic:\n");
  printf("mean     = %12.3f\n",Lm);
  printf("mean3D   = %12.3f (Lx=%g, Ly=%g, Lz=%g)\n",Lm3D,Lmean[0],Lmean[1],Lmean[2]);
  
  fprintf(fpout,"#velocity statistic:\n");
  fprintf(fpout,"# min      = %12.3g km/sec\n",vmin);
  fprintf(fpout,"# max      = %12.3g km/sec\n",vmax);
  fprintf(fpout,"# mean     = %12.3g km/sec  +- %12.3g km/sec (1-sigma)\n",vm,sigV);
  fprintf(fpout,"# mean3D   = %12.3g km/sec  +- %12.3g km/sec (1-sigma)\n",vm3D,sigV3D);
  fprintf(fpout,"#angular momentum statistic:\n");
  fprintf(fpout,"# mean     = %12.3g km/sec\n",Lm);
  fprintf(fpout,"# mean3D   = %12.3g km/sec\n",Lm3D);
  
  fprintf(stderr,"\n => writing velocity distribution '%s'... ",outdata);
  statistic(vmin,vmax,fpout);
  fprintf(stderr,"done\n");

  fclose(fpout);
}


/*===========================================================================
* do you own stuff with simulation data given now in proper units...
*===========================================================================*/
void statistic(double vmin, double vmax, FILE *fpout)
{
   double          vlog[1000];
   long unsigned   i_vel[1000];
   
   double v1log, v2log, vdlog, vAlog, v, v_part_log;
   long unsigned   i, Nbin, Ibin, i_all;
   double vx, vy, vz, x, y, z;
   partptr cur_part;
   
   Nbin  = NBINS;
#ifdef VLOG
   v1log = log10(vmin);                   
   v2log = log10(vmax);
#else
   v1log = (vmin);                   
   v2log = (vmax);
#endif
   vdlog = (v2log - v1log)/(double)(Nbin);                  
   vAlog = v1log + vdlog/2;
   
   
   for(i=0; i<1000; i++)
     {
      vlog[i]  = 0.0;
      i_vel[i] = 0;
     }
   
   for(i=0; i<Nbin; i++)
     {        
      vlog[i] = v1log + (double)i * vdlog;
     }
   
   
   for(i=0; i<io.header.no_part; i++)
     {
      cur_part   = io.fst_part + i;
      
      vx         = cur_part->mom[X] * v_fac;
      vy         = cur_part->mom[Y] * v_fac;
      vz         = cur_part->mom[Z] * v_fac;
      
      v          = sqrt(pow2(vx)+pow2(vy)+pow2(vz));
#ifdef VLOG
      v_part_log = log10(v);
#else
      v_part_log = (v);
#endif
      Ibin       = (int) ((v_part_log-v1log)/vdlog + 0.5);
      
      if(Ibin < 0 || Ibin > Nbin) 
        {
         fprintf(stderr,"\nsomething going on...\n");
         exit(0);
        }
      else
         i_vel[Ibin] += 1;      
     }
   
   for(i=0, i_all=0; i<Nbin; i++)
     {
      i_all += i_vel[i];
#ifdef VLOG
      fprintf(fpout,"%g %g %d\n", pow(10.,vlog[i]), (double)i_vel[i]/(double)io.header.no_part/vdlog, i_all);
#else
      fprintf(fpout,"%g %g %d\n", vlog[i], (double)i_vel[i]/(double)Nmax/vdlog, i_all);
#endif
     }
}
