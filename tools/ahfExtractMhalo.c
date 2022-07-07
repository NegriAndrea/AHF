#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PGAS  0
#define PDM   1
#define PSTAR 4


/*==================================================================================================
 * main:
 *
 *
 *==================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
   long   nhalos, ihalo, npart, ipart;
   
   char   outfile[2048], indata[2048];
   int    slen;
   
   struct {
      long   id;
      long   itype;
      double mass;
   } p;
   
   struct {
      double mDM;
      double mgas;
      double mstars;
      long   ncont;
   } h;
   
   FILE *fpin, *fpout;
   
   /*===================================================================
    * check consistency
    *===================================================================*/
   if(argc<2)
     {
      fprintf(stderr,"usage: %s .AHF_particles\n",*argv);
      exit(1);
     }
   
   /*===================================================================
    * prepare output filename
    *===================================================================*/
   strcpy(outfile,argv[1]);
   
   slen=strlen(outfile);
   
   outfile[slen-9]='m';
   outfile[slen-8]='h';
   outfile[slen-7]='a';
   outfile[slen-6]='l';
   outfile[slen-5]='o';
   outfile[slen-4]='\0';
   
   /*===================================================================
    * be verbose
    *===================================================================*/
   printf("================================================================================\n");
   printf("  Read *.AHF_particles file and extract DM, gas, and stellar mass for each halo\n");
   printf("================================================================================\n");
   
   if(!(fpin=fopen(argv[1],"r")))
     {
      fprintf(stderr,"FATAL: could not open input file %s\n",argv[1]);
      exit(0);
     }
   if(!(fpout=fopen(outfile,"w")))
     {
      fprintf(stderr,"FATAL: could not open output file %s\n",outfile);
      exit(0);
     }

   fprintf(stderr,"o reading: %s\n",argv[1]);
   fprintf(stderr,"o writing: %s\n",outfile);
   
   
   /*===================================================================
    * start reading and writing
    *===================================================================*/
   /* total number of halos */
   fgets(indata,2048,fpin);
   sscanf(indata,"%ld",&nhalos);
   
   /* loop over all halos in file */
   for(ihalo=0; ihalo<nhalos; ihalo++)
     {
      /* read number of particles in halo */
      fgets(indata,2048,fpin);
      sscanf(indata,"%ld",&npart);
      
      /* reset halo masses */
      h.mDM    = 0.;
      h.mgas   = 0.;
      h.mstars = 0.;
      h.ncont  = 0;
      
      /* loop over all particles in halo */
      for(ipart=0; ipart<npart; ipart++)
        {
         /* read particle information */
         fgets(indata,2048,fpin);
         sscanf(indata,"%ld %ld %lf",&(p.id), &(p.itype), &(p.mass));
         
         if(p.itype == PDM)
            h.mDM    += p.mass;
         else if(p.itype == PGAS)
            h.mgas   += p.mass;
         else if(p.itype == PSTAR)
            h.mstars += p.mass;
         else
            h.ncont  += 1;
        }
      
      /* write halo masses */
      fprintf(fpout,"%e %e %e %ld\n",h.mDM,h.mgas,h.mstars,h.ncont);
     }
}