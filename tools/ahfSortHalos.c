#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"

int main(argc,argv)
int argc;
char **argv;
{
   char           indata[1024], outdata[1024], dummyline[1024], fstring[1024], *there_is_more;
   FILE          *fpin, *fpout;
   long unsigned  nhalos, ihalo, jhalo, ncolumns, icolumn, *idx, sortcolumn;
   double       **halodata, *space, iVar, *fsort;
   
   
   
   if(argc<5)
     {
      fprintf(stderr,"Use: %s .AHF_halos .AHF_halos.sorted NCOLUMNS sortcolumn\n",*argv);
      exit(1);
     }
   
   printf("==========================================================================\n");
   printf("  Read *.AHF_halos files and sort the halos according to specific column\n");
   printf("==========================================================================\n");
   
   strcpy(indata,argv[1]);
   strcpy(outdata,argv[2]);
   ncolumns   = (long unsigned) atof(argv[3]);
   sortcolumn = (long unsigned) atof(argv[4]);
   if(sortcolumn > ncolumns)
    {
     fprintf(stderr,"Use: %s .AHF_halos .AHF_halos.sorted NCOLUMNS sortcolumn\n",*argv);
     fprintf(stderr,"     sortcolumn has to be obviously smaller than NCOLUMNS!\n");
     exit(1);
    }
   
   fprintf(stderr,"o reading from input file: %s\n",indata);
   fprintf(stderr,"o writing  to output file: %s\n\n",outdata);
	  
   if((fpin=fopen(indata,"r"))==NULL)
     {
      printf("I cannot open %s\n", indata);
      exit(1);
     }
   
   if((fpout=fopen(outdata,"w"))==NULL)
     {
      printf("I cannot open %s\n", outdata);
      exit(1);
     }
   
   
   /* read/write header line */
   fgets(dummyline,1024,fpin);
   fprintf(fpout,"%s",dummyline);

   /* count the number of halos */
   nhalos = 0;
   while((there_is_more = fgets(fstring, 1024, fpin)) != NULL)
      nhalos++;
   fprintf(stderr,"o found %10d halos\n\n",nhalos);

   /* allocate memory for halodata */
   space    = (double *)   calloc(nhalos*ncolumns , sizeof(double));
   halodata = (double **)  calloc(nhalos,           sizeof(double *));
   for(ihalo=0; ihalo<nhalos; ihalo++)
      halodata[ihalo] = space + ihalo*ncolumns; 	

   
   
   /* rewind to first halo */
   rewind(fpin);
   fgets(dummyline,1024,fpin);

   /* eventually read all halos into halodata */
   for(ihalo=0; ihalo<nhalos; ihalo++)
     {
      for(icolumn=0; icolumn<ncolumns; icolumn++) 
        {
         fscanf(fpin,"%lf",&iVar);
         halodata[ihalo][icolumn] = iVar;
        }
     }
   
   /* sort halos according to [sortcolumn] */
   fsort = (double *)        calloc(nhalos+1, sizeof(double));
   idx   = (long unsigned *) calloc(nhalos+1, sizeof(long unsigned));
   for(ihalo=0; ihalo<nhalos; ihalo++)
      fsort[ihalo] = halodata[ihalo][sortcolumn];
   indexx(nhalos, fsort-1, idx-1);
   
   /* write sorted halos */
   for(jhalo=nhalos-1; jhalo>0; jhalo--)
     {
      ihalo = idx[jhalo]-1;
      
      fprintf(fpout,"%ld   ",(long unsigned)halodata[ihalo][0]);
      for(icolumn=1; icolumn<ncolumns; icolumn++)
        {
         fprintf(fpout,"%lf   ",halodata[ihalo][icolumn]);
        }
      fprintf(fpout,"\n");
      fflush(fpout);
     }
   
   
   free(fsort);
   free(space);
   free(halodata);
   
   
   fclose(fpin);
   fclose(fpout);
}
