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
   char           outdata[1024], dummyline[1024], fstring[1024], *there_is_more;
   FILE          *fpin, *fpout;
   long unsigned  npart;
   double         nvpart, Xc, Yc, Zc, VXc, VYc, VZc, Mvir, Rvir;
   float          r, g, b;
   int            iHost, jHost, Nsat;
   double         Box, Xhost, Yhost, Zhost, Rhost;
   double         Dist, dX, dY, dZ;
   long           haloID, hostHalo, numSubStruct;
   long           ldummy;
   double         ddummy; 
   
   if(argc<7)
     {
      fprintf(stderr,"Use: %s .AHF_halos r g b iHost Box\n",*argv);
      exit(1);
     }
   
   printf("==============================================================\n");
   printf("  Read *.AHF_halos file and convert to stereo2's geom format\n");
   printf("==============================================================\n");
   
   r     = (float) atof(argv[2]);
   g     = (float) atof(argv[3]);
   b     = (float) atof(argv[4]);
   iHost = (int)   atof(argv[5]);
   Box   = (double)atof(argv[6]);
   
   strcpy(outdata,argv[1]);
   strcat(outdata,".geom");
   fprintf(stderr,"writing GEOM file: %s\n\n",outdata);
	  
   if((fpin=fopen(argv[1],"r"))==NULL)
     {
      printf("I cannot open %s\n", argv[1]);
      exit(1);
     }
   
   if((fpout=fopen(outdata,"w"))==NULL)
     {
      printf("I cannot open %s\n", outdata);
      exit(1);
     }
   
   
   /* overread header line */
   fgets(dummyline,1024,fpin);

   if(iHost >= 0)
     {
      jHost = 0;
      while((there_is_more = fgets(fstring, 1024, fpin)) != NULL && jHost < iHost)
         jHost++;
      
      sscanf(fstring,"%ld %ld %ld %lf %ld %lf %lf %lf %lf %lf %lf %lf",&haloID,&hostHalo,&numSubStruct,&Mvir,&npart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Rvir);
      fprintf(fpout,"s  %16.8g %16.8g %16.8g %16.8g  %f %f %f\n", Xc, Yc, Zc, Rvir, 1.-r,1.-g,1.-b);

      fprintf(stderr,"o getting satellites of host halo (%g,%g,%g) with mass M=%g/h Msun\n",Xc,Yc,Zc,Mvir);
      
      Xhost = Xc;
      Yhost = Yc;
      Zhost = Zc;
      Rhost = Rvir;
      
      Nsat  = 0;
      
      while((there_is_more = fgets(fstring, 1024, fpin)) != NULL)
        {
         /* extract information from last read dummyline */
         sscanf(fstring,"%ld %ld %ld %lf %ld %lf %lf %lf %lf %lf %lf %lf",&haloID,&hostHalo,&numSubStruct,&Mvir,&npart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Rvir);

         
         dX = Xc-Xhost;
         dY = Yc-Yhost;
         dZ = Zc-Zhost;
         
         if(dX >  Box/2.) dX -= Box;
         if(dY >  Box/2.) dY -= Box;
         if(dZ >  Box/2.) dZ -= Box;
         if(dX < -Box/2.) dX += Box;
         if(dY < -Box/2.) dY += Box;
         if(dZ < -Box/2.) dZ += Box;
                           
         Dist = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
                           
         /* write to geom file */
         if(Dist < Rhost)
           {
            fprintf(fpout,"s  %16.8g %16.8g %16.8g %16.8g  %f %f %f\n", Xc, Yc, Zc, Rvir, r,g,b);
            Nsat++;
           }
        }
      
      fprintf(stderr,"o number of satellites = %d\n",Nsat);
     }
   else
     {
      while((there_is_more = fgets(fstring, 1024, fpin)) != NULL)
        {
         /* extract information from last read dummyline */
         sscanf(fstring,"%ld %ld %ld %lf %ld %lf %lf %lf %lf %lf %lf %lf",&haloID,&hostHalo,&numSubStruct,&Mvir,&npart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Rvir);

         /* write to geom file */
         fprintf(fpout,"s  %16.8g %16.8g %16.8g %16.8g  %f %f %f\n", Xc, Yc, Zc, Rvir, r,g,b);
        }
     }
   fclose(fpin);
   fclose(fpout);
}
