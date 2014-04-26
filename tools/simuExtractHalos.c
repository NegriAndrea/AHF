#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libio_serial/io_serial.h"
#include "../src/libutility/utility.h"

//#define DEFCOLOR	// define your own colors for stereo2
//#define STEREO2		// write stereo2 file 
//#define ID_CHECK

//#define DEBUG_HOLE
#define Xhole        23.03291167
#define Yhole        32.58578480
#define Zhole        27.45872678
#define HoleDistance  0.002


/*===================
 * COMMON structures
 *===================*/
info_io io;

/*=========================================================================================*/
/*=========================================================================================*/
/*=========================================================================================*/
int main()
{
  double        xpart, ypart, zpart, vxpart, vypart, vzpart, mpart, *dist2, dist;
  long unsigned no_halos;
  partptr       cur_part;
  
  long         *id, nid;
  int           px, py, pz;
  
  long          ipart, jpart, *ihalo, nhalo, jhalo, j, mhalo, ifirst, ilast;
  double        x_fac, v_fac, m_fac;
  char          infile[MAXSTRING],outprefix[MAXSTRING],outfile[MAXSTRING],idfile[MAXSTRING],dummyline[MAXSTRING],PLfile[MAXSTRING];
  FILE          *fpin, *fpout, *fpid, *fpPL;
  int           single_file;
  double        r,g,b;
  int           ISEED=-123456;
  
  
  printf("====================================================================\n");
  printf(" extract halo particles from simulation binary using AHF_particles\n");
  printf("====================================================================\n");
  
  printf("Please give name of input file:                ");
  scanf("%s",infile);
  printf("%s\n",infile);
  printf("Please give name of AHF_particles:             ");
  scanf("%s",idfile);
  printf("%s\n",idfile);
  printf("Please give number of halos to extract:        ");
  scanf("%ld",&nhalo);
  printf("%ld\n",nhalo);
  
  /* a negative nhalo means that we want to extract all consecutive halos from halo #i -> halo #j */
  if(nhalo < 0)
    {
      printf("Please give id of first halo to extract:       ", j);
      scanf("%ld",&ifirst);
      printf("%ld\n",ifirst);
      printf("Please give id of last  halo to extract:       ", j);
      scanf("%ld",&ilast);
      printf("%ld\n",ilast);
      nhalo = ilast-ifirst+1;
      ihalo = (long *) calloc(nhalo, sizeof(long));
      for(j=0; j<nhalo; j++)
        ihalo[j] = ifirst+j;
    }
  /* otherwise specify the halo id's manually */
  else
    {
      ihalo = (long *) calloc(nhalo, sizeof(long));
      for(j=0; j<nhalo; j++)
        {
          printf("Please give id of %5ld.halo to extract:       ", j);
          scanf("%ld",&ihalo[j]);
          printf("%ld\n",ihalo[j]);
        }
    }
  printf("Please give prefix for output file:            ");
  scanf("%s",outprefix);
  printf("%s\n\n",outprefix);
  
  /* either dump all particles into a single file or write one file for each halo */
  printf("Writing particles to single file (1:yes/0:no): ");
  scanf("%d",&single_file);
  printf("%d\n\n",single_file);   
  
#ifdef DEFCOLOR
  printf("RGB color (0.0-1.0):");
  scanf("%lf %lf %lf",&r, &g, &b);
  printf("%g %g %g\n\n",r,g,b);   
#endif   
  /*-------------------------------------------------------
   * read simulation raw data
   *-------------------------------------------------------*/
  input(infile);
  
  /*-------------------------------------------------------
   * converion factors
   *-------------------------------------------------------*/    
#ifdef STEREO2
  x_fac  = io.header.boxsize;
  v_fac  = io.header.boxsize / io.header.t_unit / io.header.a_current;
  m_fac  = io.header.pmass;
#else
  x_fac = 1.0;
  v_fac = 1.0;
  m_fac = 1.0;
#endif
  
  /* dump all particles into the same STEREO2 file */
  if(single_file == 1)
    {
      fprintf(stderr,"---------------------------------------------------------------------\n");
      fprintf(stderr,"             amigaExtractHalos taking over now...\n");
      fprintf(stderr,"---------------------------------------------------------------------\n");
      sprintf(outfile,"%s",outprefix);
      fprintf(stderr,"\n o writing all particles to %s\n",outfile);
      if((fpout=fopen(outfile,"w"))==NULL)
        {
          printf("I cannot open %s\n", outfile);
          exit(1);
        }
    }
  
  /*----------------------------------------------------------------------
   * loop over all halos to extract
   *----------------------------------------------------------------------*/
  fprintf(stderr,"\n o extracting %5ld halos:\n",nhalo);
  for(mhalo=0; mhalo<nhalo; mhalo++)
    {
      /*----------------------------------------------------------------------
       * read in particles id's for halo "ihalo" from AHF_particles file
       *----------------------------------------------------------------------*/
      fprintf(stderr,"   o reading particle id's of halo #%5ld ... ",ihalo[mhalo]);
      
      if((fpid=fopen(idfile,"rb"))==NULL)
        {
          printf("I cannot open %s\n", idfile);
          exit(1);
        }
     
     /* the first line is the total number of halos */
     fscanf(fpid,"%ld\n",&no_halos);
      
      jhalo = -1;
      id    = NULL;
      
      while( jhalo < ihalo[mhalo] )
        {
          /* how many id's in current halo */
          fscanf(fpid,"%ld\n",&nid);
          
          /* allocate array to hold particle id's */
          id = (long *) realloc(id, nid*sizeof(long));
          
          for(j=0; j<nid; j++)
            {
              fgets(dummyline,MAXSTRING,fpid);
              sscanf(dummyline,"%ld",&id[j]);
            }
          
          /* next halo */
          jhalo++;
        }
      fprintf(stderr,"done\n");
      fclose(fpid);
      
      /*----------------------------------------------------------------------
       * write to STEREO2 file
       *----------------------------------------------------------------------*/
      /* write multiple STEREO2 files */
      if(single_file == 0)
        {
          sprintf(outfile,"%s-%ld",outprefix,ihalo[mhalo]);
          fprintf(stderr,"   o writing %10ld particles to %s ... ",nid,outfile);
          
          if((fpout=fopen(outfile,"w"))==NULL)
            {
              printf("I cannot open %s\n", outfile);
              exit(1);
            }
        }
      
      /* assign a color to the particles in this particular halo */
#ifndef DEFCOLOR
      r = (double) ran3(&ISEED);
      g = (double) ran3(&ISEED);
      b = (double) ran3(&ISEED);
#endif  
      
#ifdef ID_CHECK
      for(ipart=0; ipart<nid; ipart++)
        {
          for(jpart=0; jpart<ipart; jpart++)
            if(id[ipart] == id[jpart])
              fprintf(stderr,"dublicate particles in halo %d\n",ihalo[mhalo]);
          for(jpart=ipart+1; jpart<nid; jpart++)
            if(id[ipart] == id[jpart])
              fprintf(stderr,"dublicate particles in halo %d\n",ihalo[mhalo]);
        }
#endif
      
      for(ipart=0; ipart<nid; ipart++)
        {
          /* current particle coordinates */
          if(id[ipart] >= 0) cur_part = io.fst_part + id[ipart];
          
          xpart  = cur_part->pos[X] * x_fac;
          ypart  = cur_part->pos[Y] * x_fac;
          zpart  = cur_part->pos[Z] * x_fac;
          vxpart = cur_part->mom[X] * v_fac;
          vypart = cur_part->mom[Y] * v_fac;
          vzpart = cur_part->mom[Z] * v_fac;
#ifdef MULTIMASS
          mpart  = cur_part->weight * m_fac;
#else
          mpart  = m_fac;
#endif
          
#ifdef DEBUG_HOLE
          dist = pow2(xpart-Xhole)+pow2(ypart-Yhole)+pow2(zpart-Zhole);
          
          if(dist < pow2(HoleDistance))
#endif
#ifdef STEREO2        
            fprintf(fpout,"P  %20.10g %20.10g %20.10g    %4.3g %4.3g %4.3g  1\n",  xpart, ypart, zpart, r,g,b);
#else
          fprintf(fpout,"%14.8g %14.8g %14.8g    %14.8g %14.8g %14.8g    %16.8g\n", xpart, ypart, zpart, vxpart, vypart, vzpart, mpart);
#endif
          fflush(fpout);
        }
      
      if(single_file == 0)
        {
          fclose(fpout);
          fprintf(stderr,"done\n");
        }
      
    } /* mhalo */
  
  if(single_file == 1)
    {
      fclose(fpout);
      fprintf(stderr,"done\n");
    }
  
  
  free(io.fst_part); 
}
