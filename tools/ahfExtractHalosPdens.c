#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"

#define STEREO2


int main()
{
  double        *xpart, *ypart, *zpart, *pdens, *r, *g, *b, pdens_max, pdens_min, logdens;
  double        *idx;
  long unsigned *idx_inv;
  
  long          *id, nid;
  
  long           iwrite, ipart, *ihalo, nhalo, jhalo, j, mhalo, ifirst, ilast, no_part, itmp, jtmp;
  double         x_fac, v_fac, m_fac;
  char           c, infile[MAXSTRING],outprefix[MAXSTRING],outfile[MAXSTRING],idfile[MAXSTRING],dummyline[MAXSTRING];
  FILE           *fpin, *fpout, *fpid;
  int            i;
  int            single_file;
  
  
  printf("===================================================================\n");
  printf("  extract halo particles from simuPdens file using AHF_particles\n");
  printf("===================================================================\n");
  
  printf("Please give name of input file:                ");
  scanf("%s",infile);
  printf("%s\n",infile);
  printf("Please give total number of particles:         ");
  scanf("%ld",&no_part);
  printf("%ld\n",no_part);   
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
  
  
  
  /*-------------------------------------------------------
   * open files
   *-------------------------------------------------------*/
  if((fpin=fopen(infile,"rb"))==NULL)
    {
      printf("I cannot open %s\n", infile);
      exit(1);
    }
  
  
  /*-------------------------------------------------------
   * allocate memory
   *-------------------------------------------------------*/
  xpart   = (double *) calloc(no_part, sizeof(double));
  ypart   = (double *) calloc(no_part, sizeof(double));
  zpart   = (double *) calloc(no_part, sizeof(double));
  pdens   = (double *) calloc(no_part, sizeof(double));
  idx     = (double *) calloc(no_part, sizeof(double));
  idx_inv = (long unsigned  *) calloc(no_part, sizeof(long unsigned));
  r       = (double *) calloc(no_part, sizeof(double));
  g       = (double *) calloc(no_part, sizeof(double));
  b       = (double *) calloc(no_part, sizeof(double));
  
  /*-------------------------------------------------------
   * read in amigaPdens file (ASCII)
   *-------------------------------------------------------*/
  fprintf(stderr," o reading Pdens file ... ");
#ifdef PDENS_HEADER 
  /* assume 10 header lines (cf. amigaPdens.c) */
  for(i=0; i<10; i++)
    {
      fgets(dummyline,1024,fpin);
    }
#endif
  
  pdens_max = -1000.;
  pdens_min =  1000.;
  for(ipart=0; ipart<no_part; ipart++)
    {
      fgets(dummyline,1024,fpin);
      
#ifdef STEREO2
      /* assume STEREO2 input file */
      sscanf(dummyline,"%c %lf %lf %lf %lf %lf %lf %lf",&c, &xpart[ipart],&ypart[ipart],&zpart[ipart],&r[ipart],&g[ipart],&b[ipart],&idx[ipart]);
#else
      sscanf(dummyline,"%lf %lf %lf %lf %lf %lf %lf",&xpart[ipart],&ypart[ipart],&zpart[ipart],&pdens[ipart],&idx[ipart]);
      if(pdens[ipart] > pdens_max)
        pdens_max = pdens[ipart];
      if(pdens[ipart] < pdens_min)
        pdens_min = pdens[ipart];
#endif
    }
  fprintf(stderr,"done\n");
  
  fclose(fpin);
  
  /* particles are no longer ordered according to AMIGA id */
  indexx(no_part,idx-1,idx_inv-1); /* "-1" because indexx returns index in range [1,npart] */
  
  
  /* dump all particles into the same STEREO2 file */
  if(single_file == 1)
    {
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
  fprintf(stderr," o extracting %5ld halos:\n",nhalo);
  for(mhalo=0; mhalo<nhalo; mhalo++)
    {
      /*----------------------------------------------------------------------
       * read in particles id's for halo "ihalo" from AHF/BDM_particles file
       *----------------------------------------------------------------------*/
      fprintf(stderr,"   o reading particle id's of halo #%5ld ... ",ihalo[mhalo]);
      
      if((fpid=fopen(idfile,"rb"))==NULL)
        {
          printf("I cannot open %s\n", idfile);
          exit(1);
        }
      
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
       * write to output file
       *----------------------------------------------------------------------*/
      /* write multiple output files */
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
      
      for(ipart=0; ipart<nid; ipart++)
        {
          
          iwrite = idx_inv[id[ipart]]-1;
          
#ifdef STEREO2
          fprintf(fpout,"p %12.6g %12.6g %12.6g      %4.3g %4.3g %4.3g\n",
                  xpart[iwrite], ypart[iwrite], zpart[iwrite], r[iwrite],g[iwrite],b[iwrite]);
#else
          logdens  = log10(CurDens-MinDens+1.)/log10(MaxDens-MinDens+1.)
          r[ipart] =    pow(logdens,2.0);
          g[ipart] =    pow(logdens,4.0);
          b[ipart] = 1.-pow(logdens,0.2);
          fprintf(fpout,"p %12.6g %12.6g %12.6g %4.3g %4.3g %4.3g\n",
                  xpart[iwrite], ypart[iwrite], zpart[iwrite], r[iwrite],g[iwrite],b[iwrite]);        
#endif
        }
      
      if(single_file == 0)
        {
          fclose(fpout);
          fprintf(stderr,"done\n");
        }
      
    } /* mhalo */
  
  if(single_file == 1)
    fclose(fpout);
  fclose(fpid);
  
}

