/*==================================================================================================
 *  MergerRates:
 *
 *   - reads MergerTree's *_mtree files
 *   - counts 'mergers' in each *_mtree file defined via MERGER_RATIO below
 *
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <libgen.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"
#include "../src/libutility/utility.h"

#define M0                                    1000           // these are "number of particles" criteria!
#define Mi                                    40
#define MERGER_RATIO                          0.25           // writes output that readily allows to find mergers

#define GOTTLOEBER_CRITERION                  1              // maximum number of progenitors to check for a possible merger (Gottloeber rates)
#define FAKHOURI_CRITERION                    3000000        // maximum number of progenitors to check for a possible merger (non-Gottloeber rates)

#define HALOID_RESTRICTION

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct MTREE *MTREEptr;
typedef struct MTREE
{
  uint64_t  haloid;
  uint64_t  npart;
  uint64_t  nprog;
  uint64_t *progid;
  uint64_t *nshared;
  uint64_t *npartprog;
  uint64_t  nmerger;
  uint64_t *iprogmerger;
#ifdef HALOID_RESTRICTION
  int8_t    follow_halo;
  int8_t   *follow_progenitor;
#endif
} MTREE;


/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/
float   *zred; // table containing the redshifts for all the outputs
float   *tage; // table containing the corresponding ages of the universe
float    hubble;
char     haloidsfile[MAXSTRING];

/*-------------------------------------------------------------------------------------
 *                                   PROTOTYPES
 *-------------------------------------------------------------------------------------*/
int      merger_rates     (char Filename[MAXSTRING], int);
MTREE   *read_mtree       (char *, uint64_t *);
uint64_t count_merger     (MTREE *, uint64_t, uint64_t *);
uint64_t count_progs      (MTREE *, uint64_t);
void     write_merger     (FILE *, MTREE *, uint64_t);
void     write_progs      (FILE *, MTREE *, uint64_t);
void     read_zred        (char zredfile[MAXSTRING], int , int);
int      test_zredfile    (char zredfile[MAXSTRING]);
#ifdef HALOID_RESTRICTION
void     read_haloids        (char haloidsfile_izred[MAXSTRING], MTREE *, uint64_t);
void     write_haloids       (char haloidsfile_izred[MAXSTRING], MTREE *, uint64_t);
int      qcompareHaloIDs  (const void *, const void *);
int      bcompareMtreeIDs (const void *, const void *);
#endif


/*-------------------------------------------------------------------------------------
 *                              CHECK FOR USER STUPIDITY
 *-------------------------------------------------------------------------------------*/
#if (defined GOTTLOEBER_CRITERION && defined FAKHOURI_CRITERION)
#undef FAKHOURI_CRITERION
#endif

/*==================================================================================================
 * main
 *==================================================================================================*/
int main()
{
  int      i, nFiles, isimu, need_simuparams;
  char   **Filename;
  char     zredfile[MAXSTRING];
  double   omega0, lambda0, h;
  
  /*==================================================================*
   *                          USER INTERFACE                          *
   *==================================================================*/
  fprintf(stderr,"==========================================================\n");
  fprintf(stderr,"  count 1:%1d mergers in each of a list of *_mtree files\n",(int)(1./MERGER_RATIO+0.5));
  fprintf(stderr,"==========================================================\n");
#ifdef GOTTLOEBER_CRITERION
  fprintf(stderr,"(using Gottloeber et al. (2000) criterion to define mergers)\n");
#endif
#ifdef FAKHOURI_CRITERION
  fprintf(stderr,"(using Fakhouri et al. (2010) criterion to define mergers)\n");
#endif
#ifdef HALOID_RESTRICTION
  fprintf(stderr,"Please give name of haloidsfile:               ");
  scanf("%s", haloidsfile);
  fprintf(stderr,"%s\n",haloidsfile);
#endif
  fprintf(stderr,"\nPlease give number of _mtree files:      ");
  scanf("%d", &nFiles);
  fprintf(stderr,"%d\n",nFiles);
  
  /* allocate memory for nFiles filenames, each of size MAXSTRING */
  Filename  = (char **) calloc(nFiles, sizeof(char *));
  for(i=0; i<nFiles; i++) {
    Filename[i]  = (char *) calloc(MAXSTRING, sizeof(char));
  }
  
  /* read input filenames from stdin */
  for(i=0; i<nFiles; i++) {
    fprintf(stderr,"Please give name of %5d. *_mtree file:            ",i+1);
    scanf("%s", Filename[i]);
    fprintf(stderr,"%s\n",Filename[i]);
  }
  fprintf(stderr,"\n");

  fprintf(stderr,"\nPlease give name of file with redshifts: ");
  scanf("%s", &zredfile);
  fprintf(stderr,"%s\n",zredfile);
  need_simuparams = test_zredfile(zredfile);
  
  if(need_simuparams) {
    fprintf(stderr,"\nPlease give Omega0:                      ");
    if(scanf("%lf", &omega0) == 0)
      omega0 = 1.0;
    fprintf(stderr,"%lf\n",omega0);
    
    fprintf(stderr,"\nPlease give OmegaLambda0:                ");
    if(scanf("%lf", &lambda0) == 0)
      lambda0 = 0.0;
    fprintf(stderr,"%lf\n",lambda0);
    
    fprintf(stderr,"\nPlease give h:                           ");
    if(scanf("%lf", &h) == 0)
      h = 1.0;
    fprintf(stderr,"%lf\n",h);
    
    // make these parameters available to calc_t() and read_zred()
    simu.omega0  = omega0;
    simu.lambda0 = lambda0;
    hubble       = h;
  }
  
  
  /*======================================================================*
   *  READ FILE CONTAINING THE REDSHIFTS                                  *
   *======================================================================*/
  read_zred(zredfile, nFiles, need_simuparams);

  
  /*======================================================================*
   *  LOOP OVER ALL INPUT FILES                                           *
   *======================================================================*/
#ifndef HALOID_RESTRICTION
#pragma omp parallel for schedule(dynamic) default(none) private(i) shared(Filename,nFiles)
#endif
  for(i=0; i<nFiles; i++) {
    
    /* read the next file into memory */
    merger_rates(Filename[i], i);
    
  } // for(nFiles)
  
  
  /*==================================================================*
   *                             CLEANUP                              *
   *==================================================================*/
  fprintf(stderr,"\nCleaning up ... ");
  
  /* remove filename storage */
  for(i=0; i<nFiles; i++) {
    if(Filename[i])  free(Filename[i]);
  }
  if(Filename)  free(Filename);
  
  printf("finished\n");
  return(1);
}


/*==================================================================================================
 * merger_rates:
 *
 *  get statistics for multiple progenitors (e.g. mass ratios, etc.)
 *
 *==================================================================================================*/
int merger_rates(char Filename[MAXSTRING], int izred)
{
  FILE    *fpout_progs, *fpout_merger;
  char     Filename_progs[MAXSTRING], Filename_merger[MAXSTRING], haloidsfile_izred[MAXSTRING];
  uint64_t ihalo, halo_nmerger, halo_nprogs, nhalos, nhalos_M0;
  MTREE   *mtree;
  
  fprintf(stderr,"  o merger rates for %s: ",Filename);
  
  // open *_progs file
  //===================
  strcpy(Filename_progs, basename(Filename));
  strcat(Filename_progs, "_progs");
  fpout_progs = fopen(Filename_progs,"w");
  if(fpout_progs == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",Filename_progs);
    exit(0);
  }
  fprintf(fpout_progs,"#z(1) t[Gyr](2) dz(3) dt[Gyr](4) Nprogs(5) Nhalos(6)\n");
  fprintf(fpout_progs,"#N0(1) N1(2) dN(3) dN/N0(4)\n");
  
  // put the redshift into the file (and later next to it the total number of mergers and halos)
  fprintf(fpout_progs,"%f %f %f %f ", zred[izred], tage[izred], (zred[izred+1]-zred[izred]), (tage[izred]-tage[izred+1]));
  

  // open *_merger file
  //====================
  strcpy(Filename_merger, basename(Filename));
  strcat(Filename_merger, "_merger");
  fpout_merger = fopen(Filename_merger,"w");
  if(fpout_merger == NULL)  {
    fprintf(stderr,"could not open file %s\nexiting\n",Filename_merger);
    exit(0);
  }
  fprintf(fpout_merger,"#z(1) t[Gyr](2) dz(3) dt[Gyr](4) Nmerger(5) Nhalos_M0(6) M0=%d Mi=%d MERGER_RATIO=%lf\n",(int)M0, (int)Mi,(double)MERGER_RATIO);
  fprintf(fpout_merger,"#haloid(1)      Npart(2)      Nmerger(3)    progid(4) Nshared(5) Npartprog(6)\n");
  fprintf(fpout_merger,"# progid(1)      Npartprog(2)      Nshared(3)\n");
  
  // put the redshift into the file (and later next to it the total number of mergers and halos)
  fprintf(fpout_merger,"%f %f %f %f ", zred[izred], tage[izred], (zred[izred+1]-zred[izred]), (tage[izred]-tage[izred+1]));

  
  /*====================================================================
   *   read *_mtree file
   *   -> allocates memory for mtree[] array of structures from scratch
   *====================================================================*/
  mtree = read_mtree(Filename, &nhalos);
  
#ifdef HALOID_RESTRICTION
  if(izred == 0) {
    strcpy(haloidsfile_izred, haloidsfile);
  }
  else {
    sprintf(haloidsfile_izred, "haloids_%06d.txt",izred);
  }
  read_haloids(haloidsfile_izred, mtree, nhalos);
#endif
  
  /*=======================
   *   analyse mergers
   *=======================*/
  halo_nmerger = count_merger(mtree, nhalos, &nhalos_M0);

  // put the total number of mergers and halos next to the redshift in the output file
  fprintf(fpout_merger,"%"PRIu64" %"PRIu64"\n",halo_nmerger,nhalos_M0);
  
  write_merger(fpout_merger, mtree, nhalos);

  
  /*=======================
   *  analyse progenitors
   *=======================*/
  halo_nprogs = count_progs(mtree, nhalos);
  
  // put the total number of halos with progenitor and halos next to the redshift in the output file
  fprintf(fpout_progs,"%"PRIu64" %"PRIu64"\n",halo_nprogs,nhalos);
  
  write_progs(fpout_progs, mtree, nhalos);
  
  
#ifdef HALOID_RESTRICTION
  sprintf(haloidsfile_izred, "haloids_%06d.txt",izred+1);
  write_haloids(haloidsfile_izred, mtree, nhalos);
#endif
  
  
  /*====================================================================
   *   -> free the memory of the mtree[] array
   *====================================================================*/
  /* remove mtree[] from memory */
  for(ihalo=0; ihalo<nhalos; ihalo++) {
    if(mtree[ihalo].progid != NULL) free(mtree[ihalo].progid);
    if(mtree[ihalo].nshared != NULL) free(mtree[ihalo].nshared);
    if(mtree[ihalo].npartprog != NULL) free(mtree[ihalo].npartprog);
    if(mtree[ihalo].iprogmerger != NULL) free(mtree[ihalo].iprogmerger);
  }
  if(mtree != NULL) free(mtree);
  
  fclose(fpout_progs);
  fclose(fpout_merger);
  
  fprintf(stderr," done\n");
  
  
  return(1);
}

/*==================================================================================================
 * read_mtree:
 *
 *       simply reads in the *_mtree file and 
 *       puts it into the array of structures mtree[ihalo].XYZ
 *
 * Note: at this stage we just treat these entries as "lines" -> no connection to halos yet!
 *
 *==================================================================================================*/
MTREE *read_mtree(char *prefix, uint64_t *nhalos)
{
  uint64_t ihalo, iprog;
  char     line[MAXSTRING], mtreename[MAXSTRING];
  FILE    *fpin;
  MTREE   *mtree;
  
  sprintf(mtreename,"%s_mtree",prefix);
  if((fpin = fopen(mtreename,"r")) == NULL) {
    fprintf(stderr,"cannot open  %s\nEXIT\n",mtreename);
    exit(0);
  }
  
  // ignore first two header lines
  fgets(line,MAXSTRING,fpin);
  fgets(line,MAXSTRING,fpin);
  
  // read rest of file into mtree[]
  mtree   = NULL;
  *nhalos = 0;
  
  // read first halo line
  fgets(line,MAXSTRING,fpin);
  while(!feof(fpin)) {
    // add halo to mtree[]
    (*nhalos)++;
    mtree = (MTREEptr) realloc(mtree, (*nhalos)*sizeof(MTREE));
    
    // scan halo properties
    sscanf(line,"%"SCNi64" %"SCNi64" %"SCNi64,
           &(mtree[(*nhalos)-1].haloid),
           &(mtree[(*nhalos)-1].npart),
           &(mtree[(*nhalos)-1].nprog)  );
    
    
    // make room for progenitor properties
    mtree[(*nhalos)-1].progid            = (uint64_t *) calloc(mtree[(*nhalos)-1].nprog, sizeof(uint64_t));
    mtree[(*nhalos)-1].nshared           = (uint64_t *) calloc(mtree[(*nhalos)-1].nprog, sizeof(uint64_t));
    mtree[(*nhalos)-1].npartprog         = (uint64_t *) calloc(mtree[(*nhalos)-1].nprog, sizeof(uint64_t));
#ifdef HALOID_RESTRICTION
    mtree[(*nhalos)-1].follow_progenitor = (int8_t *)   calloc(mtree[(*nhalos)-1].nprog, sizeof(int8_t)); // automatically filled with zeros
#endif
    
    // read progenitor properties
    for(iprog=0; iprog<mtree[(*nhalos)-1].nprog; iprog++) {
      // read from file
      fgets(line,MAXSTRING,fpin);
      sscanf(line,"%"SCNi64" %"SCNi64" %"SCNi64,
             &(mtree[(*nhalos)-1].nshared[iprog]),
             &(mtree[(*nhalos)-1].progid[iprog]),
             &(mtree[(*nhalos)-1].npartprog[iprog])  );
      
    }
    
    // get next halo
    fgets(line,MAXSTRING,fpin);
  }
  fclose(fpin);
  
  fprintf(stderr,"nhalos=%"PRIu64" ... ",(*nhalos));
  
  return(mtree);
}


/*==================================================================================================
 * count_merger:
 *
 *     for each halo stored in mtree[].* we count the number of progenitors
 *     complying with our criterion for a merger
 *
 *==================================================================================================*/
uint64_t count_merger(MTREE *mtree, uint64_t nhalos, uint64_t *nhalos_M0)
{
  
  uint64_t ihalo, iprog, halo_nmerger;
  halo_nmerger = 0;
  double   frac, frac_max;
  
  *nhalos_M0 = 0;
  
  // loop over all haloes
  for(ihalo=0; ihalo<nhalos; ihalo++) {
#ifdef HALOID_RESTRICTION
    if(mtree[ihalo].follow_halo == 1) {
      mtree[ihalo].follow_progenitor[0] = 1; // halos are only written to _mtree if there is at laest one progenitor
#endif
      // reset counter for mergers for this halo to zero
      mtree[ihalo].nmerger = 0;
      
      // set array of merging progenitor ids to NULL
      mtree[ihalo].iprogmerger = NULL;
      
      // 1. check: is the actual progenitor a meaningful progenitor
      if( mtree[ihalo].npart > M0 )
       {
        // count the number of halos above our "mass" criterion
        (*nhalos_M0)++;
        
#ifdef GOTTLOEBER_CRITERION
        // is there at least 1 progenitor
        //if(mtree[ihalo].nprog > 0) // halos are only written to _mtree if there is at laest one progenitor
         {
          
          // mass ratio with the most likely progenitor
          frac_max = FRAC((double)mtree[ihalo].npartprog[0]/(double)mtree[ihalo].npart);
          
          // or search for the largest fraction
          for(iprog=1; iprog<MIN(GOTTLOEBER_CRITERION,mtree[ihalo].nprog); iprog++) {
            frac = FRAC((double)mtree[ihalo].npartprog[0]/(double)mtree[ihalo].npart);
            if(frac > frac_max)
              frac_max = frac;
          }
          
          
          if((1.-frac_max) > MERGER_RATIO) {
            // increase the number of mergers for this halo and store the merging progenitor iprog
            mtree[ihalo].nmerger++;
            mtree[ihalo].iprogmerger                         = (uint64_t *) realloc(mtree[ihalo].iprogmerger, mtree[ihalo].nmerger*sizeof(uint64_t));
            mtree[ihalo].iprogmerger[mtree[ihalo].nmerger-1] = iprog;
            
            // increase the number of total mergers
            halo_nmerger++;
          }
        }
#endif // GOTTLOEBER_CRITERION
        
#ifdef FAKHOURI_CRITERION
        
        // are there at least 2 progenitors to check?
        if(mtree[ihalo].nprog > 1) {
          
          // loop over all progenitor but the actual main progenitor
          for(iprog=1; iprog<MIN(FAKHOURI_CRITERION,mtree[ihalo].nprog); iprog++) {
            
            // 2. check: is this progenitor a meaningful progenitor
            if( mtree[ihalo].npartprog[iprog] > Mi )
             {
              
              // 3. check: physical definition of a merger
              frac = (double)mtree[ihalo].npartprog[iprog]/(double)(mtree[ihalo].npartprog[0]);
              
              if(FRAC(frac) > MERGER_RATIO){
                
                // increase the number of mergers for this halo and store the merging progenitor iprog
                mtree[ihalo].nmerger++;
                mtree[ihalo].iprogmerger                         = (uint64_t *) realloc(mtree[ihalo].iprogmerger, mtree[ihalo].nmerger*sizeof(uint64_t));
                mtree[ihalo].iprogmerger[mtree[ihalo].nmerger-1] = iprog;
                
                // increase the number of total mergers
                halo_nmerger++;
                
                // one halo is only allowed to have one merger at a time
                break; // this will leave the for(iprog)-loop
                
              } // if(MERGER_RATIO)
              
             } // if(mtree[ihalo].npartprog[iprog] > Mi)
            
          } // for(iprog)
          
        } // if(mtree[ihalo].npart > M0)
#endif // FAKHOURI_CRITERION
       } // if(nprog>0)
      
#ifdef HALOID_RESTRICTION
    } // if(follow_halo)
    else {
      mtree[ihalo].nmerger     = 0;
      mtree[ihalo].iprogmerger = NULL;
    }
#endif
    
  } // for(ihalo)
  
  
  return(halo_nmerger);
}

/*==================================================================================================
 * write_merger:
 *
 *     for each halo classified as merger stored in mtree[].* we write something to file
 *
 *==================================================================================================*/
void write_merger(FILE *fpout_merger, MTREE *mtree, uint64_t nhalos)
{
	uint64_t ihalo, imerger, iprog;

  // loop over all haloes
	for(ihalo=0; ihalo<nhalos; ihalo++) {
    
    // did this halo experience at least one merger?
    if(mtree[ihalo].nmerger > 0) {
      
      // write information about this merger into fpout_merger
      fprintf(fpout_merger,"%"PRIu64" %"PRIu64" %"PRIu64"     %"PRIu64" %"PRIu64" %"PRIu64"\n",
              mtree[ihalo].haloid,mtree[ihalo].npart,mtree[ihalo].nmerger,
              mtree[ihalo].progid[0], mtree[ihalo].nshared[0], mtree[ihalo].npartprog[0]);
      
      // loop over all the merger that we counted
      for(imerger=0; imerger<mtree[ihalo].nmerger; imerger++) {
        
        // retrieve the id of the progenitor corresponding to this merger
        iprog = mtree[ihalo].iprogmerger[imerger];

        // write information about this merger into fpout_merger
        fprintf(fpout_merger," %"PRIu64" %"PRIu64" %"PRIu64"\n",mtree[ihalo].progid[iprog],mtree[ihalo].npartprog[iprog],mtree[ihalo].nshared[iprog]);
        
      } // for(iprog)
    } // if(nmerger)
	} // for(ihalo)
}

/*==================================================================================================
 * count_progs:
 *
 *     just count the number of halos with progenitor
 *
 *==================================================================================================*/
uint64_t count_progs(MTREE *mtree, uint64_t nhalos)
{
	uint64_t ihalo;
	uint64_t nprogs;
  
  nprogs = 0;
	for(ihalo=0; ihalo<nhalos; ihalo++) {
    if(mtree[ihalo].nprog > 0) {
      nprogs++;
    }
	}
  
  return(nprogs);
}

/*==================================================================================================
 * write_progs:
 *
 *     just write the direct progneitors and the mass accretion information
 *
 *==================================================================================================*/
void write_progs(FILE *fpout, MTREE *mtree, uint64_t nhalos)
{
	uint64_t ihalo;
  int64_t  dnpart;
  double   fMacc;
  
	for(ihalo=0; ihalo<nhalos; ihalo++) {
    if(mtree[ihalo].nprog > 0) {
      dnpart = (int64_t)mtree[ihalo].npart - (int64_t)mtree[ihalo].npartprog[0];
      fMacc  = (double)dnpart/(double)mtree[ihalo].npart;
      
      fprintf(fpout,"%"PRIu64" %"PRIu64" %"PRIi64" %lf\n", mtree[ihalo].npart, mtree[ihalo].npartprog[0], dnpart, fMacc);
    }
	}
}

/*==================================================================================================
 * read_zred:
 *
 *     read a file than contains the redshifts of the _halos files used with MergerRates
 *     simultaneously calculates the corresponding age of the universe
 *
 *==================================================================================================*/
void read_zred(char zredfile[MAXSTRING], int nFiles, int need_simuparams)
{
  FILE *fpin;
  char  line[MAXSTRING];
  int   nzred;
  float fdummy;
  
  fpin = fopen(zredfile,"r");
  if(fpin == NULL) {
    fprintf(stderr,"Could not open %s.\nABORTING\n",zredfile);
    exit(0);
  }
  
  zred  = NULL;
  tage  = NULL;
  nzred = 0;

  fgets(line,MAXSTRING,fpin);
  while(!feof(fpin)) {
    if (strncmp(line,"#",1) != 0) {
      nzred++;
      zred          = (float *) realloc(zred, nzred*sizeof(float));
      tage          = (float *) realloc(tage, nzred*sizeof(float));
      
      if(need_simuparams == TRUE) {
        sscanf(line,"%f", &(zred[nzred-1]) );
        tage[nzred-1] = calc_t(1./(1.+zred[nzred-1])) * Mpc/H0/Gyr / hubble; // in [Gyr]
      }
      else {
        sscanf(line,"%f %f %f", &(zred[nzred-1]),&fdummy, &(tage[nzred-1]));
      }
    }
    
    fgets(line,MAXSTRING,fpin);
  }
  
  // some sanity check
  if(nzred != nFiles) {
    fprintf(stderr,"you did not provide the correct number of redshifts:\n");
    fprintf(stderr,"  there are %d redshifts in %s\n",nzred,zredfile);
    fprintf(stderr,"  but there are %d redshifts needed!\nABORTING\n",nFiles);
    exit(0);
  }
  
  // in any case, add one final entry to zred[] and tage[] with zero in it to allow easier calcuation of dz and dt
  zred          = (float *) realloc(zred, (nzred+1)*sizeof(float));
  tage          = (float *) realloc(tage, (nzred+1)*sizeof(float));
  zred[nzred]   = 0.0;
  tage[nzred]   = 0.0;
  
  fclose(fpin);
}

/*==================================================================================================
 * test_zredfile:
 *
 *     test whether the zredfile also contains the age of the universe
 *
 *==================================================================================================*/
int test_zredfile(char zredfile[MAXSTRING])
{
  FILE *fpin;
  char  line[MAXSTRING];
  int   need_simuparams;
  float a,b,c;
  
  fpin = fopen(zredfile,"r");
  if(fpin == NULL) {
    fprintf(stderr,"Could not open %s.\nABORTING\n",zredfile);
    exit(0);
  }
  
  fgets(line,MAXSTRING,fpin);
  while(strncmp(line,"#",1) == 0) {
    fgets(line,MAXSTRING,fpin);
  }
  
  if(sscanf(line,"%lf %lf %lf", &a, &b, &c) != 3)
    need_simuparams = TRUE;
  else
    need_simuparams = FALSE;
  
  fclose(fpin);

  return (need_simuparams);
}

#ifdef HALOID_RESTRICTION
/*==================================================================================================
 * read_haloids:
 *==================================================================================================*/
void read_haloids(char haloidsfile_izred[MAXSTRING], MTREE *mtree, uint64_t nhalos)
{
  FILE    *fp;
  char     line[MAXSTRING];
  uint64_t haloid, ihalo;
  MTREE   *itmp_mtree;
  uint64_t nlost, nfound;
  
  fp = fopen(haloidsfile_izred,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",haloidsfile_izred);
    exit(0);
  }
  
  // sort mtree[] with respects to mtree[].haloid
  qsort((void *)mtree, nhalos, sizeof(MTREE), qcompareHaloIDs);

  nlost  = 0;
  nfound = 0;
  fgets(line,MAXSTRING,fp);
  while(!feof(fp)) {
    sscanf(line,"%"SCNi64,&haloid);
    
    // find pointer to mtree[] memory where haloid is located
    itmp_mtree = (MTREE *)bsearch(&(haloid), mtree, nhalos, sizeof(MTREE), bcompareMtreeIDs);

    if(itmp_mtree == NULL) {
      nlost++;
    }
    else {
      nfound++;
      
      // position in mtree[] array
      ihalo = itmp_mtree - mtree;
      
      // flag halo to be followed
      mtree[ihalo].follow_halo = 1;
    }
    
    fgets(line,MAXSTRING,fp);
  }
  
  fprintf(stderr," (all haloids read in, found %"PRIu64" lost %"PRIu64")",nfound,nlost);

  fclose(fp);
}

/*==================================================================================================
 * write_haloids:
 *==================================================================================================*/
void write_haloids(char haloidsfile_izred[MAXSTRING], MTREE *mtree, uint64_t nhalos)
{
  FILE    *fp;
  uint64_t ihalo, iprog;
  
  fp = fopen(haloidsfile_izred,"w");
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",haloidsfile_izred);
    exit(0);
  }
  
  for(ihalo=0; ihalo<nhalos; ihalo++) {
    for(iprog=0; iprog<mtree[ihalo].nprog ;iprog++) {
      if(mtree[ihalo].follow_progenitor[iprog] == 1) {
        fprintf(fp,"%"PRIu64"\n",mtree[ihalo].progid[iprog]);
      }
    }
  }
  
  fclose(fp);
}

/*==============================================================================
 *  compare halo ids (used with qsort)
 *==============================================================================*/
int qcompareHaloIDs(const void *mtree1, const void *mtree2)
{
	uint64_t n1, n2;
  
	n1 = ((MTREE *)mtree1)->haloid;
	n2 = ((MTREE *)mtree2)->haloid;
  
	return n1 < n2 ? -1 : (n1 > n2 ? 1 : 0);
}

/*==============================================================================
 *  compare particle ids (used with bsearch)
 *==============================================================================*/
int bcompareMtreeIDs(const void *id, const void *mtree)
{
  const uint64_t haloid = ((MTREE *)mtree)->haloid;
	const uint64_t *i     = (const uint64_t *)id;
  
	return *i < haloid ? -1 : (*i > haloid ? 1 : 0);
}
#endif

