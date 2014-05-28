/*==================================================================================================
 *  MergerTree:   Merger Tree AHF_particles files
 *
 *
 *  input:    - how often to perform
 *            - 2x _particles files
 *
 *  output:   - 1x _mtree file
 *
 *
 * it is checked what halos in file2 make up the halos in file1, i.e.
 *
 *   file1   file2
 *
 *    0        0
 *    0       17
 *    0       31    -> halo #0 in file1 shares particles with halos #0,17,31 in file2
 *    1        2
 *    1       12
 *    1        4    -> halo #1 in file1 shares particles with halos #2,12,4  in file2
 *       etc.
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
#include <time.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"
#include "../src/libutility/utility.h"
#include "../src/libio/xstring.h"

/*-------------------------------------------------------------------------------------
 *                                     DEFINES
 *-------------------------------------------------------------------------------------*/
#define MINCOMMON      10             // we only cross-correlate haloes if they at least share MINCOMMON particles
#define ONLY_USE_PTYPE 1              // restrict analysis to particles of this type (1 = dark matter)
//#define MTREE_BOTH_WAYS               // make sure that every halo has only one descendant
//#define SUSSING2013                   // write _mtree in format used for Sussing Merger Trees 2013
//#define EXCLUSIVE_PARTICLES           // each particle is only allowed to belong to one object (i.e. the lowest mass one)
//#define WITH_QSORT                    // uses qsort() instead of indexx() when ordering the progenitors according to merit function

// support for AHF's MPI output
//#define READ_MPARTICLES               // support to read multiple _particle files, they must be of the latest AHF _particles file format!
#define NDIGITS                 4     // number of digits to be used for fileid


// support of ancient formats (without nhalos line, without haloids, whatever)
//#define THERE_IS_NO_NHALOS_LINE
//#define READ_HALOIDS_FROM_FILE
//#define USE_LINENUMBER_AS_HALOID      // overwrites(!) haloid as found in _particles
//#define CLUES_CROSS_CORRELATION         // fiddle with IDs for CLUES data to allow for CroCos

/*-------------------------------------------------------------------------------------
 *                                   DEPENDENCIES
 *-------------------------------------------------------------------------------------*/
#if (defined EXCLUSIVE_PARTICLES && defined WITH_OPENMP)
#define WITH_OPENMP2
#endif

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct MTREE *MTREEptr;
typedef struct MTREE
{
  uint64_t haloid[2];
  uint64_t id[2];
  uint64_t npart[2];
  uint64_t common;
  double   merit;
} MTREE;

typedef struct HALOS *HALOptr;
typedef struct HALOS
{
  uint64_t  haloid;
  uint64_t  npart;
  uint64_t *Pid;
  
  uint64_t  ncroco;
  MTREEptr  mtree;
}HALOS;

typedef struct PARTS *PARTptr;
typedef struct PARTS
{
  uint64_t  nhalos;
  uint64_t *Hid;
}PARTS;

/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/
HALOptr     halos[2];
PARTptr     parts[2];
uint64_t    nHalos[2];
uint64_t    PidMax[2]={0,0};
uint64_t    PidMin=(1<<62);

/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
int      read_particles         (char filename[MAXSTRING], int isimu);
int      particle_halo_mapping  (int  isimu);
int      cross_correlation      (char OutFile[MAXSTRING]);
int      create_mtree           (uint64_t ihalo, int  isimu0,int isimu1);
int      create_mtree_qsort     (uint64_t ihalo, int isimu0, int isimu1);
int      clean_connection       (uint64_t ihalo, int isimu0, int isimu1);
int      write_mtree            (char OutFile[MAXSTRING]);
uint64_t max_merit              (uint64_t ihalo, int isimu);
void     construct_filename     (char filename[MAXSTRING], int32_t, char infile[MAXSTRING]);
uint64_t count_halos            (char filename[MAXSTRING], int32_t);
int32_t  count_particles_files  (char filename[MAXSTRING]);
int      merit_sort             (const void *m1, const void *m2);

/*==================================================================================================
 * main:
 *
 *       simply a wrapper for successive calls to create_mtree()
 *
 *==================================================================================================*/
int main()
{
  int      i, nFiles, isimu;
  uint64_t ihalo, ipart;
  char   **HaloFile;
  char   **OutFile;
  
  /*==================================================================*
   *                          USER INTERFACE                          *
   *==================================================================*/
  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"  construct a cross-correlation between consecutive *_particles files\n");
  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"\nPlease give number of particles files (default=2):      ");
  scanf("%d", &nFiles);
  fprintf(stderr,"%d\n",nFiles);
  
  /* allocate memory for nFiles filenames, each of size MAXSTRING */
  HaloFile = (char **) calloc(nFiles, sizeof(char *));
  OutFile  = (char **) calloc(nFiles, sizeof(char *));
  for(i=0; i<nFiles; i++)
   {
    HaloFile[i] = (char *) calloc(MAXSTRING, sizeof(char));
    OutFile[i]  = (char *) calloc(MAXSTRING, sizeof(char));
   }
  
  /* read input filenames from stdin */
  for(i=0; i<nFiles; i++)
   {
    fprintf(stderr,"Please give name of %5d. *_particles file:            ",i+1);
    scanf("%s", HaloFile[i]);
    fprintf(stderr,"%s\n",HaloFile[i]);
   }
  
  /* read output filenames from stdin */
  for(i=0; i<nFiles-1; i++)
   {
    fprintf(stderr,"Please give prefix for %5d. output file:                 ",i+1);
    scanf("%s", OutFile[i]);
    fprintf(stderr,"%s\n",OutFile[i]);
   }
  fprintf(stderr,"\n");
  
  
  /*======================================================================*
   *  CREATE CROSS-CORRELATIONS                                           *
   *======================================================================*/
  /* read the first file into memory */
  fprintf(stderr,"Startup:\n");
  read_particles(HaloFile[0], 0);
#ifdef READ_HALOIDS
  read_haloids(HaloFile[0],0);
#endif
  particle_halo_mapping(0);
  fprintf(stderr,"\n");
  
  for(i=0; i<nFiles-1; i++) {
    
    /* be verbose */
    fprintf(stderr,"Correlating '%s' to '%s'\n           -> writing to '%s'\n",
            HaloFile[i],HaloFile[i+1],OutFile[i]);
    
    /* read the next file into memory */
    read_particles(HaloFile[i+1], 1);
    particle_halo_mapping(1);
    
    /* cross correlate HaloFile[i] to HaloFile[i+1] */
    cross_correlation(OutFile[i]);
    
    /* be verbose */
    fprintf(stderr,"  o making file 1 the new file 0 ...");
    
    /* remove HaloFile[0] from memory */
    for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
      if(halos[0][ihalo].Pid   != NULL) {
        free(halos[0][ihalo].Pid);
        halos[0][ihalo].Pid = NULL;
      }
      if(halos[0][ihalo].mtree != NULL) {
        free(halos[0][ihalo].mtree);
        halos[0][ihalo].mtree = NULL;
      }
    }
    for(ipart=0; ipart<PidMax[0]; ipart++) {
      if(parts[0][ipart].Hid != NULL) {
        free(parts[0][ipart].Hid);
        parts[0][ipart].Hid = NULL;
      }
    }
    if(halos[0] != NULL) {
      free(halos[0]);
      halos[0] = NULL;
    }
    if(parts[0] != NULL) {
      free(parts[0]);
      parts[0] = NULL;
    }
    
    /* make HaloFile[i+1] the new HaloFile[i] */
    nHalos[0] = nHalos[1];
    halos[0]  = halos[1];
    parts[0]  = parts[1];
    PidMax[0] = PidMax[1];
    
    /* be verbose */
    fprintf(stderr," done\n");
  } // for(nFiles)
  
  
  /*==================================================================*
   *                             CLEANUP                              *
   *==================================================================*/
  fprintf(stderr,"\nCleaning up ... ");
  
  /* remove HaloFile[0] from memory */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    if(halos[0][ihalo].Pid   != NULL) free(halos[0][ihalo].Pid);
    if(halos[0][ihalo].mtree != NULL) free(halos[0][ihalo].mtree);
  }
  for(ipart=0; ipart<PidMax[0]; ipart++) {
    if(parts[0][ipart].Hid != NULL) free(parts[0][ipart].Hid);
  }
  if(halos[0] != NULL) free(halos[0]);
  if(parts[0] != NULL) free(parts[0]);
  
  /* remove filename storage */
  for(i=0; i<nFiles; i++)
   {
    if(HaloFile[i]) free(HaloFile[i]);
    if(OutFile[i])  free(OutFile[i]);
   }
  if(HaloFile) free(HaloFile);
  if(OutFile)  free(OutFile);
  
  printf("finished\n");
  return(1);
}


/*==================================================================================================
 * read_particles:
 *
 * read the file storing the particle IDs for each halo
 *
 *      nHalos = number of halos found in file
 *      Pid    = id's of all those particles
 *
 *==================================================================================================*/
int read_particles(char filename[MAXSTRING], int isimu)
{
  FILE     *fpin, *fphids;
  char      line[MAXSTRING], infile[MAXSTRING], hidsname[MAXSTRING];
  int64_t   ihalo, jhalo;
  uint64_t  nPartInHalo, nPartInUse, ipart, jpart, Pid, Ptype, haloid, haloid_from_file;
  uint64_t  PidMin_local=(1<<62);
  uint64_t  PidMax_local=0;
  clock_t   elapsed;
  
  elapsed = clock();
  
#ifdef READ_MPARTICLES
  int32_t  nfiles, ifile;
  uint64_t nHalosInFile;
  
  fprintf(stderr,"  o reading multiple files %s ",filename);
  
  // count the number of particles files/snapshot
  nfiles = count_particles_files(filename);
  fprintf(stderr," (nfiles=%"PRIi32,nfiles);
  if(nfiles == 0) {
    fprintf(stderr,"Could not open multiple _particles file complying with the required filename convention\nABORTING\n");
    exit(0);
  }
  
  // accumulate the number of halos by summing the first line in all _particles files
  nHalos[isimu] = count_halos(filename, nfiles);
  
  fprintf(stderr," nhalos=%"PRIu64") file: ",nHalos[isimu]);
  
  // allocate memory for haloes as one block
  halos[isimu]  = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));
  if(halos[isimu] == NULL) {
    fprintf(stderr,"\nCould not allocate memory for halos[] array of isimu = %d (size = %f GB)\nAborting\n",
            isimu,(float)(nHalos[isimu]*sizeof(HALOS)/1024./1024./1024.));
    exit(0);
  }
  
  // reset halo and total particle counter
  ihalo = 0;
  for(ifile=0; ifile<nfiles; ifile++) {
    
    fprintf(stderr,"%"PRIi32,ifile);

    // open file
    construct_filename(filename, ifile, infile);
    fpin = fopen(infile,"r");

    // read total number of haloes in file
    fgets(line,MAXSTRING,fpin);
    sscanf(line,"%"SCNu64, &nHalosInFile);
    
    fprintf(stderr," (%"PRIu64") ",nHalosInFile);

    // loop over all haloes in this file
    for(jhalo=0; jhalo<nHalosInFile; jhalo++) {
      
      // get actual information
      fgets(line,MAXSTRING,fpin);
      sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid);
      
      // store haloid
      halos[isimu][ihalo].haloid = haloid;

      // loop over all particles in halo
      nPartInUse              = 0;
      halos[isimu][ihalo].Pid = NULL;
      for(ipart=0; ipart<nPartInHalo; ipart++)
       {
        /* read line containing ID and possibly some more information */
        fgets(line,MAXSTRING,fpin);
        
        /* check whether we are able to read a meaningful particle type, too */
        if(sscanf(line,"%"SCNu64" %"SCNu64, &Pid, &Ptype) == 1) {
          /* if not, set Ptype to 1 as this is the type we will use below */
          Ptype = 1;
        }
        else if(Ptype < 0 || Ptype > abs(PDMbndry)) {
          /* not a meaningful type, maybe something else has been stored? */
          Ptype = 1;
        }
        
        // here we can restrict the cross-correlation to a ceratain sub-set of all particles
#ifdef ONLY_USE_PTYPE
        if(Ptype == ONLY_USE_PTYPE)
#endif
         {
          halos[isimu][ihalo].Pid = (uint64_t *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+2)*sizeof(uint64_t));
          if(halos[isimu][ihalo].Pid == NULL) {
            fprintf(stderr,"\n read_particles: could not realloc() halos[%d][%ld].Pid for %"PRIu64" particles (Pid=%ld)\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1),halos[isimu][ihalo].Pid);
            exit(-1);
          }

          halos[isimu][ihalo].Pid[nPartInUse] = Pid;
          
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax[isimu])       PidMax[isimu] = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin)              PidMin        = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax_local)        PidMax_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin_local)        PidMin_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
          
          nPartInUse++;
         } // if(Ptype)
        
       } // for(ipart)
      
      // store number of particles in halo
      halos[isimu][ihalo].npart = nPartInUse;

      // move to next halo in overall list
      ihalo++;
      
    } // for(jhalo)
    
    // close file
    fclose(fpin);
    
  } // for(ifile)
  
#else // READ_MPARTICLES
  
  fprintf(stderr,"  o reading file %s ...",filename);
  
  fpin = fopen(filename,"r");
  if(fpin == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting!\n",filename);
    exit(0);
   }

#ifdef READ_HALOIDS_FROM_FILE
  sprintf(hidsname,"%s_hids",filename);
  fphids = fopen(hidsname,"r");
  if(fphids == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting!\n",hidsname);
    exit(0);
   }
#endif
  
  /* reset all variables */
  nHalos[isimu] = 0;
  ihalo         = -1;
  halos[isimu]  = NULL;
  
  /* get the first line from file */
  fgets(line,MAXSTRING,fpin);
  
#ifndef THERE_IS_NO_NHALOS_LINE
  /* for AHF_particles files the first line is numGoodHalos which we will happily ignore and count that number ourselves */
  if(strncmp(line,"#",1) != 0 && sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)
    fgets(line,MAXSTRING,fpin);
#endif
  
  do {
    if(strncmp(line,"#",1) != 0)
     {
      /* has a haloid been written */
      if(sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)
       {
        /* if not, just get the number of particles */
        sscanf(line,"%"SCNu64, &nPartInHalo);
        
        /* and use halo counter as id */
        haloid = ihalo+1; // +1, because ihalo has not been incremented yet!
       }
#ifdef USE_LINENUMBER_AS_HALOID
      haloid = ihalo+1; // +1, because ihalo has not been incremented yet!
#endif
#ifdef READ_HALOIDS_FROM_FILE
      fscanf(fphids,"%"SCNi64,&haloid_from_file);
      haloid = haloid_from_file;
#endif
      
      /* found yet another halo */
      ihalo++;
      nHalos[isimu] += 1;
      halos[isimu]   = (HALOptr) realloc(halos[isimu], (nHalos[isimu]+1)*sizeof(HALOS));
      
      /* store haloid */
      halos[isimu][ihalo].haloid = haloid;
      
      /* halos[][].Pid will be incrementally filled using realloc() */
      halos[isimu][ihalo].Pid   = NULL;
      halos[isimu][ihalo].mtree = NULL;
      
      /* read all their id's */
      nPartInUse = 0;
      for(ipart=0; ipart<nPartInHalo; ipart++)
       {
        /* read line containing ID and possibly some more information */
        fgets(line,MAXSTRING,fpin);
        
        /* check whether we are able to read a meaningful particle type, too */
        if(sscanf(line,"%"SCNu64" %"SCNu64, &Pid, &Ptype) == 1) {
          /* if not, set Ptype to 1 as this is the type we will use below */
          Ptype = 1;
        }
        else if(Ptype < 0 || Ptype > abs(PDMbndry)) {
          /* not a meaningful type, maybe something else has been stored? */
          Ptype = 1;
        }
        
        // here we can restrict the cross-correlation to a ceratain sub-set of all particles
#ifdef ONLY_USE_PTYPE
        if(Ptype == ONLY_USE_PTYPE)
#endif
         {
          halos[isimu][ihalo].Pid             = (uint64_t *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(uint64_t));
          if(halos[isimu][ihalo].Pid == NULL) {
            fprintf(stderr,"read_particles: could not realloc() halos[%d][%ld].Pid for %"PRIu64"particles\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1));
            exit(-1);
          }
#ifdef CLUES_CROSS_CORRELATION
          if(isimu == 0)
            Pid -= 52953088;
#endif
          halos[isimu][ihalo].Pid[nPartInUse] = Pid;
          
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax[isimu])       PidMax[isimu] = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin)              PidMin        = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax_local)        PidMax_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin_local)        PidMin_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
          
          nPartInUse++;
         }
       }
      
      /* store number of particles in halo */
      halos[isimu][ihalo].npart = nPartInUse;
     }
  } while( fgets(line,MAXSTRING,fpin) != NULL);
  
  fclose(fpin);
#ifdef READ_HALOIDS_FROM_FILE
  fclose(fphids);
#endif
  
#endif // READ_MPARTICLES
  
  elapsed = clock()-elapsed;

  fprintf(stderr," done in %4.2f sec. (nhalos = %"PRIu64", full ID range = %"PRIu64" -> %"PRIu64", local ID range = %"PRIu64" -> %"PRIu64")\n",
          (float)elapsed/CLOCKS_PER_SEC,nHalos[isimu],PidMin,PidMax[isimu],PidMin_local,PidMax_local);
  

  return(1);
}


/*==================================================================================================
 * particle_halo_mapping:
 *
 *  for each particle remember to which halo(s) it belongs
 *
 *==================================================================================================*/
int particle_halo_mapping(int isimu)
{
  int64_t  ihalo;         // the downwards for-loop is running until ihalo=-1
  uint64_t ipart, jpart, PidMax_global;
  clock_t  elapsed;
  
  elapsed = clock();
  
  PidMax_global = MAX(PidMax[0],PidMax[1]);
  fprintf(stderr,"  o creating particle<->halo mapping for file %d (PidMax=%"PRIu64", PidMax_global=%"PRIu64") ... ",isimu,PidMax[isimu],PidMax_global);

  parts[isimu] = (PARTptr) calloc((PidMax_global+1), sizeof(PARTS)); // +1 because we are accessing an array like [PidMax_global] !
  if(parts[isimu] == NULL) {
    fprintf(stderr,"\nCould not allocate memory for parts[] array of isimu = %d (size = %f GB)\nAborting\n",
            isimu,(float)((PidMax_global+1)*sizeof(PARTS)/1024./1024./1024.));
    exit(0);
  }
  
  /* recording every halo it belongs to: running from low to high mass objects to allow for unique assignment! */
#ifdef WITH_OPENMP2
#pragma omp parallel for schedule(dynamic) private(ihalo,jpart,ipart) shared(nHalos,halos,parts,isimu)
#endif
  for(ihalo=nHalos[isimu]-1; ihalo>=0; ihalo--)
   {
    for(jpart=0; jpart<halos[isimu][ihalo].npart; jpart++)
     {
      ipart = halos[isimu][ihalo].Pid[jpart];

#ifdef EXCLUSIVE_PARTICLES
      if(parts[isimu][ipart].nhalos == 0)
#endif
       {
        parts[isimu][ipart].nhalos++;
        parts[isimu][ipart].Hid = (uint64_t *) realloc(parts[isimu][ipart].Hid, parts[isimu][ipart].nhalos*sizeof(uint64_t)); // valgrind says "16 bytes lost"
        if(parts[isimu][ipart].Hid == NULL) {
          fprintf(stderr,"particle_halo_mapping(): could not realloc() memory for Hid array\n");
          exit(0);
        }
        
        parts[isimu][ipart].Hid[parts[isimu][ipart].nhalos-1] = ihalo;
       }
      
     }
   }
  
  elapsed = clock()-elapsed;
  fprintf(stderr,"done in %4.2f sec.\n",(float)elapsed/CLOCKS_PER_SEC);
  
  return(1);
}

/*==================================================================================================
 * cross_correlation:
 *
 *  for each halo at isimu=0 figure out how many particles are in common with khalo at isimu=1
 *
 *==================================================================================================*/
int cross_correlation(char OutFile[MAXSTRING])
{
  uint64_t  ihalo;
  clock_t  elapsed;
  
  
  /*---------------------------------------------------------
   * backwards correlation
   *---------------------------------------------------------*/
  elapsed = clock();
  fprintf(stderr,"  o generating cross-correlation 0->1 for %"PRIu64" haloes ...",nHalos[0]);
#ifdef WITH_OPENMP
#  pragma omp parallel for schedule (dynamic)	shared(nHalos) private(ihalo)
#endif
  /* cross-correlation simu0->simu1 */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    create_mtree(ihalo, 0, 1);
  }
  elapsed = clock()-elapsed;
  fprintf(stderr," done in %4.2f sec.\n", (double)elapsed/CLOCKS_PER_SEC);
  
  
  
  
#ifdef MTREE_BOTH_WAYS
  
  /*---------------------------------------------------------
   * forward correlation
   *---------------------------------------------------------*/
  elapsed = clock();
  fprintf(stderr,"  o generating cross-correlation 1->0 for %"PRIu64" haloes ...",nHalos[1]);
#ifdef WITH_OPENMP
#  pragma omp parallel for schedule (dynamic)	shared(nHalos) private(ihalo)
#endif
  /* cross-correlation simu0<-simu1 */
  for(ihalo=0; ihalo<nHalos[1]; ihalo++) {
    create_mtree(ihalo, 1, 0);
  }
  elapsed = clock()-elapsed;
  fprintf(stderr," done in %4.2f sec.\n", (double)elapsed/CLOCKS_PER_SEC);
  
  
  
  elapsed = clock();
  fprintf(stderr,"  o removing network connections ...");
#ifdef WITH_OPENMP2
#  pragma omp parallel for schedule (dynamic)	shared(nHalos) private(ihalo)
#endif
  /* clean connections simu0->simu1 */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    clean_connection(ihalo, 0, 1);
  }
  elapsed = clock()-elapsed;
  fprintf(stderr," done in %4.2f sec.\n", (double)elapsed/CLOCKS_PER_SEC);
  
#endif // MTREE_BOTH_WAYS
  
  
  
  write_mtree(OutFile);
  
  return(1);
}

/*==================================================================================================
 * clean_connection
 *==================================================================================================*/
int clean_connection(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t jhalo, icroco, ncroco_new;
  MTREE    *mtree;
#ifdef DEBUG
  uint64_t idesc;
#endif
  
  /* count number of new crocos */
  ncroco_new = 0;
  mtree      = NULL;
  
  /* loop over all cross-correlated haloes */
  for(icroco=0; icroco<halos[isimu0][ihalo].ncroco; icroco++) {
    jhalo = halos[isimu0][ihalo].mtree[icroco].id[1];
    
    /* check whether the present halo is the most likely descendant of this progenitor */
    if(max_merit(jhalo, isimu1) == ihalo) {
      // keep jhalo in mtree-list
      ncroco_new++;
      mtree = (MTREEptr) realloc(mtree, (ncroco_new+1)*sizeof(MTREE));
      
      mtree[ncroco_new-1].id[0]     = halos[isimu0][ihalo].mtree[icroco].id[0];
      mtree[ncroco_new-1].haloid[0] = halos[isimu0][ihalo].mtree[icroco].haloid[0];
      mtree[ncroco_new-1].npart[0]  = halos[isimu0][ihalo].mtree[icroco].npart[0];
      mtree[ncroco_new-1].common    = halos[isimu0][ihalo].mtree[icroco].common;
      mtree[ncroco_new-1].id[1]     = halos[isimu0][ihalo].mtree[icroco].id[1];
      mtree[ncroco_new-1].haloid[1] = halos[isimu0][ihalo].mtree[icroco].haloid[1];
      mtree[ncroco_new-1].npart[1]  = halos[isimu0][ihalo].mtree[icroco].npart[1];
      
#ifdef DEBUG
      fprintf(stderr,"icroco=%ld (of %ld) for ihalo=%ld: jhalo=%ld is     a real progenitor of ihalo=%ld (jhalo has %ld descendants)\n",
              icroco,halos[isimu0][ihalo].ncroco,ihalo,
              halos[isimu1][jhalo].haloid,halos[isimu0][ihalo].haloid,
              halos[isimu1][jhalo].ncroco);
#endif
    }
    else {
      // remove jhalo from mtree-list and hence do not add it to the new mtree[] list
#ifdef DEBUG
      fprintf(stderr,"icroco=%ld (of %ld) for ihalo=%ld: jhalo=%ld is NOT a real progenitor of ihalo=%ld (jhalo has %ld descendants)\n",
              icroco,halos[isimu0][ihalo].ncroco,ihalo,
              halos[isimu1][jhalo].haloid,halos[isimu0][ihalo].haloid,
              halos[isimu1][jhalo].ncroco);
      for(idesc=0; idesc<halos[isimu1][jhalo].ncroco; idesc++) {
        fprintf(stderr,"    idesc=%ld haloidesc=%ld\n",idesc,halos[isimu1][jhalo].mtree[idesc].haloid[1]);
      }
#endif
    }
  } // for(icroco)
  
  /* replace halos[isimu0][ihalo].mtree[] with new structure array */
  if(halos[isimu0][ihalo].mtree != NULL) {
    free(halos[isimu0][ihalo].mtree);
    halos[isimu0][ihalo].mtree = NULL;
  }
  halos[isimu0][ihalo].ncroco = ncroco_new;
  halos[isimu0][ihalo].mtree  = mtree;
}

/*==================================================================================================
 * max_merit
 *==================================================================================================*/
uint64_t max_merit(uint64_t jhalo, int isimu)
{
  uint64_t ihalo;
  
  /* mtree[] is ordered by merit and hence we only need to check the first entry */
  if(halos[isimu][jhalo].ncroco > 0) {
    return(halos[isimu][jhalo].mtree[0].id[1]);
  }
  else {
#ifdef DEBUG
    fprintf(stderr,"jhalo=%ld in isimu=%d does not point to anywhere!?\n",jhalo,isimu);
#endif
    return(0);
  }
}

#ifndef WITH_QSORT
/*==================================================================================================
 * create_mtree_index
 *==================================================================================================*/
int create_mtree(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t  jhalo, khalo, ipart, jpart, ncroco, icroco;
  int64_t   jcroco;
  
  uint64_t      *common;
  MTREE         *mtree;
  double        *merit;
  long unsigned *idx;
  
  // reset the actual mtree[] pointer
  if(halos[isimu0][ihalo].mtree != NULL)
    free(halos[isimu0][ihalo].mtree);
  halos[isimu0][ihalo].mtree = NULL;
  
  /* temporary array pointers */
  common = NULL;
  mtree  = NULL;
  merit  = NULL;
  idx    = NULL;
  
  /* common[] records how many particles ihalo(isimu0) has in common with khalo(isimu1) */
  common = (uint64_t *) calloc(nHalos[isimu1], sizeof(uint64_t));
  
  for(jpart=0; jpart<halos[isimu0][ihalo].npart; jpart++) {
    ipart = halos[isimu0][ihalo].Pid[jpart];
    
    //fprintf(stderr,"jpart=%"PRIu64" ipart=%"PRIu64"\n",jpart,ipart);
    
    /* ipart belongs to nhalos halos in isimu1 */
    for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
      khalo          = parts[isimu1][ipart].Hid[jhalo];
      common[khalo] += 1;
    }
  }
  
  /* determine number of credible cross-correlations */
  ncroco = 0;
  for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
    if(common[khalo] > MINCOMMON)
      ncroco++;
  }
  halos[isimu0][ihalo].ncroco = ncroco;
  
  /* does not make sense to continue if there are no cross-correlations */
  if(ncroco > 0) {
    
    /* allocate memory for cross-correlations */
    halos[isimu0][ihalo].mtree  = (MTREEptr) calloc(ncroco, sizeof(MTREE));
    mtree  = (MTREEptr)        calloc(ncroco, sizeof(MTREE));
    idx    = (long unsigned *) calloc(ncroco, sizeof(long unsigned));
    merit  = (double *)        calloc(ncroco, sizeof(double));
    
    /* store cross-correlations temporarily for sorting by indexx() */
    icroco = 0;
    for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
      if(common[khalo] > MINCOMMON){
        mtree[icroco].id[0]     = ihalo;
        mtree[icroco].haloid[0] = halos[isimu0][ihalo].haloid;
        mtree[icroco].npart[0]  = halos[isimu0][ihalo].npart;
        mtree[icroco].common    = common[khalo];
        mtree[icroco].id[1]     = khalo;
        mtree[icroco].haloid[1] = halos[isimu1][khalo].haloid;
        mtree[icroco].npart[1]  = halos[isimu1][khalo].npart;
        mtree[icroco].merit     = pow2((double)common[khalo])/((double)halos[isimu0][ihalo].npart*(double)halos[isimu1][khalo].npart);
        
        merit[icroco] = mtree[icroco].merit;
        
        icroco++;
      }
    }
    
    /* order by merit function */
    indexx((long unsigned)ncroco, merit-1, idx-1);
    
    /* store in halos[isimu0][*].mtree structure (descending order!) */
    for(jcroco=0; jcroco<ncroco; jcroco++) {
      icroco = idx[ncroco-1-jcroco]-1;
      
      /* store mtree[] inside halos[][] structure in correct order */
      halos[isimu0][ihalo].mtree[jcroco].id[0]    = mtree[icroco].id[0];
      halos[isimu0][ihalo].mtree[jcroco].haloid[0]= mtree[icroco].haloid[0];
      halos[isimu0][ihalo].mtree[jcroco].npart[0] = mtree[icroco].npart[0];
      halos[isimu0][ihalo].mtree[jcroco].common   = mtree[icroco].common;
      halos[isimu0][ihalo].mtree[jcroco].id[1]    = mtree[icroco].id[1];
      halos[isimu0][ihalo].mtree[jcroco].haloid[1]= mtree[icroco].haloid[1];
      halos[isimu0][ihalo].mtree[jcroco].npart[1] = mtree[icroco].npart[1];
    }
    
    /* free temporary structures */
    if(mtree) {
      free(mtree);
      mtree = NULL;
    }
    if(idx) {
      free(idx);
      idx = NULL;
    }
    if(merit) {
      free(merit);
      merit = NULL;
    }
  } // if(ncroco)
  
  /* free temporary structures */
  if(common) {
    free(common);
    common = NULL;
  }
}
#else // WITH_QSORT
/*==================================================================================================
 * create_mtree_qsort
 *==================================================================================================*/
int create_mtree(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t  jhalo, khalo, ipart, jpart, ncroco, icroco;
  int64_t   jcroco;
  
  uint64_t      *common;
  
  // reset the actual mtree[] pointer
  if(halos[isimu0][ihalo].mtree != NULL)
    free(halos[isimu0][ihalo].mtree);
  halos[isimu0][ihalo].mtree = NULL;

  /* temporary array pointers */
  common = NULL;

  /* common[] records how many particles ihalo(isimu0) has in common with khalo(isimu1) */
  common = (uint64_t *) calloc(nHalos[isimu1], sizeof(uint64_t));
  
  for(jpart=0; jpart<halos[isimu0][ihalo].npart; jpart++) {
    ipart = halos[isimu0][ihalo].Pid[jpart];
    
    /* ipart belongs to nhalos halos in isimu1 */
    for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
      khalo          = parts[isimu1][ipart].Hid[jhalo];
      common[khalo] += 1;
    }
  }
  
  /* determine number of credible cross-correlations */
  ncroco = 0;
  for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
    if(common[khalo] > MINCOMMON)
      ncroco++;
  }
  halos[isimu0][ihalo].ncroco = ncroco;
  
  /* does not make sense to continue if there are no cross-correlations */
  if(ncroco > 0) {
    /* allocate memory for cross-correlations */
    halos[isimu0][ihalo].mtree  = (MTREEptr) calloc(ncroco, sizeof(MTREE));
    
    /* store cross-correlations temporarily for sorting by indexx() */
    icroco = 0;
    for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
      if(common[khalo] > MINCOMMON){
        halos[isimu0][ihalo].mtree[icroco].id[0]     = ihalo;
        halos[isimu0][ihalo].mtree[icroco].haloid[0] = halos[isimu0][ihalo].haloid;
        halos[isimu0][ihalo].mtree[icroco].npart[0]  = halos[isimu0][ihalo].npart;
        halos[isimu0][ihalo].mtree[icroco].common    = common[khalo];
        halos[isimu0][ihalo].mtree[icroco].id[1]     = khalo;
        halos[isimu0][ihalo].mtree[icroco].haloid[1] = halos[isimu1][khalo].haloid;
        halos[isimu0][ihalo].mtree[icroco].npart[1]  = halos[isimu1][khalo].npart;
        halos[isimu0][ihalo].mtree[icroco].merit     = pow2((double)common[khalo])/((double)halos[isimu0][ihalo].npart*(double)halos[isimu1][khalo].npart);
        
        icroco++;
      }
    }

    /* qsort particles according to merit function */
    qsort((void *)halos[isimu0][ihalo].mtree, ncroco, sizeof(MTREE), merit_sort);
    
  } else { // if(ncroco)
    halos[isimu0][ihalo].mtree = NULL;
  }
  

  /* free temporary structures */
  if(common) {
    free(common);
    common = NULL;
  }
}
#endif // WITH_QSORT

/*==================================================================================================
 * write_mtree:
 *==================================================================================================*/
int write_mtree(char OutFile[MAXSTRING])
{
  uint64_t  ihalo;
  int64_t   icroco;
  FILE *fpout, *fpout_idx;
  char outname[MAXSTRING], outname_idx[MAXSTRING];
  clock_t   elapsed;
  
  elapsed = clock();
  fprintf(stderr,"  o writing cross-correlation for %"PRIu64" haloes ...",nHalos[0]);
  
  sprintf(outname,"%s_mtree",OutFile);
  strcpy(outname_idx, outname);
  strcat(outname_idx, "_idx");
  
  fpout = fopen(outname,"w");
  if(fpout == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",outname);
    exit(0);
  }
  
  fpout_idx = fopen(outname_idx,"w");
  if(fpout_idx == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",outname_idx);
    exit(0);
  }
  
  
#ifdef SUSSING2013
  fprintf(fpout,"%"PRIu64"\n",nHalos[0]);
#else // SUSSING2013
  fprintf(fpout,"#   HaloID(1)   HaloPart(2)  NumProgenitors(3)\n");
  fprintf(fpout,"#      SharedPart(1)    HaloID(2)   HaloPart(3)\n");
  fprintf(fpout_idx,"# HaloID(1) HaloID(2)\n");
#endif // SUSSING2013
  fflush(fpout);
  fflush(fpout_idx);
  
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    
    if(halos[0][ihalo].ncroco > 0) {
      //      this is the old format where the haloid corresponds to the linenumber
      //      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].id[1]);
      
      // this is the case where we use the haloid as found in *_particles
      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].mtree[0].haloid[1]);
      fflush(fpout_idx);
      
      
#ifdef SUSSING2013
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].ncroco);
#else // SUSSING2013
      fprintf(fpout,"%"PRIu64"  %"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].npart,
              halos[0][ihalo].ncroco);
#endif // SUSSING2013
      fflush(fpout);
      
      for(icroco=0; icroco<halos[0][ihalo].ncroco; icroco++) {
#ifdef SUSSING2013
        fprintf(fpout,"%"PRIu64"\n",
                halos[0][ihalo].mtree[icroco].haloid[1]);
#else // SUSSING2013
        fprintf(fpout,"  %"PRIu64"  %"PRIu64"  %"PRIu64"\n",
                halos[0][ihalo].mtree[icroco].common,
                halos[0][ihalo].mtree[icroco].haloid[1],
                halos[0][ihalo].mtree[icroco].npart[1]);
#endif // SUSSING2013
        fflush(fpout);
      }
    }
#ifdef SUSSING2013
    else {
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].ncroco);
    }
#endif // SUSSING2013
  }
  
  /* close files */
  fclose(fpout);
  fclose(fpout_idx);
  
  elapsed = clock()-elapsed;
  fprintf(stderr," done in %4.2f sec.\n",(float)elapsed/CLOCKS_PER_SEC);
  return(1);
}

/*==================================================================================================
 * construct_filename
 *==================================================================================================*/
void construct_filename(char filename[MAXSTRING], int32_t ifile, char infile[MAXSTRING])
{
  char  *base, *path;
  char  prefix[MAXSTRING], suffix[MAXSTRING];
  char *d, *p;
  int   id, ip;
  
  path   = xdirname(filename);
  base   = xbasename(filename);
  
  
  d      = strstr(base,"0000");
  id     = d-base;
  strncpy(prefix,base,id);
  prefix[id] = '\0';
  
  p      = strstr(base,"_particles");
  ip     = p-base+11;
  strncpy(suffix,base+id+NDIGITS,ip);
  suffix[ip] = '\0';
  
  sprintf(infile,"%s/%s%04"PRIi32"%s",path,prefix,ifile,suffix);
  
#ifdef DEBUG3
  fprintf(stderr,"\npath  =%s",path);
  fprintf(stderr,"\nbase  =%s",base);
  fprintf(stderr,"\nprefix=%s",prefix);
  fprintf(stderr,"\nsuffix=%s",suffix);
  fprintf(stderr,"\ninfile=%s",infile);
#endif
  
  free(path);
  free(base);
}

/*==================================================================================================
 * count_particles_files
 *==================================================================================================*/
int32_t count_particles_files(char filename[MAXSTRING])
{
  char infile[MAXSTRING];
  int32_t nfiles;
  FILE *fp;
  
  nfiles = 0;
  construct_filename(filename, nfiles, infile);
  while ((fp=fopen(infile,"r")) != NULL) {
    fclose(fp);
    nfiles++;
    construct_filename(filename, nfiles, infile);
  }
  
  if(nfiles > pow(10,NDIGITS)-1) {
    fprintf(stderr,"There are nfiles = %"PRIi32" but you are only using NDIGITS = %d to construct the filenames\nABORTING\n",
            nfiles,NDIGITS);
    exit(0);
  }
  
  return (nfiles);
}

/*==================================================================================================
 * count_particles_files
 *==================================================================================================*/
uint64_t count_halos(char filename[MAXSTRING], int32_t nfiles)
{
  char infile[MAXSTRING], line[MAXSTRING];
  int32_t ifile;
  FILE *fp;
  uint64_t nh, ntmp;
  int i;
  
  nh = 0;
  for(ifile=0; ifile<nfiles; ifile++) {
    construct_filename(filename, ifile, infile);
    fp = fopen(infile,"r");
    fgets(line,MAXSTRING,fp);
    sscanf(line,"%"SCNi64,&ntmp);
    fclose(fp);
    nh += ntmp;
  }
  
  return(nh);
}

/*==================================================================================================
 * merit_sort
 *==================================================================================================*/
int merit_sort(const void *mtree1, const void *mtree2)
{
  double merit1, merit2;
  
  merit1 = ((MTREEptr)mtree1)->merit;
  merit2 = ((MTREEptr)mtree2)->merit;
  
  if(merit1 > merit2)
    return(-1);
  else if(merit1 < merit2)
    return(+1);
  else
    return(0);
  
  
}

