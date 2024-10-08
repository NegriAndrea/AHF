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
#define ONLY_USE_PTYPE                // restrict analysis to certain particle species (define them yourself in read_particles()!)
                                      // note: when restricting to stars (or gas) it will only work with -DUSE_PIDMAP!
#define MINCOMMON      10             // we only cross-correlate haloes if they at least share MINCOMMON particles
#define MTREE_BOTH_WAYS               // make sure that every halo has only one descendant (this is relevant for SAM models)
//#define USE_PIDMAP                    // the Pids are not used as array indices anymore, e.g. one can now also use star particles for the tree building
//#define EXCLUSIVE_PARTICLES           // each particle is only allowed to belong to one object (i.e. the lowest mass one): NOT fully tested yet!!!
//#define WITH_QSORT                    // uses qsort() instead of indexx() when ordering the progenitors according to merit function: NOT fully implemented yet!!!

// flags that only work together with SNAPSKIPPING
//#define SNAPSKIPPING                         // whenever a connection [0]->[1] is not considered credible, the halo will be copied and considered in the connection [1]->[2] (recursively)
#define SNAPSKIPPING_UNCREDIBLEMASSRATIO 2   // do not allow for mass jumps Mdesc = UNCREDIBLEMASSRATIO * Mprog
//#define SNAPSKIPPING_CONSIDERALLPROGENITORS  // consider all progenitors for UNCREDIBLEMASSRATIO test, in order of merit function
//#define DEBUG_SNAPSKIPPING


// support for AHF's MPI output
//#define READ_MPARTICLES               // support to read multiple _particle files, they must be of the latest AHF _particles file format!
#define NDIGITS                 4     // number of digits to be used for fileid


// support of ancient formats (without nhalos line, without haloids, whatever)
//#define THERE_IS_NO_NHALOS_LINE
//#define READ_HALOIDS_FROM_FILE
//#define USE_LINENUMBER_AS_HALOID      // overwrites(!) haloid as found in _particles

// only relevant for read_particles_bin()
#define SWAPBYTES 0

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
  uint64_t  npart;       // CAREFUL: we restrict the tree building to ONLY_USE_PTYPE and hence this npart is not necessarily what is given in _halos as the total number of particles!
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

#ifdef USE_PIDMAP
uint64_t *PidMap[2];
uint64_t  NPids[2];
void create_PidMap(int isimu);
int pid_cmp(const void *id1, const void *id2);
#endif

/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/
HALOptr     halos[2];
PARTptr     parts[2];
uint64_t    nHalos[2];
uint64_t    PidMax[2]={0,0}, PidMax_global=0;
uint64_t    PidMin=((uint64_t)1<<62);

/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
int      read_particles         (char filename[MAXSTRING], int isimu);
int      read_particles_bin     (char filename[MAXSTRING], int isimu);
int      particle_halo_mapping  (int  isimu);
int      cross_correlation      (char OutFile[MAXSTRING]);
int      write_mtree            (char OutFile[MAXSTRING]);
void     clean_connection       (uint64_t ihalo, int isimu0, int isimu1);
void     create_mtree           (uint64_t ihalo, int  isimu0,int isimu1);
void     create_mtree_qsort     (uint64_t ihalo, int isimu0, int isimu1);
uint64_t max_merit              (uint64_t ihalo, int isimu);
void     construct_filename     (char filename[MAXSTRING], int32_t, char infile[MAXSTRING]);
uint64_t count_halos            (char filename[MAXSTRING], int32_t);
int32_t  count_particles_files  (char filename[MAXSTRING]);
int      merit_sort             (const void *m1, const void *m2);
#ifdef SNAPSKIPPING
void     check_connections      ();
#endif

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
#ifdef USE_PIDMAP
  create_PidMap(0);
#endif
  particle_halo_mapping(0);
  fprintf(stderr,"\n");
  
  for(i=0; i<nFiles-1; i++) {
    
    /* be verbose */
    fprintf(stderr,"Correlating [%d/%d] '%s' to '%s'\n           -> writing to '%s'\n", i, nFiles-2, HaloFile[i], HaloFile[i+1],OutFile[i]);

    /* read the next file into memory */
    read_particles(HaloFile[i+1], 1);
#ifdef USE_PIDMAP
    create_PidMap(1);
#endif
    particle_halo_mapping(1);
    
#ifdef MTREE_BOTH_WAYS
    /* capture the case where PidMax[1] > PidMax[0] */
    if(PidMax[1] > PidMax[0]) {
      parts[0] = (PARTptr) realloc(parts[0], (PidMax[1]+1)*sizeof(PARTS)); // +1 because we are accessing an array like [PidMax_global] !
    }
#endif
    
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
    PidMax[0] = PidMax[1];
    
#ifdef MEMCOPY_1_TO_0 // Weiguang had problems running MergerTree, but this memcpy() apparently fixed the SegFault

    // manually copy everything from halos[1] and parts[1] over to halos[0] and parts[0]
    halos[0]  = (HALOptr) calloc(nHalos[0], sizeof(HALOS));
    for(ihalo=0; ihalo<nHalos[1]; ihalo++) {
        halos[0][ihalo].Pid = (uint64_t *) calloc(halos[1][ihalo].npart, sizeof(uint64_t));
        memcpy(halos[0][ihalo].Pid, halos[1][ihalo].Pid, halos[1][ihalo].npart * sizeof(halos[1][ihalo].Pid));
        halos[0][ihalo].haloid = halos[1][ihalo].haloid;
        halos[0][ihalo].npart = halos[1][ihalo].npart;
        halos[0][ihalo].ncroco = halos[1][ihalo].ncroco;
        halos[0][ihalo].mtree = (MTREEptr) calloc(halos[1][ihalo].ncroco, sizeof(MTREE));
        memcpy(halos[0][ihalo].mtree, halos[1][ihalo].mtree, halos[1][ihalo].ncroco * sizeof(halos[1][ihalo].mtree));
    }
    parts[0] = (PARTptr) calloc(PidMax[1]+1, sizeof(PARTS));
    for(ipart=0; ipart<PidMax[1]+1; ipart++) {
        parts[0][ipart].nhalos  = parts[1][ipart].nhalos;
        parts[0][ipart].Hid = (uint64_t *) calloc(parts[1][ipart].nhalos, sizeof(uint64_t));
        memcpy(parts[0][ipart].Hid, parts[1][ipart].Hid, parts[1][ipart].nhalos * sizeof(parts[1][ipart].Hid));
    }
    // and eventually wipe halos[1] and parts[1] from memory
    for(ihalo=0; ihalo<nHalos[1]; ihalo++) {
      if(halos[1][ihalo].Pid   != NULL) free(halos[1][ihalo].Pid);
      if(halos[1][ihalo].mtree != NULL) free(halos[1][ihalo].mtree);
    }
    for(ipart=0; ipart<PidMax[1]; ipart++) {
      if(parts[1][ipart].Hid != NULL) free(parts[1][ipart].Hid);
    }
    if(halos[1] != NULL) free(halos[1]);
    if(parts[1] != NULL) free(parts[1]);
#else
    // my strong belief still is that doing it this way should not cause a problem as pointers are *pointers* and nothing else...
    halos[0]  = halos[1];
    parts[0]  = parts[1];
#endif
    
#ifdef USE_PIDMAP
    free(PidMap[0]);
    PidMap[0] = PidMap[1];
    NPids[0]  = NPids[1];
#endif
    
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
  
#ifdef USE_PIDMAP
  free(PidMap[0]);
#endif
  
  printf("finished\n");
  return(1);
}


/*==================================================================================================
 * read_particles_bin:
 *
 * read AHF_particles_bin files (only existent for CLUES-WMAP3 simulation!?)
 *
 *==================================================================================================*/
int read_particles_bin(char filename[MAXSTRING], int isimu)
{
  FILE     *fpin, *fphids;
  char      line[MAXSTRING], infile[MAXSTRING], hidsname[MAXSTRING];
  int64_t   ihalo, jhalo;
  uint64_t  nPartInHalo, nPartInUse, ipart, jpart, Pid, Ptype, haloid, haloid_from_file;
  uint64_t  PidMin_local=((uint64_t)1<<62);
  uint64_t  PidMax_local=0;
  clock_t   elapsed;
  
  char c;
  int  i32, nh, ih;
  long i64;
  
  elapsed = clock();
  
  fprintf(stderr,"  o reading file %s ...",filename);
  
  fpin = fopen(filename,"rb");
  if(fpin == NULL)
  {
    fprintf(stderr,"could not open file %s\nexiting!\n",filename);
    exit(0);
  }
  
  /* reset all variables */
  nHalos[isimu] = 0;
  halos[isimu]  = NULL;
  
  // rubbish ?
  ReadInt(fpin,&i32,SWAPBYTES);
  ReadInt(fpin,&i32,SWAPBYTES);
  
  // nuber of haloes
  ReadInt(fpin,&i32,SWAPBYTES);
  nh = i32;
  
  // rubbish ?
  fread(&c,sizeof(char),1,fpin);
  ReadInt(fpin,&i32,SWAPBYTES);
  ReadInt(fpin,&i32,SWAPBYTES);
  ReadLong(fpin,&i64,SWAPBYTES);
  
  // loop over halos and their particle ids
  for(ihalo=0; ihalo<nh; ihalo++) {
    
    // use halo counter as haloid
    haloid = ihalo;
    
    // nparticles in halo
    ReadLong(fpin,&i64,SWAPBYTES);
    nPartInHalo = i64;
    
    /* found yet another halo */
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
      // read particle id
      ReadInt(fpin,&i32,SWAPBYTES);
      Pid = i32;
      
      // simulation only contains DM particles
      Ptype = 1;
      
      // here we can restrict the cross-correlation to a ceratain sub-set of all particles
#ifdef ONLY_USE_PTYPE
      if(Ptype == 1) // restrict to 1=dm
//      if(Ptype == 1 || Ptype == 4) // restrict to 1=dm and 4=stars
#endif
      {
        halos[isimu][ihalo].Pid             = (uint64_t *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(uint64_t));
        if(halos[isimu][ihalo].Pid == NULL) {
          fprintf(stderr,"read_particles: could not realloc() halos[%d][%ld].Pid for %"PRIu64"particles\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1));
          exit(-1);
        }
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
  } // for()
  
  fclose(fpin);
  
  elapsed = clock()-elapsed;
  
  fprintf(stderr," done in %4.2f sec. (nhalos = %"PRIu64", full ID range = %"PRIu64" -> %"PRIu64", local ID range = %"PRIu64" -> %"PRIu64")\n",
          (float)elapsed/CLOCKS_PER_SEC,nHalos[isimu],PidMin,PidMax[isimu],PidMin_local,PidMax_local);
  
  
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
  uint64_t  PidMin_local=((uint64_t)1<<62);
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
        if(Ptype == 1) // restrict to 1=dm
//        if(Ptype == 1 || Ptype == 4) // restrict to 1=dm and 4=stars
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
        if(Ptype == 1) // restrict to 1=dm
//        if(Ptype == 1 || Ptype == 4) // restrict to 1=dm and 4=stars
#endif
        {
          halos[isimu][ihalo].Pid = (uint64_t *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(uint64_t));
          if(halos[isimu][ihalo].Pid == NULL) {
            fprintf(stderr,"read_particles: could not realloc() halos[%d][%ld].Pid for %"PRIu64"particles\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1));
            exit(-1);
          }
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
  uint64_t ipart, jpart, *position;
  clock_t  elapsed;
  
  elapsed = clock();
  
  PidMax_global = MAX(PidMax[0],PidMax_global);
  PidMax_global = MAX(PidMax[1],PidMax_global);
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
#ifdef USE_PIDMAP
      //  find halos[isimu][ihalo].Pid[jpart] in PidMap[isimu] to get the 'ipart' to be used with parts[isimu][ipart]
      position = (uint64_t *)bsearch(&(halos[isimu][ihalo].Pid[jpart]), PidMap[isimu], NPids[isimu], sizeof(uint64_t), pid_cmp);
      if(position != NULL) {
        ipart = (uint64_t) (position - PidMap[isimu]);
        
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
      else {
        // this particle must come from some previous simulation and does not exist here...
      }
#else // USE_PIDMAP
      // simply use the Pid as the 'ipart'
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
      } // if(nhalos == 0)
      
#endif // USE_PIDMAP
      
      
    } // for(jpart)
  } // for(ihalo)
  
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
  
  
#ifdef SNAPSKIPPING
  elapsed = clock();
  fprintf(stderr,"  o removing missing/uncredible connections ...");
  check_connections();
  fprintf(stderr,"  o removing missing/uncredible connections -> done in %4.2f sec.\n", (double)elapsed/CLOCKS_PER_SEC);
#endif
  
  
  write_mtree(OutFile);
  
  return(1);
}

/*==================================================================================================
 * clean_connection
 *==================================================================================================*/
void clean_connection(uint64_t ihalo, int isimu0, int isimu1)
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
      mtree[ncroco_new-1].merit     = halos[isimu0][ihalo].mtree[icroco].merit;

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

/*==================================================================================================
 * create_mtree
 *==================================================================================================*/
void create_mtree(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t  jhalo, khalo, ipart, jpart, ncroco, icroco;
  int64_t   jcroco;
  
  uint64_t      *common;
  MTREE         *mtree;
  double        *merit;
  long unsigned *idx;
#ifdef USE_PIDMAP
  uint64_t *position;
#endif
  
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
#ifdef USE_PIDMAP
    //  find halos[isimu0][ihalo].Pid[jpart] in PidMap[isimu1] to get the actual ipart to be used with parts[isimu1][ipart]
    position = (uint64_t *) bsearch(&(halos[isimu0][ihalo].Pid[jpart]), PidMap[isimu1], NPids[isimu1], sizeof(uint64_t), pid_cmp);
    
    // when using the PidMap[] we need to make this additional check as this map for each isimu only contains the actual present particles
    if(position != NULL) {
      ipart = (uint64_t) (position - PidMap[isimu1]);
      
      /* ipart belongs to nhalos halos in isimu1 */
      for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
        khalo          = parts[isimu1][ipart].Hid[jhalo];
        common[khalo] += 1;
      }
    }
    else {
      // Pid[jpart] is not in [isimu1] and hence not shared with any halo from [isimu1]
    }
#else
    // simply use Pid as 'ipart'
    ipart = halos[isimu0][ihalo].Pid[jpart];
    
    // here we do not need that check as parts[][] contains PidMax_global entries, irrespective of isimu
    // (and .nhalos will be 0 for non-existent particles...)
    
    /* ipart belongs to nhalos halos in isimu1 */
    for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
      khalo          = parts[isimu1][ipart].Hid[jhalo];
      common[khalo] += 1;
    }
#endif
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
        
        // this is the applied merit function
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
      halos[isimu0][ihalo].mtree[jcroco].merit    = mtree[icroco].merit;
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

/*==================================================================================================
 * create_mtree_qsort
 *==================================================================================================*/
void create_mtree_qsort(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t  jhalo, khalo, ipart, jpart, ncroco, icroco;
  int64_t   jcroco;
  uint64_t *common;
#ifdef USE_PIDMAP
  uint64_t *position;
#endif
  
  // reset the actual mtree[] pointer
  if(halos[isimu0][ihalo].mtree != NULL)
    free(halos[isimu0][ihalo].mtree);
  halos[isimu0][ihalo].mtree = NULL;
  
  /* temporary array pointers */
  common = NULL;
  
  /* common[] records how many particles ihalo(isimu0) has in common with khalo(isimu1) */
  common = (uint64_t *) calloc(nHalos[isimu1], sizeof(uint64_t));
  
  for(jpart=0; jpart<halos[isimu0][ihalo].npart; jpart++) {
#ifdef USE_PIDMAP
    //  find halos[isimu0][ihalo].Pid[jpart] in PidMap[isimu1] to get the actual ipart to be used with parts[isimu1][ipart]
    position = (uint64_t *) bsearch(&(halos[isimu0][ihalo].Pid[jpart]), PidMap[isimu1], NPids[isimu1], sizeof(uint64_t), pid_cmp);
    
    // when using the PidMap[] we need to make this additional check as this map for each isimu only contains the actual present particles
    if(position != NULL) {
      ipart = (uint64_t) (position - PidMap[isimu1]);
      
      /* ipart belongs to nhalos halos in isimu1 */
      for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
        khalo          = parts[isimu1][ipart].Hid[jhalo];
        common[khalo] += 1;
      }
    }
    else {
      // Pid[jpart] is not in [isimu1] and hence not shared with any halo from [isimu1]
    }
#else
    // simply use Pid as 'ipart'
    ipart = halos[isimu0][ihalo].Pid[jpart];
    
    // here we do not need that check as parts[][] contains PidMax_global entries, irrespective of isimu
    // (and .nhalos will be 0 for non-existent particles...)
    
    /* ipart belongs to nhalos halos in isimu1 */
    for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
      khalo          = parts[isimu1][ipart].Hid[jhalo];
      common[khalo] += 1;
    }
#endif
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

/*==================================================================================================
 * write_mtree:
 *==================================================================================================*/
int write_mtree(char OutFile[MAXSTRING])
{
  uint64_t  ihalo, nHalos0_good;
  int64_t   icroco;
  FILE *fpout, *fpout_idx, *fpout_croco;
  char outname[MAXSTRING], outname_idx[MAXSTRING], outname_croco[MAXSTRING];
  clock_t   elapsed;
  
  elapsed = clock();
  
  sprintf(outname,"%s_mtree",OutFile);
  sprintf(outname_croco,"%s_croco",OutFile);
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
  
  fpout_croco = fopen(outname_croco,"w");
  if(fpout_croco == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",outname_croco);
    exit(0);
  }
  
  
  
  // count the number of halos with a progenitor
  nHalos0_good = (uint64_t)0;
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    if(halos[0][ihalo].ncroco > 0) {
      nHalos0_good++;
    }
  }
  fprintf(stderr,"  o writing cross-correlation for %"PRIu64" haloes (out of %"PRIu64" in total) ...",nHalos0_good,nHalos[0]);
  fprintf(fpout,"%"PRIu64"\n",nHalos0_good);
  fprintf(fpout_croco,"#   HaloID(1)   HaloPart(2)  NumProgenitors(3)\n");
  fprintf(fpout_croco,"#      SharedPart(1)    HaloID(2)   HaloPart(3) merit(4)\n");
  fprintf(fpout_idx,"# HaloID(1) HaloID(2)\n");
  fflush(fpout);
  fflush(fpout_idx);
  fflush(fpout_croco);
  
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    
    if(halos[0][ihalo].ncroco > 0) {
      //      this is the old format where the haloid corresponds to the linenumber
      //      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].id[1]);
      
      // this is the case where we use the haloid as found in *_particles
      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].mtree[0].haloid[1]);
      fflush(fpout_idx);
      
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].ncroco);
      
      fprintf(fpout_croco,"%"PRIu64"  %"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].npart,
              halos[0][ihalo].ncroco);
      fflush(fpout);
      fflush(fpout_croco);
      
      for(icroco=0; icroco<halos[0][ihalo].ncroco; icroco++) {
        fprintf(fpout,"%"PRIu64"\n",
                halos[0][ihalo].mtree[icroco].haloid[1]);
        
        fprintf(fpout_croco,"  %"PRIu64"  %"PRIu64"  %"PRIu64" %lf\n",
                halos[0][ihalo].mtree[icroco].common,
                halos[0][ihalo].mtree[icroco].haloid[1],
                halos[0][ihalo].mtree[icroco].npart[1],
                halos[0][ihalo].mtree[icroco].merit);
        fflush(fpout);
        fflush(fpout_croco);
      }
    } // ncroco > 0
  } // for(ihalo)
  
  /* close files */
  fclose(fpout);
  fclose(fpout_idx);
  fclose(fpout_croco);
  
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

#ifdef USE_PIDMAP
/*==================================================================================================
 * pid_cmp
 *==================================================================================================*/
int pid_cmp(const void *id1, const void *id2)
{
  uint64_t *pid1, *pid2;
  
  pid1 = (uint64_t *)id1;
  pid2 = (uint64_t *)id2;
  
  if(*pid1 < *pid2)
    return(-1);
  else if(*pid1 > *pid2)
    return(+1);
  else
    return(0);
}

/*==================================================================================================
 * create_PidMap
 *
 * PidMap[isimu] shall be used like this:
 *
 *     PidMap[isimu][ipart] = Pid;
 *     ipart                =  (uint64_t) bsearch(&(Pid), PidMap[isimu], NPids[isimu], sizeof(uint64_t), pid_cmp);
 *
 *       -> PidMap[isimu] maps the unique 'Pid' onto a unique 'ipart = [0, NPids[simu]-1]'
 *==================================================================================================*/
void create_PidMap(int isimu)
{
  uint64_t  ihalo, ipart;
  uint64_t *tmp_PidMap, tmp_NPids;
  clock_t   elapsed;
  elapsed = clock();
  
  fprintf(stderr,"  o creating PidMap for simu=%d ... ",isimu);
  
  tmp_NPids   = (uint64_t)0;
  tmp_PidMap  = NULL;
  
  for(ihalo=(uint64_t)0; ihalo<nHalos[isimu]; ihalo++) {
    for(ipart=(uint64_t)0; ipart<halos[isimu][ihalo].npart; ipart++) {
      tmp_PidMap            = (uint64_t *)realloc(tmp_PidMap, (tmp_NPids+1)*sizeof(uint64_t));
      tmp_PidMap[tmp_NPids] = halos[isimu][ihalo].Pid[ipart];
      tmp_NPids++;
    }
  }
  
  // make PidMap[] searchable via bsearch()
  qsort((void *)tmp_PidMap, tmp_NPids, sizeof(uint64_t), pid_cmp);
  
  // remove duplicates from sorted list
  PidMap[isimu] = NULL;
  NPids[isimu]  = (uint64_t)0;
  for(ipart=0; ipart<tmp_NPids-1; ipart++) {
    while(tmp_PidMap[ipart] == tmp_PidMap[ipart+1]) {
      ipart++;
    }
    PidMap[isimu]               = (uint64_t *)realloc(PidMap[isimu] , (NPids[isimu]+1)*sizeof(uint64_t));
    PidMap[isimu][NPids[isimu]] = tmp_PidMap[ipart];
    NPids[isimu]++;
  }
  PidMap[isimu]               = (uint64_t *)realloc(PidMap[isimu] , (NPids[isimu]+1)*sizeof(uint64_t));
  PidMap[isimu][NPids[isimu]] = tmp_PidMap[ipart];
  NPids[isimu]++;
  free(tmp_PidMap);
  
  PidMax[isimu] = NPids[isimu];
  
  elapsed = clock()-elapsed;
  
  fprintf(stderr," done in %4.2f sec. (removed %"PRIu64" duplicates, keeping %"PRIu64" unique Pids)\n", (float)elapsed/CLOCKS_PER_SEC,(tmp_NPids-NPids[isimu]),NPids[isimu]);
}
#endif


#ifdef SNAPSKIPPING
/*==================================================================================================
 * check for uncredible progenitors:
 *
 *   if a halo at [0] does not have a credible progenitor at [1] it will be added to the
 *   list of halos at [1] and considered when checking the connection [1]->[2]
 *==================================================================================================*/
int connectionrejection(HALOS halo, MTREE progmtree)
{
  double Mratio;
  
  //  Mratio = (double)halo.npart/(double)halo.mtree[0].npart[1];
  Mratio = (double)halo.npart/(double)progmtree.npart[1];
  if(Mratio<1) Mratio=1.0/Mratio;
  
#ifdef SNAPSKIPPING_UNCREDIBLEMASSRATIO
  if(Mratio > SNAPSKIPPING_UNCREDIBLEMASSRATIO)
    return(1);
  else
    return(0);
#else
  return(0);
#endif
}

void check_connections()
{
  MTREE    mtree_tmp;
  uint64_t ihalo, ipart, nHalos1, icroco, jcroco;
  int      rejected, missing;
  
  PidMax_global = MAX(PidMax[0],PidMax_global);
  PidMax_global = MAX(PidMax[1],PidMax_global);
  nHalos1       = nHalos[1];
  
  fprintf(stderr,"(PidMax_global=%"PRIu64")",PidMax_global);
  
  // loop over all halos at [0]
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    
    // we assume that this halo is *not* yet missing/rejected
    rejected = 0;
    missing  = 0;
    
    // if a progenitor exists check for credibility...
    if(halos[0][ihalo].ncroco > 0) {
      
#ifdef SNAPSKIPPING_CONSIDERALLPROGENITORS
      // loop over all progenitors, checking the connectionrejection condition for halos[0][ihalo]
      for(icroco=0; icroco<halos[0][ihalo].ncroco; icroco++) {
        if(connectionrejection(halos[0][ihalo], halos[0][ihalo].mtree[icroco]) == FALSE) {
          // the first one not giving 'connectionrejection(halos[0][ihalo]) == TRUE' should be the new top-of-the-list
          break;
        }
      }
      
      // if all progenitors failed the test we cannot assign a credible progenitor
      if(icroco == halos[0][ihalo].ncroco) {
#ifdef DEBUG_SNAPSKIPPING
        fprintf(stderr,"   rejected connection:\n");
        fprintf(stderr,"      [0]: haloid=%12"PRIu64" npart=%12"PRIu64" ncroco=%12"PRIu64" mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64" mtree.common=%12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].npart, halos[0][ihalo].ncroco, halos[0][ihalo].mtree[0].haloid[0], halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].npart[0], halos[0][ihalo].mtree[0].common);
        fprintf(stderr,"      [1]:                                                               mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", halos[0][ihalo].mtree[0].haloid[1], halos[0][ihalo].mtree[0].id[1], halos[0][ihalo].mtree[0].npart[1]);
#endif
        rejected = 1;
        
        // flag halos[0][ihalo] to not be written (only halos with ncroco>0 will be written to file...)
        halos[0][ihalo].ncroco = 0;
      }
      // but the first one that complies with our criterion (i.e. icroco) will now become the top-of-the-list
      else {
#ifdef DEBUG_SNAPSKIPPING
        if(icroco>0) {
          fprintf(stderr,"   using minor branch connection:\n");
          fprintf(stderr,"      [0]: haloid=%12"PRIu64" npart=%12"PRIu64" ncroco=%12"PRIu64" mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].npart, halos[0][ihalo].ncroco, halos[0][ihalo].mtree[0].haloid[0], halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].npart[0]);
          fprintf(stderr,"      [1]: old =>                  merit=%lf                                mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", halos[0][ihalo].mtree[0].merit, halos[0][ihalo].mtree[0].haloid[1], halos[0][ihalo].mtree[0].id[1], halos[0][ihalo].mtree[0].npart[1]);
          fprintf(stderr,"      [1]: new => icroco=%"PRIu64" merit=%lf                       mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", icroco,halos[0][ihalo].mtree[icroco].merit,halos[0][ihalo].mtree[icroco].haloid[1], halos[0][ihalo].mtree[icroco].id[1], halos[0][ihalo].mtree[icroco].npart[1]);
        }
#endif
        
        // re-shuffle the mtree[] list up to this point
        mtree_tmp.haloid[0] = halos[0][ihalo].mtree[icroco].haloid[0];
        mtree_tmp.haloid[1] = halos[0][ihalo].mtree[icroco].haloid[1];
        mtree_tmp.id[0]     = halos[0][ihalo].mtree[icroco].id[0];
        mtree_tmp.id[1]     = halos[0][ihalo].mtree[icroco].id[1];
        mtree_tmp.npart[0]  = halos[0][ihalo].mtree[icroco].npart[0];
        mtree_tmp.npart[1]  = halos[0][ihalo].mtree[icroco].npart[1];
        mtree_tmp.common    = halos[0][ihalo].mtree[icroco].common;
        mtree_tmp.merit     = halos[0][ihalo].mtree[icroco].merit;
        for(jcroco=1; jcroco<=icroco; jcroco++) {
          halos[0][ihalo].mtree[jcroco].haloid[0] = halos[0][ihalo].mtree[jcroco-1].haloid[0];
          halos[0][ihalo].mtree[jcroco].haloid[1] = halos[0][ihalo].mtree[jcroco-1].haloid[1];
          halos[0][ihalo].mtree[jcroco].id[0]     = halos[0][ihalo].mtree[jcroco-1].id[0];
          halos[0][ihalo].mtree[jcroco].id[1]     = halos[0][ihalo].mtree[jcroco-1].id[1];
          halos[0][ihalo].mtree[jcroco].npart[0]  = halos[0][ihalo].mtree[jcroco-1].npart[0];
          halos[0][ihalo].mtree[jcroco].npart[1]  = halos[0][ihalo].mtree[jcroco-1].npart[1];
          halos[0][ihalo].mtree[jcroco].common    = halos[0][ihalo].mtree[jcroco-1].common;
          halos[0][ihalo].mtree[jcroco].merit     = halos[0][ihalo].mtree[jcroco-1].merit;
        }
        halos[0][ihalo].mtree[0].haloid[0] = mtree_tmp.haloid[0];
        halos[0][ihalo].mtree[0].haloid[1] = mtree_tmp.haloid[1];
        halos[0][ihalo].mtree[0].id[0]     = mtree_tmp.id[0];
        halos[0][ihalo].mtree[0].id[1]     = mtree_tmp.id[1];
        halos[0][ihalo].mtree[0].npart[0]  = mtree_tmp.npart[0];
        halos[0][ihalo].mtree[0].npart[1]  = mtree_tmp.npart[1];
        halos[0][ihalo].mtree[0].common    = mtree_tmp.common;
        halos[0][ihalo].mtree[0].merit     = mtree_tmp.merit;
      }
      
#else // SNAPSKIPPING_CONSIDERALLPROGENITORS
      
      if(connectionrejection(halos[0][ihalo], halos[0][ihalo].mtree[0]) == TRUE) {
#ifdef DEBUG_SNAPSKIPPING
        fprintf(stderr,"   rejected connection:\n");
        fprintf(stderr,"      [0]: haloid=%12"PRIu64" npart=%12"PRIu64" ncroco=%12"PRIu64" mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64" mtree.common=%12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].npart, halos[0][ihalo].ncroco, halos[0][ihalo].mtree[0].haloid[0], halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].npart[0], halos[0][ihalo].mtree[0].common);
        fprintf(stderr,"      [1]:                                                               mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", halos[0][ihalo].mtree[0].haloid[1], halos[0][ihalo].mtree[0].id[1], halos[0][ihalo].mtree[0].npart[1]);
#endif
        rejected = 1;
        
        // flag halos[0][ihalo] to not be written (only halos with ncroco>0 will be written to file...)
        halos[0][ihalo].ncroco = 0;
      }
#endif // SNAPSKIPPING_CONSIDERALLPROGENITORS
    }
    // if a progenitor has not been found, treat as rejected and keep searching...
    else {
#ifdef DEBUG_SNAPSKIPPING
      fprintf(stderr,"   missing connection:\n");
      fprintf(stderr,"      [0]: haloid=%12"PRIu64" npart=%12"PRIu64" ncroco=%12"PRIu64" mtree=%12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].npart, halos[0][ihalo].ncroco, halos[0][ihalo].mtree);
#endif
      missing = 1;
    }
    
    
    // now deal with the missing/rejected conncetion...
    if(missing || rejected) {
      
      // only try to follow the missing/rejected halo, if it has enough particles
      if(halos[0][ihalo].npart > MINCOMMON) {
        
        // copy halos[0][ihalo] over to halos[1] (we copy as halos[0][ihalo] memory will be free'd!)
#ifdef DEBUG_SNAPSKIPPING
        fprintf(stderr,"           copying information:\n");
        fprintf(stderr,"             [0]:  haloid=%"PRIu64" npart=%"PRIu64" ncroco=%"PRIu64"\n",halos[0][ihalo].haloid,halos[0][ihalo].npart,halos[0][ihalo].ncroco);
        fprintf(stderr,"             [1]:  nHalos=%"PRIu64" -> ",nHalos[1]);
#endif
        halos[1] = (HALOptr) realloc(halos[1], (nHalos[1]+1)*sizeof(HALOS));
        halos[1][nHalos[1]].haloid = halos[0][ihalo].haloid;
        halos[1][nHalos[1]].npart  = halos[0][ihalo].npart;
        halos[1][nHalos[1]].ncroco = halos[0][ihalo].ncroco;
        halos[1][nHalos[1]].mtree  = NULL;  // we do not have any merger tree information for this one (yet)
        halos[1][nHalos[1]].Pid    = (uint64_t *) calloc(halos[0][ihalo].npart, sizeof(uint64_t));
        for(ipart=0; ipart<halos[1][nHalos[1]].npart; ipart++) {
          halos[1][nHalos[1]].Pid[ipart] = halos[0][ihalo].Pid[ipart];
          //if(halos[1][nHalos[1]].Pid[ipart] > PidMax[1]) PidMax[1] = halos[1][nHalos[1]].Pid[ipart]; // not needed as particle_halo_mapping takes MAX(PidMax[0],PidMax[1])
        }
        
        // increment number of halos at [1]
        nHalos[1]++;
#ifdef DEBUG_SNAPSKIPPING
        fprintf(stderr,"%"PRIu64"\n",nHalos[1]);
#endif
        
      }// if(MINCOMMON)
    } // if (missing || rejected)
  } // for(ihalo)
  
  // have we added additional halos to [1]?
  if(nHalos[1] > nHalos1) {
#ifdef USE_PIDMAP
    // update the PidMap for [1]
    free(PidMap[1]);
    create_PidMap(1);
#endif
    
    // we need to update the particle_halo_maping for [1]
    for(ipart=0; ipart<PidMax_global+1; ipart++) {   // free() old memory first
      if(parts[1][ipart].Hid != NULL){
        free(parts[1][ipart].Hid);
        parts[1][ipart].Hid = NULL;
      }
    }
    free(parts[1]);
    parts[1] = NULL;
    
    fprintf(stderr,"\n");
    particle_halo_mapping(1);                       // now call for re-creation of that mapping
  }
  
  return;
}
#endif // SNAPSKIPPING

