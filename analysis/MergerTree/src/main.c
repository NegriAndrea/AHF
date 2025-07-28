/*==================================================================================================
 *  MergerTree:   Merger Tree AHF_particles files
 *
 *
 *  input:    - how often to perform (N)
 *            - (N+1) _particles files
 *
 *  output:   - N _mtree file
 *
 *
 *==================================================================================================*/
/*-------------------------------------------------------------------------------------
 *             standard C includes as well as 
 *             all relevant AHF includes,
 *             define.h, and
 *             tdef.h, and
 *             all own prototypes
 *-------------------------------------------------------------------------------------*/
#include "include.h"

/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/
HALOptr     halos[2];
PARTptr     parts[2];
uint64_t    nHalos[2];
uint64_t    PidMax[2]={0,0}, PidMax_global=0;
uint64_t    PidMin=((uint64_t)1<<62);

uint64_t   *PidMap[2]; // only used with "#define USE_PIDMAP"
uint64_t    NPids[2];

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
  
  // dump the #define's to stderr
  dump_defines();
  
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
  
  // this cannot really be parallelized, primarily because of the SNAPSKIPPING
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
    
#else // MEMCOPY_1_TO_0
    
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
