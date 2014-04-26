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

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"
#include "../src/libutility/utility.h"

#include <mpi.h>
#include <omp.h>

#define MINCOMMON       10            // we only cross-correlate haloes if they at least share MINCOMMON particles
#define MAX_PARENT_HALO 5	      // maximum number of haloes to which a particle can belong to
#define ONLY_USE_PTYPE 1              // restrict analysis to particles of this type (1 = dark matter)
//#define EXCLUSIVE_PARTICLES           // each particle is only allowed to belong to one object (i.e. the lowest mass one)
//#define MTREE_BOTH_WAYS               // make sure that every halo has only one descendant
//#define SUSSING2013                   // write _mtree in format used for Sussing Merger Trees 2013
//#define USE_LINENUMBER_AS_HALOID      // do not use the haloid as found in _particles
//#define MERGER_RATIO   0.25           // writes output that readily allows to find mergers

//#define DEBUG_MPI

#define NUM_OMP_THREADS 2

#define IMPROVE_LB	// As the workload on each task depends on both the total number of particles and 
			// the total number of halos, this flag enables a L/B algorithm that balances the
			// data using the product N_part * LB_PART_FAC * sqrt(1 + N_halos * LB_HALO_FAC) where
			// the two factors are introduced below	
#ifdef IMPROVE_LB
#define LB_PART_FAC 1.0
#define LB_HALO_FAC 1.0
#endif

#define COUNTER 5000
#define ORDER				// Order halos in merger tree by merit function


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
} MTREE;

typedef struct HALOS *HALOptr;
typedef struct HALOS
{
  uint64_t  haloid;
  uint64_t  npart;
  uint64_t *Pid;	// Pid holds the global particle ID
  uint64_t *Pidord;	// Pid holds the global particle ID ordered using qsort
  uint64_t *Pindex;	// Index is the local particle position in struct parts[i][index]
  uint64_t  ncroco;	// Local number of cross correlation
  uint64_t  global_ncroco;	// Global number - sums all the haloes swapped so far
  MTREEptr  mtree;
}HALOS;

typedef struct PARTS *PARTptr;
typedef struct PARTS
{
  uint64_t  nhalos;
  uint64_t  Pid;				// This now stores the particle id corresponding to the indexed particle
  uint64_t  Hid[MAX_PARENT_HALO];		// Hid is the global halo id 
  uint64_t  Hindex[MAX_PARENT_HALO];		// Hid is the local halo position in halo[isimu][index]
}PARTS;

/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/

HALOptr halos[2];
HALOptr halos_tmp[2];
PARTptr parts[2];
PARTptr parts_buffer[2];	// on each task we need a buffer structure to save the incoming swapped particles


uint64_t    nHalos[2];
uint64_t    nHalosTmp[2];
uint64_t    nPart[2];		// total number of particles on task, use this instead of PidMax
uint64_t    nPartTmp[2];		

size_t      totHaloSize;	// when reading-in the halo keep track of its size including dynamically allocated mem
size_t      totHaloSizeTmp;	

int NReadTask;		// NReadTask is the number of tasks reading the files
int TotTask; 
int LocTask; 
int SendTask; 
int RecvTask;
int DeltaTask;
int filesPerTask; 
int extraFilesPerTask;
MPI_Status status;

#ifdef MERGER_RATIO
uint64_t nlines;  // number of lines in *_mtree file
#endif

/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
int      read_particles         (char filename[MAXSTRING], int isimu, int ifile);
int      particle_halo_mapping  (int  isimu);
int      halo_particle_mapping  (int  isimu);
int      cross_correlation      (int  iloop); 
int      create_mtree           (uint64_t ihalo, int isimu0, int isimu1, int iloop);
int      clean_connection       (uint64_t ihalo, int isimu0, int isimu1);
int      write_mtree            (int isimu0, char OutFile[MAXSTRING]);
uint64_t max_merit              (uint64_t ihalo, int isimu);

int	 assign_input_files_to_tasks	(char *partList, char *tempDir, char ***locPartFile, int nFiles);
int 	 cmpfunc			(const void * a, const void * b); 
int	 add_halos			(int ifile, int isimu);
int      alloc_halos			(int isimu);
int      free_halos			(int isimu);
int 	 load_balance			(int isimu); 
void 	 intersection 			(int isimu0, int isimu1, uint64_t ihalo, uint64_t khalo, uint64_t *common);
uint64_t constructHaloId(uint64_t ihalo);
#ifdef MERGER_RATIO
int  	 assign_progenitors	        (char OutFile[MAXSTRING]);
void	 read_mtree             	(char *infile);
#endif

#ifdef DEBUG_MPI
void 	check_parts			(int isimu);
void 	check_halos			(int isimu);
#endif


/*==================================================================================================
 * main:
 *
 *       simply a wrapper for successive calls to create_mtree()
 *
 *==================================================================================================*/
int main(int argv, char **argc)
{
  int    ifile, jfile, jchunk, count, nFiles;   // nFiles is the number of snapshots (composed of several chunks)
  char   *tempDir;		// Temporary files are stored here - to be cleaned at the end of the run
  char   *partList;		// This file contains all the files to be submitted to task 0 
  char   *prefixOut;

  parts_buffer[1]=NULL;

  /* local variables */
  uint64_t buffer_npart;	// Buffer to recieve the new total number of particles
  uint64_t buffer_nhalo;	// Buffer to recieve the new total number of halos
  time_t    elapsed = (time_t)0, total = (time_t) 0;

  char   ***locPartFile=NULL;		// Each task stores _particle urls here
  char   **locOutFile=NULL;		// Each task will dump to this file

  count = 1 ;
  
  /* set some global variables common to all MPI tasks */
  NReadTask = atoi(argc[count++]); 
  nFiles = atoi(argc[count++]);	  
  partList = argc[count++];	   
  prefixOut = argc[count++];
  tempDir = argc[count++];

#ifdef DEBUG_MPI
  /* check that the global variables have been read correctly */
  for(ifile=0; ifile<5; ifile++)
    fprintf(stderr, "argc[%d]=%s\n", ifile+1, argc[ifile+1]);
#endif

  /* general variables are set, now init MPI parallel part */
  MPI_Init(&argv, &argc);
  MPI_Comm_rank(MPI_COMM_WORLD, &LocTask); 
  MPI_Comm_size(MPI_COMM_WORLD, &TotTask); 

  /* RecvTask is recieving from LocTask and SendTask is sending to LocTask */
  SendTask = (LocTask+TotTask-1) % (TotTask); 
  RecvTask = (LocTask+1) % (TotTask);
  DeltaTask = TotTask - NReadTask; 

  /* This is useful if you have more file-chunks than tasks */
  if(NReadTask > TotTask)
  {
     filesPerTask = (int) NReadTask / TotTask;  
     extraFilesPerTask = (int) NReadTask % TotTask;  

     if(LocTask < extraFilesPerTask)
       filesPerTask++; 
  }
  else
  {
     filesPerTask = 1;
  }

  /* allocate memory for locOutFiles, each of size MAXSTRING */
  locOutFile  = (char **) calloc(nFiles, sizeof(char *));

  for(ifile=0; ifile<nFiles; ifile++)
    locOutFile[ifile]  = (char *) calloc(MAXSTRING, sizeof(char));

#ifdef DEBUG_MPI
  fprintf(stderr, "Task=%d is sending to and Task=%d recieving from LocTask=%d\n", SendTask, RecvTask, LocTask);
#endif

  if(LocTask == 0)
     fprintf(stderr, "MergerTree MPI mode, reading %d files per task on %d tasks on a total of %d cpus.\n", 
		nFiles, NReadTask, TotTask);
  
  /* read the first file into memory */
  if(LocTask == 0)
    fprintf(stderr,"\nStartup:\n");

  /* TODO read_particles routine should allow for a single task to read more files */
  if(LocTask < NReadTask)
  {
     /* alloc memory for local urls storage */
     locPartFile = (char ***) calloc(filesPerTask, sizeof(char **));

     for(ifile=0; ifile<filesPerTask; ifile++)
       locPartFile[ifile] = (char **) calloc(nFiles, sizeof(char *));

     /* now assign a list of urls to every task and alloc local memory */
     assign_input_files_to_tasks(partList, tempDir, locPartFile, nFiles);
  
     /* read particles now reads into the tmp file which is reallocated in the main halos struct */
     for(ifile=0; ifile<filesPerTask; ifile++)
     {
        read_particles(locPartFile[ifile][0], 0, ifile);
        add_halos(ifile, 0);
     }
  }

    /*  If TotTask > NReadTask then redistribute particles across different tasks */
    if(DeltaTask > 0) 
      load_balance(0);  

    /* map the particles from the halos after load balancing */
    particle_halo_mapping(0);

  if(LocTask == 0)
    fprintf(stderr,"\n");

  /* now loop over all files */
  for(ifile=0; ifile<nFiles-1; ifile++)
  {
    sprintf(locOutFile[ifile], "%s_%03d-%03d.%04d", prefixOut, ifile, ifile+1, LocTask); 

    /* every task reads the next file into memory - the next file should be a chunk of _particle files
       at a different redshift */
  if(LocTask < NReadTask)
    for(jfile=0; jfile<filesPerTask; jfile++) /* this works if each task reads in more than one file */
    {
       /* be verbose */
       if(LocTask == 0 && jfile == 0)
         fprintf(stderr,"Correlating '%s' to '%s'\n           -> writing to '%s'\n",
           locPartFile[jfile][ifile],locPartFile[jfile][1+ifile],locOutFile[ifile]);
#ifdef DEBUG_MPI
  if(LocTask == 0)
       fprintf(stderr, "Loop=%d), jfile=%d, filename=%s\n", ifile, jfile, locPartFile[jfile][1+ifile]);
#endif
       read_particles(locPartFile[jfile][1+ifile], 1, jfile);
       add_halos(jfile, 1);
    }

    /* if NTask > N files to be read then scatter the data through the tasks */
    if(DeltaTask > 0)   
       load_balance(1);  

    MPI_Barrier(MPI_COMM_WORLD);  

    particle_halo_mapping(1);

    /* now swap the files across all the tasks and look for correlations */
    for(jchunk=0; jchunk<TotTask; jchunk++)
    {
      elapsed -= time(NULL);

      /* cross correlate locHaloFile[i] to locHaloFile[i+1] */
      cross_correlation(jchunk);

      elapsed += time(NULL);
      total += elapsed;

      if(LocTask == 0)
        fprintf(stderr, "Cross correlation completed step %d/%d of file %d in %ld sec, total elapsed time is %ld s.\n", 
	        jchunk+1, TotTask, jfile, elapsed, total);

      elapsed = (time_t) 0;

#ifdef DEBUG_MPI 
      fprintf(stderr, "Total elapsed time on task=%d is %ld sec.\n", LocTask, total);
#endif

      /* wait for all the processes to end the correlation */
      MPI_Barrier(MPI_COMM_WORLD);  

#ifdef DEBUG_MPI
	fprintf(stderr, "Task=%d is sending %"PRIu64" particles and %"PRIu64" halos to task=%d\n",
	 		LocTask, nPart[1], nHalos[1], RecvTask);
#endif

	/* swap the number of particles and halos across tasks */
	MPI_Sendrecv(&nPart[1], sizeof(uint64_t), MPI_BYTE, RecvTask, 0,
                 &buffer_npart, sizeof(uint64_t), MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

	MPI_Sendrecv(&nHalos[1], sizeof(uint64_t), MPI_BYTE, RecvTask, 0,
                  &buffer_nhalo, sizeof(uint64_t), MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

#ifdef DEBUG_MPI
	fprintf(stderr, "Task=%d has recieved %"PRIu64" particles and %"PRIu64" halos from task=%d\n",
	 		LocTask, buffer_npart, buffer_nhalo, SendTask);
#endif
 
	/* allocate on LocTask a buffer to recieve parts data from SendTask */
       parts_buffer[1] = (PARTS *) calloc(buffer_npart , sizeof(PARTS));

	  /* now swap the particles across all tasks - halos will be reconstructed locally later on */	
	  MPI_Sendrecv(parts[1], nPart[1] * sizeof(PARTS), MPI_BYTE, RecvTask, 0, 
	    parts_buffer[1], buffer_npart * sizeof(PARTS), MPI_BYTE, SendTask, 0, MPI_COMM_WORLD, &status);

	/* Free the old simu1 halo structures before reconstructing the new ones from
           the swapped particles. Do this before resetting nHaloes[1] from the buffer value */
        free_halos(1);
	
	/* update the local number of particles after Sendrecv */
	nPart[1] = buffer_npart;
	nHalos[1] = buffer_nhalo;

	/* free the old pointer and set it equal to the new one */
	if(parts[1] != NULL)
	  free(parts[1]);

	parts[1] = parts_buffer[1];

	/* now reallocate the halo structs and map back the particles into their respective halos */
        alloc_halos(1);
        halo_particle_mapping(1);

        MPI_Barrier(MPI_COMM_WORLD);  

#ifdef MERGER_RATIO
	// FIXME not tested for the MPI version
      assign_progenitors(locOutFile[ifile]); //dumps information about progenitors
#endif
    } /* End swapping data */

	write_mtree(0, locOutFile[ifile]);  
  
    /* be verbose */
    if(LocTask == 0)
      fprintf(stderr,"  o making file 1 the new file 0 ...");

     /* moves the mtree properties from file 0 to file 1 */ 
     // FIXME not implemented for the MPI version
     //remove_network_connections();

    /* remove halo[0] structs from memory */
    free_halos(0);

     if(parts[0] != NULL) {
        free(parts[0]);
        parts[0] = NULL;
    }

    /* make HaloFile[i+1] the new HaloFile[i] */
    nHalos[0] = nHalos[1];
    nPart[0]  = nPart[1];
    halos[0]  = halos[1];
    parts[0]  = parts[1];

    /* be verbose */
    if(LocTask == 0)
    fprintf(stderr," done\n");
  } /* for(nFiles) */
  


  /*==================================================================*
   *                             CLEANUP                              *
   *==================================================================*/
    if(LocTask == 0)
      fprintf(stderr,"\nCleaning up ...\n");

  free_halos(0);
  if(parts[0] != NULL) free(parts[0]);

  /* free input filename storage */

  if(LocTask < NReadTask)
  {
    for(jfile=0; jfile<filesPerTask; jfile++)
      for(ifile=0; ifile<nFiles; ifile++)
        if(locPartFile[jfile][ifile]) free(locPartFile[jfile][ifile]);
         if(locPartFile) free(locPartFile);
  }

  /* free output filename storage*/
  for(ifile=0; ifile<nFiles; ifile++)
    if(locOutFile[ifile])  free(locOutFile[ifile]);

  if(locOutFile)  free(locOutFile);

  if(LocTask == 0)
    printf("finished\n");

  MPI_Finalize();

  return(1);
}


/*==================================================================================================
 * load_balance:
 * 
 * when using more tasks than input files (parts of a given catalogue) we need to redistribute
 * halos among them. 
 * At this point only the halos structures have been read-in; the inverse
 * mapping to the particles will be done after the halos have been scattered.
 * Since tasks that read in the halos use dynamically allocated memory (the for the Pid and Pindex 
 * relative to the particles belonging to the halo) we use the MPI_Pack/MPI_Unpack functions.
 *
 *==================================================================================================*/
int load_balance(int isimu)
{
  int irecv=0, itask=0, nTaskPerSend=0, extraTask=0, buffer_pos_recv=0, ihalo=0, jhalo=0, ipart=0, jpart=0; 
  int *nRecvTasks=NULL, *sendTasks=NULL, **recvTasks=NULL, *nHalosBuffer=NULL;
  uint64_t locNPart=0, **haloBuffer=NULL;

  size_t sizeUInt64=0, sizeSendBuffer=0, sizeTempBuffer=0, sizeLocHalo=0; 
  size_t *sizeBuffer=NULL, *haloSendBuffer=NULL; 

  HALOptr temp_halo[2];
  int *buffer_position=NULL;

#ifdef IMPROVE_LB
  uint64_t *partPerTask=NULL, minNPartTmp=0, minNPart=0;
  double nHaloFac=0, nPartFac=0, nPartFacTmp=0;
#endif

  sizeUInt64 = sizeof(uint64_t);

  /* recv tasks keeps track of all the tasks to which a single (reading) tasks 
   * has to broadcast data, send task tells the recieving task who is sending.
   * Each reading tasks can send to many tasks, but each recieving task
   * only gets data from a single one. */
  if(LocTask < NReadTask)
  { 
    recvTasks = (int **) calloc(NReadTask, sizeof(int *));
    nRecvTasks = (int *) calloc(NReadTask, sizeof(int));
  }
  else
  {
    sendTasks = (int *) calloc(DeltaTask, sizeof(int));
  }

  nTaskPerSend = (int) DeltaTask / NReadTask;
  extraTask = DeltaTask % NReadTask;

  if(LocTask == 0)
    fprintf(stderr, "Load balancing from %d ReadTasks on %d TotalTasks, each SendTask has %d RecvTasks\n.",
             NReadTask, TotTask, nTaskPerSend);

  /* we figure out who is sending to how many tasks */ 
  if(LocTask < NReadTask)
  {
     if(LocTask < extraTask)
     {
	nRecvTasks[LocTask] = nTaskPerSend + 1;
        recvTasks[LocTask] = (int *) calloc(nTaskPerSend+1, sizeof(int));
     }
     else if (LocTask >= extraTask && nTaskPerSend > 0)
     {
	nRecvTasks[LocTask] = nTaskPerSend;
        recvTasks[LocTask] = (int *) calloc(nTaskPerSend, sizeof(int));
     }
     else 
     {
	nRecvTasks[LocTask] = 0;
        recvTasks[LocTask] = NULL;
     }

     for(itask=0; itask<nRecvTasks[LocTask]; itask++)
     {
       recvTasks[LocTask][itask] = LocTask + (itask + 1) * NReadTask;
     }
   } 
   else /* the recieving tasks have to figure out who is sending */ 
   {
      sendTasks[LocTask-NReadTask] = LocTask % NReadTask;
   }

  /* Figure out how much every task has to send */

  if(LocTask < NReadTask && nRecvTasks[LocTask]>0)
  {
    haloBuffer = (uint64_t **) calloc(nRecvTasks[LocTask], sizeof(uint64_t *)); 
    sizeBuffer = (size_t *) calloc(nRecvTasks[LocTask], sizeof(size_t));
    buffer_position = (int *) calloc(nRecvTasks[LocTask], sizeof(int));
    nHalosBuffer = (int *) calloc(nRecvTasks[LocTask]+1, sizeof(int));

    haloBuffer[irecv] = (uint64_t *) calloc(1, sizeUInt64);
 
    temp_halo[isimu] = (HALOptr) calloc(1, sizeof(HALOS)); 

#ifdef IMPROVE_LB
    partPerTask = (uint64_t *) calloc(nRecvTasks[LocTask]+1, sizeUInt64); 
#endif

    /* Loop over haloes and assign them to the different tasks, packing them into separate buffers */
    for(ihalo=1; ihalo<=nHalos[isimu]; ihalo++)
    {
#ifdef IMPROVE_LB
      jhalo = ihalo-1;
#else
      jhalo = nHalos[isimu] - ihalo;
#endif

      irecv = (jhalo) % (nRecvTasks[LocTask]+1);
      locNPart = halos[isimu][jhalo].npart;
      sizeLocHalo = (size_t) locNPart * sizeUInt64;
     
#ifdef IMPROVE_LB
      minNPart = partPerTask[0];
      nHaloFac = sqrt(nHalosBuffer[0]);
      nPartFac = nHaloFac * (double) minNPart;
      
      for(itask=1; itask<nRecvTasks[LocTask]+1; itask++)
      {
         minNPartTmp = partPerTask[itask];
 	 nPartFacTmp = (double) minNPartTmp * LB_PART_FAC * sqrt(nHalosBuffer[itask] * LB_HALO_FAC);

         if(nPartFacTmp < nPartFac) 
         {
	    irecv = itask;
	    nPartFac = nPartFacTmp;
         }
      }

         partPerTask[irecv] += locNPart;

#endif


      /* Each halo holds two uint64_t and two arrays Pid and Pindex */
      sizeTempBuffer = 2 * sizeUInt64 * (halos[isimu][jhalo].npart + 1);

      if(irecv == nRecvTasks[LocTask]) /* keep the halo locally */
      {
         temp_halo[isimu] = (HALOptr) realloc(temp_halo[isimu], (nHalosBuffer[irecv] + 1) * sizeof(HALOS));
         temp_halo[isimu][nHalosBuffer[irecv]].npart = locNPart;
         temp_halo[isimu][nHalosBuffer[irecv]].haloid = halos[isimu][jhalo].haloid;

         temp_halo[isimu][nHalosBuffer[irecv]].Pid = (uint64_t *) malloc(sizeLocHalo);

         memcpy(temp_halo[isimu][nHalosBuffer[irecv]].Pid, halos[isimu][jhalo].Pid, sizeLocHalo);   
      }
      else
      {
         /* now allocate and initialize the buffers to which will be sent to the next task */
         sizeBuffer[irecv] += sizeTempBuffer;
         haloBuffer[irecv]  = (uint64_t *) realloc(haloBuffer[irecv], sizeBuffer[irecv]);

         /* Now pack the halo and particles information into a single buffer */
         MPI_Pack(&halos[isimu][jhalo].haloid, sizeUInt64, MPI_BYTE, haloBuffer[irecv],
          sizeBuffer[irecv], &buffer_position[irecv], MPI_COMM_WORLD);

         MPI_Pack(&halos[isimu][jhalo].npart, sizeUInt64, MPI_BYTE, haloBuffer[irecv],
          sizeBuffer[irecv], &buffer_position[irecv], MPI_COMM_WORLD);

         MPI_Pack(&halos[isimu][jhalo].Pid[0], sizeLocHalo, MPI_BYTE, haloBuffer[irecv], 
          sizeBuffer[irecv], &buffer_position[irecv], MPI_COMM_WORLD);

         nPart[isimu] -= locNPart;
      }

      nHalosBuffer[irecv]++;

      /* reallocate the local halos struct at each step to free up some memory */
#ifndef IMPROVE_LB
      free(halos[isimu][jhalo].Pid);
      free(halos[isimu][jhalo].Pindex);
      halos[isimu] = (HALOptr) realloc(halos[isimu], (jhalo+1) * sizeof(HALOS));
#endif
   } /* for ihalo */

   nHalos[isimu] = nHalosBuffer[nRecvTasks[LocTask]];
   nPart[isimu] = 0;

   free(halos[isimu]);
   halos[isimu] = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));

   ipart=0;
   for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
   { 
      halos[isimu][ihalo].npart = temp_halo[isimu][ihalo].npart;
      halos[isimu][ihalo].haloid = temp_halo[isimu][ihalo].haloid;

      locNPart = halos[isimu][ihalo].npart;
      sizeLocHalo = (size_t) locNPart * sizeUInt64;
      nPart[isimu] += locNPart;

      halos[isimu][ihalo].Pid = (uint64_t *) malloc (sizeLocHalo);
      halos[isimu][ihalo].Pindex = (uint64_t *) malloc (sizeLocHalo);

      memcpy(halos[isimu][ihalo].Pid, temp_halo[isimu][ihalo].Pid, sizeLocHalo);   

      for(jpart=0; jpart<locNPart; jpart++)
      {
        halos[isimu][ihalo].Pindex[jpart] = ipart;
        ipart++;
      }
      
      free(temp_halo[isimu][ihalo].Pid);
   }

#ifdef DEBUG_MPI
   fprintf(stderr, "Task=%d is left with %"PRIu64" and %"PRIu64" particles.\n",
	LocTask, nHalos[isimu], nPart[isimu]);
#endif

   free(temp_halo[isimu]);

#ifdef DEBUG_LOG 
   if(LocTask < NReadTask)
   {
      for(khalo=0; khalo<nHalos[isimu]; khalo++)
	dump_log_halo(isimu, khalo);
   }
#endif

     if(nHalos[isimu] == 0)
     {
       fprintf(stderr, "WARNING! No haloes left on Task=%d. Are you using too many MPI tasks?\n", LocTask);
     }
     else
     {
#ifdef DEBUG_MPI
	for(itask=0; itask<nRecvTasks[LocTask]; itask++)
 	{
	  totHaloSize -= sizeBuffer[itask];
	    fprintf(stderr, "Task=%d is sending %zd bytes and %d halos to task %d\n",
	            LocTask, sizeBuffer[itask], nHalosBuffer[itask], recvTasks[LocTask][itask]);
	}

	    fprintf(stderr, "Task=%d has done packing data, local data size is now=%zd, with %"PRIu64" halos.\n", 
		    LocTask, totHaloSize, nHalos[isimu]);
#endif
     }
  } /* if(LocTask < NReadTask)*/

  /* Send all the buffers */
  if(LocTask < NReadTask && nRecvTasks[LocTask]>0)
  {
    for(irecv=0; irecv<nRecvTasks[LocTask]; irecv++)
    {
      MPI_Send(&sizeBuffer[irecv], sizeof(size_t), MPI_BYTE, recvTasks[LocTask][irecv], 0, MPI_COMM_WORLD);
      MPI_Send(&nHalosBuffer[irecv], 1, MPI_INT, recvTasks[LocTask][irecv], 0, MPI_COMM_WORLD);
      MPI_Send(&haloBuffer[irecv][0], sizeBuffer[irecv], MPI_BYTE, recvTasks[LocTask][irecv], 0, MPI_COMM_WORLD);

#ifdef DEBUG_MPI
    fprintf(stderr, "Sending %zd bytes of messages from task=%d to task=%d.\n", 
	sizeBuffer[irecv], LocTask, recvTasks[LocTask][irecv]);
#endif

    }
  }
  else if(LocTask >= NReadTask) /* the task is a recieving one */
  {
      MPI_Recv(&sizeSendBuffer, sizeof(size_t), MPI_BYTE, sendTasks[LocTask-NReadTask], 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&nHalos[isimu], 1, MPI_INT, sendTasks[LocTask-NReadTask], 0, MPI_COMM_WORLD, &status);
   
      /* Use this buffer to store the incoming packed halos */
      haloSendBuffer = (uint64_t *) malloc(sizeSendBuffer);

#ifdef DEBUG_MPI
    fprintf(stderr, "\nRecieved %zd bytes of messages and %"PRIu64" halos on task=%d from task=%d.\n", 
	sizeSendBuffer, nHalos[isimu], LocTask, sendTasks[LocTask-NReadTask]);
#endif

      MPI_Recv(haloSendBuffer, sizeSendBuffer, MPI_BYTE, sendTasks[LocTask-NReadTask], 0, MPI_COMM_WORLD, &status);
  }

  /* Unpack the recieved buffers on each task */
  if(LocTask >= NReadTask)
  {
    halos[isimu]   = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));
    buffer_pos_recv = 0;
    nPart[isimu] = 0;
    ipart = 0;

    for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
    {
       /* Now unpack halos particles information from the buffers into the halos structures */
       MPI_Unpack(&haloSendBuffer[0], sizeSendBuffer, &buffer_pos_recv, &halos[isimu][ihalo].haloid, 
	sizeUInt64, MPI_BYTE, MPI_COMM_WORLD);
	
       MPI_Unpack(&haloSendBuffer[0], sizeSendBuffer, &buffer_pos_recv, &halos[isimu][ihalo].npart, 
 	sizeUInt64, MPI_BYTE, MPI_COMM_WORLD);

       locNPart = halos[isimu][ihalo].npart;
       sizeLocHalo = (size_t) locNPart * sizeUInt64;

       halos[isimu][ihalo].Pid = (uint64_t *) malloc(sizeLocHalo);
       halos[isimu][ihalo].Pindex = (uint64_t *) malloc(sizeLocHalo);

       MPI_Unpack(&haloSendBuffer[0], sizeSendBuffer, &buffer_pos_recv, &halos[isimu][ihalo].Pid[0], 
        sizeLocHalo, MPI_BYTE, MPI_COMM_WORLD);

      for(jpart=0; jpart<locNPart; jpart++)
      {
        halos[isimu][ihalo].Pindex[jpart] = ipart;
        ipart++;
      }      
       nPart[isimu] += locNPart;
    }

   } /* if(LocTask >= NReadTask)*/

#ifdef DEBUG_MPI
	size_t totHaloSizeTest=0;

	for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
	  totHaloSizeTest += 2 * (1 + halos[isimu][ihalo].npart) * sizeUInt64;

	fprintf(stderr, "\nDone load balancing.\nTask=%d has %"PRIu64" halos and %"PRIu64" particles, size=%zd.\n", 
		LocTask, nHalos[isimu], nPart[isimu], totHaloSizeTest);
#endif

  /* Now free local buffers */
  if(LocTask < NReadTask && nRecvTasks[LocTask]>0)
  {
     for(ihalo=0; ihalo<nRecvTasks[LocTask]; ihalo++)
       free(haloBuffer[ihalo]);
     
     free(haloBuffer);
     free(sizeBuffer);
     free(haloSendBuffer);
     free(nHalosBuffer);
  }

  return(1);
}



/*==================================================================================================
 * assign_input_files_to_tasks:
 *
 * assigns to the local task the list of files which will have to be read.
 * For the moment, each task reads in one file per redshift.
 *
 *      partList  = a list with all the particle files which will be read by task 0
 *      tempDir	  = a temporary folder where all the lists of files corresponding to each task are stored
 *      NReadTask = number of processors actually reading (must be smaller or equal to the number of tasks) 
 *
 *==================================================================================================*/
int assign_input_files_to_tasks(char *partList, char *tempDir, char ***locPartFile, int nFiles)
{
  int ifile, jfile;
  FILE *locPartListFile=NULL;
  char command[MAXSTRING];
  char locPartList[MAXSTRING];
  char dummy[MAXSTRING-1];
  char task[6];
  
  if(LocTask == 0)
    fprintf(stderr, "\n  o assigning files to each task from %s\n", partList);

//  locPartList[MAXSTRING] = (char *) calloc(filesPerTask, sizeof(char));

  if(LocTask < NReadTask)
  {
    
    for(ifile=0; ifile<filesPerTask; ifile++)
    {
      sprintf(task, "%04d", LocTask + ifile * TotTask); 
      sprintf(locPartList, "%s.%s.tmp", tempDir, task);
      sprintf(command, "sed s/0000/%s/ <%s >%s", task, partList, locPartList);
      system(command);

      locPartListFile = fopen(locPartList, "r");

      for(jfile=0; jfile<nFiles; jfile++)
      {
        fgets(dummy, MAXSTRING-1, locPartListFile);
        locPartFile[ifile][jfile] = (char*) calloc(strlen(dummy)+5, sizeof(char));
        strcpy(locPartFile[ifile][jfile], dummy);
        locPartFile[ifile][jfile][strlen(dummy)-1]='\0';
#ifdef DEBUG_MPI
        fprintf(stderr, "Task=%d will read from file %s\n", LocTask, locPartFile[ifile][jfile]);
#endif
      }

      /* close the temp file generated with the list of the files to be read */
      fclose(locPartListFile);
      locPartListFile = NULL;
    } /* ifile, loop on the file chunks that each task is loading in */
  } /* else task remains idle */

  if(LocTask == 0)
    fprintf(stderr, " done.\n");

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
int read_particles(char filename[MAXSTRING], int isimu, int ifile)
{
  FILE     *fpin;
  char      line[MAXSTRING];
  int64_t   ihalo;
  uint64_t  nPartInHalo, nPartInUse, ipart, Pid, Pindex, Ptype, haloid;
  time_t    elapsed = (time_t)0;

  elapsed -= time(NULL);

#ifdef DEBUG_MPI
    fprintf(stderr,"  o reading file %s on task %d ...\n",filename, LocTask);
#else
  if(LocTask == 0)
    fprintf(stderr,"  o reading file %s on task %d ...",filename, LocTask);
#endif

  fpin = fopen(filename,"r");
  if(fpin == NULL)
   {
    fprintf(stderr,"could not open file %s on task %d \nexiting!\n",filename, LocTask);
    exit(0);
   }
  
  /* reset all variables */
  nHalosTmp[isimu] = 0;
  nPartTmp[isimu]  = 0;
  ihalo         = -1;
  halos_tmp[isimu]  = NULL;

  /* this array keeps track of the size of each halo's dynamically allocated particle informations, Pid and Pindex */
  totHaloSizeTmp = 0;

  /* get the first line from file */
  fgets(line,MAXSTRING,fpin);
 
  /* for AHF_particles files the first line is numGoodHalos which we can happily ignore */
  if(strncmp(line,"#",1) != 0 && sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)

    fgets(line,MAXSTRING,fpin);  
  
      /* local numbering of the particle */
  if(ifile == 0)
      Pindex = 0;
  else
      Pindex = nPart[isimu];

  do {
    if(strncmp(line,"#",1) != 0)
     {
      /* has a haloid been written */
      if(sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)
       {
        /* if not, just get the number of particles */
        sscanf(line,"%"SCNu64, &nPartInHalo);
        
        /* construct a meaningful haloid using linenumber and MPI_rank */
        haloid = constructHaloId((uint64_t)(ihalo+1)); // +1, because ihalo has not been incremented yet!

       }
#ifdef USE_LINENUMBER_AS_HALOID
      haloid = ihalo+1; // +1, because ihalo has not been incremented yet!
#endif
      
      /* found yet another halo */
      ihalo++;
      nHalosTmp[isimu] += 1;
      halos_tmp[isimu]   = (HALOptr) realloc(halos_tmp[isimu], nHalosTmp[isimu]*sizeof(HALOS));

      /* store haloid */
      halos_tmp[isimu][ihalo].haloid = haloid;
      
      /* halos_tmp[][].Pid will be incrementally filled using realloc() */
      halos_tmp[isimu][ihalo].Pid = NULL;
      halos_tmp[isimu][ihalo].Pindex = NULL;

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

         else if(Ptype > abs(PDMbndry)) {
          /* not a meaningful type, maybe something else has been stored? */
          Ptype = 1;
         }
        
        // here we can restrict the cross-correlation to a ceratain sub-set of all particles
#ifdef ONLY_USE_PTYPE
        if(Ptype == ONLY_USE_PTYPE)
#endif
         {
          halos_tmp[isimu][ihalo].Pid    = (uint64_t *) realloc(halos_tmp[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(uint64_t));
          halos_tmp[isimu][ihalo].Pindex = (uint64_t *) realloc(halos_tmp[isimu][ihalo].Pindex, (nPartInUse+1)*sizeof(uint64_t));
          if(halos_tmp[isimu][ihalo].Pid == NULL) {
            fprintf(stderr,"read_particles: could not realloc() halos_tmp[%d][%ld].Pid for %"PRIu64"particles\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1));
            exit(-1);
          }

          halos_tmp[isimu][ihalo].Pid[nPartInUse] = Pid;
          halos_tmp[isimu][ihalo].Pindex[nPartInUse] = Pindex;
          nPartInUse++;
	  Pindex++;
         }
       }
      
      /* store number of particles in halo */
      halos_tmp[isimu][ihalo].npart = nPartInUse;  

      /* add to the total number of particles on task, plus the haloid and haloindex */    
      nPartTmp[isimu] += nPartInUse;
      totHaloSizeTmp += 2 * (nPartInUse + 1) * sizeof(uint64_t);
     }

  } while( fgets(line,MAXSTRING,fpin) != NULL);
 
  fclose(fpin);
  
  elapsed += time(NULL);

  if(LocTask == 0)
   fprintf(stderr,"\n done in %ld sec, temp num halos=%"PRIu64", total temp halo size=%zd kb.\n", 
		elapsed, nHalosTmp[isimu], totHaloSize/1024); 

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
  int64_t  ihalo;         	// the downwards for-loop is running until ihalo=-1
  uint64_t ipart, jpart, npart, nparts;
  time_t   elapsed = (time_t)0;

  npart = nPart[isimu];
  parts[isimu] = (PARTptr) calloc(npart, sizeof(PARTS));

  elapsed -= time(NULL);

#ifdef DEBUG_MPI
  fprintf(stderr,"  o creating particle<->halo mapping for file %d (PidMax=%"PRIu64", halos=%"PRIu64") on task %d ...",
		isimu, npart, nHalos[isimu], LocTask);
  fprintf(stderr, "\nAllocated %"PRIu64" parts on task %d\n", npart, LocTask);
#else
  if(LocTask == 0)
    fprintf(stderr,"  o creating particle<->halo mapping for file %d (n part=%"PRIu64") ...",isimu,npart);
#endif

  /* recording every halo it belongs to: running from low to high mass objects to allow for unique assignment! */
  for(ihalo=nHalos[isimu]-1; ihalo>=0; ihalo--)
   {
      /* Order particles here */
      nparts = halos[isimu][ihalo].npart;
      halos[isimu][ihalo].Pidord = calloc(nparts, sizeof(uint64_t));
      memcpy(halos[isimu][ihalo].Pidord, halos[isimu][ihalo].Pid, nparts * sizeof(uint64_t));   
      qsort(halos[isimu][ihalo].Pidord, nparts, sizeof(uint64_t), &cmpfunc);

    for(jpart=0; jpart<halos[isimu][ihalo].npart; jpart++)
     {
      ipart = halos[isimu][ihalo].Pindex[jpart];
#ifdef EXCLUSIVE_PARTICLES
      if(parts[isimu][ipart].nhalos == 0)
#endif
       {
           parts[isimu][ipart].nhalos++;
	 if(parts[isimu][ipart].nhalos<MAX_PARENT_HALO)
	  {
            parts[isimu][ipart].Hindex[parts[isimu][ipart].nhalos-1] = ihalo;
            parts[isimu][ipart].Hid[parts[isimu][ipart].nhalos-1] = halos[isimu][ihalo].haloid;
            parts[isimu][ipart].Pid = halos[isimu][ihalo].Pid[jpart];
	  } 
  	 else
	  {
            parts[isimu][ipart].nhalos--;
	  }
       }
     }
   }
  
  elapsed += time(NULL);

  if(LocTask == 0)
    fprintf(stderr,"\n done in %ld sec.\n",elapsed);
  return(1);
}


/*==================================================================================================
 * halo_particle_mapping:
 *
 *  Reconstruct the halo - particle map
 *  Since we are swapping particles only (not halos) between the different snapshots, 
 *  we need to reconstruct the halo structs, using the haloid data in the part structs.
 *
 *==================================================================================================*/
int halo_particle_mapping(int isimu)
{
  uint64_t ipart, locn, ihalo, jhalo, nparts;
  time_t   elapsed = (time_t)0;

  elapsed -= time(NULL);


  if(LocTask == 0)
    fprintf(stderr,"  o creating halo<->particle mapping for simu %d on task=%d...\n",isimu, LocTask);

  /* Loop over particles */
  for(ipart=0; ipart<nPart[isimu]; ipart++)
   {
    /* Now for each particle loop over the haloes it belongs to */
    for(jhalo=0; jhalo<parts[isimu][ipart].nhalos; jhalo++)
     {
         ihalo = parts[isimu][ipart].Hindex[jhalo]; // Hindex contains the local halo number
	 halos[isimu][ihalo].npart++;
         locn = halos[isimu][ihalo].npart;

         halos[isimu][ihalo].haloid = parts[isimu][ipart].Hid[jhalo]; // Hid contains the absolute halo ID
         halos[isimu][ihalo].Pid = (uint64_t *) realloc(halos[isimu][ihalo].Pid, locn*sizeof(uint64_t));
         halos[isimu][ihalo].Pindex = (uint64_t *) realloc(halos[isimu][ihalo].Pindex, locn*sizeof(uint64_t));

         halos[isimu][ihalo].Pindex[locn-1] = ipart;
         halos[isimu][ihalo].Pid[locn-1] = parts[isimu][ipart].Pid;
     }
   }

   /* now reconstruct the local ordered paricles list */
   for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
   {      
      nparts = halos[isimu][ihalo].npart;
      halos[isimu][ihalo].Pidord = calloc(nparts, sizeof(uint64_t));
      memcpy(halos[isimu][ihalo].Pidord, halos[isimu][ihalo].Pid, nparts * sizeof(uint64_t));   
      qsort(halos[isimu][ihalo].Pidord, nparts, sizeof(uint64_t), &cmpfunc);
  }

  elapsed += time(NULL);

  if(LocTask == 0)
   fprintf(stderr,"\n done in %ld sec.\n",elapsed);

  return(1);
}


/*==================================================================================================
 * cross_correlation:
 *
 *  for each halo at isimu=0 figure out how many particles are in common with khalo at isimu=1
 *
 *==================================================================================================*/
int cross_correlation(int iloop)
{
  uint64_t  ihalo, count=COUNTER;
  time_t   elapsed = (time_t)0;
  
  /*---------------------------------------------------------
   * backwards correlation
   *---------------------------------------------------------*/
    
  if(LocTask == 0 && iloop == 0)
    fprintf(stderr,"  o generating cross-correlation 0->1 for %"PRIu64" haloes on %d tasks...",
		nHalos[0], TotTask);

#ifdef DEBUG_MPI
   fprintf(stderr, "\ncreating mtree for %"PRIu64" haloes on task %d\n", nHalos[0], LocTask);
#endif 

  elapsed -= time(NULL);
  /* cross-correlation simu0->simu1. When creating the m_tree now we need to 
   * take into account that simu1 is split in TotTask files, and we now are
   * creating the m_tree for the iLoop-th one
   *  */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
        create_mtree(ihalo, 0, 1, iloop);

        if(ihalo == count) 
        {
#ifdef DEBUG_MPI	
           fprintf(stderr, "on task=%d, jpart=%"PRIu64"/%"PRIu64" done.\n", LocTask, ihalo, nHalos[0]);
#else
           fprintf(stderr, ".");
#endif
           count += 500;
        }
 }

  elapsed += time(NULL);

  if(LocTask == 0)
    fprintf(stderr,"\n done in %ld sec.\n",elapsed);

  return(1);
}


/*==================================================================================================
 * create_mtree
 *==================================================================================================*/
int create_mtree(uint64_t ihalo, int isimu0, int isimu1, int iloop)
{
  uint64_t  khalo, ncroco, jcroco;
  uint64_t  *common=NULL;

  /* At each redshift you need to allocate and initialize the mtree struct only for the first chunk 
   * of the total halo catalogue, then realloc the other connections for the haloes first read in by
   * the other tasks */
  if(iloop == 0)
  {
    halos[isimu0][ihalo].mtree  = (MTREEptr) calloc(1,sizeof(MTREE)); 
    halos[isimu0][ihalo].global_ncroco = 0;  
  }
  
  /* common[] records how many particles ihalo(isimu0) has in common with khalo(isimu1) */
  common = (uint64_t *) calloc(nHalos[isimu1], sizeof(uint64_t));


/* use this custom linear algorithm to check for the matching particles */
  omp_set_num_threads(NUM_OMP_THREADS);
#pragma omp parallel for private(khalo) shared(nHalos, common, ihalo, isimu0, isimu1)
  for(khalo=0; khalo<nHalos[isimu1]; khalo++) 
    intersection(isimu0, isimu1, ihalo, khalo, common);

  /* determine number of credible cross-correlations */
  ncroco = 0;
  for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
    if(common[khalo] > MINCOMMON)
      ncroco++;
  }

  /* Update the number of connections found in total */
  halos[isimu0][ihalo].global_ncroco += ncroco; 

  /* does not make sense to continue if there are no cross-correlations */
  if(ncroco > 0) 
  {
    /* allocate memory for cross-correlations */
    halos[isimu0][ihalo].mtree  = 
	(MTREEptr) realloc(halos[isimu0][ihalo].mtree, halos[isimu0][ihalo].global_ncroco * sizeof(MTREE)); 
  
    /* The halo[][].mtree for iloop>0 may already contain some connections to haloes 
     * so we need to initialize the variables */
    jcroco = halos[isimu0][ihalo].global_ncroco - ncroco; 

    /* Store mtree inside halos[][] structure in incorrect order. 
     * The correct order according to merit function will be implemented when all haloes are stored */
    for(khalo=0; khalo<nHalos[isimu1]; khalo++) 
    {
      if(common[khalo] > MINCOMMON)
      {
        halos[isimu0][ihalo].mtree[jcroco].id[0]    = ihalo;
        halos[isimu0][ihalo].mtree[jcroco].haloid[0]= halos[isimu0][ihalo].haloid;
        halos[isimu0][ihalo].mtree[jcroco].npart[0] = halos[isimu0][ihalo].npart;
        halos[isimu0][ihalo].mtree[jcroco].common   = common[khalo];
        halos[isimu0][ihalo].mtree[jcroco].id[1]    = khalo;
        halos[isimu0][ihalo].mtree[jcroco].haloid[1]= halos[isimu1][khalo].haloid;
        halos[isimu0][ihalo].mtree[jcroco].npart[1] = halos[isimu1][khalo].npart;
        jcroco++;
      }
    }  
  } /* if ncroco>0 */
  
    /* free temporary structures */
    if(common) free(common);

  return(1);
}


/*==================================================================================================
 * write_mtree
 *==================================================================================================*/
int write_mtree(int isimu0, char OutFile[MAXSTRING])
{
  uint64_t  ihalo, jhalo, icroco;
  FILE *fpout, *fpout_idx;
  char outname[MAXSTRING], outname_idx[MAXSTRING];
  time_t   elapsed = (time_t)0;

#ifdef ORDER
  MTREE    *mtree;
  
  uint64_t npart[2], jcroco, ncroco, common;
  double *merit;
  long unsigned *idx;

  npart[0] = 0;
  npart[1] = 0;

    /* Order haloes in mtree by merit function */
    for(ihalo=0; ihalo<nHalos[isimu0]; ihalo++) 
    {
      ncroco = halos[isimu0][ihalo].global_ncroco;

      mtree  = (MTREEptr)        calloc(ncroco, sizeof(MTREE));
      merit = (double *) calloc(ncroco, sizeof(double));
      idx   = (long unsigned *) calloc(ncroco, sizeof(long unsigned));

      if(ncroco > 0)
      {
	for(icroco=0; icroco<ncroco; icroco++)
        {
          common = halos[isimu0][ihalo].mtree[icroco].common;
          npart[0] = halos[isimu0][ihalo].mtree[icroco].npart[0]; 
          npart[1] = halos[isimu0][ihalo].mtree[icroco].npart[1]; 

          mtree[icroco].id[0]     = halos[isimu0][ihalo].mtree[icroco].id[0];
          mtree[icroco].haloid[0] = halos[isimu0][ihalo].mtree[icroco].haloid[0];
          mtree[icroco].npart[0]  = halos[isimu0][ihalo].mtree[icroco].npart[0];
          mtree[icroco].common    = halos[isimu0][ihalo].mtree[icroco].common;
          mtree[icroco].id[1]     = halos[isimu0][ihalo].mtree[icroco].id[1];
          mtree[icroco].haloid[1] = halos[isimu0][ihalo].mtree[icroco].haloid[1];
          mtree[icroco].npart[1]  = halos[isimu0][ihalo].mtree[icroco].npart[1];

          merit[icroco] = pow2((double)common)/((double)npart[0]*(double)npart[1]);
	} 

	    /* order by merit function */
	    indexx((long unsigned)ncroco, merit-1, idx-1);
	
	/* Free the mtree of halos and realloc it from scratch so we can store the tree in the right order */
	free(halos[isimu0][ihalo].mtree);
	halos[isimu0][ihalo].mtree = (MTREEptr) calloc(ncroco, sizeof(MTREE));

	/* Now gather again the */
        for(jcroco=0; jcroco<ncroco; jcroco++)
        {
          icroco = idx[ncroco-1-jcroco]-1;

	  halos[isimu0][ihalo].mtree[jcroco].id[0] = mtree[icroco].id[0];
	  halos[isimu0][ihalo].mtree[jcroco].haloid[0] = mtree[icroco].haloid[0];
	  halos[isimu0][ihalo].mtree[jcroco].npart[0] = mtree[icroco].npart[0];
	  halos[isimu0][ihalo].mtree[jcroco].common = mtree[icroco].common;
	  halos[isimu0][ihalo].mtree[jcroco].id[1] = mtree[icroco].id[1];
	  halos[isimu0][ihalo].mtree[jcroco].haloid[1] = mtree[icroco].haloid[1];
	  halos[isimu0][ihalo].mtree[jcroco].npart[1] = mtree[icroco].npart[1];
        }
      }
    }
#endif // ORDER

  elapsed -= time(NULL);
  fprintf(stderr,"  o writing cross-correlation for %"PRIu64" haloes on task %d...\n",nHalos[0], LocTask);
  
  sprintf(outname,"%s_mtree",OutFile);
  strcpy(outname_idx, outname);
  strcat(outname_idx, "_idx");
  
  fpout = fopen(outname,"w");
  if(fpout == NULL) {
    fprintf(stderr,"could not open file %s on task %d,\nexiting\n",outname, LocTask);
    exit(0);
  }
  
  fpout_idx = fopen(outname_idx,"w");
  if(fpout_idx == NULL) {
    fprintf(stderr,"could not open file %s on task %d,\nexiting\n",outname_idx, LocTask);
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

    if(DeltaTask > 0) /* in this case you had to load balance and the largest halos will be in the last positions */
      jhalo = nHalos[0] - ihalo -1;
    else
      jhalo = ihalo;

    if(halos[isimu0][jhalo].global_ncroco > 0) {
      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[isimu0][jhalo].mtree[0].haloid[0], halos[isimu0][jhalo].mtree[0].haloid[1]);
      fflush(fpout_idx);
      
#ifdef SUSSING2013
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[isimu0][jhalo].haloid,
              halos[isimu0][jhalo].global_ncroco);
#else // SUSSING2013
      fprintf(fpout,"%"PRIu64"  %"PRIu64"  %"PRIu64"\n",
              halos[isimu0][jhalo].haloid,
              halos[isimu0][jhalo].npart,
              halos[isimu0][jhalo].global_ncroco);
#endif // SUSSING2013
      fflush(fpout);
      
      for(icroco=0; icroco<halos[isimu0][jhalo].global_ncroco; icroco++) {
#ifdef SUSSING2013
        fprintf(fpout,"%"PRIu64"\n",
                halos[isimu0][jhalo].mtree[icroco].haloid[1]);
#else // SUSSING2013
        fprintf(fpout,"  %"PRIu64"  %"PRIu64"  %"PRIu64"\n",
                halos[isimu0][jhalo].mtree[icroco].common,
                halos[isimu0][jhalo].mtree[icroco].haloid[1],
                halos[isimu0][jhalo].mtree[icroco].npart[1]);
#endif // SUSSING2013
        fflush(fpout);
      }
    }
#ifdef SUSSING2013
    else {
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[isimu0][jhalo].haloid,
              halos[isimu0][jhalo].ncroco);      
    }
#endif // SUSSING2013
  }
  
  /* close files */
  fclose(fpout);
  fclose(fpout_idx);

  elapsed += time(NULL);

   if(LocTask == 0)
    fprintf(stderr,"\n done in %ld sec.\n",elapsed);
  return(1);
}


/*==================================================================================================
 * remove_network_connections
 *==================================================================================================*/
int remove_network_connections(void)
{
  unsigned int ihalo;
  time_t   elapsed = (time_t)0;

  elapsed -= time(NULL);

  if(LocTask == 0)
    fprintf(stderr,"  o removing network connections");

  /* clean connections simu0->simu1 */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    clean_connection(ihalo, 0, 1);
  }
  elapsed += time(NULL);

  if(LocTask == 0)
   fprintf(stderr,"\n done in %ld sec.\n",elapsed);

  return(1);
}


/*==================================================================================================
 * clean_connection //FIXME
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
  for(icroco=0; icroco<halos[isimu0][ihalo].global_ncroco; icroco++) {
    jhalo = halos[isimu0][ihalo].mtree[icroco].id[1];
    
    /* check whether the present halo is the most likely descendant of this progenitor */
    if(max_merit(jhalo, isimu1) == ihalo) {
      // keep jhalo in mtree-list
      ncroco_new++;
      mtree = (MTREEptr) realloc(mtree, ncroco_new*sizeof(MTREE));
      
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
      // remove jhalo from mtree-list and henc do not add it to the new mtree[] list
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
  if(halos[isimu0][ihalo].mtree) free(halos[isimu0][ihalo].mtree);
  halos[isimu0][ihalo].mtree  = mtree;

  return(1);
}

/*==================================================================================================
 * max_merit
 *==================================================================================================*/
uint64_t max_merit(uint64_t jhalo, int isimu)
{
  /* mtree[] is ordered by merit and hence we only need to check the first entry */
  if(halos[isimu][jhalo].global_ncroco > 0) {
    return(halos[isimu][jhalo].mtree[0].id[1]);
  }
  else {
#ifdef DEBUG
    fprintf(stderr,"jhalo=%ld in isimu=%d does not point to anywhere!?\n",jhalo,isimu);
#endif
    return(0);
  }
}

  /* Reallocate halos after having swapped particles among tasks */
int alloc_halos(int isimu)
{
  uint64_t i;

  halos[isimu] = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));

  for(i=0; i<nHalos[isimu]; i++)
  {
    halos[isimu][i].Pid = (uint64_t *) calloc(1, sizeof(uint64_t));
    halos[isimu][i].Pindex = (uint64_t *) calloc(1, sizeof(uint64_t));
  }

  /* Init halo particle number to zero */
  for(i=0; i<nHalos[isimu]; i++)
    halos[isimu][i].npart = 0;

 return(1);
}



int add_halos(int ifile, int isimu)
{
  uint64_t nparts, ihalo, min_halo;
  size_t sizeHalo;
  
  if(ifile == 0) /* init the halo struct when reading the first file */ 
  {  
     min_halo = 0;
     nPart[isimu] = nPartTmp[isimu];
     nHalos[isimu] = nHalosTmp[isimu];
     totHaloSize = totHaloSizeTmp;
     alloc_halos(isimu);
  } else { /* realloc and make room for the new halos */
     min_halo = nHalos[isimu];
     nPart[isimu] += nPartTmp[isimu];
     nHalos[isimu] += nHalosTmp[isimu];
     totHaloSize += totHaloSizeTmp;
     halos[isimu] = (HALOptr) realloc(halos[isimu], ( (nHalos[isimu]+1)) * sizeof(HALOS));
  }

  /* copy the halos into the old structure and free up memory */

  for(ihalo=min_halo; ihalo<nHalos[isimu]; ihalo++)
  {
     nparts = halos_tmp[isimu][ihalo-min_halo].npart;
     halos[isimu][ihalo].npart = halos_tmp[isimu][ihalo-min_halo].npart;
     halos[isimu][ihalo].haloid = halos_tmp[isimu][ihalo-min_halo].haloid;
     sizeHalo = nparts * sizeof(uint64_t);

     halos[isimu][ihalo].Pid = (uint64_t *) malloc(sizeHalo);
     halos[isimu][ihalo].Pindex = (uint64_t *) malloc(sizeHalo);

     memcpy(halos[isimu][ihalo].Pid, halos_tmp[isimu][ihalo-min_halo].Pid, sizeHalo);   
     memcpy(halos[isimu][ihalo].Pindex, halos_tmp[isimu][ihalo-min_halo].Pindex, sizeHalo);

     /* free temp memory as the particles are copied */
     free(halos_tmp[isimu][ihalo-min_halo].Pindex);
     free(halos_tmp[isimu][ihalo-min_halo].Pid);
  }

  /* free the temporary structures */	
  free(halos_tmp[isimu]);
  nHalosTmp[isimu]=0;
  nPartTmp[isimu]=0;
  totHaloSizeTmp=0;

  return(1);
} 

  /* Clean up the halos */
int free_halos(int isimu)
{
  uint64_t i;

  for(i=0; i<nHalos[isimu]; i++)
  {
    if(halos[isimu][i].Pid != NULL)
      free(halos[isimu][i].Pid);

    if(halos[isimu][i].Pindex != NULL)
      free(halos[isimu][i].Pindex);

    if(halos[isimu][i].mtree != NULL)
      free(halos[isimu][i].mtree);
  }

  free(halos[isimu]);

 return(1);
}


/* custiom linear time algorithm */
void intersection(int isimu0, int isimu1, uint64_t ihalo, uint64_t khalo, uint64_t *common)
{
	uint64_t ipart=0, kpart=0, apart=0, bpart=0;

	while (ipart < halos[isimu0][ihalo].npart && kpart < halos[isimu1][khalo].npart)
	{
		apart = halos[isimu0][ihalo].Pidord[ipart];
		bpart = halos[isimu1][khalo].Pidord[kpart];

		if (apart == bpart)
		{
			common[khalo]++;
			ipart++;
			kpart++;
		}
		else if (apart < bpart)
		{
			ipart++;
		}
		else
		{
			kpart++;
		}
	}
    
}

/* This function is called by qsort() and bsearch() needed for the matching of the particles */
int cmpfunc(const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

/* construct a unique haloid based upon line number (=ihalo) and MPI_rank */
uint64_t constructHaloId(uint64_t ihalo)
{
  uint64_t haloid;
  
  // we allow for 64-48=16 bits for the MPI_rank (up to 65536 MPI tasks)
  haloid = (uint64_t) LocTask;
  haloid = haloid << 48;

  // simply overlay line number now
  haloid = haloid | ihalo;
  
  return(haloid);
}

/* Debug routines to check the content of the structures sent across tasks */
#ifdef DEBUG_MPI
void check_parts(int isimu)
{
  int i=0, ntot;
  ntot = nPart[isimu];
	
    for(i=0; i<ntot; i++)
      fprintf(stderr,"Local particle %d on task=%d has ID %"PRIu64"\n", i, LocTask, parts[isimu][i].Pid);

}

void check_halos(int isimu)
{
  unsigned int i=0, j=0, ntot;
  ntot = nHalos[isimu];

    for(i=0; i<ntot; i++)
    {
      fprintf(stderr,"Local halo %d on task=%d has ID %"PRIu64"\n", i, LocTask, halos[isimu][i].haloid);
        for(j=0; j<halos[isimu][i].npart; j++)
           fprintf(stdout, "\t\tpart_id=%"PRIu64"\n", halos[isimu][i].Pid[j]);
    }
}


#ifdef DEBUG_LOG
void dump_log_halo(int isimu, int ihalo)
{
  uint64_t i=0, ntot, haloid;
  char log_name[200];
  FILE *log_halo;

  ntot = halos[isimu][ihalo].npart;
  haloid = halos[isimu][ihalo].haloid;

  sprintf(log_name, "/home/carlesi/MERGER_TREE/TEST/haloID_%03d.task_%02d.halo", ihalo, LocTask);
  log_halo = fopen(log_name, "w");

	//fprintf(stderr, "Printing halo %03d on task %d to %s\n", ihalo, LocTask, log_name);

	fprintf(log_halo, "# NPART=%"PRIu64"\t ID=%"PRIu64"\n", 
		halos[isimu][ihalo].npart,
		halos[isimu][ihalo].haloid);

    for(i=0; i<ntot; i++)
    {
	fprintf(log_halo, "%"PRIu64"\t%"PRIu64"\n", 
		halos[isimu][ihalo].Pid[i],
		halos[isimu][ihalo].Pindex[i]);
    }
 
  fclose(log_halo);
}
#endif
#endif

#ifdef MERGER_RATIO
/*==================================================================================================
 * assign_progenitors:
 *
 *  get statistics for multiple progenitors (e.g. mass ratios, etc.)
 *
 *  update 10/10/2007:
 *       each subhalo shares the most particles with its host :-(
 *       -> tried to fix this issue..
 *
 *==================================================================================================*/
int assign_progenitors(char OutFile[MAXSTRING])
{
  FILE   *fpin, *fpout, *fpout_merger;
  char    OutFile_mtree[MAXSTRING], OutFile_idx[MAXSTRING], OutFile_merger[MAXSTRING], line[MAXSTRING];
  long unsigned   *id1, *npart1, *common, *id2, *npart2, *idx, iline;
  double  xcommon, xnpart1, xnpart2, *ratio;
  long    prev_id1, nprog, iprog;
  
  fprintf(stderr,"  o assigning progenitors ...");
  
  strcpy(OutFile_mtree, OutFile);
  strcat(OutFile_mtree, "_mtree");
  strcpy(OutFile_idx, OutFile_mtree);
  strcat(OutFile_idx, "_progs");
  strcpy(OutFile_merger, OutFile_mtree);
  strcat(OutFile_merger, "_merger");
  
  /* open output file */
  fpout = fopen(OutFile_idx,"w");
  if(fpout == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting\n",OutFile_idx);
    exit(0);
   }
  fprintf(fpout,"#id(1) Np(2) iprog1(3) Np1(4) ncommon1(5) iprog2(6) Np2(7) ncommon2(8) iprog3(9) Np3(10) ncommon3(11)\n");
  
  /* open output file */
  fpout_merger = fopen(OutFile_merger,"w");
  if(fpout_merger == NULL)
   {
    fprintf(stderr,"could not open file %s\nexiting\n",OutFile_merger);
    exit(0);
   }
  fprintf(fpout_merger,"#id(1) iprog1(2) iprog2(3) common2/ncommon1(4)\n");

  /* read *_mtree file */
  read_mtree(OutFile);
  
  prev_id1  = mtree[0].id[0];  
  nprog     = 0;
  id1    = (long *) calloc(1, sizeof(long));
  id2    = (long *) calloc(1, sizeof(long));
  npart1 = (long *) calloc(1, sizeof(long));
  npart2 = (long *) calloc(1, sizeof(long));
  common = (long *) calloc(1, sizeof(long));
  
  for(iline=0; iline<nlines; iline++)
   {
    if(mtree[iline].id[0] == prev_id1)
     {
      /* make room for one more progenitor */
      nprog++;
      id1    = (long *) realloc(id1,    (nprog)*sizeof(long));
      id2    = (long *) realloc(id2,    (nprog)*sizeof(long));
      npart1 = (long *) realloc(npart1, (nprog)*sizeof(long));
      npart2 = (long *) realloc(npart2, (nprog)*sizeof(long));
      common = (long *) realloc(common, (nprog)*sizeof(long));
      
      /* copy information from mtree[] */
      id1[nprog-1]    = mtree[iline].id[0];
      id2[nprog-1]    = mtree[iline].id[1];
      npart1[nprog-1] = mtree[iline].npart[0];
      npart2[nprog-1] = mtree[iline].npart[1];
      common[nprog-1] = mtree[iline].common;
     }
    else
     {
      ratio = (double *)        calloc(nprog+1, sizeof(double));
      idx   = (long unsigned *) calloc(nprog+1, sizeof(long unsigned));
      
      for(iprog=0; iprog<nprog; iprog++)
       {
        /* calculate progenitor criterion */
        xcommon          = (double)common[iprog];
        xnpart1          = (double)npart1[iprog];
        xnpart2          = (double)npart2[iprog];
        ratio[iprog]     = pow2(xcommon)/(xnpart1*xnpart2);
       }
      /* sort all progenitor by merit function */
      indexx(nprog, ratio-1, idx-1);
      
      /*-----------------------------------------------------
       * dump information to file
       *-----------------------------------------------------*/
      if(nprog > 2)
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1],
                id2[idx[nprog-2]-1],npart2[idx[nprog-2]-1],common[idx[nprog-2]-1],
                id2[idx[nprog-3]-1],npart2[idx[nprog-3]-1],common[idx[nprog-3]-1]);
        
        /* iprog2 is most credible second progenitor */
        if(id2[idx[nprog-1]-1]<id2[idx[nprog-2]-1] && id2[idx[nprog-2]-1]<id2[idx[nprog-3]-1])
         {
          if((double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-2]-1],
                    (double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1]);
         }
        /* iprog3 is most credible second progenitor */
        else if(id2[idx[nprog-1]-1]>id2[idx[nprog-2]-1] && id2[idx[nprog-2]-1]<id2[idx[nprog-3]-1])
         {
          if((double)common[idx[nprog-3]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-3]-1],
                    (double)common[idx[nprog-3]-1]/(double)common[idx[nprog-1]-1]);
         }
        /* nothing else to do as only one real progenitor exists */
       }
      
      else if (nprog > 1)
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       %10ld %10ld %10ld       -1 -1 -1\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1],
                id2[idx[nprog-2]-1],npart2[idx[nprog-2]-1],common[idx[nprog-2]-1]);            

        /* this indicates that the halo itself is a subhalo */
        if(id2[idx[nprog-1]-1]>id2[idx[nprog-2]-1])
         {
          /* nothing to do as only one real progenitor exists */
         }
        else
         {
          /* iprog2 is most credible second progenitor */
          if((double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1] > MERGER_RATIO)
            fprintf(fpout_merger,"%10ld %10ld %10ld %16.8lf\n",
                    prev_id1,id2[idx[nprog-1]-1],id2[idx[nprog-2]-1],
                    (double)common[idx[nprog-2]-1]/(double)common[idx[nprog-1]-1]);
         }
       }
      else
       {
        fprintf(fpout,"%10ld %10ld       %10ld %10ld %10ld       -1 -1 -1       -1 -1 -1\n",
                prev_id1,npart1[0],
                id2[idx[nprog-1]-1],npart2[idx[nprog-1]-1],common[idx[nprog-1]-1]);
        
        /* nothing else to do as there are no multiple progenitors */
       }
      
      
      free(ratio);
      free(idx);
      
      /* start a new progenitor list */
      free(id1);
      free(id2);
      free(npart1);
      free(npart2);
      free(common);
      
      nprog  = 1;
      id1    = (long *) calloc(nprog, sizeof(long));
      id2    = (long *) calloc(nprog, sizeof(long));
      npart1 = (long *) calloc(nprog, sizeof(long));
      npart2 = (long *) calloc(nprog, sizeof(long));
      common = (long *) calloc(nprog, sizeof(long));
      
      id1[nprog-1]    = mtree[iline].id[0];
      id2[nprog-1]    = mtree[iline].id[1];
      npart1[nprog-1] = mtree[iline].npart[0];
      npart2[nprog-1] = mtree[iline].npart[1];
      common[nprog-1] = mtree[iline].common;
      
      prev_id1        = id1[nprog-1];
      
     }
   }
  
  if(id1 != NULL) free(id1);
  if(id2 != NULL) free(id2);
  if(npart1 != NULL) free(npart1);
  if(npart2 != NULL) free(npart2);
  if(common != NULL) free(common);
  
  fclose(fpout);
  fclose(fpout_merger);
  
   if(LocTask == 0)
    fprintf(stderr," done\n");
  
  return(1);
}

/*==================================================================================================
 * read_mtree:
 *
 *       simply reads in the *_mtree file and 
 *       puts it into the array of structures mtree[iline].XYZ
 *
 * Note: at this stage we just treat these entries as "lines" -> no connection to halos yet!
 *
 *==================================================================================================*/
void read_mtree(char *prefix)
{
  long unsigned iline;
  char          line[MAXSTRING], outname[MAXSTRING];
  FILE         *fpin;
  
  sprintf(outname,"%s_mtree",prefix);
  if((fpin = fopen(outname,"r")) == NULL)
   {
    fprintf(stderr,"cannot open  %s\nEXIT\n",outname);
    exit(0);
   }
  
  // count number of lines
  nlines = 0;
  while(!feof(fpin))
   {
    nlines++;
    fgets(line,MAXSTRING,fpin);
   }
  
  //fprintf(stderr,"  o found %ld lines in:  %s\n",nlines,infile);
  
  // allocate memory
  mtree = (MTREEptr) realloc((MTREEptr)mtree, nlines*sizeof(MTREE));
  
  // actually read the file
  rewind(fpin);
  for(iline=0; iline<nlines; iline++)
   {
    // read next line from file
    fgets(line,MAXSTRING,fpin);
    sscanf(line,"%"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64,
           &(mtree[iline].id[0]), 
           &(mtree[iline].npart[0]), 
           &(mtree[iline].common), 
           &(mtree[iline].id[1]), 
           &(mtree[iline].npart[1]));
    //fprintf(stderr,"iline=%ld id[0]=%ld\n",iline,mtree[iline].id[0]);
   }
  
  fclose(fpin);
}

#endif // MERGER_RATIO
