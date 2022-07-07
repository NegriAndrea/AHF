/* $Id: io_mgadget.c,v 1.20 2007/12/10 13:13:08 knolli Exp $ */

/**
 * \file io_gadget.c
 *
 * Provides functions for reading and writing Gadget files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "../define.h"    // needed for all those #ifdef features here
#include "../common.h"    // needed for global_io.params->reader (but f->size_mycomm should equally work!?)

#include "io_mgadget.h"
#include "io_gadget.h"
#include "io_util.h"

/* includes needed for adjustments of the file descriptors */
#ifdef CHECK_RLIMIT_NOFILE
#include <errno.h>
#include <sys/resource.h>
#include <unistd.h>
#define NFILES_EXTRA 50
#endif



/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_mgadget_t
io_mgadget_open_withoutNCPUREADING_EQ_NFILES(io_logging_t log,
                char *fname,
                io_file_swap_t swapped,
                io_file_mode_t mode,
                uint32_t reader)
{
	int32_t i;
	io_mgadget_t f;
	char **fnames;
#ifdef CHECK_RLIMIT_NOFILE
  struct rlimit limit;
#endif
#if (defined WITH_MPI && defined BCASTHEADER)
  int32_t *slen;
#endif
  
	/* XXX THIS IS CURRENTLY ONLY FOR READING! */
	if (mode != IO_FILE_READ)
		return NULL;
  
	/* Get memory for the structure */
	f = (io_mgadget_t)malloc(sizeof(io_mgadget_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_mgdaget_open(): sizeof(io_mgadget_struct_t)");
		return NULL;
	}
  
	/* Start filling the structure */
  
	/* Okay, we are a Multiple Gadget file */
	f->ftype = IO_FILE_MGADGET;
  
	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	MPI_Comm_split(MPI_COMM_WORLD, 1, f->rank, &(f->mycomm));
	MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	MPI_Comm_size(f->mycomm, &(f->size_mycomm));
#	endif
  
	/* Split the filename in path and stem */
	f->path = NULL;
	f->stem = NULL;
	if (  io_util_split_pathfname(fname, &(f->path), &(f->stem)) == NULL) {
		io_logging_fatal(log, "io_mgdaget_open(): Could not split %s in path and filename.", fname);
    if(f->path) free(f->path);
    if(f->stem) free(f->stem);
		free(f);
		return NULL;
	}
	io_logging_msg(log, INT32_C(1), "Will look in %s for %s", f->path, f->stem);
  
#if (defined WITH_MPI && defined BCASTHEADER)
  if(f->rank_mycomm == 0)
#endif
   {
    /* Get the filenames */
    f->numfiles = io_util_findfiles(f->path, f->stem, "%i", "", &fnames);
    
    if (f->numfiles <= 0) {
      io_logging_fatal(log, "io_mgdaget_open(): Could not open anything starting with %s in %s.", f->stem, f->path);
      if(f->path)  free(f->stem);
      if(f->stem)  free(f->path);
      free(f);
      return NULL;
    }
   }
  
#if (defined WITH_MPI && defined BCASTHEADER)
  /* first round of broadcasts (needed to array allocations) */
  MPI_Bcast(&(f->numfiles), sizeof(int32_t), MPI_BYTE, 0, f->mycomm);
  
  /* determine the length of each string from root MPI task */
  slen = (int32_t *) malloc(f->numfiles*sizeof(int32_t));
  if(f->rank_mycomm == 0) {
    for(i=0; i<f->numfiles; i++) {
      slen[i] = strlen(fnames[i]) + 1; // incl. the terminating character
    }
  }
  MPI_Bcast(slen, f->numfiles*sizeof(int32_t), MPI_BYTE, 0, f->mycomm);
  
  /* allocate memory for fnames[][] on all MPI tasks but the root one */
  if(f->rank_mycomm != 0) {
    fnames = (char **)malloc(f->numfiles*sizeof(char *));
    for(i=0; i<f->numfiles; i++) {
      fnames[i] = (char *)malloc(slen[i]*sizeof(char));
    }
  }
  
  /* second round of broadcasts (actual filenames) */
  for(i=0; i<f->numfiles; i++) {
    MPI_Bcast(fnames[i], slen[i]*sizeof(char), MPI_BYTE, 0, f->mycomm);
  }
#endif /* defined WITH_MPI && defined BCASTHEADER */
  
#if (defined WITH_MPI && defined NCPUREADING_EQ_NFILES)
  if(f->numfiles != global_io.params->reader) {
		io_logging_fatal(log, "io_mgdaget_open(): You are using NCPUREADING_EQ_NFILES and hence the number of files %d should be the number of MPI tasks reading %d (size_mycomm=%d).",
                     f->numfiles, global_io.params->reader, f->size_mycomm);
		if(f->path)  free(f->stem);
		if(f->stem)  free(f->path);
		free(f);
    return NULL;
  }
#endif
  
#ifdef DEBUG_MGADGET
  io_logging_msg(log,6,"\nio_mgdaget_open(): found the following %d files:\n",f->numfiles);
  for (i=0; i<f->numfiles; i++){
    io_logging_msg(log,6,"%s",fnames[i]);
  }
#endif
  
#ifdef CHECK_RLIMIT_NOFILE
  /* check whether we are able to use that many readers */
  if(getrlimit(RLIMIT_NOFILE, &limit) == 0)
   {
#ifdef WITH_MPI
    fprintf(stderr,"\nPRESENT FILE DESCRIPTOR LIMITS (rank=%d of %d):\n",f->rank,f->size);
#else
    fprintf(stderr,"\nPRESENT FILE DESCRIPTOR LIMITS :\n");
#endif
    fprintf(stderr,"   rlim_cur = %llu\n",limit.rlim_cur);
    fprintf(stderr,"   rlim_max = %llu\n",limit.rlim_max);
   }
  else
   {
#ifdef WITH_MPI
    fprintf(stderr,"getrlimit() failed with errno=%d (rank %d of %d)\n", errno,f->rank,f->size);
#else
    fprintf(stderr,"getrlimit() failed with errno=%d\n", errno);
#endif
    exit(1);
   }
  if(limit.rlim_max < f->numfiles+NFILES_EXTRA)
   {
#ifdef WITH_MPI
    fprintf(stderr,"\nPROBLEM (rank=%d of %d):\n",f->rank,f->size);
#else
    fprintf(stderr,"\nPROBLEM:\n");
#endif
    fprintf(stderr,"   only the super-user can increase beyond the maximum limit:\n");
    fprintf(stderr,"     numfiles = %d\n",f->numfiles);
    fprintf(stderr,"     rlim_max = %llu (incl. 3 for stdin, stdout, and stderr)\n",limit.rlim_max);
    fprintf(stderr,"   ABORTING\n");
    exit(1);
   }
  if (limit.rlim_cur < f->numfiles+NFILES_EXTRA)
   {
#ifdef WITH_MPI
    fprintf(stderr,"\nincreasing number of file descriptors on rank %d (of %d) to %d\n",f->rank,f->size,f->numfiles+NFILES_EXTRA);
#else
    fprintf(stderr,"\nincreasing number of file descriptors to %d\n",f->numfiles+NFILES_EXTRA);
#endif
    limit.rlim_cur = (rlim_t) (f->numfiles+NFILES_EXTRA);
    if (setrlimit(RLIMIT_NOFILE, &limit) != 0)
     {
#ifdef WITH_MPI
      fprintf(stderr,"setrlimit() failed with errno=%d (rank %d of %d)\n", errno,f->rank,f->size);
#else
      fprintf(stderr,"setrlimit() failed with errno=%d\n", errno);
#endif
      exit(1);
     }
    
    // check whether it worked or not
    if(getrlimit(RLIMIT_NOFILE, &limit) == 0)
     {
#ifdef WITH_MPI
      fprintf(stderr,"\nNEW LIMITS (rank=%d of %d):\n",f->rank,f->size);
#else
      fprintf(stderr,"\nNEW LIMITS:\n");
#endif
      fprintf(stderr,"   rlim_cur = %llu\n",limit.rlim_cur);
      fprintf(stderr,"   rlim_max = %llu\n",limit.rlim_max);
     }
    else
     {
#ifdef WITH_MPI
      fprintf(stderr,"getrlimit() failed with errno=%d (rank %d of %d)\n", errno,f->rank,f->size);
#else
      fprintf(stderr,"getrlimit() failed with errno=%d\n", errno);
#endif
      exit(1);
     }
   }
#endif
  
	/* Glue the files into the MGadget structure */
	f->files = (io_gadget_t *)malloc( sizeof(io_gadget_t)*(f->numfiles));
	if (f->files == NULL) {
    io_logging_memfatal(log,  "io_mgadget structure (2)");
		for (i=0; i<f->numfiles; i++)
			free(fnames[i]);
		free(fnames);
		free(f->stem);
		free(f->path);
		free(f);
		return NULL;
	}
  
#if (defined WITH_MPI && defined BCASTHEADER)
  /* only rank==0 will open all those files */
  if(f->rank_mycomm == 0) {
#endif
    for (i=0; i<f->numfiles; i++)
     {
#ifdef DEBUG_MGADGET
      io_logging_msg(log,6,"\ntrying to open file %s (%d of %d) for reading ... ", fnames[i],i+1,f->numfiles);
#endif
      
#		ifdef WITH_MPI
      /* TODO
       * THIS	IS JUST A NASTY HACK TO PREVENT io_gadget.c FROM REDOING
       * THE MPI-SPLIT. CALLED FROM HERE IT IS ONLY SUPPOSED TO ACT AS
       * A DUMMY INTERFACE TO A GADGET FILE.
       * TODO
       */
      (f->files)[i] = io_gadget_open(log, fnames[i], swapped, mode, f->size+1);
#		else
      (f->files)[i] = io_gadget_open(log, fnames[i], swapped, mode, reader);
#		endif
      
      if ((f->files)[i] == NULL) {
        int32_t j;
        for (j=i; i<f->numfiles; j++)
          free(fnames[j]);
        free(fnames);
        while (i>0) {
          i--;
          io_gadget_close(log, &((f->files)[i]));
        }
        free(f->stem);
        free(f->path);
        free(f);
        return NULL;
      }
      
#ifdef FOPENCLOSE
      //fprintf(stderr,"FOPENCLOSE: opened %s ... ",((f->files)[i])->fname);
      fclose(((f->files)[i])->file);
      // DO NOT SET THE POINTER TO NULL AS file!=NULL INDICATES FOR FOLLOWING ROUTINES THAT FILE COULD BE OPENED
      //fprintf(stderr,"and closed file (temporarily)\n");
#endif
     }
#if (defined WITH_MPI && defined BCASTHEADER)
  }
  /* we nevertheless need to generate all the file structures ready to receive the broadcast */
  else {
    for (i=0; i<f->numfiles; i++) {
      ((f->files)[i]) = (io_gadget_t)calloc(sizeof(io_gadget_struct_t),1);
      ((f->files)[i])->fname = NULL;
    }
  }
#endif
  
  
#if (defined WITH_MPI && defined BCASTHEADER)
  // BROADCAST (f->files)[i] values to all other MPI tasks!
  for (i=0; i<f->numfiles; i++) {
    MPI_Bcast((f->files)[i], sizeof(io_gadget_struct_t), MPI_BYTE, 0, f->mycomm);
    
    /* only rank == 0 properly allocated the memory for the filename and copied the data there */
    if(f->rank_mycomm != 0) {
      ((f->files)[i])->fname = (char *)malloc(sizeof(char) * (strlen(fnames[i]) + 1));
      strncpy(((f->files)[i])->fname, fnames[i], strlen(fnames[i])+1);
    }
  }
#endif
  
  /* the whole char *fnames[] is not needed anymore */
  for (i=0; i<f->numfiles; i++) {
    free(fnames[i]);
  }
	free(fnames);
  
	/* Set initial values */
	f->no_part     = UINT64_C(0);
	f->multimass   = INT8_C(0);
	f->mmass       = 1e40;
	f->minweight   = 1e40;
	f->maxweight   = 0.0;
	f->sumweight   = 0.0;
	f->no_species  = INT32_C(0);
	f->minpos[0]   = f->minpos[1] = f->minpos[2] = 1e40;
	f->maxpos[0]   = f->maxpos[1] = f->maxpos[2] = -1e40;
	f->posscale    = 1.0;
	f->weightscale = 1.0;
  
	return f;
}

extern io_mgadget_t
io_mgadget_open(io_logging_t log,
                char *fname,
                io_file_swap_t swapped,
                io_file_mode_t mode,
                uint32_t reader)
{
	int32_t i;
	io_mgadget_t f;
	char **fnames;
#ifdef CHECK_RLIMIT_NOFILE
  struct rlimit limit;
#endif
#if (defined WITH_MPI && (defined BCASTHEADER || defined NCPUREADING_EQ_NFILES))
  int32_t *slen;
#endif

	/* XXX THIS IS CURRENTLY ONLY FOR READING! */
	if (mode != IO_FILE_READ)
		return NULL;

	/* Get memory for the structure */
	f = (io_mgadget_t)malloc(sizeof(io_mgadget_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_mgdaget_open(): sizeof(io_mgadget_struct_t)");
		return NULL;
	}

	/* Start filling the structure */

	/* Okay, we are a Multiple Gadget file */
	f->ftype = IO_FILE_MGADGET;

	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	MPI_Comm_split(MPI_COMM_WORLD, 1, f->rank, &(f->mycomm));
	MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	MPI_Comm_size(f->mycomm, &(f->size_mycomm));
#	endif

	/* Split the filename in path and stem */
	f->path = NULL;
	f->stem = NULL;
	if (  io_util_split_pathfname(fname, &(f->path), &(f->stem)) == NULL) {
		io_logging_fatal(log, "io_mgdaget_open(): Could not split %s in path and filename.", fname);
    if(f->path) free(f->path);
    if(f->stem) free(f->stem);
		free(f);
		return NULL;
	}
	io_logging_msg(log, INT32_C(1), "Will look in %s for %s", f->path, f->stem);

#if (defined WITH_MPI && (defined BCASTHEADER || defined NCPUREADING_EQ_NFILES)) // even for NCPUREADING_EQ_NFILES only the master determines 'numfiles'
  if(f->rank_mycomm == 0)
#endif
   {
    /* Get the filenames */
    f->numfiles = io_util_findfiles(f->path, f->stem, "%i", "", &fnames);
    
    if (f->numfiles <= 0) {
      io_logging_fatal(log, "io_mgdaget_open(): Could not open anything starting with %s in %s.", f->stem, f->path);
      if(f->path)  free(f->stem);
      if(f->stem)  free(f->path);
      free(f);
      return NULL;
    }
   }
  
  
#if (defined WITH_MPI && (defined BCASTHEADER || defined NCPUREADING_EQ_FILES))
  /* first round of broadcasts from root=0 (needed to array allocations) */
  MPI_Bcast(&(f->numfiles), sizeof(int32_t), MPI_BYTE, 0, f->mycomm);

  /* determine the length of each string from root MPI task */
  slen = (int32_t *) malloc(f->numfiles*sizeof(int32_t));
  if(f->rank_mycomm == 0) {
    for(i=0; i<f->numfiles; i++) {
      slen[i] = strlen(fnames[i]) + 1; // incl. the terminating character
    }
  }
  MPI_Bcast(slen, f->numfiles*sizeof(int32_t), MPI_BYTE, 0, f->mycomm);
  
  /* allocate memory for fnames[][] on all MPI tasks but the root one */
  if(f->rank_mycomm != 0) {
    fnames = (char **)malloc(f->numfiles*sizeof(char *));
    for(i=0; i<f->numfiles; i++) {
      fnames[i] = (char *)malloc(slen[i]*sizeof(char));
    }
  }
  
  /* second round of broadcasts from root=0 (actual filenames) */
  for(i=0; i<f->numfiles; i++) {
    MPI_Bcast(fnames[i], slen[i]*sizeof(char), MPI_BYTE, 0, f->mycomm);
  }  
#endif /* (defined WITH_MPI && (defined BCASTHEADER || defined NCPUREADING_EQ_FILES)) */
  
  
  
#if (defined WITH_MPI && defined NCPUREADING_EQ_NFILES)
  if(f->numfiles != global_io.params->reader) {
		io_logging_fatal(log, "io_mgdaget_open(): You are using NCPUREADING_EQ_NFILES and hence the number of files %d should be the number of MPI tasks reading %d (size_mycomm=%d).",
                     f->numfiles, global_io.params->reader, f->size_mycomm);
		if(f->path)  free(f->stem);
		if(f->stem)  free(f->path);
		free(f);    
    return NULL;
  }
#endif /* (defined WITH_MPI && defined NCPUREADING_EQ_NFILES) */

  
  
#ifdef DEBUG_MGADGET
  io_logging_msg(log,6,"\nio_mgdaget_open(): found the following %d files:\n",f->numfiles);
  for (i=0; i<f->numfiles; i++){
    io_logging_msg(log,6,"%s",fnames[i]);
  }
#endif

  
  
  
  
  
  
  
  
#ifdef CHECK_RLIMIT_NOFILE
  /* check whether we are able to use that many readers */
  if(getrlimit(RLIMIT_NOFILE, &limit) == 0)
   {
#ifdef WITH_MPI
    fprintf(stderr,"\nPRESENT FILE DESCRIPTOR LIMITS (rank=%d of %d):\n",f->rank,f->size);
#else
    fprintf(stderr,"\nPRESENT FILE DESCRIPTOR LIMITS :\n");
#endif
    fprintf(stderr,"   rlim_cur = %llu\n",limit.rlim_cur);
    fprintf(stderr,"   rlim_max = %llu\n",limit.rlim_max);
   }
  else
   {
#ifdef WITH_MPI
    fprintf(stderr,"getrlimit() failed with errno=%d (rank %d of %d)\n", errno,f->rank,f->size);
#else
    fprintf(stderr,"getrlimit() failed with errno=%d\n", errno);
#endif
    exit(1);
   }
  if(limit.rlim_max < f->numfiles+NFILES_EXTRA)
   {
#ifdef WITH_MPI
    fprintf(stderr,"\nPROBLEM (rank=%d of %d):\n",f->rank,f->size);
#else
    fprintf(stderr,"\nPROBLEM:\n");
#endif
    fprintf(stderr,"   only the super-user can increase beyond the maximum limit:\n");
    fprintf(stderr,"     numfiles = %d\n",f->numfiles);
    fprintf(stderr,"     rlim_max = %llu (incl. 3 for stdin, stdout, and stderr)\n",limit.rlim_max);
    fprintf(stderr,"   ABORTING\n");
    exit(1);
   }
  if (limit.rlim_cur < f->numfiles+NFILES_EXTRA)
   {
#ifdef WITH_MPI
    fprintf(stderr,"\nincreasing number of file descriptors on rank %d (of %d) to %d\n",f->rank,f->size,f->numfiles+NFILES_EXTRA);
#else
    fprintf(stderr,"\nincreasing number of file descriptors to %d\n",f->numfiles+NFILES_EXTRA);
#endif
    limit.rlim_cur = (rlim_t) (f->numfiles+NFILES_EXTRA);
    if (setrlimit(RLIMIT_NOFILE, &limit) != 0)
     {
#ifdef WITH_MPI
      fprintf(stderr,"setrlimit() failed with errno=%d (rank %d of %d)\n", errno,f->rank,f->size);
#else
      fprintf(stderr,"setrlimit() failed with errno=%d\n", errno);
#endif
      exit(1);
     }
    
    // check whether it worked or not
    if(getrlimit(RLIMIT_NOFILE, &limit) == 0)
     {
#ifdef WITH_MPI
      fprintf(stderr,"\nNEW LIMITS (rank=%d of %d):\n",f->rank,f->size);
#else
      fprintf(stderr,"\nNEW LIMITS:\n");
#endif
      fprintf(stderr,"   rlim_cur = %llu\n",limit.rlim_cur);
      fprintf(stderr,"   rlim_max = %llu\n",limit.rlim_max);
     }
    else
     {
#ifdef WITH_MPI
      fprintf(stderr,"getrlimit() failed with errno=%d (rank %d of %d)\n", errno,f->rank,f->size);
#else
      fprintf(stderr,"getrlimit() failed with errno=%d\n", errno);
#endif
      exit(1);
     }
   }
#endif /* CHECK_RLIMIT_NOFILE */
  
  
  
  
  
  
  
  
  
  
  
  
	/* Glue the files into the MGadget structure */
	f->files = (io_gadget_t *)malloc( sizeof(io_gadget_t)*(f->numfiles));
	if (f->files == NULL) {
    io_logging_memfatal(log,  "io_mgadget structure (2)");
		for (i=0; i<f->numfiles; i++)
			free(fnames[i]);
		free(fnames);
		free(f->stem);
		free(f->path);
		free(f);
		return NULL;
	}
  
  
#if (defined WITH_MPI && defined NCPUREADING_EQ_NFILES)

  /* generate all the file structures and copy file names into them */
  if(f->rank_mycomm != 0) {
    for (i=0; i<f->numfiles; i++) {
      ((f->files)[i]) = (io_gadget_t)calloc(sizeof(io_gadget_struct_t),1);
      ((f->files)[i])->fname = (char *)malloc(sizeof(char) * (strlen(fnames[i]) + 1));
      strncpy(((f->files)[i])->fname, fnames[i], strlen(fnames[i])+1);
    }
  }

  /* loop over all filenames... */
  for (i=0; i<f->numfiles; i++) {
    /* ...but only rank==i will open its own file */
    if(i == f->rank_mycomm) {
#ifdef DEBUG_MGADGET
      io_logging_msg(log,6,"\nrank %d trying to open file %s (%d of %d) for reading ... ",f->rank_mycomm,fnames[i],i+1,f->numfiles);
#endif
      
#		ifdef WITH_MPI
      /* TODO
       * THIS	IS JUST A NASTY HACK TO PREVENT io_gadget.c FROM REDOING
       * THE MPI-SPLIT. CALLED FROM HERE IT IS ONLY SUPPOSED TO ACT AS
       * A DUMMY INTERFACE TO A GADGET FILE.
       * TODO
       */
      (f->files)[i] = io_gadget_open(log, fnames[i], swapped, mode, f->size+1);
#		else
      (f->files)[i] = io_gadget_open(log, fnames[i], swapped, mode, reader);
#		endif
      
      if ((f->files)[i] == NULL) {
        int32_t j;
        for (j=i; i<f->numfiles; j++)
          free(fnames[j]);
        free(fnames);
        while (i>0) {
          i--;
          io_gadget_close(log, &((f->files)[i]));
        }
        free(f->stem);
        free(f->path);
        free(f);
        return NULL;
      }
      
#ifdef FOPENCLOSE
      //fprintf(stderr,"FOPENCLOSE: opened %s ... ",((f->files)[i])->fname);
      fclose(((f->files)[i])->file);
      // DO NOT SET THE POINTER TO NULL AS file!=NULL INDICATES FOR FOLLOWING ROUTINES THAT FILE COULD BE OPENED
      //fprintf(stderr,"and closed file (temporarily)\n");
#endif
     }
  }
  
  
#else /* (defined WITH_MPI && defined NCPUREADING_EQ_NFILES) */
  
  
#if (defined WITH_MPI && defined BCASTHEADER)
  /* only rank==0 will open all those files */
  if(f->rank_mycomm == 0) {
#endif
    for (i=0; i<f->numfiles; i++)
     {
#ifdef DEBUG_MGADGET
      io_logging_msg(log,6,"\ntrying to open file %s (%d of %d) for reading ... ", fnames[i],i+1,f->numfiles);
#endif
      
#		ifdef WITH_MPI
      /* TODO
       * THIS	IS JUST A NASTY HACK TO PREVENT io_gadget.c FROM REDOING
       * THE MPI-SPLIT. CALLED FROM HERE IT IS ONLY SUPPOSED TO ACT AS
       * A DUMMY INTERFACE TO A GADGET FILE.
       * TODO
       */
      (f->files)[i] = io_gadget_open(log, fnames[i], swapped, mode, f->size+1);
#		else
      (f->files)[i] = io_gadget_open(log, fnames[i], swapped, mode, reader);
#		endif
      
      if ((f->files)[i] == NULL) {
        int32_t j;
        for (j=i; i<f->numfiles; j++)
          free(fnames[j]);
        free(fnames);
        while (i>0) {
          i--;
          io_gadget_close(log, &((f->files)[i]));
        }
        free(f->stem);
        free(f->path);
        free(f);
        return NULL;
      }
      
#ifdef FOPENCLOSE
      //fprintf(stderr,"FOPENCLOSE: opened %s ... ",((f->files)[i])->fname);
      fclose(((f->files)[i])->file);
      // DO NOT SET THE POINTER TO NULL AS file!=NULL INDICATES FOR FOLLOWING ROUTINES THAT FILE COULD BE OPENED
      //fprintf(stderr,"and closed file (temporarily)\n");
#endif
     }
#if (defined WITH_MPI && defined BCASTHEADER)
  }
  /* we nevertheless need to generate all the file structures ready to receive the broadcast */
  else {
    for (i=0; i<f->numfiles; i++) {
      ((f->files)[i]) = (io_gadget_t)calloc(sizeof(io_gadget_struct_t),1);
      ((f->files)[i])->fname = NULL;
    }
  }
#endif /* (defined WITH_MPI && defined BCASTHEADER) */
  
#endif /* (defined WITH_MPI && defined NCPUREADING_EQ_NFILES) */
 

  
  
  
  
  // synchronisation of all headers across MPI tasks:
  //==================================================
  
#if (defined WITH_MPI && defined BCASTHEADER)
  // BROADCAST (f->files)[i] values to all other MPI tasks, root=0!
  for (i=0; i<f->numfiles; i++) {
    MPI_Bcast((f->files)[i], sizeof(io_gadget_struct_t), MPI_BYTE, 0, f->mycomm);

    /* only rank == 0 properly allocated the memory for the filename and copied the data there */
    if(f->rank_mycomm != 0) {
      ((f->files)[i])->fname = (char *)malloc(sizeof(char) * (strlen(fnames[i]) + 1));
      strncpy(((f->files)[i])->fname, fnames[i], strlen(fnames[i])+1);
    }
  }
#endif /* (defined WITH_MPI && defined BCASTHEADER) */
    
#if (defined WITH_MPI && defined NCPUREADING_EQ_NFILES)
  // BROADCAST (f->files)[i] values to all other MPI tasks, root=i!
  for (i=0; i<f->numfiles; i++) {
    MPI_Bcast((f->files)[i], sizeof(io_gadget_struct_t), MPI_BYTE, i, f->mycomm);
    
    /* only rank == i properly allocated the memory for the filename and copied the data there */
    if(f->rank_mycomm != i) {
      ((f->files)[i])->fname = (char *)malloc(sizeof(char) * (strlen(fnames[i]) + 1));
      strncpy(((f->files)[i])->fname, fnames[i], strlen(fnames[i])+1);
    }
  }
#endif /* (defined WITH_MPI && defined NCPUREADING_EQ_NFILES) */
  
  
  
  
  
  
  
  
  
  
  
  
  /* the whole char *fnames[] is not needed anymore */
  for (i=0; i<f->numfiles; i++) {
    free(fnames[i]);
  }
	free(fnames);
  
	/* Set initial values */
	f->no_part     = UINT64_C(0);
	f->multimass   = INT8_C(0);
	f->mmass       = 1e40;
	f->minweight   = 1e40;
	f->maxweight   = 0.0;
	f->sumweight   = 0.0;
	f->no_species  = INT32_C(0);
	f->minpos[0]   = f->minpos[1] = f->minpos[2] = 1e40;
	f->maxpos[0]   = f->maxpos[1] = f->maxpos[2] = -1e40;
	f->posscale    = 1.0;
	f->weightscale = 1.0;

	return f;
}

extern void
io_mgadget_close(io_logging_t log,
                 io_mgadget_t *f)
{
	int32_t i;

	/* Catch NULLs */
	if (f == NULL || *f == NULL)
		return;

	/* Put header to the file if necessary */
	/* XXX Only relevant for writing */

	/* free() all related structures and fclose() the files (in case the file pointer is not NULL) */
	for (i=0; i<(*f)->numfiles; i++)
		io_gadget_close(log, &(((*f)->files)[i]));
  
	if((*f)->files) free((*f)->files);
  if((*f)->stem)  free((*f)->stem);
  if((*f)->path)  free((*f)->path);
	if(*f)          free(*f);
	*f = NULL;

	return;
}

extern void
io_mgadget_init(io_logging_t log,
                io_mgadget_t f)
{
	int32_t i;

	if (f == NULL)
		return;

	if (f->files[0]->mode != IO_FILE_READ) {
		io_logging_warn(log, INT32_C(1),
		                "%s (first file of %" PRIi32 ") is not opened "
		                "for reading. Will do nothing.",
		                f->files[0]->fname, f->numfiles);
		return;
	}

	/* Check for multimass and sum up the particle count */
	f->no_part = UINT64_C(0);
	f->multimass = INT8_C(0);
  
#if (defined WITH_MPI && defined BCASTHEADER)
  if(f->rank_mycomm == 0) {
#endif
    for (i=0; i<f->numfiles; i++) {
      io_gadget_init(log, (f->files)[i]);
      f->multimass |= f->files[i]->multimass;
      f->no_part += f->files[i]->no_part;
    }
#if (defined WITH_MPI && defined BCASTHEADER)
  }
  else {
    /* allocate memory on this MPI task to hold the header for each file */
    for (i=0; i<f->numfiles; i++) {
      ((f->files)[i])->header = (io_gadget_header_t)malloc((size_t)GADGET_HEADER_SIZE+1);
    }
  }
    
  /* Broadcast the header of each file from rank==0 to all other MPI tasks */
  for (i=0; i<f->numfiles; i++) {
    MPI_Bcast(((f->files)[i])->header,((size_t)GADGET_HEADER_SIZE+1),MPI_BYTE,0,f->mycomm);

    /* io_gadget_init() also set some more relevant values in f */
    MPI_Bcast(&(((f->files)[i])->multimass),           sizeof(int8_t),   MPI_BYTE, 0, f->mycomm);
    MPI_Bcast(&(((f->files)[i])->no_part),             sizeof(uint64_t), MPI_BYTE, 0, f->mycomm);
    MPI_Bcast(&(((f->files)[i])->no_part_with_mass),   sizeof(uint64_t), MPI_BYTE, 0, f->mycomm);
    MPI_Bcast(&(f->multimass),                         sizeof(int8_t),   MPI_BYTE, 0, f->mycomm);
    MPI_Bcast(&(f->no_part),                           sizeof(uint64_t), MPI_BYTE, 0, f->mycomm);
  }
  
#endif

	return;
}

extern uint64_t
io_mgadget_readpart(io_logging_t log,
                    io_mgadget_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg)
{
	uint64_t part_read, tmp;
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight;

	/*
	 * First read the particles unscaled.
	 */
  f->minpos[0] = 1e40;
  f->minpos[1] = 1e40;
  f->minpos[2] = 1e40;
  f->maxpos[0] = -1e40;
  f->maxpos[1] = -1e40;
  f->maxpos[2] = -1e40;

	part_read = io_mgadget_readpart_raw(log, f, pskip, pread, strg);
	if (part_read != pread) {
		return UINT64_C(0);
	}

	/* And do the scaling */
#ifdef WITH_MPI
	io_gadget_scale_global(log, f->mycomm,  f->maxpos, f->minpos, &(f->mmass));
#endif
	tmp = io_gadget_scale_particles(log, f->maxpos, f->minpos,
	                                &(f->files[0]->header->boxsize),
	                                f->files[0]->header->expansion,
	                                f->posscale, f->mmass,
	                                part_read, strg);

	if (tmp != part_read) {
		return tmp;
	}

	/* Wow, we are done! */
	return part_read;
}


extern uint64_t
io_mgadget_readpart_raw(io_logging_t log,
                        io_mgadget_t f,
                        uint64_t pskip,
                        uint64_t pread,
                        io_file_strg_struct_t strg)
{
	long tmp;
	uint64_t partread, pread_file, pread_done, partinfile;
	uint64_t pskip_file, pskip_done;
	bool something_to_read = false;
	int32_t i;

	/* See if there is anything to do */
	if ( (f == NULL) || (f->files == NULL) )
		return UINT64_C(0);

	/* Initialize accounting of skipping and reading */
	pskip_done = pread_done = UINT64_C(0);
  
	/* Read the particles from the different files */
	for (i=0; i<f->numfiles; i++) {
		/*
		 * First figure out how many particles are in the file (use the
		 * tmporary long variable for that and copy it over to the 64
		 * bit integer afterwards.
		 */
		if (io_gadget_get(log, f->files[i], IO_FILE_GET_NOPART, &tmp) != true) {
			io_logging_fatal(log, "Could not get number of particles from %s", f->files[i]->fname);
		}
		partinfile = (uint64_t)tmp;

		/* Then do some arithmetic and set the skipping and reading
		 * numbers correctly for the file 
		 */
		pskip_file = (pskip_done<pskip) ? pskip-pskip_done : UINT64_C(0);
		pread_file = (pread_done<pread) ? pread-pread_done : UINT64_C(0);
		if (pskip_file > partinfile) {
			pskip_file = partinfile;
			pread_file = UINT64_C(0);
		} else {
			if (pread_file > partinfile - pskip_file)
				pread_file = partinfile - pskip_file;
		}

#if (defined WITH_MPI && defined NCPUREADING_EQ_NFILES)
    
    /* only read f->files[rank] */
    if(i == f->rank_mycomm) {
      pread_file = partinfile;
      pskip_file = UINT64_C(0);

      /* Now read the particles */
      partread = io_gadget_readpart_raw(log, f->files[i], pskip_file, pread_file, strg);
      
      if (pread_file != partread) {
        io_logging_fatal(log, "Something went wrong. Wanted to read %" PRIu64 " particles, but got %" PRIu64 ". Aborting.",
                         pread_file, partread);
        return pread_done + partread;
      }
      /* Update the skipping arithmetic */
      pskip_done += pskip_file;
      pread_done += pread_file;
      
      /* Move the particle pointers */
      strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride*partread);
      strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride*partread);
      strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride*partread);
      strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride*partread);
      strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride*partread);
      strg.momz.val = (void *)(((char *)strg.momz.val)  + strg.momz.stride*partread);
      if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride*partread);
      if (strg.id.val     != NULL) strg.id.val     = (void *)(((char *)strg.id.val)     + strg.id.stride    *partread);
      if (strg.u.val      != NULL) strg.u.val      = (void *)(((char *)strg.u.val)      + strg.u.stride     *partread);
#ifdef METALHACK
      if (strg.z.val      != NULL) strg.z.val      = (void *)(((char *)strg.z.val)      + strg.z.stride     *partread);
      if (strg.age.val    != NULL) strg.age.val    = (void *)(((char *)strg.age.val)    + strg.age.stride   *partread);
#endif
      
      /* Update the important things needed for scaling */
      f->sumweight = f->files[i]->sumweight;  // no accumulation of sumweight as we do an MPI_SUM reduction below
      if (isless(f->files[i]->mmass, f->mmass))       		 f->mmass     = f->files[i]->mmass;
      if (isless(f->files[i]->minweight, f->minweight))		 f->minweight = f->files[i]->minweight;
      if (isgreater(f->files[i]->maxweight, f->maxweight)) f->maxweight = f->files[i]->maxweight;
      if (isgreater(f->files[i]->maxpos[0], f->maxpos[0])) f->maxpos[0] = f->files[i]->maxpos[0];
      if (isgreater(f->files[i]->maxpos[1], f->maxpos[1])) f->maxpos[1] = f->files[i]->maxpos[1];
      if (isgreater(f->files[i]->maxpos[2], f->maxpos[2])) f->maxpos[2] = f->files[i]->maxpos[2];
      if (isless(f->files[i]->minpos[0], f->minpos[0]))    f->minpos[0] = f->files[i]->minpos[0];
      if (isless(f->files[i]->minpos[1], f->minpos[1]))    f->minpos[1] = f->files[i]->minpos[1];
      if (isless(f->files[i]->minpos[2], f->minpos[2]))    f->minpos[2] = f->files[i]->minpos[2];
      /* TODO Something indeed needs be done with no_species... */
    }
    
#else /* (defined WITH_MPI && defined NCPUREADING_EQ_NFILES) */
    
		/* Now read the particles */
		partread = io_gadget_readpart_raw(log, f->files[i], pskip_file, pread_file, strg);
		if (pread_file != partread) {
			io_logging_fatal(log, "Something went wrong. Wanted to read %" PRIu64 " particles, but got %" PRIu64 ". Aborting.",
                       pread_file, partread);
			return pread_done + partread;
		}

		/* Update the skipping arithmetic */
		pskip_done += pskip_file;
		pread_done += pread_file;

		/* Move the particle pointers */
		strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride*partread);
		strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride*partread);
		strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride*partread);
		strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride*partread);
		strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride*partread);
		strg.momz.val = (void *)(((char *)strg.momz.val)  + strg.momz.stride*partread);
		if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride*partread);
		if (strg.id.val != NULL)     strg.id.val = (void *)(((char *)strg.id.val) + strg.id.stride*partread);
		if (strg.u.val != NULL)      strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride*partread);
#ifdef METALHACK
		if (strg.z.val != NULL)      strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride*partread);
		if (strg.age.val != NULL)    strg.age.val = (void *)(((char *)strg.age.val) + strg.age.stride*partread);
#endif

		/* Update the important things needed for scaling */
		f->sumweight += f->files[i]->sumweight;
		if (isless(f->files[i]->mmass, f->mmass))       		 f->mmass = f->files[i]->mmass;
		if (isless(f->files[i]->minweight, f->minweight))		 f->minweight = f->files[i]->minweight;
		if (isgreater(f->files[i]->maxweight, f->maxweight)) f->maxweight = f->files[i]->maxweight;
		if (isgreater(f->files[i]->maxpos[0], f->maxpos[0])) f->maxpos[0] = f->files[i]->maxpos[0];
		if (isgreater(f->files[i]->maxpos[1], f->maxpos[1])) f->maxpos[1] = f->files[i]->maxpos[1];
		if (isgreater(f->files[i]->maxpos[2], f->maxpos[2])) f->maxpos[2] = f->files[i]->maxpos[2];
		if (isless(f->files[i]->minpos[0], f->minpos[0]))    f->minpos[0] = f->files[i]->minpos[0];
		if (isless(f->files[i]->minpos[1], f->minpos[1]))    f->minpos[1] = f->files[i]->minpos[1];
		if (isless(f->files[i]->minpos[2], f->minpos[2]))    f->minpos[2] = f->files[i]->minpos[2];
		/* TODO Something indeed needs be done with no_species... */

#endif /* (defined WITH_MPI && defined NCPUREADING_EQ_NFILES) */
	} /* f->numfiles */

#if (defined WITH_MPI && defined NCPUREADING_EQ_NFILES)
 {
  /* Reduce all values */
  double buffer[3];
  
  MPI_Allreduce(&(f->mmass),     buffer, 1, MPI_DOUBLE, MPI_MIN, f->mycomm);
  f->mmass = buffer[0];
  
  MPI_Allreduce(&(f->minweight), buffer, 1, MPI_DOUBLE, MPI_MIN, f->mycomm);
  f->minweight = buffer[0];
  
  MPI_Allreduce(&(f->maxweight), buffer, 1, MPI_DOUBLE, MPI_MAX, f->mycomm);
  f->maxweight = buffer[0];
  
  MPI_Allreduce(&(f->sumweight), buffer, 1, MPI_DOUBLE, MPI_SUM, f->mycomm);
  f->sumweight = buffer[0];

  /* note, this reduction is in fact repeated by io_gadget_scale_global() but what the heck... */
  MPI_Allreduce(&(f->maxpos),    buffer, 3, MPI_DOUBLE, MPI_MAX, f->mycomm);
  f->maxpos[0] = buffer[0];
  f->maxpos[1] = buffer[1];
  f->maxpos[2] = buffer[2];

  MPI_Allreduce(&(f->minpos),    buffer, 3, MPI_DOUBLE, MPI_MIN, f->mycomm);
  f->minpos[0] = buffer[0];
  f->minpos[1] = buffer[1];
  f->minpos[2] = buffer[2];
  
  /* indicate that all files are ready for closure */
	for (i=0; i<f->numfiles; i++) {
    (f->files[i])->file = NULL;
  }
}
#endif
  
	return pread_done;
}

extern bool
io_mgadget_get(io_logging_t log,
               io_mgadget_t f,
               io_file_get_t what,
               void *res)
{
	if ( (f == NULL) || (f->files[0]->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_NOPART_IN_FILE:
		case IO_FILE_GET_NOPART:
			*((long *)res) = (long)(f->no_part);
			break;
		case IO_FILE_GET_NOVPART:
			if (f->no_part > UINT64_C(0)) {
				if (    isgreater(f->sumweight, 0.0)
				     && isgreater(f->mmass, 0.0))
					*((double *)res) = f->sumweight / f->mmass;
				else
					*((double *)res) = (double)f->no_part;
			} else {
					io_logging_warn(log, INT32_C(0),
					                "Cannot calculate novpart yet. "
					                "You first need to read the "
					                "particles.");
					return false;
			}
			break;
		case IO_FILE_GET_NOSPECIES:
			if (isgreater(f->files[0]->header->massarr[1], 0.0))
				*((int *)res) = 1;
			else {
				*((int *)res) = 1;
				io_logging_warn(log, INT32_C(1),
				                "Not implemented yet.");
				return false;
			}
			break;
		case IO_FILE_GET_BOXSIZE:
			*((double *)res) =   f->files[0]->header->boxsize
			                   * f->posscale;
			break;
		case IO_FILE_GET_PMASS:
			*((double *)res) =   f->mmass * f->weightscale;
			break;
		case IO_FILE_GET_ZINITIAL:
			io_logging_warn(log, INT32_C(1),
			                "zinitial is not set in a Gadget file, "
			                "using current redshift");
		case IO_FILE_GET_Z:
			*((double *)res) = f->files[0]->header->redshift;
			break;
		case IO_FILE_GET_AINITIAL:
			io_logging_warn(log, INT32_C(1),
			                "ainitial is not set in a Gadget file, "
			                "using current expansion");
		case IO_FILE_GET_A:
			*((double *)res) = f->files[0]->header->expansion;
			break;
		case IO_FILE_GET_OMEGA0:
			*((double *)res) = f->files[0]->header->omega0;
			break;
		case IO_FILE_GET_OMEGAL:
			*((double *)res) = f->files[0]->header->omegalambda;
			break;
		case IO_FILE_GET_H:
			*((double *)res) = f->files[0]->header->hubbleparameter;
			break;
		case IO_FILE_GET_DOUBLE:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't store the use of "
			                "double precision. Assuming it is not "
			                "double precision.");
			*((int *)res) = 0;
			break;
		case IO_FILE_GET_MMASS:
			if (isgreater(f->files[0]->header->massarr[1], 0.0))
				*((int *)res) = 0;
			else
				*((int *)res) = 1;
			break;
		case IO_FILE_GET_NOTSTEP:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't store the step number. "
			                "Setting to 0.");
			*((int32_t *)res) = 0;
			break;
		case IO_FILE_GET_TSTEP:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't store the timestep. "
			                "Setting to 0.0");
			*((double *)res) = 0.0;
			break;
		case IO_FILE_GET_HEADERSTR:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't have a header string. "
			                "Using a dummy one.");
			*((char **)res) = "No header string.";
			break;
		case IO_FILE_GET_MINWEIGHT:
			if (isgreater(f->mmass, 0.0))
				*((double *)res) = f->minweight / f->mmass;
			else {
				io_logging_warn(log, INT32_C(1),
				                "Don't know minweight yet, setting to "
				                "0.0.");
				*((double *)res) = 0.0;
			}
			break;
		case IO_FILE_GET_MAXWEIGHT:
			if (isgreater(f->mmass, 0.0))
				*((double *)res) = f->maxweight / f->mmass; 
			else {
				io_logging_warn(log, INT32_C(1),
				                "Don't know maxweight yet, setting to "
				                "0.0.");
				*((double *)res) = 0.0;
			}
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern bool
io_mgadget_set(io_logging_t log,
               io_mgadget_t f,
               io_file_get_t what,
               void *res)
{
	int32_t i;

	if (f == NULL)
		return false;
	for (i=0; i<f->numfiles; i++)
		if (f->files[i]->header == NULL)
			return false;

	switch (what) {
		case IO_FILE_GET_BOXSIZE:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->boxsize = *((double *)res);
			}
			break;
		case IO_FILE_GET_PMASS:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->massarr[1] = *((double *)res);
			}
			break;
		case IO_FILE_GET_Z:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->redshift = *((double *)res);
			}
			break;
		case IO_FILE_GET_A:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->expansion = *((double *)res);
			}
			break;
		case IO_FILE_GET_OMEGA0:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->omega0 = *((double *)res);
			}
			break;
		case IO_FILE_GET_OMEGAL:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->omegalambda = *((double *)res);
			}
			break;
		case IO_FILE_GET_H:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->hubbleparameter = *((double *)res);
			}
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern void
io_mgadget_log(io_logging_t log, io_mgadget_t f)
{
	int32_t i;

	io_logging_msg(log, INT32_C(5),
	               "Fileobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Filetype:            %s",
	               io_file_typestr(f->ftype));
	io_logging_msg(log, INT32_C(5),
	               "  Path:                %s",
	               f->path);
	io_logging_msg(log, INT32_C(5),
	               "  Stem:                %s",
	               f->stem);
	io_logging_msg(log, INT32_C(5),
	               "  Number of files:     %" PRIi32,
	               f->numfiles);
	io_logging_msg(log, INT32_C(5),
	               "  Number of particles: %" PRIu64,
	               f->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass         :  %" PRIi8,
	               f->multimass);
	io_logging_msg(log, INT32_C(5),
	               "  Mmass             :  %g",
	               f->mmass);
	io_logging_msg(log, INT32_C(5),
	               "  Minimal Weight    :  %g",
	               f->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  Maximal Weight    :  %g",
	               f->maxweight);
	io_logging_msg(log, INT32_C(5),
	               "  Sum of weights    :  %g",
	               f->sumweight);

  
  for (i=0; i<f->numfiles; i++) {
		io_logging_msg(log, INT32_C(5),
		               "  ---> File %"PRIi32 ":",
		               i);
		io_gadget_log(log, (f->files)[i]);
	}

	return;
}

extern void
io_mgadget_resetscale(io_logging_t log,
                      io_mgadget_t f,
                      double posscale,
                      double weightscale) {
	int32_t i;

	if (f == NULL)
		return;

	for (i=0; i<f->numfiles; i++) {
		io_gadget_resetscale(log, f->files[i], posscale, weightscale);
	}
	f->posscale = posscale;
	f->weightscale = weightscale;

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
