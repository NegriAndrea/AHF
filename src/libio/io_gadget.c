/* $Id: io_gadget.c,v 1.31 2008/07/31 08:22:59 knolli Exp $ */

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
#ifdef WITH_MPI
#	include <mpi.h>
#endif
#include <limits.h>

#include "io_gadget.h"
#include "io_gadget_header.h"
#include "io_util.h"

#include "../define.h"
#include "../common.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/

/**
 * \brief Helper function to open the file
 *
 * This function is supposed to get inlined anyway. Makes
 * io_gadget_open more readable.
 *
 * \param log   The logging object.
 * \param f     The Gadget file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_gadget_t
local_openopen(io_logging_t log, io_gadget_t f, io_file_mode_t mode);

 /**
  * \brief Helper funtion to set the swap-state and write some
  *        messages
  *
  * \param log      The logging object.
  * \param f        The Gadget file object sofar.
  * \param swapped  The swap state.
  *
  * \return Nothing.
  */
inline static void
local_openswapped(io_logging_t log,
                  io_gadget_t f,
                  io_file_swap_t swapped);

/**
 * \brief Try to find out swapping status.
 *
 * \param log  The logging object.
 * \param f    The file object.
  *
  * \return Returns the file object or NULL in case of an error.
 */
inline static io_gadget_t
local_opengetswap(io_logging_t log, io_gadget_t f);

/**
 * \brief Tries to figure out which Gadget file version is to be used.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Nothing.
 */
inline static void
local_openversion(io_logging_t log, io_gadget_t f);

/**
 * \brief Check before writing particles and place file pointer at the
 *        right spot.
 *
 * \param log
 * \param f
 * \param pskip
 * \param bytes
 *
 * \return 
 */
inline static int32_t
local_write_common(io_logging_t log,
                   io_gadget_t f,
                   uint64_t pskip,
                   int32_t bytes);

static uint64_t
local_get_block_pos(io_logging_t log,
                    io_gadget_t f,
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg);

static uint64_t
local_get_block_vel(io_logging_t log,
                    io_gadget_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg);

static uint64_t
local_get_block_id(io_logging_t log,
                   io_gadget_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg);

static uint64_t
local_get_block_mass(io_logging_t log,
                     io_gadget_t f,
                     uint64_t pskip,
                     uint64_t pread,
                     io_file_strg_struct_t strg);

static uint64_t
local_get_block_u(io_logging_t log,
                  io_gadget_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg);

void local_find_block(io_gadget_t, char *);

#ifdef METALHACK
static uint64_t
local_get_block_z(io_logging_t log,
                  io_gadget_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg);

static uint64_t
local_get_block_age(io_logging_t log,
                    io_gadget_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg);

static void
local_skip_GADGET1_blocks(FILE *, int, int);
#endif



/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_gadget_t
io_gadget_open(io_logging_t log,
               char *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t reader)
{
	io_gadget_t f;

	/* Get memory for the structure */
	f = (io_gadget_t)malloc(sizeof(io_gadget_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_gadget structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Store the filename */
	f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
	if (f->fname == NULL) {
		io_logging_memfatal(log, "filename of GadgetFile");
		free(f);
		return NULL;
	}
	strncpy(f->fname, fname, strlen(fname)+1);

	/* Okay, we are a Gadget file */
	f->ftype = IO_FILE_GADGET;

	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	if (f->size >= reader) {
		/* TODO 
		 * THIS IS JUST A QUICK HACK TO PREVENT MGADGET FILES TO
		 * AGAIN TRY TO SPLIT THE COMMUNICATOR, THAT IS ALREADY DONE
		 * IN mgadget.c 
		 * TODO 
		 */
		MPI_Comm_split(MPI_COMM_WORLD, 1, f->rank, &(f->mycomm));
		MPI_Comm_size(f->mycomm, &(f->size_mycomm));
		MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	} else {
		f->mycomm = MPI_COMM_NULL;
		f->size_mycomm = -1;
		f->rank_mycomm = -1;
	}

#	endif

	/* Try to open the file and set mode */
	if (local_openopen(log, f, mode) == NULL) {
		free(f->fname);
		free(f);
		return NULL;
	}
    
	/* Set swapping */
	local_openswapped(log, f, swapped);
	if (    (f->mode == IO_FILE_READ)
	     && (f->swapped == IO_FILE_UNKOWN_SWAPPING)) {
		if (local_opengetswap(log, f) != f) {
			io_logging_fatal(log, "Cannot open this file.");
			free(f->fname);
			free(f);
			return NULL;
		}
	}

	/* Identify Gadget format */
	local_openversion(log, f);

	/* Nothing for the header for now */
	f->header = NULL;

	/* Set some dummy values */
	f->no_part           = UINT64_C(0);
	f->no_part_with_mass = UINT64_C(0);
	f->multimass         = INT8_C(0);
	f->mmass             = 1e40;
	f->minweight         = 1e40;
	f->maxweight         = 0.0;
	f->sumweight         = 0.0;
	f->no_species        = INT32_C(0);
	f->posscale          = 1.0;
	f->weightscale       = 1.0;

	return f;
}

extern void
io_gadget_close(io_logging_t log,
                io_gadget_t *f)
{
	/* Catch NULLs */
	if (f == NULL || *f == NULL)
		return;

	/* Put header to the file if necessary */
	if (    ((*f)->mode == IO_FILE_WRITE)
	     && ((*f)->header != NULL)) {
		io_gadget_header_write(log, (*f)->header, *f);
	}

	/* Close */
	if ((*f)->header != NULL)
		io_gadget_header_del(log, &((*f)->header));
	if ((*f)->fname != NULL)
		free((*f)->fname);
#	ifdef WITH_MPI
	if ((*f)->mycomm != MPI_COMM_NULL)
		MPI_Comm_free(&((*f)->mycomm));
#	endif

	/* Actually close the file */
	if ((*f)->file != NULL)
		fclose((*f)->file);

	/* Cleaning */
	free(*f);
	*f = NULL;

	return;
}

extern void
io_gadget_init(io_logging_t log,
               io_gadget_t f)
{
	if (f == NULL)
		return;

	if (f->header != NULL) {
		io_logging_warn(log, INT32_C(1),
		                "Already have the header information! Rereading.");
		io_gadget_header_del(log, &(f->header));
	}

	if (f->mode != IO_FILE_READ) {
		io_logging_warn(log, INT32_C(1),
		                "%s is not opened for reading. "
		                "Will do nothing.",
		                f->fname);
		return;
	}

	io_logging_msg(log, INT32_C(5),
	               "Starting to initialize file object from %s",
	               f->fname);
	f->header = io_gadget_header_get(log, f);
	io_logging_msg(log, INT32_C(5),
	               "Done with initializing file object from %s",
	               f->fname);

	/* Check for multimass file and also sum up the particle count */
	f->multimass = 0;
	{ 
		int i;
		for (i=0; i<6; i++) {
			f->no_part += (uint64_t)(f->header->np[i]);
			if (    islessequal(f->header->massarr[i], 0.0)
			     && (f->header->np[i]>0)) {
				f->multimass |= (1<<i);
				io_logging_msg(log, INT32_C(5),
				               "Particle type %d requires reading "
				               "of a mass array.", i);
				f->no_part_with_mass += (uint64_t)(f->header->np[i]);
			}
		}
	}

	return;
}

extern uint64_t
io_gadget_readpart(io_logging_t log,
                   io_gadget_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg)
{
	uint64_t particles_read, tmp;
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight;

	/* 
	 * First read the particles unscaled. This will set important
	 * scaling information
	 */
	particles_read = io_gadget_readpart_raw(log, f, pskip, pread, strg);

	if (particles_read != pread) {
		return UINT64_C(0);
	}

	/* And do the scaling */
#ifdef WITH_MPI
	io_gadget_scale_global(log, f->mycomm,  f->maxpos, f->minpos, &(f->mmass));
#endif
	tmp = io_gadget_scale_particles(log, f->maxpos, f->minpos,
	                                &(f->header->boxsize),
	                                f->header->expansion,
	                                f->posscale, f->mmass,
	                                particles_read, strg);
	if (tmp != particles_read) {
		return tmp;
	}

	/* Wow, we are done! */
	return particles_read;
}

/* No we define a bunch of macros to make life easier */
#define SKIP  {io_util_readuint32(f->file, &blocksize,  f->swapped);}
#define SKIP2 {io_util_readuint32(f->file, &blocksize2, f->swapped);}
#define CHECK_BLOCK {\
	if (blocksize != blocksize2) {\
		io_logging_fatal(log,\
		                 "The block boundaries (beginning: %" PRIu32\
		                 " end: %" PRIu32 ") are not identical. "\
		                 "Corrupt file?", blocksize, blocksize2);\
		fflush(NULL);\
		exit(999);\
		return UINT64_C(0);\
	} else {\
		io_logging_msg(log, INT32_C(5),\
		               "Block claimed correctly to be %f MB long.", \
		               (float)(blocksize/1024./1024.));\
	}\
}
#define VERIFY_BLOCK(name) {\
	if (f->ver == 2) {\
		SKIP;\
		io_util_readstring(f->file, str, (size_t)4);\
		io_util_readuint32(f->file, &nextblocksize, f->swapped);\
		if (strncmp(name, str, 4) != 0) {\
			io_logging_fatal(log,\
			                 "Wrong block! Expected: %s  Found: %s",\
			                 name, str);\
		} else {\
			io_logging_msg(log, INT32_C(1),\
			               "Arrived at block %s, size of it will be %" \
			               PRIu32, str, nextblocksize);\
		}\
		SKIP2;\
		CHECK_BLOCK;\
	}\
}
extern uint64_t
io_gadget_readpart_raw(io_logging_t log,
                       io_gadget_t f,
                       uint64_t pskip,
                       uint64_t pread,
                       io_file_strg_struct_t strg)
{
	long skipsize;
	char str[5];
	uint32_t blocksize, blocksize2;
	uint64_t funcrtn;
	uint32_t nextblocksize;
  
#ifdef FOPENCLOSE
  //fprintf(stderr,"FOPENCLOSE: io_gadget_readpart_raw() opening %s ... ",f->fname);
  f->file = fopen(f->fname,IO_FILE_MODE_READ);
  if(f->file == NULL) {
    io_logging_fatal(log,"io_gadget_readpart_raw(): could not open file %s for reading",f->fname);
    return UINT64_C(0);
   }
#endif
  
	/* Check if we actually have to do something */
	if ( (f == NULL) || (f->header == NULL) )
		return UINT64_C(0);

	/* Make sure that we are at the beginning of the file */
	rewind(f->file);

	/* Now go to the first block of stuff, skipping the header and in
	 * case of a Gadget-2 file, the descriptive block */
	skipsize = 2*sizeof(int)+GADGET_HEADER_SIZE;
	if (f->ver == 2)
		skipsize += (sizeof(int) * 3 + 4*sizeof(char));
	fseek(f->file, skipsize, SEEK_SET);

	/* Positions */
#ifdef GADGET_MAGNETICUM
  local_find_block(f, "POS ");
#else
	VERIFY_BLOCK("POS ");
#endif
	funcrtn = local_get_block_pos(log, f, &pskip, &pread, strg);
	if (funcrtn != pread) {
		io_logging_fatal(log, 
		                 "Expected to read %"PRIu64
		                 " particle positions, but only got %"PRIu64
		                 ".  Aborting.");
		return UINT64_C(0);
	}

	/* Velocities */
	VERIFY_BLOCK("VEL ");
	funcrtn = local_get_block_vel(log, f, pskip, pread, strg);
	if (funcrtn != pread) {
		io_logging_fatal(log, 
		                 "Expected to read %"PRIu64
		                 " particle velocities, but only got %"PRIu64
		                 ".  Aborting.");
		return UINT64_C(0);
	}

	/* Identities */
	VERIFY_BLOCK("ID  ");
	funcrtn = local_get_block_id(log, f, pskip, pread, strg);
	if (funcrtn != pread) {
		io_logging_fatal(log, 
		                 "Expected to read %"PRIu64
		                 " particle identities, but only got %"PRIu64
		                 ".  Aborting.");
		return UINT64_C(0);
	}

	/* Mass (might not exist!) */
	if (f->multimass != 0) {
		VERIFY_BLOCK("MASS");
	}
	funcrtn = local_get_block_mass(log, f, pskip, pread, strg);
	if (funcrtn != pread) {
		io_logging_fatal(log, 
		                 "Expected to obtain %"PRIu64
		                 " particle masses, but only got %"PRIu64
		                 ".  Aborting.");
		return UINT64_C(0);
	}

	/* Gas energies */
	if (f->header->np[0] > 0) {
		VERIFY_BLOCK("U   ");
	}
	funcrtn = local_get_block_u(log, f, pskip, pread, strg);
	if (funcrtn != pread) {
		io_logging_fatal(log, 
		                 "Expected to read %"PRIu64
		                 " particle identities, but only got %"PRIu64
		                 ".  Aborting.");
		return UINT64_C(0);
	}

#	ifdef METALHACK
	local_get_block_age(log, f, pskip, pread, strg);
	local_get_block_z(log, f, pskip, pread, strg);
#	endif

#ifdef FOPENCLOSE
  fclose(f->file);
  f->file = NULL; // we are not going to read any more from the file and hence can indicate to close it by setting f->file=NULL
  //fprintf(stderr,"and closed again (for good!)\n");
#endif

	/* Return the number of particles read */
	return pread;
}
/* Getting rid of the macros */
#undef VERIFY_BLOCK


extern uint64_t
io_gadget_writepart(io_logging_t log,
                    io_gadget_t f,
                    uint64_t pskip,
                    uint64_t pwrite,
                    io_file_strg_struct_t strg)
{
	return 0;
}

extern uint64_t
io_gadget_writepart_ord(io_logging_t log,
                        io_gadget_t f,
                        uint64_t pskip,
                        uint64_t pwrite,
                        void *nxt_part,
                        io_file_strg_struct_t strg)
{
	return UINT64_C(0);
}

extern bool
io_gadget_get(io_logging_t log,
              io_gadget_t f,
              io_file_get_t what,
              void *res)
{
	if ( (f == NULL) || (f->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_NOPART_IN_FILE:
		case IO_FILE_GET_NOPART:
			*((long *)res) = (long)f->no_part;
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
			if (f->multimass) {
				if (f->no_species >= 1)
					*((int *)res) = f->no_species;
				else {
					io_logging_warn(log, INT32_C(0),
					                "Don't know the number of particle "
					                "species yet. You first need to "
					                "read the particles.");
					return false;
				}
			} else {
				*((int *)res) = 1;
			}
			break;
		case IO_FILE_GET_BOXSIZE:
			*((double *)res) =   f->header->boxsize * f->posscale;
			break;
		case IO_FILE_GET_PMASS:
			if (f->multimass) {
				*((double *)res) = f->mmass * f->weightscale;
			} else {
				*((double *)res) =   f->header->massarr[1] * f->weightscale;
			}
			break;
		case IO_FILE_GET_ZINITIAL:
			io_logging_warn(log, INT32_C(1),
			                "zinitial is not set in a Gadget file, "
			                "using current redshift");
		case IO_FILE_GET_Z:
			*((double *)res) = f->header->redshift;
			break;
		case IO_FILE_GET_AINITIAL:
			io_logging_warn(log, INT32_C(1),
			                "ainitial is not set in a Gadget file, "
			                "using current expansion.");
		case IO_FILE_GET_A:
			*((double *)res) = f->header->expansion;
			break;
		case IO_FILE_GET_OMEGA0:
			*((double *)res) = f->header->omega0;
			break;
		case IO_FILE_GET_OMEGAL:
			*((double *)res) = f->header->omegalambda;
			break;
		case IO_FILE_GET_H:
			*((double *)res) = f->header->hubbleparameter;
			break;
		case IO_FILE_GET_DOUBLE:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't store the use of "
			                "double precision. Assuming it is not "
			                "double precision.");
			*((int *)res) = 0;
			break;
		case IO_FILE_GET_MMASS:
			if (isgreater(f->header->massarr[1], 0.0))
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
io_gadget_set(io_logging_t log,
              io_gadget_t f,
              io_file_get_t what,
              void *res)
{
	if ( (f == NULL) || (f->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_BOXSIZE:
			f->header->boxsize = *((double *)res);
			break;
		case IO_FILE_GET_PMASS:
			f->header->massarr[1] = *((double *)res);
			break;
		case IO_FILE_GET_Z:
			f->header->redshift = *((double *)res);
			break;
		case IO_FILE_GET_A:
			f->header->expansion = *((double *)res);
			break;
		case IO_FILE_GET_OMEGA0:
			f->header->omega0 = *((double *)res);
			break;
		case IO_FILE_GET_OMEGAL:
			f->header->omegalambda = *((double *)res);
			break;
		case IO_FILE_GET_H:
			f->header->hubbleparameter = *((double *)res);
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern void
io_gadget_log(io_logging_t log, io_gadget_t f)
{
	io_logging_msg(log, INT32_C(5),
	               "Fileobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Filetype:             %s",
	               io_file_typestr(f->ftype));
	io_logging_msg(log, INT32_C(5),
	               "  Filename:             %s",
	               f->fname);
	io_logging_msg(log, INT32_C(5),
	               "  Mode:                 %" PRIi8,
	               f->mode);
	io_logging_msg(log, INT32_C(5),
	               "  Swapping:             %" PRIi8,
	               f->swapped);
	io_logging_msg(log, INT32_C(5),
	               "  File version:         %" PRIi8,
	               f->ver);
	io_logging_msg(log, INT32_C(5),
	               "  Header size:          %" PRIi8,
	               f->ver);
	io_logging_msg(log, INT32_C(5),
	               "  No. particles:        %" PRIu64,
	               f->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  No. particles w/mass: %" PRIu64,
	               f->no_part_with_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass:            %" PRIi8,
	               f->multimass);
   {
      int i;
      for (i=0; i<6; i++)
         io_logging_msg(log, INT32_C(5),
                        "      Part type %d:      %" PRIi8, i,
                        ((int8_t)(f->multimass) & (1<<i)) >> i);
   }
	io_logging_msg(log, INT32_C(5),
	               "  MMass (Halo parts):   %g",
	               f->mmass);
	io_logging_msg(log, INT32_C(5),
	               "  Minimal Weight:       %g",
	               f->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  Maximal Weight:       %g",
	               f->maxweight);
	io_logging_msg(log, INT32_C(5),
	               "  Sum of all weights:   %g",
	               f->sumweight);
	io_logging_msg(log, INT32_C(5),
	               "  No. of species:       %" PRIi32,
	               f->no_species);
	io_logging_msg(log, INT32_C(5),
	               "  Position scale:       %g",
	               f->posscale);
	io_logging_msg(log, INT32_C(5),
	               "  Weight scale:         %g",
	               f->weightscale);
	io_gadget_header_log(log, f->header);

	return;
}

extern void
io_gadget_resetscale(io_logging_t log,
                     io_gadget_t f,
                     double posscale,
                     double weightscale) {
	if (f == NULL)
		return;

	io_logging_msg(log, INT32_C(8),
	               "Old posscale: %g   New posscale: %g",
	               f->posscale, posscale);
	io_logging_msg(log, INT32_C(8),
	               "Old weightscale: %g   New weightscale: %g",
	               f->weightscale, weightscale);
	f->posscale = posscale;
	f->weightscale = weightscale;

	return;
}

extern uint64_t
io_gadget_scale_particles(io_logging_t log,
                          double maxpos[],
                          double minpos[],
                          double *boxsize,
                          double expansion,
                          double posscale,
                          double mmass,
                          uint64_t particles_read,
                          io_file_strg_struct_t strg)
{
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight, scale_u;
	uint64_t i;

	/* Now we can do the scaling */
	box[0] = fabs(maxpos[0] - minpos[0]);
	box[1] = fabs(maxpos[1] - minpos[1]);
	box[2] = fabs(maxpos[2] - minpos[2]);
	if (isgreater(box[0], *boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "x-Separation of particles exceeds boxsize "
		                "(%g > %g), resetting boxsize.",
		                box[0], *boxsize);
		*boxsize = box[0];
	}
	if (isgreater(box[1], *boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "y-Separation of particles exceeds boxsize "
		                "(%g > %g), resetting boxsize.",
		                box[1], *boxsize);
		*boxsize = box[1];
	}
	if (isgreater(box[2], *boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "z-Separation of particles exceeds boxsize "
		                "(%g > %g), resetting boxsize.",
		                box[2], *boxsize);
		*boxsize = box[2];
	}
	io_logging_msg(log, INT32_C(4),
	               "Extreme positions: xmin = %g  xmax = %g",
	               minpos[0], maxpos[0]);
	io_logging_msg(log, INT32_C(4),
	               "                   ymin = %g  ymax = %g",
	               minpos[1], maxpos[1]);
	io_logging_msg(log, INT32_C(4),
	               "                   zmin = %g  zmax = %g",
	               minpos[2], maxpos[2]);
	shift[0] = (isless(minpos[0], 0.0) ? -minpos[0] : 0.0);
	shift[1] = (isless(minpos[1], 0.0) ? -minpos[1] : 0.0);
	shift[2] = (isless(minpos[2], 0.0) ? -minpos[2] : 0.0);
	io_logging_msg(log, INT32_C(4),
	               "Applying shift: (%g, %g, %g)",
	               shift[0], shift[1], shift[2]);

  
  
	/*================================================================
   *             conversion to AHF internal units
   *
   *    NOTE: the thermal energy is scaled to (km/sec)^2 !!!!!
   *================================================================*/
	scale_pos    = 1.0/(*boxsize);
#ifdef NO_EXPANSION
	scale_mom    = 1.0 / (*boxsize * posscale * 100.);
#else
	scale_mom    = sqrt(expansion) * expansion / (*boxsize * posscale * 100.);
#endif
	scale_weight = 1.0/mmass;
  scale_u      = 1.0; // already in (km/sec)^2 !?
  
  
    
	io_logging_msg(log, INT32_C(3),  "Scaling by:  positions:  %g", scale_pos);
	io_logging_msg(log, INT32_C(3),  "             velocities: %g", scale_mom);
	io_logging_msg(log, INT32_C(3),  "             weights:    %g", scale_weight);
	io_logging_msg(log, INT32_C(3),  "             u:          %g", scale_u);

  
  /* keep track of the applied shifts and scales */
  simu.pos_shift[0] = shift[0];
  simu.pos_shift[1] = shift[1];
  simu.pos_shift[2] = shift[2];
  simu.pos_scale    = scale_pos;
  

  
  
	/* Define the actual scaling calls type independent */
#	define SCALE_CALL(type) {\
		*((type *)strg.posx.val) += (type)(shift[0]); \
		*((type *)strg.posx.val) *= (type)(scale_pos); \
		*((type *)strg.posy.val) += (type)(shift[1]); \
		*((type *)strg.posy.val) *= (type)(scale_pos); \
		*((type *)strg.posz.val) += (type)(shift[2]); \
		*((type *)strg.posz.val) *= (type)(scale_pos); \
		*((type *)strg.momx.val) *= (type)(scale_mom); \
		*((type *)strg.momy.val) *= (type)(scale_mom); \
		*((type *)strg.momz.val) *= (type)(scale_mom); \
		if (strg.weight.val != NULL) *((type *)strg.weight.val) *= (type)(scale_weight); \
    if (strg.u.val != NULL) {if(*((type *)strg.u.val) >= 0.) *((type *)strg.u.val) *= (type)(scale_u);} \
		strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride); \
		strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride); \
		strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride); \
		strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride); \
		strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride); \
		strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride); \
		if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride); \
    if (strg.u.val != NULL) strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride); \
	}

	/* Do the scaling depending on the storage type */
	if (strg.bytes_float == sizeof(float)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(float);
		}
	} else if (strg.bytes_float == sizeof(double)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(double);
		}
	} else if (strg.bytes_float == sizeof(long double)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(long double);
		}
	} else {
		io_logging_fatal(log,
		                 "Don't know which floating point types "
		                 "has %" PRIi32 " bytes. Aborting read.",
		                 strg.bytes_float);
		return UINT64_C(0);
	}

	/* Clean the macro away */
#	undef SCALE_CALL

	/* And we are done */
	return particles_read;
}

#ifdef WITH_MPI
extern void
io_gadget_scale_global(io_logging_t log,
                       MPI_Comm comm,
                       double *maxpos,
                       double *minpos,
                       double *mmass)
{
	int size, rank;
	double buffer[3];

	io_logging_msg(log, INT32_C(5),
	               "Updating local scale values to global values.");
	MPI_Allreduce((void *)maxpos, (void *)buffer, 3,
	              MPI_DOUBLE, MPI_MAX, comm);
	io_logging_msg(log, INT32_C(5),
	               "local : maxpos[0] = %g \t"
	               "maxpos[1] = %g \t"
	               "maxpos[2] = %g",
	               maxpos[0], maxpos[1], maxpos[2]);
	maxpos[0] = buffer[0];
	maxpos[1] = buffer[1];
	maxpos[2] = buffer[2];
	io_logging_msg(log, INT32_C(5),
	               "global: maxpos[0] = %g \t"
	               "maxpos[1] = %g \t"
	               "maxpos[2] = %g",
	               maxpos[0], maxpos[1], maxpos[2]);

	MPI_Allreduce((void *)minpos, (void *)buffer, 3,
	              MPI_DOUBLE, MPI_MIN, comm);
	io_logging_msg(log, INT32_C(5),
	               "local : minpos[0] = %g \t"
	               "minpos[1] = %g \t"
	               "minpos[2] = %g",
	               minpos[0], minpos[1], minpos[2]);
	minpos[0] = buffer[0];
	minpos[1] = buffer[1];
	minpos[2] = buffer[2];
	io_logging_msg(log, INT32_C(5),
	               "global: minpos[0] = %g \t"
	               "minpos[1] = %g \t"
	               "minpos[2] = %g",
	               minpos[0], minpos[1], minpos[2]);

	MPI_Allreduce((void *)mmass, (void *)buffer, 1,
	              MPI_DOUBLE, MPI_MIN, comm);
	io_logging_msg(log, INT32_C(5), "local : mmass = %g", *mmass);
	*mmass = buffer[0];
	io_logging_msg(log, INT32_C(5), "global: mmass = %g", *mmass);

	return;
}
#endif /* WITH_MPI */


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/

inline static io_gadget_t
local_openopen(io_logging_t log, io_gadget_t f, io_file_mode_t mode)
{
	if (mode == IO_FILE_READ) {
		f->file = fopen(f->fname, IO_FILE_MODE_READ);
		if (f->file == NULL) {
			io_logging_fatal(log, "Could not open '%s' for reading.", f->fname);
			return NULL;
		}
	} else {
		f->file = fopen(f->fname, IO_FILE_MODE_WRITE);
		if (f->file == NULL) {
			io_logging_fatal(log, "Could not open '%s' for writing.", f->fname);
			return NULL;
		}
	}

	f->mode = mode;

	return f;
}

inline static void
local_openswapped(io_logging_t log,
                  io_gadget_t f,
                  io_file_swap_t swapped)
{
	if (f->mode == IO_FILE_READ) {
		switch(swapped) {
			case IO_FILE_ISNOT_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming unswapped file");
				break;
			case IO_FILE_IS_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming swapped file");
				break;
			case IO_FILE_UNKOWN_SWAPPING:
			default:
				io_logging_msg(log, INT32_C(3),
				               "Will try to find out swap status");
		}
	}

	/* Now set the swapping */
	f->swapped = swapped;

	return;
}

inline static io_gadget_t
local_opengetswap(io_logging_t log, io_gadget_t f)
{
	uint32_t bbound1, bbound2, tmp;

	/* Get the first block size from the start of the file */
	rewind(f->file);
	fread((void *)(&tmp), sizeof(int32_t), 1, f->file);

	/* Do a byteswap on a copy of that block size */
	bbound1 = tmp;
	io_util_sexchange(&bbound1, sizeof(uint32_t));

	io_logging_msg(log, INT32_C(0),
	               "Boundary (file): %" PRIu32
	               " Boundary (swapped): %" PRIu32,
	               tmp, bbound1);

	/* If the bytewapped block size is smaller, then it is probable that
	 * the file is byteswapped.
	 */
	if (bbound1 > tmp) {
		bbound1 = tmp;
		f->swapped = IO_FILE_ISNOT_SWAPPED;
		io_logging_msg(log, INT32_C(0),
		               "Trying nonswapped...");
	} else {
		f->swapped = IO_FILE_IS_SWAPPED;
		io_logging_msg(log, INT32_C(0),
		               "Trying swapped...");
	}

	io_logging_msg(log, INT32_C(0),
	               "Will skip %" PRIu32 " bytes...", bbound1);

	/* Test the assumption by skipping the block and reading the end
	 * block delimiter */
	fseek(f->file, (long)bbound1, SEEK_CUR);
	if (io_util_readuint32(f->file, &bbound2, f->swapped) == 0) {
		io_logging_fatal(log,
		                 "Could not read the second block delimiter. "
		                 "Corrupt file?");
		return NULL;
	}

	io_logging_msg(log, INT32_C(0),
	               "Second boundary (file): %" PRIu32,
	               bbound2);

	if (bbound1 == bbound2) {
		if (f->swapped == IO_FILE_IS_SWAPPED) {
			io_logging_msg(log, INT32_C(2),
			               "Identified a byte swapped file.");
		} else {
			io_logging_msg(log, INT32_C(2),
			               "Identified a not byte swapped file.");
		}
	} else {
		/* Swap assumption failed, trying it the other way around */
		io_util_sexchange(&bbound1, sizeof(uint32_t));
		fseek(f->file, (long)(sizeof(uint32_t)+bbound1), SEEK_SET);

		if (f->swapped == IO_FILE_IS_SWAPPED)
			f->swapped = IO_FILE_ISNOT_SWAPPED;
		else
			f->swapped = IO_FILE_IS_SWAPPED;

		if (io_util_readuint32(f->file, &bbound2, f->swapped) == 0) {
			io_logging_fatal(log,
			                 "Could not read the second block delimiter. "
			                 "Corrupt file?");
			return NULL;
		}

		/* See if it worked now */
		if (bbound1 == bbound2) {
			if (f->swapped == IO_FILE_IS_SWAPPED) {
				io_logging_msg(log, INT32_C(2),
				               "Identified a byte swapped file (2. try).");
			} else {
				io_logging_msg(log, INT32_C(2),
				               "Identified a not byte swapped file "
				               "(2. try).");
			}
		}
		else {
			/* Bad, cannot read this file */
			io_logging_fatal(log, "Cannot identify swapping status!");
			return NULL;
		}
	}

	return f;
}

inline static void
local_openversion(io_logging_t log, io_gadget_t f)
{
	uint32_t bbound;
	char tmp[5];

	/* We always want to write Gadget2 files */
	if (f->mode == IO_FILE_WRITE) {
		f->ver = 2;
		return;
	}

	/* Read the first block boundary and the first 4 bytes of the block */
	fseek(f->file, 0L, SEEK_SET); // this is equivalent to frewind(f->file)
	io_util_readuint32(f->file, &bbound, f->swapped);
	io_util_readstring(f->file, tmp, 4);

	/* If it is a Gadget 2 file, only HEAD and an integer will be
	 * stored in the block.
	 */
	if ( (bbound == 4+sizeof(int)) && (strncmp("HEAD", tmp, 4) == 0) ) {
		f->ver = 2;
		io_logging_msg(log, INT32_C(2), "Found a Gadget2 file.");
	} else {
		f->ver = 1;
		io_logging_msg(log, INT32_C(2), "Assuming a Gadget1 file.");
	}

	return;
}

inline static int32_t
local_write_common(io_logging_t log,
                   io_gadget_t f,
                   uint64_t pskip,
                   int32_t bytes)
{
	return 0;
}

#define CHECK_FLOATBYTES(bfile, bstore) {\
	if (bfile == sizeof(float)) { \
		io_logging_msg(log, INT32_C(1), \
		               "Obviously the file uses float for " \
		               "floating point values (%" PRIi32 " bytes).", \
		               bfile); \
	} else if (bfile == sizeof(double)) { \
		io_logging_msg(log, INT32_C(1), \
		               "Obviously the file uses double for " \
		               "floating point values (%" PRIi32 " bytes).", \
		               bfile); \
	} else { \
		io_logging_fatal(log, \
		                 "No clue what kind of floating point uses " \
		                 "%" PRIi32 " bytes. Aborting reading.", \
		                 bfile); \
		return UINT64_C(0); \
	}\
	if (bfile < bstore) { \
		io_logging_msg(log, INT32_C(1), \
		               "The floating point values in the file have " \
		               "less precision than the particle storage " \
		               "(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
		               "No problem, will upcast.", \
		                bfile, bstore); \
	} else if (bfile > bstore) { \
		io_logging_warn(log, INT32_C(1), \
		                "The floating point values in the file have " \
		                "a higher precision than the particle storage " \
		                "(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
		                "Will downcast, but precision might be lost.", \
		                bfile, bstore); \
	} \
}
static uint64_t
local_get_block_pos(io_logging_t log,
                    io_gadget_t f,
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg)
{
	uint32_t blocksize, blocksize2;
	int32_t partsize;
	uint32_t bytes_file;
	double fposx, fposy, fposz;
	float dummy;
	uint64_t i;

	/* Figure out how many bytes are used for float storage */
	SKIP;
	if (blocksize >  (unsigned)INT_MAX) {
		io_logging_warn(log, INT32_C(1),
		                "Exceeding capacity of Gadget files. "
		                "We will try to deal with this, but in "
		                "principle you should not exceed block sizes "
		                "larger than %i bytes.  Blocks larger than "
		                "%u bytes will fail completely.  Consider "
		                "splitting the snapshot over more files.",
		                INT_MAX, UINT_MAX);
	}
	bytes_file = blocksize / (3*f->no_part);
#ifdef AHOBBS_GADGET_FILE_FIX
  bytes_file = sizeof(float);
#endif
	CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " particle positions "
	               "with %" PRIi32 " bytes per float (%f MB total) "
	               "are stored.",
	               f->no_part, bytes_file,
	               (float)(blocksize/1024./1024.));

	/* Set the number of particles to loop over correctly */
	io_logging_msg(log, INT32_C(3),
	               "Asked to read %" PRIu64 " and to skip %" PRIu64
	               " particles. Checking those numbers.",
	               *pread, *pskip);
	if (*pskip > f->no_part) {
		*pskip = f->no_part;
		io_logging_msg(log, INT32_C(3),
		               "Cannot skip more than there is, will now "
		               "only skip %" PRIu64 " particles.",
		               *pskip);
	}
	if ( f->no_part - *pskip < *pread ) {
		*pread = f->no_part - *pskip;
		io_logging_msg(log, INT32_C(3),
		               "Cannot read more than there is left after "
		               "skipping. Will only read %" PRIu64
		               "particles.", *pread);
	}

	/* Set the particle size */
	partsize = 3*bytes_file;

	/* Go to the first particle we want to read */
	fseek(f->file, partsize*(*pskip), SEEK_CUR);

	/* Set extreme position detectors */
	f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e40;
	f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e40;

	/* Loop over the particle positions */
	for (i=0; i<*pread; i++) {
		/* STEP 1:  Read the particle data form the file */
		if (bytes_file == sizeof(float)) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposx = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposy = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposz = (double)dummy;
		} else if (bytes_file == sizeof(double)){
			io_util_readdouble(f->file, &fposx, f->swapped);
			io_util_readdouble(f->file, &fposy, f->swapped);
			io_util_readdouble(f->file, &fposz, f->swapped);
		}
		/* STEP 2:  Detect extreme positions */
		if (isless(fposx, f->minpos[0]))
			f->minpos[0] = fposx;
		if (isless(fposy, f->minpos[1]))
			f->minpos[1] = fposy;
		if (isless(fposz, f->minpos[2]))
			f->minpos[2] = fposz;
		if (isgreater(fposx, f->maxpos[0]))
			f->maxpos[0] = fposx;
		if (isgreater(fposy, f->maxpos[1]))
			f->maxpos[1] = fposy;
		if (isgreater(fposz, f->maxpos[2]))
			f->maxpos[2] = fposz;
		/* STEP 3:  Store the particle in the array */
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.posx.val) = (float)fposx;
			*((float *)strg.posy.val) = (float)fposy;
			*((float *)strg.posz.val) = (float)fposz;
		} else {
			*((double *)strg.posx.val) = fposx;
			*((double *)strg.posy.val) = fposy;
			*((double *)strg.posz.val) = fposz;
		}
		/* STEP 4:  Increment the pointers to the next particle */
		strg.posx.val = (void *)(((char *)strg.posx.val)
		                          + strg.posx.stride);
		strg.posy.val = (void *)(((char *)strg.posy.val)
		                          + strg.posy.stride);
		strg.posz.val = (void *)(((char *)strg.posz.val)
		                          + strg.posz.stride);
	} /* End of particle position loop */

	/* Go to the end of the particle position block */
	fseek(f->file, partsize*(f->no_part - (*pread + *pskip)), SEEK_CUR);
	SKIP2;
	CHECK_BLOCK;

	/* And return the number of read particles for error checking */
	return *pread;
}

static uint64_t
local_get_block_vel(io_logging_t log,
                    io_gadget_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg)
{
	uint32_t blocksize, blocksize2;
	int32_t partsize;
	uint32_t bytes_file;
	double fmomx, fmomy, fmomz;
	float dummy;
	uint64_t i;

	/* Figure out how many bytes are used for float storage */
	SKIP;
	bytes_file = blocksize / (3*f->no_part);
#ifdef AHOBBS_GADGET_FILE_FIX
  bytes_file = sizeof(float);
#endif
	CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " particle velocities "
	               "with %" PRIi32 " bytes per float (%f MB total) "
	               "are stored.",
	               f->no_part, bytes_file,
	               (float)(blocksize/1024./1024.));

	/* Set the particle size */
	partsize = 3*bytes_file;

	/* Go to the first particle we want to read */
	fseek(f->file, partsize*pskip, SEEK_CUR);

	/* Loop over the particle velocities */
	for (i=0; i<pread; i++) {
		/* STEP 1:  Read the particle data form the file */
		if (bytes_file == sizeof(float)) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomx = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomy = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomz = (double)dummy;
		} else {
			io_util_readdouble(f->file, &fmomx, f->swapped);
			io_util_readdouble(f->file, &fmomy, f->swapped);
			io_util_readdouble(f->file, &fmomz, f->swapped);
		}
		/* STEP 2:  Store the particle in the array */
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.momx.val) = (float)fmomx;
			*((float *)strg.momy.val) = (float)fmomy;
			*((float *)strg.momz.val) = (float)fmomz;
		} else {
			*((double *)strg.momx.val) = fmomx;
			*((double *)strg.momy.val) = fmomy;
			*((double *)strg.momz.val) = fmomz;
		}
		/* STEP 3:  Increment the pointers to the next particle */
		strg.momx.val = (void *)(((char *)strg.momx.val)
		                         + strg.momx.stride);
		strg.momy.val = (void *)(((char *)strg.momy.val)
		                         + strg.momx.stride);
		strg.momz.val = (void *)(((char *)strg.momz.val)
		                         + strg.momx.stride);
	} /* End of particle velocity loop */

	/* Go to the end of the particle velocity block */
	fseek(f->file, partsize*(f->no_part - (pread + pskip)), SEEK_CUR);
	SKIP2;
	CHECK_BLOCK;

	/* And return the number of read particles for error checking */
	return pread;
}

static uint64_t
local_get_block_id(io_logging_t log,
                   io_gadget_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg)
{
	uint32_t blocksize, blocksize2;
	int32_t partsize;
	uint32_t bytes_int_file;
	uint64_t fid;
	uint64_t i;
	uint32_t dummy_int;

	/* Figure out how many bytes are used for int storage */
	SKIP;
	bytes_int_file = blocksize / f->no_part;
#ifdef AHOBBS_GADGET_FILE_FIX
  bytes_int_file = sizeof(uint32_t);
#endif
	if (    (bytes_int_file != sizeof(uint32_t))
	     && (bytes_int_file != sizeof(uint64_t))) {
		io_logging_fatal(log,
		                 "Can't handle reading of integers "
		                 "with %" PRIi32 " bytes. Aborting.",
		                 bytes_int_file);
		return pread;
	}
	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " particle IDs "
	               "with %" PRIi32 " bytes per int (%f MB total) "
	               "are stored.",
	               f->no_part, bytes_int_file,
	               (float)(blocksize/1024./1024.));
	if (bytes_int_file < strg.bytes_int) {
		io_logging_warn(log, INT32_C(1),
		                "File uses %" PRIi32 " bytes per integer for "
		                "the IDs, the particle storage has %" PRIi32
		                " bytes available. Will upcast.",
		                bytes_int_file, strg.bytes_int);
	}
	if (bytes_int_file > strg.bytes_int) {
		io_logging_warn(log, INT32_C(1),
		                "File uses %" PRIi32 " bytes per integer for "
		                "the IDs, the particle storage has only %" PRIi32
		                " bytes available. Will downcast, be aware that "
		                "this might lead to bogus values...",
		                bytes_int_file, strg.bytes_int);
	}

	/* Reset partsize to the bytes user for integer storage value */
	partsize = bytes_int_file;

	/* Go to the first particle we want to read */
	fseek(f->file, partsize*pskip, SEEK_CUR);

	/* See if we have to read the IDs */
	if (strg.id.val == NULL) {
		io_logging_warn(log, INT32_C(1),
	    	            "Discarding IDs as no storage for the IDs "
		                "has been specified.");
		fseek(f->file, partsize*pread, SEEK_CUR);
	} else {
		/* Loop over the particle IDs */
		for (i=0; i<pread; i++) {
			/* STEP 1:  Read the ID from the file */
			if (bytes_int_file == sizeof(uint32_t)) {
				io_util_readuint32(f->file, &dummy_int, f->swapped);
				fid = (uint64_t)dummy_int;
			} else  {
				io_util_readuint64(f->file, &fid, f->swapped);
			}
			/* STEP 2:  Store the ID in the array */
			if (strg.bytes_int == 4) {
				*((uint32_t *)strg.id.val) = (uint32_t)fid;
			} else {
				*((uint64_t *)strg.id.val) = fid;
			}
			/* STEP 3:  Increment the pointers to the next particle */
			strg.id.val = (void *)(((char *)strg.id.val) + strg.id.stride);
		} /* End of particle ID loop */
	} /* End of catch NULL-Id */

	/* Go to the end of the particle ID block */
	fseek(f->file, partsize*(f->no_part - (pread + pskip)), SEEK_CUR);
	SKIP2;
	CHECK_BLOCK;

	/* And return the number of read particles for error checking */
	return pread;
}

static uint64_t
local_get_block_mass(io_logging_t log,
                     io_gadget_t f,
                     uint64_t pskip,
                     uint64_t pread,
                     io_file_strg_struct_t strg)
{
	uint32_t blocksize=0, blocksize2=0;
	int32_t partsize;
	uint32_t bytes_file;
	double fweight, oldfweight;
	uint32_t curprtt;
	float dummy;
	uint64_t i, j;


	/* We are going to need that a few times for book-keeping */
#	define BOOK_KEEPING {\
		if (   isgreater(fweight, oldfweight) \
		    || isless(fweight, oldfweight)) { \
			f->no_species++; \
			oldfweight = fweight; \
			if (isless(fweight,f->minweight)) \
				f->minweight = fweight; \
			if (isgreater(fweight, f->maxweight)) \
				f->maxweight = fweight; \
			if (    (i >= f->header->np[0]) \
			     && (i < f->header->np[0]+f->header->np[1]) \
			     && isless(fweight, f->mmass) ) \
				f->mmass = fweight; \
		} \
	}

	/* First see if there is a mass block at all, set the partsize */
	if (f->multimass != 0) {
		/* Figure out how many bytes are used for float storage */
		SKIP;
		bytes_file = blocksize / (f->no_part_with_mass);
		io_logging_msg(log, INT32_C(1), "blocksize = %i",
		               (int)blocksize);
#ifdef AHOBBS_GADGET_FILE_FIX
    bytes_file = sizeof(float);
#endif
		CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
		partsize = bytes_file;
		io_logging_msg(log, INT32_C(1),
		               "A total of %" PRIu64 " particle weights "
		               "with %" PRIi32 " bytes per float (%f MB total) "
		               "are stored.",
		               f->no_part_with_mass, bytes_file,
		               (float)(blocksize/1024./1024.));
	
	} else {
		partsize = strg.bytes_float;
		io_logging_msg(log, INT32_C(1),
		               "Will construct mass information solely from"
		               " the header.");
	}

	/* Check if we can store the masses, if not they will be
	 * discarded */
	if (strg.weight.val == NULL) {
		io_logging_warn(log, INT32_C(1),
	    	            "Discarding masses as no storage for the "
		                "masses has been specified.");
	}

	/* Initialize some things */
	f->sumweight = 0.0;
	oldfweight = 0.0;
	f->no_species = 0;

	/* Set the particle type to the first type that actually occurs */
	j = 0;
	curprtt = 0;
	while (f->header->np[curprtt] == 0) {
		curprtt++;
	}

	/* Now we really read the particles */
	for (i=0; i<f->no_part; i++) {
		/* STEP 1:  Update the current particle type */
		if ((long)j >= f->header->np[curprtt]) {
			do {
				curprtt++;
			} while (f->header->np[curprtt] == 0);
			io_logging_msg(log, INT32_C(0),
			               "Detected new particle species! "
			               "i=%" PRIu64 " j=%" PRIu64
			               " curprtt=%" PRIi32 " np[curprtt]=%" PRIu64,
			               i, j, curprtt, f->header->np[curprtt]);
			j = 1;
		} else {
			j++;
		}
		/* STEP 2:  Get the particle weight */
		if (f->multimass & (1<<curprtt)) {
			/* Okay, for this particle type it is in the file */
			if (bytes_file == sizeof(float)) {
				io_util_readfloat(f->file, &dummy, f->swapped);
				fweight = (double)dummy;
			} else {
				io_util_readdouble(f->file, &fweight, f->swapped);
			}
		} else {
			/* For this particle type it is in the header */
			fweight = f->header->massarr[curprtt];
		}
		/* STEP 3:  Do the book-keeping */
		BOOK_KEEPING;
		f->sumweight += fweight;
		/* STEP 4:  Store the particle if it is requested */
		if ( (pskip <= i) && (i-pskip < pread) ) {
			if (strg.weight.val != NULL) {
				/* Okay, there is something to write to */
				if (strg.bytes_float == 4) {
					*((float *)strg.weight.val) = (float)fweight;
				} else {
					*((double *)strg.weight.val) = fweight;
				}
				strg.weight.val = (void *)(((char *)strg.weight.val)
				                           + strg.weight.stride);
			} else {
				/* WE DISCARD IT! */
				;
			}
		}
	} /* End of particle weight loop */

	/* If necessary, verify that the reading went okay */
	if (f->multimass != 0) {
		SKIP2;
		CHECK_BLOCK;
	}
#	undef BOOK_KEEPING

	/* And return the number of read particles for error checking */
	return pread;
}

static uint64_t
local_get_block_u(io_logging_t log,
                  io_gadget_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg)
{
	uint32_t blocksize, blocksize2;
	int32_t partsize;
	uint32_t bytes_file;
	double fweight, oldfweight;
	uint32_t curprtt;
	uint64_t i, psum;
	float dummy;
	double fu;
	int ptype;

	/* See if there is a gas block */
	if (f->header->np[0] > 0) {
		/* Figure out how many bytes are used for float storage */
		SKIP;
		bytes_file = blocksize / (f->header->np[0]);
		CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
		io_logging_msg(log, INT32_C(1),
		               "A total of %" PRIu64 " gas particle energies "
		               "with %" PRIi32 " bytes per float (%f MB total) "
		               "are stored.",
		               f->header->np[0], bytes_file,
		               (float)(blocksize/1024./1024.));
	}
  else {
    /* There are no gas particles in the file */
    bytes_file = 0;
   }

	/* Set the particle size */
	partsize = bytes_file;

	/* Only go to the gas block if required and if it is actually there */
	if (f->header->np[0] > 0 && pskip <= (uint64_t)(f->header->np[0]) && strg.u.val != NULL) {
		/* Go to the first particle we want to read */
		fseek(f->file, partsize*pskip, SEEK_CUR);

		/* Loop over gas particles */
		for (i=0; i<f->header->np[0]-pskip && i<pread; i++) {
			/* STEP 1:  Read the energy from the file */
			if (bytes_file == sizeof(float)) {
				io_util_readfloat(f->file, &dummy, f->swapped);
				fu = (double)dummy;
			} else  {
				io_util_readdouble(f->file, &fu, f->swapped);
			}
			/* STEP 2:  Store the energy in the array */
			if (strg.bytes_float == sizeof(float)) {
				*((float *)strg.u.val) = (float)fu;
			} else {
				*((double *)strg.u.val) = fu;
			}
			/* STEP 3:  Increment the pointers to the next particle */
			strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride);
		}

		/* Skip to the end of the energy block */
		fseek(f->file, partsize*(f->header->np[0]-pskip-i), SEEK_CUR);
	} else {
		/* i carries the information of how many energies got read,
		 * since we read none, set it to 0 */
		i = 0;
		/* Skip to the end of the energy block */
		fseek(f->file, partsize*(f->header->np[0]), SEEK_CUR);
	}

	/* If there was an energy block, finish reading it. */
	if (f->header->np[0] > 0) {
		SKIP2;
		CHECK_BLOCK;
	}

	/* Set the energy to the (negative) particle type for the rest*/
	/* Set the sum of all particles already looped over (either skipped,
	 * or actually read) */
	psum = f->header->np[0];
	/* The sum incorporates all particles up to and including this type */
	ptype = 0;
	/* Now we do something ugly and use the energy to store the type of
	 * the particles.  However, we will not do that, if there is no
	 * energy storage provided in the particle structure. */
	if (strg.u.val != NULL) {
		for (; i<pread; i++) {
			/* Figure out if we are still at the right particle type */
			while (i+pskip>=psum && ptype<5) {
				ptype++;
				psum += f->header->np[ptype];
			}
			fu = (double)(-ptype);   // check PGAS, PDM, PSTAR
			/* STEP 1:  Store the energy in the array */
			if (strg.bytes_float == sizeof(float)) {
				*((float *)strg.u.val) = (float)fu;
			} else {
					*((double *)strg.u.val) = fu;
			}
			/* STEP 2:  Increment the pointers to the next particle */
			strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride);
		}
	}

	/* And return the number of read particles for error checking */
	return pread;
}

#ifdef GADGET_MAGNETICUM
#include "../define.h"
#define GET_BLOCK {\
SKIP;\
io_util_readstring(f->file, str, (size_t)4);\
io_util_readuint32(f->file, &nextblocksize, f->swapped);\
io_logging_msg(log, INT32_C(1),\
"Arrived at block %s, size of it will be %" \
PRIi32, str, nextblocksize);\
SKIP2;\
CHECK_BLOCK;\
}
void local_find_block(io_gadget_t f, char *blockname)
{
	uint32_t blocksize, blocksize2, nextblocksize;
	int32_t partsize, bytes_file;
	char str[5];
	int tries;
  
  str[0] = 'A';
  tries = 0;
  nextblocksize = 0;
  while ( (strncmp(str, blockname, 4) != 0) && (tries < 10)) {
    fprintf(stderr,"skipping block %s of size %"PRIu32"\n",str,nextblocksize);
    fseek(f->file, nextblocksize, SEEK_CUR);
    GET_BLOCK;
    tries++;
  }
  fprintf(stderr,"found block %s of size %"PRIu32"\n",str,nextblocksize);
}
#endif


#ifdef METALHACK
#include "../define.h"
#define GET_BLOCK {\
	SKIP;\
	io_util_readstring(f->file, str, (size_t)4);\
	io_util_readuint32(f->file, &nextblocksize, f->swapped);\
	io_logging_msg(log, INT32_C(1),\
	               "Arrived at block %s, size of it will be %" \
	               PRIi32, str, nextblocksize);\
	SKIP2;\
	CHECK_BLOCK;\
}

static uint64_t
local_get_block_z(io_logging_t log,
                  io_gadget_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg)
{
	uint32_t blocksize, blocksize2, nextblocksize;
	int32_t partsize, bytes_file;
	char str[5];
	int tries;
	uint64_t i, metalzread, starsoffset, metalzskip;
	double fz;
	float dummy;
	
	/* if this is not a Gadget 2 file, skip until METAL block */
  /* (assuming the orderting POS,VEL,ID,MASS,UGAS,RHO,NE,NH,HSML,SFR,AGE,Z) */
	if (f->ver != 2) {
    local_skip_GADGET1_blocks(f->file, f->swapped, 11); // we skip the first 11 blocks
	}

  /* we have a Gadget 2 file and hence use the HEAD to find Z */
  else {
    /* Go to the metal block */
    str[0] = 'A';
    tries = 0;
    nextblocksize = 0;
    while ( (strncmp(str, "Z   ", 4) != 0) && (tries < 100)) {
#ifdef GADGET_MAGNETICUM
      fprintf(stderr,"skipping block %s of size %"PRIu32"\n",str,nextblocksize);
#endif
      fseek(f->file, nextblocksize, SEEK_CUR);
      GET_BLOCK;
      tries++;
    }
#ifdef GADGET_MAGNETICUM
    fprintf(stderr,"found block %s of size %"PRIu32"\n",str,nextblocksize);
#endif
    if (tries >= 150) {
      METALDIE;
    }
  }

	/* Start with the block */
	SKIP;
	bytes_file = blocksize / (f->header->np[0]+f->header->np[4]);
	CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " metallicities "
	               "with %" PRIi32 " bytes per float (%f MB total) "
	               "are stored.",
	               f->header->np[0]+f->header->np[4], bytes_file,
	               (float)(blocksize/1024./1024.));

	/* Set the particle size */
	partsize = bytes_file;

	/* Find the right spot */
	metalzskip = UINT64_C(0);
	metalzread = UINT64_C(0);
	starsoffset =   f->header->np[0] + f->header->np[1]
	              + f->header->np[2] + f->header->np[3];
	if (pskip < f->header->np[0])
		metalzskip = pskip; // Skip the first gas particles
	else if (pskip < starsoffset)
		metalzskip = f->header->np[0]; // Skip all gas
	else if (pskip < starsoffset + f->header->np[4])
		metalzskip = f->header->np[0] + pskip - starsoffset; // Skip all gas and a few stars
	else
		metalzskip = f->header->np[0] + f->header->np[4]; // Skip all gas and all stars

	/* Now we skip to right position in the file */
	fseek(f->file, metalzskip * partsize, SEEK_CUR);

	/* Read the gas metalz */
	for (i=0; i<pread && i+pskip<f->header->np[0]; i++) {
		/* STEP 1:  Read the metal from the file */
		if (bytes_file == sizeof(float)) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fz = (double)dummy;
		} else  {
			io_util_readdouble(f->file, &fz, f->swapped);
		}
		/* STEP 2:  Store the metal in the array */
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.z.val) = (float)fz;
		} else {
			*((double *)strg.z.val) = fz;
		}
		/* STEP 3:  Increment the pointers to the next particle */
		strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
		metalzread++;
	}

	/* Conjure up metalz for the nonmetal particles */
	for (; i<pread && i+pskip<starsoffset; i++) {
		fz = 0.0;
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.z.val) = (float)fz;
		} else {
			*((double *)strg.z.val) = fz;
		}
		strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
	}

	/* Read the star metalz */
	for (; i<pread && i+pskip<starsoffset+f->header->np[4]; i++) {
		/* STEP 1:  Read the metal from the file */
		if (bytes_file == sizeof(float)) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fz = (double)dummy;
		} else  {
			io_util_readdouble(f->file, &fz, f->swapped);
		}
		/* STEP 2:  Store the metal in the array */
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.z.val) = (float)fz;
		} else {
			*((double *)strg.z.val) = fz;
		}
		/* STEP 3:  Increment the pointers to the next particle */
		strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
		metalzread++;
	}

	/* Conjure up metalz for the nonmetal particles */
	for (; i<pread; i++) {
		fz = 0.0;
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.z.val) = (float)fz;
		} else {
			*((double *)strg.z.val) = fz;
		}
		strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
	}

	/* Go to the end of the block */
	fseek(f->file,
	      partsize*(f->header->np[0]+f->header->np[4]-metalzread-metalzskip),
	      SEEK_CUR);

	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " metals were read.",
	               metalzread);

	/* Finish the block */
	SKIP2;
	CHECK_BLOCK;

	/* And return the number of read particles for error checking */
	return pread;
}

static uint64_t
local_get_block_age(io_logging_t log,
                    io_gadget_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg)
{
	uint32_t blocksize, blocksize2, nextblocksize;
	int32_t partsize, bytes_file;
	char str[5];
	int tries;
	uint64_t i, agesread, agesskip, starsoffset;
	double fage;
	float dummy;
	
	/* if this is not a Gadget 2 file, skip until AGE block */
  /* (assuming the orderting POS,VEL,ID,MASS,UGAS,RHO,NE,NH,HSML,SFR,AGE,METAL) */
	if (f->ver != 2) {
    local_skip_GADGET1_blocks(f->file, f->swapped, 10); // we skip the first 10 blocks
	}

  /* we have a Gadget 2 file and hence use the HEAD to find AGE */
  else {
    /* Go to the AGE block */
    str[0] = 'X';
    tries = 0;
    nextblocksize = 0;
    while ( (strncmp(str, "AGE ", 4) != 0) && (tries < 100)) {
#ifdef GADGET_MAGNETICUM
      fprintf(stderr,"skipping block %s of size %"PRIu32"\n",str,nextblocksize);
#endif
      fseek(f->file, nextblocksize, SEEK_CUR);
      GET_BLOCK;
      tries++;
    }
#ifdef GADGET_MAGNETICUM
    fprintf(stderr,"found block %s of size %"PRIu32"\n",str,nextblocksize);
#endif
    if (tries >= 150) {
      METALDIE;
    }
  }

	/* Start with the block */
	SKIP;
	bytes_file = blocksize / (f->header->np[4]);
	CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " star "
	               "with %" PRIi32 " bytes per float (%f MB total) "
	               "are stored.",
	               f->header->np[4], bytes_file,
	               (float)(blocksize/1024./1024.));

	/* Set the particle size */
	partsize = bytes_file;

	/* Find the right spot in the file */
	agesskip = UINT64_C(0);
	agesread = UINT64_C(0);
	starsoffset =   f->header->np[0] + f->header->np[1]
	              + f->header->np[2] + f->header->np[3];
	if (pskip < starsoffset)
		agesskip = 0;
	else if (pskip < starsoffset + f->header->np[4])
		agesskip = pskip - starsoffset;  // Skip the first few stars
	else
		agesskip = f->header->np[4]; // Skip all stars

	fseek(f->file, partsize*(agesskip), SEEK_CUR);

	/* Set ages of particles type 0-3*/
	for (i=0; i<pread && i+pskip<starsoffset; i++) {
		fage = 0.0;
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.age.val) = (float)fage;
		} else {
			*((double *)strg.age.val) = fage;
		}
		strg.age.val = (void *)(((char *)strg.age.val) + strg.age.stride);
	}

	/* Read the star ages */
	for (; i<pread && i+pskip<starsoffset+f->header->np[4]; i++) {
		/* STEP 1:  Read the age from the file */
		if (bytes_file == sizeof(float)) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fage = (double)dummy;
		} else  {
			io_util_readdouble(f->file, &fage, f->swapped);
		}
		/* STEP 2:  Store the age in the array */
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.age.val) = (float)fage;
		} else {
			*((double *)strg.age.val) = fage;
		}
		/* STEP 3:  Increment the pointers to the next particle */
		strg.age.val = (void *)(((char *)strg.age.val) + strg.age.stride);
		agesread++;
	}

	/* Set ages of particle type 5 */
	for (; i<pread; i++) {
		fage = 0.0;
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.age.val) = (float)fage;
		} else {
			*((double *)strg.age.val) = fage;
		}
		strg.age.val = (void *)(((char *)strg.age.val) + strg.age.stride);
	}

	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " ages were read.",
	               agesread);

	/* Go to the end of the block */
	fseek(f->file,
	      partsize*(f->header->np[4]-agesread-agesskip),
	      SEEK_CUR);

	/* Finish the block */
	SKIP2;
	CHECK_BLOCK;

	/* And return the number of read particles for error checking */
	return pread;
}
#undef GET_BLOCK
#endif

#undef CHECK_FLOATBYTES
#undef SKIP
#undef SKIP2
#undef CHECK_BLOCK

#ifdef METALHACK
void local_skip_GADGET1_blocks(FILE *fpgadget, int swapped, int nblocks)
{
  double       ddummy;
  int          idummy;
  int          nskipblocks;
  unsigned int uidummy;
  unsigned int blocksize1, blocksize2;
  int          iblock;
  char         unused[GADGET_HEADER_FILLHEADER+1];
  
  /* set the number of blocks to be skipped (excluding the header!) */
  nskipblocks = nblocks;
  
  /* rewind file back to the beginning */
  rewind(fpgadget);
  
#ifdef VERBOSE_GADGET1_HEADER
  /* skip header */
  io_util_readuint32(fpgadget, &blocksize1, swapped); // GADGET1-SKIP

	io_util_readint32(fpgadget, &idummy, swapped);   //fprintf(stderr,"np0=%d\n",idummy);
	io_util_readint32(fpgadget, &idummy, swapped);   //fprintf(stderr,"np1=%d\n",idummy);
	io_util_readint32(fpgadget, &idummy, swapped);   //fprintf(stderr,"np2=%d\n",idummy);
	io_util_readint32(fpgadget, &idummy, swapped);   //fprintf(stderr,"np3=%d\n",idummy);
	io_util_readint32(fpgadget, &idummy, swapped);   //fprintf(stderr,"np4=%d\n",idummy);
	io_util_readint32(fpgadget, &idummy, swapped);   //fprintf(stderr,"np5=%d\n",idummy);
	io_util_readdouble(fpgadget, &ddummy, swapped);
	io_util_readdouble(fpgadget, &ddummy, swapped);
	io_util_readdouble(fpgadget, &ddummy, swapped);
	io_util_readdouble(fpgadget, &ddummy, swapped);
	io_util_readdouble(fpgadget, &ddummy, swapped);
	io_util_readdouble(fpgadget, &ddummy, swapped);
	io_util_readdouble(fpgadget, &ddummy, swapped);  //fprintf(stderr,"expansion=%g\n",ddummy);
	io_util_readdouble(fpgadget, &ddummy, swapped);  //fprintf(stderr,"redshift=%g\n",ddummy);
	io_util_readint32(fpgadget, &idummy, swapped);
	io_util_readint32(fpgadget, &idummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readint32(fpgadget, &idummy, swapped);
	io_util_readint32(fpgadget, &idummy, swapped);
	io_util_readdouble(fpgadget, &ddummy, swapped);  //fprintf(stderr,"boxsize=%g\n",ddummy);
	io_util_readdouble(fpgadget, &ddummy, swapped);  //fprintf(stderr,"omega=%g\n",ddummy);
	io_util_readdouble(fpgadget, &ddummy, swapped);  //fprintf(stderr,"omegalambda=%g\n",ddummy);
	io_util_readdouble(fpgadget, &ddummy, swapped);  //fprintf(stderr,"hubble=%g\n",ddummy);
	io_util_readint32(fpgadget, &idummy, swapped);
	io_util_readint32(fpgadget, &idummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readuint32(fpgadget, &uidummy, swapped);
	io_util_readint32(fpgadget, &idummy, swapped);
	io_util_readstring(fpgadget, unused, GADGET_HEADER_FILLHEADER);
  
  io_util_readuint32(fpgadget, &blocksize2, swapped); // GADGET1-SKIP
  
  //fprintf(stderr,"skipped HEADER of size %u vs. %u bytes\n",blocksize1,blocksize2);
  
#else /* VERBOSE_GADGET1_HEADER */
  nskipblocks++;
#endif
  
  /* skip GADGET1 blocks */
  for(iblock=0; iblock<nskipblocks; iblock++)
   {
    io_util_readuint32(fpgadget, &blocksize1, swapped); // GADGET1-SKIP
    fseek(fpgadget, blocksize1, SEEK_CUR);              // block-SKIP
    io_util_readuint32(fpgadget, &blocksize2, swapped); // GADGET1-SKIP
    
    //fprintf(stderr,"skipped block #%d (out of %d in total) of size %u vs. %d MB\n",iblock,nblocks,blocksize1/1024/1024,blocksize2/1024/1024);
    
    /* the file appears to be corrupted */
    if(blocksize1 != blocksize2)
     {
      fprintf(stderr,"We are already trying to help you with your ancient GADGET1 file, but enough is enough!\n");
      exit(-1);
     }
   }
}
#endif
