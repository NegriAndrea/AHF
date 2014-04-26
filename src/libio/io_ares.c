/* $Id: io_ares.c,v 1.6 2008/07/30 11:48:07 knolli Exp $ */

/**
 * \file io_ares.c
 *
 * Provides functions for reading and writing AMIGA restart (ARES) files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdbool.h>

#include "io_ares.h"
#include "io_ares_header.h"
#include "io_util.h"


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
 * io_ares_open more readable.
 *
 * \param log   The logging object.
 * \param f     The ARES file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_ares_t
local_openopen(io_logging_t log, io_ares_t f, io_file_mode_t mode);

/**
 * \brief Helper funtion to set the swap-state and write some
 *        messages
 *
 * \param log      The logging object.
 * \param f        The ARES file object sofar.
 * \param swapped  The swap state.
 *
 * \return Nothing.
 */
inline static void
local_openswapped(io_logging_t log,
                  io_ares_t f,
                  io_file_swap_t swapped);

/**
 * \brief Check before writing particles and place the file pointer at
 *        the right place.
 *
 * \param log          The logging object.
 * \param f            The file object.
 * \param *weight      If NULL no particle weight will be written.
 * \param *u           If NULL no particle energy will be written.
 * \param pskip        Number of particles to skip in the file.
 * \param bytes_float  Number of bytes per floating point value.
 * \param bytes_int    Number of bytes per integer value.
 *
 * \return Returns true if everything went fine, false otherwise.
 */
inline static bool
local_write_common(io_logging_t log,
                   io_ares_t f,
                   void *weight,
                   void *u,
                   uint64_t pskip,
                   int32_t bytes_float,
                   int32_t bytes_int);


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_ares_t
io_ares_open(io_logging_t log,
             char *fname,
             io_file_swap_t swapped,
             io_file_mode_t mode,
             uint32_t reader)
{
	io_ares_t f;

	/* Get memory for the structure */
	f = (io_ares_t)malloc(sizeof(io_ares_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_ares structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Store the filename */
	f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
	if (f->fname == NULL) {
		io_logging_memfatal(log, "filename of ARES file");
		free(f);
		return NULL;
	}
	strncpy(f->fname, fname, strlen(fname)+1);

	/* Okay, we are an ARES file */
	f->ftype = IO_FILE_ARES;

	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	if (reader > 0) {
		MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
		MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
		MPI_Comm_split(MPI_COMM_WORLD, 1,
		               f->rank, &(f->mycomm));
		MPI_Comm_size(f->mycomm, &(f->size_mycomm));
		MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	} else {
		f->rank = f->size = f->size_mycomm = f->rank_mycomm = -1;
		f->mycomm = MPI_COMM_NULL;
	}
#	endif

	/* Try to open the file */
	if (local_openopen(log, f, mode) == NULL)
		return NULL;

	/* Set the swapping and be a bit verbose about it */
	local_openswapped(log, f, swapped);

	/* Nothing for the header for now */
	f->header = NULL;

	return f;
}

extern void
io_ares_close(io_logging_t log,
              io_ares_t *f)
{
	/* Catch NULLs */
	if (f == NULL || *f == NULL)
		return;

	/* Put the header to the file if necessary */
	if (    ((*f)->mode == IO_FILE_WRITE)
	     && ((*f)->header != NULL)) {
		io_ares_header_write(log, (*f)->header, *f);
	}

	/* Close */
	if ((*f)->file != NULL)
		fclose((*f)->file);
	
	/* Clean up */
	if ((*f)->header != NULL)
		io_ares_header_del(log, &((*f)->header));
	if ((*f)->fname != NULL)
		free((*f)->fname);
#	ifdef WITH_MPI
	if ((*f)->mycomm != MPI_COMM_NULL)
		MPI_Comm_free(&((*f)->mycomm));
#	endif
	free(*f);
	*f = NULL;

	return;
}

extern void
io_ares_init(io_logging_t log,
             io_ares_t f)
{
	if (f->mode == IO_FILE_READ) {
		io_logging_msg(log, INT32_C(5),
		               "Starting to initialize file object from %s",
		               f->fname);
		f->header = io_ares_header_get(log, f);
		io_logging_msg(log, INT32_C(5),
		               "Done with initializing file object from %s",
		               f->fname);
	} else {
		io_logging_warn(log, INT32_C(1),
		                "%s is not opened for reading. "
		                "Will do nothing.",
		                f->fname);
	}

	return;
}

/* Will be used to advance to the next particle. */
#define INCR(a,b) \
	(a) = (void *)(((char *)(a)) + (b))

extern uint64_t
io_ares_readpart(io_logging_t log,
                 io_ares_t f,
                 uint64_t pskip,
                 uint64_t pread,
                 io_file_strg_struct_t strg)
{
	uint64_t i;
	uint32_t partsize;
	uint32_t dummy_id;
	float dummy;
	double fposx, fposy, fposz;
	double fmomx, fmomy, fmomz;
	double fweight, fu;
#ifdef METALHACK
	double fz;
	double fage;
#endif
	uint64_t fid = UINT64_C(0);

	/* Check if we actually need to do something */
	if (    (f == NULL)
		 || ( f->header == NULL)
	     || (f->header->no_part_in_file == 0)
	     || (((uint64_t)(f->header->no_part_in_file)) < pskip) )
		return UINT64_C(0);

	/* Show what is going to happen */
#define show(bf, bs, what) {\
		if (bf > bs) {\
			io_logging_warn(log, INT32_C(3), \
			                "   " what " (read, downcast)");\
		} else if (bf == bs) {\
			io_logging_msg(log, INT32_C(3), "   " what " (read)");\
		} else {\
			io_logging_warn(log, INT32_C(3),\
			                "   " what " (read, upcast)");\
		}\
	}
	io_logging_msg(log, INT32_C(3), "Will do the following things:");
	show(f->header->bytes_float, strg.bytes_float, "posx ");
	show(f->header->bytes_float, strg.bytes_float, "posy ");
	show(f->header->bytes_float, strg.bytes_float, "posz ");
	show(f->header->bytes_float, strg.bytes_float, "momx ");
	show(f->header->bytes_float, strg.bytes_float, "momy ");
	show(f->header->bytes_float, strg.bytes_float, "momz ");
	if ( (f->header->has_mass != UINT64_C(0))) {
		if (strg.weight.val != NULL) {
			show(f->header->bytes_float, strg.bytes_float, "weight");
		} else {
			io_logging_warn(log, INT32_C(3), "   weight (read, ignored)");
		}
	} else {
		if (strg.weight.val != NULL)
			io_logging_warn(log, INT32_C(3), "   weight (set to 1.0)");
	}
	if (strg.id.val == NULL) {
		io_logging_warn(log, INT32_C(3),
		                "   id     (read, ignored)");
	} else {
		show(f->header->bytes_int, strg.bytes_int, "ID");
	}
	if (f->header->has_gas > UINT64_C(0)) {
		if (strg.u.val != NULL) {
			show(f->header->bytes_float, strg.bytes_float, "u");
		} else {
			io_logging_warn(log, INT32_C(3), "   u (read, ignored)");
		}
	} else {
		if (strg.u.val != NULL) {
			io_logging_warn(log, INT32_C(3), "   i (set to 0.0)");
		}
	}
#undef show

	/* Set the number of particles to loop over correctly */
	if ( (uint64_t)(f->header->no_part_in_file) - pskip < pread)
		pread = (uint64_t)(f->header->no_part_in_file) - pskip;

	io_logging_msg(log, INT32_C(3),
	               "Starting to read %" PRIu64 " particles, "
	               "skipping %" PRIu64,
	               pread, pskip);

	/* Figure out how many bytes a particle has in the file */
	partsize = f->header->bytes_float * 6;
	if (f->header->has_mass > UINT64_C(0))
		partsize +=  f->header->bytes_float;
	if (f->header->has_gas > UINT64_C(0))
		partsize += f->header->bytes_float;
	partsize += f->header->bytes_int;

	/* Position the file at the right spot */
	fseek(f->file,
	      (long)(pskip * partsize) + (long)ARES_HEADER_SIZE,
	      SEEK_SET);

	/* Loop over all particles to be read */
	for (i=0; i<pread; i++) {
		/************************************************\
		 * Step 1: Read the particle data from the file *
		\************************************************/
		if (f->header->bytes_float == 4) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposx = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomx = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposy = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomy = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposz = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomz = (double)dummy;
		} else {
			io_util_readdouble(f->file, &fposx, f->swapped);
			io_util_readdouble(f->file, &fmomx, f->swapped);
			io_util_readdouble(f->file, &fposy, f->swapped);
			io_util_readdouble(f->file, &fmomy, f->swapped);
			io_util_readdouble(f->file, &fposz, f->swapped);
			io_util_readdouble(f->file, &fmomz, f->swapped);
		}
		if (f->header->has_mass > UINT64_C(0)) {
			if (f->header->bytes_float == 4) {
				io_util_readfloat(f->file, &dummy, f->swapped);
				fweight = (double)dummy;
			} else {
				io_util_readdouble(f->file, &fweight, f->swapped);
			}
			if (isless(fweight, f->header->minweight))
				f->header->minweight = fweight;
			else if (isgreater(fweight, f->header->maxweight))
				f->header->maxweight = fweight;
		} else {
			fweight = 1.0;
		}
		if (f->header->has_gas > UINT64_C(0)) {
			if (f->header->bytes_float == 4) {
				io_util_readfloat(f->file, &dummy, f->swapped);
				fu = (double)dummy;
			} else {
				io_util_readdouble(f->file, &fu, f->swapped);
			}
		} else {
			fu = 0.0;
		}
#ifdef METALHACK
		if (f->header->bytes_float == 4) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fz = (double)dummy;
		} else {
			io_util_readdouble(f->file, &fz, f->swapped);
		}
		if (f->header->bytes_float == 4) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fage = (double)dummy;
		} else {
			io_util_readdouble(f->file, &fage, f->swapped);
		}
#endif
		if (f->header->bytes_int == 4) {
			io_util_readuint32(f->file, &dummy_id, f->swapped);
			fid = (uint64_t)dummy_id;
		} else if (f->header->bytes_int == 8) {
			io_util_readuint64(f->file, &fid, f->swapped);
		}

		/************************************************\
		 * Step 2: Store the particle in the array      *
		\************************************************/
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.posx.val) = (float)fmod(fposx + 1.0, 1.0);
			*((float *)strg.posy.val) = (float)fmod(fposy + 1.0, 1.0);
			*((float *)strg.posz.val) = (float)fmod(fposz + 1.0, 1.0);
			*((float *)strg.momx.val) = (float)fmomx;
			*((float *)strg.momy.val) = (float)fmomy;
			*((float *)strg.momz.val) = (float)fmomz;
			if (strg.weight.val != NULL)
				*((float *)strg.weight.val) = (float)fweight;
			if (strg.u.val != NULL)
				*((float *)strg.u.val) = (float)fu;
#ifdef METALHACK
			*((float *)strg.z.val) = (float)fz;
			*((float *)strg.age.val) = (float)fage;
#endif
		} else {
			*((double *)strg.posx.val) = fmod(fposx + 1.0, 1.0);
			*((double *)strg.posy.val) = fmod(fposy + 1.0, 1.0);
			*((double *)strg.posz.val) = fmod(fposz + 1.0, 1.0);
			*((double *)strg.momx.val) = fmomx;
			*((double *)strg.momy.val) = fmomy;
			*((double *)strg.momz.val) = fmomz;
			if (strg.weight.val != NULL)
				*((double *)strg.weight.val) = fweight;
			if (strg.u.val != NULL)
				*((double *)strg.u.val) = fu;
#ifdef METALHACK
			*((double *)strg.z.val) = fz;
			*((double *)strg.age.val) = fage;
#endif
		}
		if (strg.bytes_int == sizeof(uint32_t)) {
			if (strg.id.val != NULL)
				*((uint32_t *)strg.id.val) = (uint32_t)fid;
		} else {
			if (strg.id.val != NULL)
				*((uint64_t *)strg.id.val) = fid;
		}

		/********************************************************\
		 * Step 3: Increment the pointers to the next particle  *
		\********************************************************/
		INCR(strg.posx.val, strg.posx.stride);
		INCR(strg.momx.val, strg.momx.stride);
		INCR(strg.posy.val, strg.posy.stride);
		INCR(strg.momy.val, strg.momy.stride);
		INCR(strg.posz.val, strg.posz.stride);
		INCR(strg.momz.val, strg.momz.stride);
		if (strg.weight.val != NULL)
			INCR(strg.weight.val, strg.weight.stride);
		if (strg.u.val != NULL)
			INCR(strg.u.val, strg.u.stride);
#ifdef METALHACK
		INCR(strg.z.val, strg.z.stride);
		INCR(strg.age.val, strg.age.stride);
#endif
		if (strg.id.val != NULL)
			INCR(strg.id.val, strg.id.stride);

	} /* particle loop ends here*/

	io_logging_msg(log, INT32_C(3),
	               "Done with reading %" PRIu64 " particles from %s",
	               i, f->fname);

	return i;
}

extern uint64_t
io_ares_writepart(io_logging_t log,
                  io_ares_t f,
                  uint64_t pskip,
                  uint64_t pwrite,
                  io_file_strg_struct_t strg)
{
	uint64_t no_part = 0;

	if (!local_write_common(log, f, strg.weight.val, strg.u.val,
	                        pskip, strg.bytes_float,
	                        strg.bytes_int))
		return UINT64_C(0);

	/* Start putting the particles to the file */
	for (no_part=0; no_part<pwrite; no_part++) {
		/* Write the current particle */
		fwrite(strg.posx.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momx.val, strg.bytes_float, 1, f->file);
		fwrite(strg.posy.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momy.val, strg.bytes_float, 1, f->file);
		fwrite(strg.posz.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momz.val, strg.bytes_float, 1, f->file);
		if (strg.weight.val != NULL)
			fwrite(strg.weight.val, strg.bytes_float, 1, f->file);
		if (strg.u.val != NULL)
			fwrite(strg.u.val, strg.bytes_float, 1, f->file);
#ifdef METALHACK
		fwrite(strg.z.val, strg.bytes_float, 1, f->file);
		fwrite(strg.age.val, strg.bytes_float, 1, f->file);
#endif
		if (strg.id.val != NULL)
			fwrite(strg.id.val, strg.bytes_int, 1, f->file);

		/* Advance to next particle */
		INCR(strg.posx.val, strg.posx.stride);
		INCR(strg.momx.val, strg.momx.stride);
		INCR(strg.posy.val, strg.posy.stride);
		INCR(strg.momy.val, strg.momy.stride);
		INCR(strg.posz.val, strg.posz.stride);
		INCR(strg.momz.val, strg.momz.stride);
		if (strg.weight.val != NULL)
			INCR(strg.weight.val, strg.weight.stride);
		if (strg.u.val != NULL)
			INCR(strg.u.val, strg.u.stride);
#ifdef METALHACK
		INCR(strg.z.val, strg.z.stride);
		INCR(strg.age.val, strg.age.stride);
#endif
		if (strg.id.val != NULL)
			INCR(strg.id.val, strg.id.stride);
	}

	/* Update the header */
	f->header->no_part_in_file = (long)(pskip + no_part);
	f->header->bytes_int = strg.bytes_int;
	f->header->bytes_float = strg.bytes_float;
	if (strg.weight.val != NULL)
		f->header->has_mass = UINT64_C(1);
	else
		f->header->has_mass = UINT64_C(0);
	if (strg.u.val != NULL)
		f->header->has_gas = UINT64_C(1);
	else
		f->header->has_gas = UINT64_C(0);

	/* Update logfile */
	io_logging_msg(log, INT32_C(3),
	               "Wrote %" PRIu64 " particles to file",
	               no_part);

	/* Done */
	return no_part;
}

extern uint64_t
io_ares_writepart_ord(io_logging_t log,
                      io_ares_t f,
                      uint64_t pskip,
                      uint64_t pwrite,
                      void *nxt_part,
                      io_file_strg_struct_t strg)
{
	uint64_t no_part = 0;
	ptrdiff_t stride;

	if (!local_write_common(log, f, strg.weight.val, strg.u.val,
	                        pskip, strg.bytes_float,
	                        strg.bytes_int))
		return UINT64_C(0);

	/* Start putting the particles to the file */
	do {
		/* Write the current particle */
		fwrite(strg.posx.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momx.val, strg.bytes_float, 1, f->file);
		fwrite(strg.posy.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momy.val, strg.bytes_float, 1, f->file);
		fwrite(strg.posz.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momz.val, strg.bytes_float, 1, f->file);
		if (strg.weight.val != NULL)
			fwrite(strg.weight.val, strg.bytes_float, 1, f->file);
		if (strg.u.val != NULL)
			fwrite(strg.u.val, strg.bytes_float, 1, f->file);
#ifdef METALHACK
		fwrite(strg.z.val, strg.bytes_float, 1, f->file);
		fwrite(strg.age.val, strg.bytes_float, 1, f->file);
#endif
		if (strg.id.val != NULL)
			fwrite(strg.id.val, strg.bytes_int, 1, f->file);
		no_part++;

		/* Advance to the next particle */
		stride = (char *)*((char **)nxt_part) - (char *)nxt_part;
		nxt_part = (void *)((char *)nxt_part + stride);
		INCR(strg.posx.val, stride);
		INCR(strg.momx.val, stride);
		INCR(strg.posy.val, stride);
		INCR(strg.momy.val, stride);
		INCR(strg.posz.val, stride);
		INCR(strg.momz.val, stride);
		if (strg.weight.val != NULL)
			INCR(strg.weight.val, stride);
		if (strg.u.val != NULL)
			INCR(strg.u.val, stride);
#ifdef METALHACK
		INCR(strg.z.val, stride);
		INCR(strg.age.val, stride);
#endif
		if (strg.id.val != NULL)
			INCR(strg.id.val, stride);
	} while (    (*((char **)nxt_part) != NULL)
	          && (no_part < pwrite) );

	/* Update the header */
	f->header->no_part_in_file = (long)(pskip + no_part);
	f->header->bytes_int = strg.bytes_int;
	f->header->bytes_float = strg.bytes_float;
	if (strg.weight.val != NULL)
		f->header->has_mass = UINT64_C(1);
	else
		f->header->has_mass = UINT64_C(0);
	if (strg.u.val != NULL)
		f->header->has_gas = UINT64_C(1);
	else
		f->header->has_gas = UINT64_C(0);

	/* Update logfile */
	io_logging_msg(log, INT32_C(3),
	               "Wrote %" PRIu64 " particles to file",
	               no_part);

	/* Done */
	return no_part;
}

#undef INCR

extern bool
io_ares_get(io_logging_t log,
            io_ares_t f,
            io_file_get_t what,
            void *res)
{
	if ( (f == NULL) || (f->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_NOPART_IN_FILE:
			*((long *)res) = (long)(f->header->no_part_in_file);
			break;
		case IO_FILE_GET_NOPART:
			*((long *)res) = (long)(f->header->no_part);
			break;
		case IO_FILE_GET_NOVPART:
			*((double *)res) = f->header->no_vpart;
			break;
		case IO_FILE_GET_NOSPECIES:
			*((int *)res) = f->header->no_species;
			break;
		case IO_FILE_GET_BOXSIZE:
			*((double *)res) = f->header->boxsize;
			break;
		case IO_FILE_GET_PMASS:
			*((double *)res) = f->header->pmass;
			break;
		case IO_FILE_GET_ZINITIAL:
			*((double *)res) = 1./(f->header->a_initial) - 1.;
			break;
		case IO_FILE_GET_Z:
			*((double *)res) = 1./(f->header->a_current) - 1.;
			break;
		case IO_FILE_GET_AINITIAL:
			*((double *)res) = f->header->a_initial;
			break;
		case IO_FILE_GET_A:
			*((double *)res) = f->header->a_current;
			break;
		case IO_FILE_GET_OMEGA0:
			*((double *)res) = f->header->omega0;
			break;
		case IO_FILE_GET_OMEGAL:
			*((double *)res) = f->header->lambda0;
			break;
		case IO_FILE_GET_H:
			io_logging_fatal(log,
			                 "The Hubble parameter is not available "
			                 "for ARES files.");
			return false;
		case IO_FILE_GET_DOUBLE:
			*((int *)res) = (f->header->bytes_float == 4) ? 0 : 1;
			break;
		case IO_FILE_GET_MMASS:
			*((int *)res) = (f->header->has_mass == 0) ? 0 : 1;
			break;
		case IO_FILE_GET_NOTSTEP:
			*((int32_t *)res) = INT32_C(1000);
			break;
		case IO_FILE_GET_TSTEP:
			*((double *)res) = f->header->timestep;
			break;
		case IO_FILE_GET_HEADERSTR:
			*((char **)res) = &(f->header->id[0]);
			break;
		case IO_FILE_GET_MINWEIGHT:
			*((double *)res) = f->header->minweight;
			break;
		case IO_FILE_GET_MAXWEIGHT:
			*((double *)res) = f->header->maxweight;
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern bool
io_ares_set(io_logging_t log,
            io_ares_t f,
            io_file_get_t what,
            void *res)
{
	if (f == NULL)
		return false;

	if (f->header == NULL) {
		io_logging_warn(log, INT32_C(3),
		                "File does not have a header yet, "
		                "creating one.");
		f->header = io_ares_header_new(log);
		if (f->header == NULL)
			return false;
	}

	switch (what) {
		case IO_FILE_GET_BOXSIZE:
			f->header->boxsize = *((double *)res);
			break;
		case IO_FILE_GET_PMASS:
			f->header->pmass = *((double *)res);
			break;
		case IO_FILE_GET_Z:
			f->header->a_current = 1./(*((double *)res) + 1.);
			f->header->a_initial = 0.1/(*((double *)res) + 1.);
			break;
		case IO_FILE_GET_A:
			f->header->a_current = *((double *)res);
			f->header->a_initial = *((double *)res)/10.0;
			break;
		case IO_FILE_GET_OMEGA0:
			f->header->omega0 = *((double *)res);
			break;
		case IO_FILE_GET_OMEGAL:
			f->header->lambda0 = *((double *)res);
			break;
		case IO_FILE_GET_H:
			io_logging_fatal(log,
			                 "The Hubble parameter is not available "
			                 "for ARES files.");
			return false;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern void
io_ares_log(io_logging_t log, io_ares_t f)
{
	io_logging_msg(log, INT32_C(5),
	               "Fileobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Filetype:          %s",
	               io_file_typestr(f->ftype));
	io_logging_msg(log, INT32_C(5),
	               "  Filename:          %s",
	               f->fname);
	io_logging_msg(log, INT32_C(5),
	               "  Mode:              %" PRIi8,
	               f->mode);
	io_logging_msg(log, INT32_C(5),
	               "  Swapping:          %" PRIi8,
	               f->swapped);
	io_ares_header_log(log, f->header);

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
inline static io_ares_t
local_openopen(io_logging_t log, io_ares_t f, io_file_mode_t mode)
{
	if (mode == IO_FILE_READ) {
		f->file = fopen(f->fname, IO_FILE_MODE_READ);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for reading.",
			                 f->fname);
			free(f);
			return NULL;
		}
	} else {
		f->file = fopen(f->fname, IO_FILE_MODE_WRITE);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for writing.",
			                 f->fname);
			free(f);
			return NULL;
		}
	}
	
	f->mode = mode;

	return f;
}

inline static void
local_openswapped(io_logging_t log,
                  io_ares_t f,
                  io_file_swap_t swapped)
{
	if (f->mode == IO_FILE_READ) {
		switch (swapped) {
			case IO_FILE_ISNOT_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming unswapped file.");
				break;
			case IO_FILE_IS_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming swapped file.");
				break;
			case IO_FILE_UNKOWN_SWAPPING:
			default:
				io_logging_msg(log, INT32_C(3),
				               "Will find out swap status!");
		}
	}
	if (f->mode == IO_FILE_WRITE) {
		io_logging_msg(log, INT32_C(3),
		               "I don't care about swapping when writing.");
	}

	/* Now set the swapping */
	f->swapped = swapped;

	return;
}

inline static bool
local_write_common(io_logging_t log,
                   io_ares_t f,
                   void *weight,
                   void *u,
                   uint64_t pskip,
                   int32_t bytes_float,
                   int32_t bytes_int)
{
	int32_t partsize;

	/* File sanity checks */
	if (f == NULL) {
		io_logging_warn(log, INT32_C(1),
		                "File object does not exist. Not writing.");
		return false;
	}
	if (f->mode != IO_FILE_WRITE) {
		io_logging_warn(log, INT32_C(1),
		                "File not opened for writing. Not writing.");
		return false;
	}
	if (f->file == NULL) {
		io_logging_warn(log, INT32_C(1),
		                "File claims to be opened, but it is not. "
		                "Not writing.");
		return false;
	}
	if (f->header == NULL) {
		io_logging_warn(log, INT32_C(3),
		                "File does not have a header yet, "
		                "creating one.");
		f->header = io_ares_header_new(log);
		if (f->header == NULL)
			return false;
	}

	/* See if we have to write NULLs to the header part */
	if (fseek(f->file, 0L, SEEK_END) != 0) {
		io_logging_fatal(log, "Could not seek in the file!");
		return false;
	}
	if ((uint64_t)ftell(f->file) < ARES_HEADER_SIZE) {
		int8_t nix = 0;
		int32_t i;

		io_logging_msg(log, INT32_C(3),
		               "Need to write a dummy header");
		rewind(f->file);
		for (i=0; i<ARES_HEADER_SIZE; i++)
			fwrite((void *)&nix,
			       sizeof(int8_t),
			       1,
			       f->file);
		io_logging_msg(log, INT32_C(3),
		               "Wrote a dummy header");
	}

	/* Put the file pointer to the right place */
	partsize = bytes_float * 6 + bytes_int;
	if (weight != NULL)
		partsize += bytes_float;
	if (u != NULL)
		partsize += bytes_float;
	fseek(f->file, ARES_HEADER_SIZE + (pskip*partsize), SEEK_SET);

	return true;
}
