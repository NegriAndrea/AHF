/**
 * \file io_ascii.c
 *
 * Provides functions for reading and writing ASCII files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "io_ascii.h"
#include "io_ascii_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/
/** The line buffer when reading from the file */
#define LOCAL_ASCII_BUFFERSIZE 512


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/

/**
 * \brief Helper function to open the file
 *
 * This function is supposed to get inlined anyway. Makes
 * io_ascii_open more readable.
 *
 * \param log   The logging object.
 * \param f     The ASCII file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_ascii_t
local_openopen(io_logging_t log, io_ascii_t f, io_file_mode_t mode);


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_ascii_t
io_ascii_open(io_logging_t log,
              char *fname,
              io_file_mode_t mode,
              uint32_t reader)
{
	io_ascii_t f;

	/* Get memory for the structure */
	f = (io_ascii_t)malloc(sizeof(io_ascii_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_ascii structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Store the filename */
	f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
	if (f->fname == NULL) {
		io_logging_memfatal(log, "filename of ASCII file");
		free(f);
		return NULL;
	}
	strncpy(f->fname, fname, strlen(fname)+1);

	/* Okay, we are an ASCII file */
	f->ftype = IO_FILE_ASCII;

	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	MPI_Comm_split(MPI_COMM_WORLD, 1,
	               f->rank, &(f->mycomm));
	MPI_Comm_size(f->mycomm, &(f->size_mycomm));
	MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
#	endif

	/* Try to open the file */
	if (local_openopen(log, f, mode) == NULL)
		return NULL;

	/* Nothing for the header for now */
	f->header = NULL;

	/* Initialise the rest to safe parameters */
	f->minweight = 1e40;
	f->maxweight = 0.0;

	return f;
}

extern void
io_ascii_close(io_logging_t log,
               io_ascii_t *f)
{
	/* Catch NULLs */
	if (f == NULL || *f == NULL)
		return;

	/* Put the header to the file if necessary */
	if (    ((*f)->mode == IO_FILE_WRITE)
	     && ((*f)->header != NULL)) {
		io_ascii_header_write(log, (*f)->header, *f);
	}

	/* Close */
	if ((*f)->file != NULL)
		fclose((*f)->file);
	
	/* Clean up */
	if ((*f)->header != NULL)
		io_ascii_header_del(log, &((*f)->header));
	if ((*f)->fname != NULL)
		free((*f)->fname);
#	ifdef WITH_MPI
	MPI_Comm_free(&((*f)->mycomm));
#	endif
	free(*f);
	*f = NULL;

	return;
}

extern void
io_ascii_init(io_logging_t log,
              io_ascii_t f)
{
	if (f->mode == IO_FILE_READ) {
		io_logging_msg(log, INT32_C(5),
		               "Starting to initialize file object from %s",
		               f->fname);
		f->header = io_ascii_header_get(log, f);
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
io_ascii_readpart(io_logging_t log,
                  io_ascii_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg)
{
	uint64_t i;
	float dummy;
	double fposx, fposy, fposz;
	double fmomx, fmomy, fmomz;
	double fweight;
	double x_fac, v_fac, m_fac;
	char buffer[LOCAL_ASCII_BUFFERSIZE];

	/* Check if we actually need to do something */
	if (    (f == NULL)
		 || ( f->header == NULL)
	     || (f->header->no_part <= 0)
	     || (((uint64_t)(f->header->no_part)) < pskip) )
		return UINT64_C(0);

	/* Show what is going to happen */
	if ( (f->header->multi_mass == 1)) {
		if (strg.weight.val == NULL) {
			io_logging_warn(log, INT32_C(3), "   weight (read, ignored)");
		}
	} else {
		if (strg.weight.val != NULL)
			io_logging_warn(log, INT32_C(3), "   weight (set to 1.0)");
	}
	if (strg.id.val != NULL)
		io_logging_warn(log, INT32_C(3),
		                "   id     (set to position in file)");
	if (strg.u.val != NULL)
		io_logging_warn(log, INT32_C(3), "   u      (set to 0.0)");

	/* Set the number of particles to loop over correctly */
	if ( (uint64_t)(f->header->no_part) - pskip < pread)
		pread = (uint64_t)(f->header->no_part) - pskip;

	io_logging_msg(log, INT32_C(3),
	               "Starting to read %" PRIu64 " particles, "
	               "skipping %" PRIu64,
	               pread, pskip);

	/* Position the file at the right spot */
	rewind(f->file);
	do {
		fgets(buffer, LOCAL_ASCII_BUFFERSIZE, f->file);
	} while (buffer[0] == '#');
	for (i=0; i<pskip; i++)
		fgets(buffer, LOCAL_ASCII_BUFFERSIZE, f->file);
	
	/* Set the scaling factors */
	x_fac = 1./f->header->boxsize;
	v_fac = f->header->a_current/100.*x_fac;
	m_fac = 1./f->header->pmass;

	/*
	 * Needed so that min and max weight are correct (not touched in the
	 * reading when single mass, hence set here)
	 */
	if (f->header->multi_mass == 0) {
		f->minweight = f->maxweight = 1.0;
	}

	/* Loop over all particles to be read */
	for (i=0; i<pread; i++) {
		/**************************************************\
		 * Step 1: Read the particle data from the buffer *
		\**************************************************/
		if (f->header->multi_mass == 1) {
			sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf",
			       &fposx, &fposy, &fposz,
			       &fmomx, &fmomy, &fmomz,
			       &fweight);
			fweight *= m_fac;
			if (isless(fweight, f->minweight))
				f->minweight = fweight;
			else if (isgreater(fweight, f->maxweight))
				f->maxweight = fweight;
		} else {
			sscanf(buffer, "%lf %lf %lf %lf %lf %lf",
			       &fposx, &fposy, &fposz,
			       &fmomx, &fmomy, &fmomz);
			fweight = 1.;
		}

		/************************************************\
		 * Step 2: Store the particle in the array      *
		\************************************************/
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.posx.val) = (float)(fposx*x_fac);
			*((float *)strg.posy.val) = (float)(fposy*x_fac);
			*((float *)strg.posz.val) = (float)(fposz*x_fac);
			*((float *)strg.momx.val) = (float)(fmomx*v_fac);
			*((float *)strg.momy.val) = (float)(fmomy*v_fac);
			*((float *)strg.momz.val) = (float)(fmomz*v_fac);
			if (strg.weight.val != NULL)
				*((float *)strg.weight.val) = (float)fweight;
		} else {
			*((double *)strg.posx.val) = fposx*x_fac;
			*((double *)strg.posy.val) = fposy*x_fac;
			*((double *)strg.posz.val) = fposz*x_fac;
			*((double *)strg.momx.val) = fmomx*v_fac;
			*((double *)strg.momy.val) = fmomy*v_fac;
			*((double *)strg.momz.val) = fmomz*v_fac;
			if (strg.weight.val != NULL)
				*((double *)strg.weight.val) = fweight;
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

		/********************************************************\
		 * Step 4: Get the next line                            *
		\********************************************************/
		fgets(buffer, LOCAL_ASCII_BUFFERSIZE, f->file);
	} /* particle loop ends here*/

	/* Loop again if required to set the ID */
	if (strg.id.val != NULL) {
		if (strg.bytes_int == sizeof(uint64_t)) {
			for (i=0; i<pread; i++) {
				*((uint64_t *)strg.id.val) = pskip + i;
				INCR(strg.id.val, strg.id.stride);
			}
		} else {
			for (i=0; i<pread; i++) {
				*((uint32_t *)strg.id.val) = (uint32_t)(pskip + i);
				INCR(strg.id.val, strg.id.stride);
			}
		}
	}
	/* Loop again if required to set the interal energy to 0.0 */
	if (strg.u.val != NULL) {
		if (strg.bytes_int == sizeof(float)) {
			for (i=0; i<pread; i++) {
				*((float *)strg.u.val) = 0.0;
				INCR(strg.u.val, strg.u.stride);
			}
		} else {
			for (i=0; i<pread; i++) {
				*((double *)strg.u.val) = 0.0;
				INCR(strg.u.val, strg.u.stride);
			}
		}
	}

	io_logging_msg(log, INT32_C(3),
	               "Done with reading %" PRIu64 " particles from %s",
	               i, f->fname);

	return i;
}

extern uint64_t
io_ascii_writepart(io_logging_t log,
                   io_ascii_t f,
                   uint64_t pskip,
                   uint64_t pwrite,
                   io_file_strg_struct_t strg)
{
	/* THIS IS NOT IMPLEMENTED YET */

	/* Done */
	return pwrite;
}

extern uint64_t
io_ascii_writepart_ord(io_logging_t log,
                       io_ascii_t f,
                       uint64_t pskip,
                       uint64_t pwrite,
                       void *nxt_part,
                       io_file_strg_struct_t strg)
{
	/* THIS IS NOT IMPLEMENTED YET */

	/* Done */
	return pwrite;
}

#undef INCR

extern bool
io_ascii_get(io_logging_t log,
             io_ascii_t f,
             io_file_get_t what,
             void *res)
{
	if ( (f == NULL) || (f->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_NOPART_IN_FILE:
		case IO_FILE_GET_NOPART:
			*((long *)res) = (long)(f->header->no_part);
			break;
		case IO_FILE_GET_NOVPART:
			*((double *)res) = f->header->total_mass/f->header->pmass;
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
			                 "for ASCII files.");
			return false;
		case IO_FILE_GET_DOUBLE:
#ifdef DOUBLE
			*((int *)res) = 1;
#else
			*((int *)res) = 0;
#endif
			break;
		case IO_FILE_GET_MMASS:
			*((int *)res) = f->header->multi_mass;
			break;
		case IO_FILE_GET_NOTSTEP:
			*((int32_t *)res) = (int32_t)(f->header->no_timestep);
			break;
		case IO_FILE_GET_TSTEP:
			*((double *)res) = 0.0;
			break;
		case IO_FILE_GET_HEADERSTR:
			*((char **)res) = &(f->header->header[0]);
			break;
		case IO_FILE_GET_MINWEIGHT:
			*((double *)res) = f->minweight;
			break;
		case IO_FILE_GET_MAXWEIGHT:
			*((double *)res) = f->maxweight;
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern bool
io_ascii_set(io_logging_t log,
             io_ascii_t f,
             io_file_get_t what,
             void *res)
{
	if (f == NULL)
		return false;

	if (f->header == NULL) {
		io_logging_warn(log, INT32_C(3),
		                "File does not have a header yet, "
		                "creating one.");
		f->header = io_ascii_header_new(log);
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
			                 "for ASCII files.");
			return false;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern void
io_ascii_log(io_logging_t log, io_ascii_t f)
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
	               "  minweight:         %g",
	               f->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  maxweight:         %g",
	               f->maxweight);
	io_ascii_header_log(log, f->header);

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
inline static io_ascii_t
local_openopen(io_logging_t log, io_ascii_t f, io_file_mode_t mode)
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
