/**
 * \file io_cubep3m.c
 *
 * Provides functions for reading and writing CUBEP3M files.
 */

/**********************************************************************
 *    Includes                                                        *
 **********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <assert.h>
#ifdef WITH_MPI
#  include <mpi.h>
#endif

#include "io_cubep3m.h"
#include "io_cubep3m_header.h"
#include "io_util.h"


/**********************************************************************
 *    Local defines, structure definitions and typedefs               *
 **********************************************************************/
#define SKIP(f) {fseek(f, 4L, SEEK_CUR); }

#define H0 100.               // not known as param.h is not included

// ptype values as used by AHF
// NOTE: this has to be identical to the definitions found in param.h
#define PGAS     (0.0)              /* identifier for gas particles; has to be exactly 0.0!!!! */
#define PDM      (-1.0)             /* identifier for dm particles; whatever negative value */
#define PSTAR    (-4.0)             /* identifier for star particles; whatever negative value */
#define PDMbndry (-5.0)             /* identifier for star particles; whatever negative value */

/**********************************************************************
 *    Prototypes of local functions                                   *
 **********************************************************************/

/**
 * \brief Helper function to open the file
 *
 * This function is supposed to get inlined anyway. Makes
 * io_cubep3m_open more readable.
 *
 * \param log   The logging object.
 * \param f     The CUBEP3M file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_cubep3m_t
local_openopen(io_logging_t log, io_cubep3m_t f, io_file_mode_t mode);


/**
 * \brief  Tries to identify the CubePM file (checks if chunked or not).
 *
 * \param[in,out]  log
 *                    The logging module.
 * \param[in,out]  f
 *                    The file to work with.
 * \param[in]      swapped
 *                    The swapping state.  The function will try to
 *                    autodetect the swapping but only set the swap state
 *                    from the file if \c swapped is
 *                    #IO_FILE_UNKOWN_SWAPPING otherwise the provided value
 *                    is passed through.
 *
 * \return  Returns nothing.
 */
inline static void
local_open_identify_file(io_logging_t   log,
                         io_cubep3m_t   f,
                         io_file_swap_t swapped);


/**
 * \brief  Given a magic number, find out if the file it was read from is
 *         swapped or not.
 *
 * \param[in]  magicNumber
 *                The magic number read from the file.  It is supposed to be
 *                a 1, but might be swapped.
 *
 * \return  Returns the swap state of the file or an unknown swapping if the
 *          status could not be determined.
 */
inline static io_file_swap_t
local_getSwapStateFile(int32_t magicNumber);

inline static void
local_print_offset_information(io_logging_t log, io_cubep3m_t f);


/**********************************************************************
 *    Implementation of global functions                              *
 **********************************************************************/
extern io_cubep3m_t
io_cubep3m_open(io_logging_t   log,
                char           *fname,
                io_file_swap_t swapped,
                io_file_mode_t mode,
                uint32_t       reader)
{
	io_cubep3m_t f;

	/* Get memory for the structure */
	f = (io_cubep3m_t)malloc(sizeof(io_cubep3m_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log, "io_cubep3m structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Store the filename (must be done BEFORE reading the header) */
	f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
	if (f->fname == NULL) {
		io_logging_memfatal(log, "filename of CUBEP3Mfile");
		free(f);
		return NULL;
	}
	strncpy(f->fname, fname, strlen(fname) + 1);

	f->ftype   = IO_FILE_CUBEP3M;
	f->header  = NULL;
	f->no_part = UINT64_C(0);

#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	if (f->size >= reader) {
		MPI_Comm_split(MPI_COMM_WORLD, 1,
		               f->rank, &(f->mycomm));
		MPI_Comm_size(f->mycomm, &(f->size_mycomm));
		MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	} else {
		f->mycomm      = MPI_COMM_NULL;
		f->size_mycomm = -1;
		f->rank_mycomm = -1;
	}
#endif

	if (local_openopen(log, f, mode) == NULL) {
		free(f->fname);
		free(f);
		return NULL;
	}

	if (f->mode = IO_FILE_READ) {
		local_open_identify_file(log, f, swapped);
	} else {
		f->isChunked = false;
		assert(swapped != IO_FILE_UNKOWN_SWAPPING);
		f->swapped   = swapped;
	}

	return f;
} /* io_cubep3m_open */

extern void
io_cubep3m_close(io_logging_t log,
                 io_cubep3m_t *f)
{
	/* Catch NULLs */
	if ((f == NULL) || (*f == NULL))
		return;

	/* Put header to the file if necessary */
	if (((*f)->mode == IO_FILE_WRITE)
	    && ((*f)->header != NULL)) {
		io_cubep3m_header_write(log, (*f)->header, *f);
	}

	/* Close */
	if ((*f)->header != NULL)
		io_cubep3m_header_del(log, &((*f)->header));
	if ((*f)->fname != NULL)
		free((*f)->fname);
#ifdef WITH_MPI
	if ((*f)->mycomm != MPI_COMM_NULL)
		MPI_Comm_free(&((*f)->mycomm));
#endif

	/* Actually close the file */
	if ((*f)->file != NULL)
		fclose((*f)->file);

	/* Cleaning */
	free(*f);
	*f = NULL;

	return;
}

extern void
io_cubep3m_init(io_logging_t log,
                io_cubep3m_t f)
{
	if (f == NULL)
		return;

	if (f->header != NULL) {
		io_logging_warn(log, INT32_C(1),
		                "Already have the header information! Rereading.");
		io_cubep3m_header_del(log, &(f->header));
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
	f->header = io_cubep3m_header_get(log, f);
	io_logging_msg(log, INT32_C(5),
	               "Done with initializing file object from %s",
	               f->fname);

	if (f->isChunked) {
		f->no_part = f->header->np_in_chunk_file;
		local_print_offset_information(log, f);
	} else {
		f->no_part = f->header->np_local;
	}

	return;
} /* io_cubep3m_init */

extern uint64_t
io_cubep3m_readpart(io_logging_t          log,
                    io_cubep3m_t          f,
                    uint64_t              pskip,
                    uint64_t              pread,
                    io_file_strg_struct_t strg)
{
	uint64_t particles_read, tmp;

	/*
	 * First read the particles unscaled.  Not quite true, the particles
	 * will be offsetted according to their position in the file grid
	 * and periodically wrapped to be in [0..ngrid[.
	 */
	particles_read = io_cubep3m_readpart_raw(log, f, pskip, pread, strg);

	if (particles_read != pread)
		return UINT64_C(0);

	/* And do the scaling */
	tmp = io_cubep3m_scale_particles(log,
	                                 f->header->boxsize,
	                                 f->header->a,
	                                 f->header->lunit,
	                                 f->header->vunit,
	                                 particles_read,
	                                 strg);
	if (tmp != particles_read) {
		return tmp;
	}

	/* Wow, we are done! */
	return particles_read;
} /* io_cubep3m_readpart */

extern uint64_t
io_cubep3m_readpart_raw(io_logging_t          log,
                        io_cubep3m_t          f,
                        uint64_t              pskip,
                        uint64_t              pread,
                        io_file_strg_struct_t strg)
{
	uint64_t i, pid;
	uint32_t bytes_file = 4;  // Hard-coding to use only float
	float    fposx, fposy, fposz;
	float    fmomx, fmomy, fmomz;
	double   boxsizeFileInternal;

	/* Check if we actually have to do something */
	if ((f == NULL) || (f->header == NULL))
		return UINT64_C(0);

#ifdef CUBEP3M_WITH_PIDS
  if (strg.bytes_int == sizeof(uint32_t)) {
    io_logging_fatal(log, "You are planning to read 64bit particle ids, but the internal id storage is only designed for 32bit ids. Aborting in %s.", __func__);
  }
#endif
  
	/* Set extreme position detectors */
	f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e40;
	f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e40;

	/* Make sure that we are at the beginning of the file */
	rewind(f->file);

	/*******************************************************************
	 *  Skip HEADER                                                    *
	 *******************************************************************/
	if (f->isChunked)
		fseek(f->file, CUBEP3M_HEADER_CHUNK_SIZE, SEEK_CUR);
	else
		fseek(f->file, CUBEP3M_HEADER_SIZE, SEEK_CUR);


	/*******************************************************************
	 *  Parallel I/O: what fraction should the current CPU read?       *
	 *******************************************************************/

	/* Set the number of particles to loop over correctly */
	io_logging_msg(log, INT32_C(3),
	               "Asked to read %" PRIu64 " and to skip %" PRIu64
	               " particles. Checking those numbers.",
	               pread, pskip);
	if (pskip > f->no_part) {
		pskip = f->no_part;
		io_logging_msg(log, INT32_C(3),
		               "Cannot skip more than there is, will now "
		               "only skip %" PRIu64 " particles.",
		               pskip);
	}
	if (f->no_part - pskip < pread) {
		pread = f->no_part - pskip;
		io_logging_msg(log, INT32_C(3),
		               "Cannot read more than there is left after "
		               "skipping. Will only read %" PRIu64
		               "particles.", pread);
	}

	/*******************************************************************
	 *  SKIP to the first particle to be read                          *
	 *******************************************************************/
#ifdef CUBEP3M_WITH_PIDS
	// 7 variables: x,y,z,vx,vy,vz,pid
  fseek(f->file, (bytes_file * 6 + sizeof(uint64_t)) * pskip, SEEK_CUR); // pids are using uint64_t as storage
#else /* CUBEP3M_WITH_PIDS */
	// 6 variables: x,y,z,vx,vy,vz
	fseek(f->file, bytes_file * 6 * pskip, SEEK_CUR);
#endif /* CUBEP3M_WITH_PIDS */


	/*******************************************************************
	 *  READ particles                                                 *
	 *******************************************************************/

	boxsizeFileInternal = (double)(f->header->ngrid);
	for (i = pskip; i < pskip + pread; i++) {
		/* Read positions */
		io_util_readfloat(f->file, &fposx, f->swapped);
		io_util_readfloat(f->file, &fposy, f->swapped);
		io_util_readfloat(f->file, &fposz, f->swapped);

		/* Read velocities */
		io_util_readfloat(f->file, &fmomx, f->swapped);
		io_util_readfloat(f->file, &fmomy, f->swapped);
		io_util_readfloat(f->file, &fmomz, f->swapped);

		/* Adding the box-offset and clip to [0..box] */
		fposx = (float)fmod(fposx + f->header->offset[0] + (double)(f->header->ngrid), boxsizeFileInternal);
		fposy = (float)fmod(fposy + f->header->offset[1] + (double)(f->header->ngrid), boxsizeFileInternal);
		fposz = (float)fmod(fposz + f->header->offset[2] + (double)(f->header->ngrid), boxsizeFileInternal);

		/* Store valuable values */
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.posx.val) = (float)fposx;
			*((float *)strg.posy.val) = (float)fposy;
			*((float *)strg.posz.val) = (float)fposz;
			*((float *)strg.momx.val) = (float)fmomx;
			*((float *)strg.momy.val) = (float)fmomy;
			*((float *)strg.momz.val) = (float)fmomz;
		} else {
			*((double *)strg.posx.val) = (double)fposx;
			*((double *)strg.posy.val) = (double)fposy;
			*((double *)strg.posz.val) = (double)fposz;
			*((double *)strg.momx.val) = (double)fmomx;
			*((double *)strg.momy.val) = (double)fmomy;
			*((double *)strg.momz.val) = (double)fmomz;
		}
#ifdef CUBEP3M_WITH_PIDS
		/* Read the ID */
		if (strg.id.val != NULL) {
      io_util_readuint64(f->file, &pid, f->swapped);
      *((uint64_t *)strg.id.val) = pid;
		}
    else {
      io_logging_fatal(log, "You are planning to read particle ids, but there is no memory allocated for it. Aborting in %s.", __func__);
    }
#else /* CUBEP3M_WITH_PIDS */
		/* Invent an ID (NOT unique for multiple files) */
		if (strg.id.val != NULL) {
			if (strg.bytes_int == sizeof(uint32_t))
				*((uint32_t *)strg.id.val) = (uint32_t)i;
			else
				*((uint64_t *)strg.id.val) = i;
		}
#endif /* CUBEP3M_WITH_PIDS */

		/* Increment the pointers to the next particle */
		strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride);
		strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride);
		strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride);
		strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride);
		strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride);
		strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride);
		if (strg.id.val != NULL) strg.id.val = (void *)(((char *)strg.id.val) + strg.id.stride);

		/* Detect extreme positions */
		if (isless(fposx, f->minpos[0]))			f->minpos[0] = fposx;
		if (isless(fposy, f->minpos[1]))			f->minpos[1] = fposy;
		if (isless(fposz, f->minpos[2]))			f->minpos[2] = fposz;
		if (isgreater(fposx, f->maxpos[0]))		f->maxpos[0] = fposx;
		if (isgreater(fposy, f->maxpos[1]))		f->maxpos[1] = fposy;
		if (isgreater(fposz, f->maxpos[2]))		f->maxpos[2] = fposz;
	}

	/* Return the number of particles read */
	return pread;
} /* io_cubep3m_readpart_raw */

extern bool
io_cubep3m_get(io_logging_t  log,
               io_cubep3m_t  f,
               io_file_get_t what,
               void          *res)
{
	if ((f == NULL) || (f->header == NULL))
		return false;

	switch (what) {
	case IO_FILE_GET_NOPART_IN_FILE:
		*((long *)res) = (long)f->no_part;
		break;
	case IO_FILE_GET_NOPART:
		if (f->isChunked)
			*((long *)res) = (long)f->header->nptotal;
		else
			*((long *)res) = (long)f->header->tot_np_in_chunk;
		break;
	case IO_FILE_GET_NOVPART:
		if (f->isChunked)
			*((double *)res) = (double)f->header->tot_np_in_chunk;
		else
			*((double *)res) = (double)f->header->nptotal;
		break;
	case IO_FILE_GET_NOSPECIES:
		*((int *)res) = (int)1;
		break;
	case IO_FILE_GET_BOXSIZE:
		*((double *)res) = f->header->boxsize;
		break;
	case IO_FILE_GET_PMASS:
		*((double *)res) = f->header->munit * f->header->mass_p;
		break;
	case IO_FILE_GET_ZINITIAL:
		io_logging_warn(log, INT32_C(1),
		                "zinitial is not set in a CUBEP3M file, "
		                "using current redshift");
	case IO_FILE_GET_Z:
		*((double *)res) = 1. / f->header->a - 1.;
		break;
	case IO_FILE_GET_AINITIAL:
		io_logging_warn(log, INT32_C(1),
		                "ainitial is not set in a CUBEP3M file, "
		                "using current expansion.");
	case IO_FILE_GET_A:
		*((double *)res) = f->header->a;
		break;
	case IO_FILE_GET_OMEGA0:
		*((double *)res) = f->header->omega0;
		break;
	case IO_FILE_GET_OMEGAL:
		*((double *)res) = f->header->lambda0;
		break;
	case IO_FILE_GET_H:
		io_logging_warn(log, INT32_C(1),
		                "CUBEP3M files don't store the Hubble parameter."
		                "Setting to 100.0");
		*((double *)res) = H0;
		break;
	case IO_FILE_GET_DOUBLE:
		io_logging_warn(log, INT32_C(1),
		                "CUBEP3M files don't store the use of "
		                "double precision. Assuming it is not "
		                "double precision.");
		*((int *)res) = 0;
		break;
	case IO_FILE_GET_MMASS:
		*((int *)res) = 0;
		break;
	case IO_FILE_GET_NOTSTEP:
		*((int32_t *)res) = f->header->nts;
		break;
	case IO_FILE_GET_TSTEP:
		io_logging_warn(log, INT32_C(1),
		                "CUBEP3M files don't store the timestep. "
		                "Setting to 0.0.");
		*((double *)res) = 0.0;
		break;
	case IO_FILE_GET_HEADERSTR:
		io_logging_warn(log, INT32_C(1),
		                "CUBEP3M files don't have a header string. "
		                "Using a dummy one.");
		*((char **)res) = "No header string.";
		break;
	case IO_FILE_GET_MINWEIGHT:
	case IO_FILE_GET_MAXWEIGHT:
		*((double *)res) = 1.0;
		break;
	default:
		io_logging_fatal(log, "Requesting something unkown in %s.",
		                 __func__);
		return false;
	} /* switch */

	return true;
} /* io_cubep3m_get */

extern void
io_cubep3m_log(io_logging_t log, io_cubep3m_t f)
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
	               "  isChunked:            %s",
	               f->isChunked ? "true" : "false");
	io_logging_msg(log, INT32_C(5),
	               "  No. particles:        %" PRIu64,
	               f->no_part);
	io_cubep3m_header_log(log, f->header);

	return;
} /* io_cubep3m_log */

extern uint64_t
io_cubep3m_scale_particles(io_logging_t          log,
                           double                boxsize,
                           double                expansion,
                           double                lunit,
                           double                vunit,
                           uint64_t              particles_read,
                           io_file_strg_struct_t strg)
{
	uint64_t i;
	double   scale_pos, scale_mom;

	/*******************************************************************\
	 * Scaling Units                                                   *
	 *******************************************************************/
	scale_pos = lunit / boxsize;
	scale_mom = vunit * expansion / (H0 * boxsize);

	io_logging_msg(log, INT32_C(3),
	               "Scaling by:  positions:  %g", scale_pos);
	io_logging_msg(log, INT32_C(3),
	               "             velocities: %g", scale_mom);

	/* Define the actual scaling calls type independent */
#define SCALE_CALL(type) {                                                \
		*((type *)strg.posx.val)  = (type)fmod((*((type *)strg.posx.val)) \
		                                       * scale_pos, 1.0);         \
		*((type *)strg.posy.val)  = (type)fmod((*((type *)strg.posy.val)) \
		                                       * scale_pos, 1.0);         \
		*((type *)strg.posz.val)  = (type)fmod((*((type *)strg.posz.val)) \
		                                       * scale_pos, 1.0);         \
		*((type *)strg.momx.val) *= (type)(scale_mom);                    \
		*((type *)strg.momy.val) *= (type)(scale_mom);                    \
		*((type *)strg.momz.val) *= (type)(scale_mom);                    \
}

#define MOVE_TO_NEXT {                                                       \
		strg.posx.val                                                        \
		    = (void *)(((char *)strg.posx.val) + strg.posx.stride);          \
		strg.posy.val                                                        \
		    = (void *)(((char *)strg.posy.val) + strg.posy.stride);          \
		strg.posz.val                                                        \
		    = (void *)(((char *)strg.posz.val) + strg.posz.stride);          \
		strg.momx.val                                                        \
		    = (void *)(((char *)strg.momx.val) + strg.momx.stride);          \
		strg.momy.val                                                        \
		    = (void *)(((char *)strg.momy.val) + strg.momy.stride);          \
		strg.momz.val                                                        \
		    = (void *)(((char *)strg.momz.val) + strg.momz.stride);          \
		if (strg.weight.val != NULL)                                         \
			strg.weight.val                                                  \
			    = (void *)(((char *)strg.weight.val) + strg.weight.stride);  \
}

	/* Do the scaling depending on the storage type */
	if (strg.bytes_float == sizeof(float)) {
		for (i = 0; i < particles_read; i++) {
			SCALE_CALL(float);
			MOVE_TO_NEXT;
		}
	} else if (strg.bytes_float == sizeof(double)) {
		for (i = 0; i < particles_read; i++) {
			SCALE_CALL(double);
			MOVE_TO_NEXT;
		}
	} else {
		io_logging_fatal(log,
		                 "Don't know which floating point types "
		                 "has %" PRIi32 " bytes. Aborting read.",
		                 strg.bytes_float);
		return UINT64_C(0);
	}

	/* Clean the macro away */
#undef SCALE_CALL
#undef MOVE_TO_NEXT

	/* And we are done */
	return particles_read;
} /* io_cubep3m_scale_particles */

/**********************************************************************\
 *    Implementation of local functions                               *
 \**********************************************************************/

inline static io_cubep3m_t
local_openopen(io_logging_t log, io_cubep3m_t f, io_file_mode_t mode)
{
	if (mode == IO_FILE_READ) {
		f->file = fopen(f->fname, IO_FILE_MODE_READ);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for reading.",
			                 f->fname);
			return NULL;
		}
	} else {
		f->file = fopen(f->fname, IO_FILE_MODE_WRITE);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for writing.",
			                 f->fname);
			return NULL;
		}
	}

	f->mode = mode;

	return f;
}

inline static void
local_open_identify_file(io_logging_t   log,
                         io_cubep3m_t   f,
                         io_file_swap_t swapped)
{
	struct magic {
		int32_t number;
		char    string[17];
	} magic;

	rewind(f->file);
	fread(&(magic.number), sizeof(int32_t), 1, f->file);
	fread(&(magic.string), sizeof(char), 16, f->file);
	magic.string[16] = '\0';
	io_logging_msg(log, INT32_C(3), "Magic number: %" PRIi32, magic.number);
	io_logging_msg(log, INT32_C(3), "Magic string: %s", magic.string);
	if (strncmp(magic.string, "CHUNKED CUBEPM  ", 16) == 0) {
		io_file_swap_t swapStateFile;
		io_logging_msg(log, INT32_C(3), "File seems to be chunked.");
		f->isChunked  = true;
		swapStateFile = local_getSwapStateFile(magic.number);
		if (swapped == IO_FILE_UNKOWN_SWAPPING) {
			f->swapped = swapStateFile;
		} else {
			f->swapped = swapped;
			if (swapped != swapStateFile) {
				io_logging_warn(log, INT32_C(1),
				                "Detected swap status is not the provided "
				                "one.  Will trust the user.");
			}
		}
	}
  else {
		f->isChunked  = false;
  }
	rewind(f->file);
}

inline static io_file_swap_t
local_getSwapStateFile(int32_t magicNumber)
{
	io_file_swap_t swapStateFile;

	if (magicNumber != 1) {
		io_util_sexchange(&(magicNumber), sizeof(int32_t));
		if (magicNumber != 1)
			swapStateFile = IO_FILE_UNKOWN_SWAPPING;
		else
			swapStateFile = IO_FILE_IS_SWAPPED;
	} else {
		swapStateFile = IO_FILE_ISNOT_SWAPPED;
	}

	return swapStateFile;
}

inline static void
local_print_offset_information(io_logging_t log, io_cubep3m_t f)
{
	double xfac = f->header->lunit;

	io_logging_msg(log, INT32_C(1),
	               "Offset information:",
	               f->fname);
	io_logging_msg(log, INT32_C(1),
	               "  Full data start: %15.12e %15.12e %15.12e",
	               f->header->chunk_start_full_data[0] * xfac,
	               f->header->chunk_start_full_data[1] * xfac,
	               f->header->chunk_start_full_data[2] * xfac);
	io_logging_msg(log, INT32_C(1),
	               "  Full data end  : %15.12e %15.12e %15.12e",
	               f->header->chunk_end_full_data[0] * xfac,
	               f->header->chunk_end_full_data[1] * xfac,
	               f->header->chunk_end_full_data[2] * xfac);
	io_logging_msg(log, INT32_C(1),
	               "  Real data start: %15.12e %15.12e %15.12e",
	               f->header->chunk_start_real_data[0] * xfac,
	               f->header->chunk_start_real_data[1] * xfac,
	               f->header->chunk_start_real_data[2] * xfac);
	io_logging_msg(log, INT32_C(1),
	               "  Real data end  : %15.12e %15.12e %15.12e",
	               f->header->chunk_end_real_data[0] * xfac,
	               f->header->chunk_end_real_data[1] * xfac,
	               f->header->chunk_end_real_data[2] * xfac);
	io_logging_msg(log, INT32_C(
	                   1),
	               "  Catalogue must be shifted by:  %15.12e %15.12e %15.12e",
	               f->header->chunk_offset[0] * xfac,
	               f->header->chunk_offset[1] * xfac,
	               f->header->chunk_offset[2] * xfac);
	io_logging_msg(log, INT32_C(1),
	               "Offset information:",
	               f->fname);
} /* local_print_offset_information */
