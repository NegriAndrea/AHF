/**
 * \file io_mcubep3m.c
 *
 * Provides functions for reading and writing Multiple CubeP3M files.
 */


/**********************************************************************
 *    Includes                                                        *
 **********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "io_mcubep3m.h"
#include "io_cubep3m.h"
#include "io_util.h"


/**********************************************************************
 *    Local defines, structure definitions and typedefs               *
 **********************************************************************/


/**********************************************************************
 *    Prototypes of local functions                                   *
 **********************************************************************/


/**********************************************************************
 *    Implementation of global functions                              *
 **********************************************************************/
extern io_mcubep3m_t
io_mcubep3m_open(io_logging_t   log,
                 char           *fname,
                 io_file_swap_t swapped,
                 io_file_mode_t mode,
                 uint32_t       reader)
{
	int32_t       i;
	io_mcubep3m_t f;
	char          **fnames;

	/* XXX THIS IS CURRENTLY ONLY FOR READING! */
	if (mode != IO_FILE_READ)
		return NULL;

	/* Get memory for the structure */
	f = (io_mcubep3m_t)malloc(sizeof(io_mcubep3m_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log, "io_mcubep3m structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Okay, we are a Multiple Gadget file */
	f->ftype = IO_FILE_MCUBEP3M;

	/* And we can just copy in the parallel information */
#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	MPI_Comm_split(MPI_COMM_WORLD, 1,
	               f->rank, &(f->mycomm));
	MPI_Comm_size(f->mycomm, &(f->size_mycomm));
	MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
#endif

	/* Split the filename in path and stem */
	f->path = NULL;
	f->stem = NULL;
	if (io_util_split_pathfname(fname, &(f->path), &(f->stem))
	    == NULL) {
		io_logging_fatal(log,
		                 "Could not split %s in path and filename.",
		                 fname);
		free(f);
		return NULL;
	}
	io_logging_msg(log, INT32_C(1),
	               "Will look in %s for %s.",
	               f->path, f->stem);

	/* Get the filenames */
	f->numfiles = io_util_findfiles(f->path, f->stem, "%i", ".dat",
	                                &fnames);
	if (f->numfiles <= 0) {
		io_logging_fatal(log,
		                 "Could not open anything starting with %s "
		                 "in %s.",
		                 f->stem, f->path);
		free(f->stem);
		free(f->path);
		free(f);
		return NULL;
	}

	/* Glue the files into the MCubeP3M structure */
	f->files = (io_cubep3m_t *)malloc(sizeof(io_cubep3m_t) * (f->numfiles));
	if (f->files == NULL) {
		for (i = 0; i < f->numfiles; i++)
			free(fnames[i]);
		free(fnames);
		free(f->stem);
		free(f->path);
		free(f);
		return NULL;
	}
	for (i = 0; i < f->numfiles; i++) {
#ifdef WITH_MPI
		f->files[i] = io_cubep3m_open(log, fnames[i], swapped,
		                              mode, f->size + 1);
#else
		f->files[i] = io_cubep3m_open(log, fnames[i], swapped, mode,
		                              reader);
#endif
		if ((f->files)[i] == NULL) {
			int32_t j;
			for (j = i; i < f->numfiles; j++)
				free(fnames[j]);
			free(fnames);
			while (i > 0) {
				i--;
				io_cubep3m_close(log, &((f->files)[i]));
			}
			free(f->stem);
			free(f->path);
			free(f);
			return NULL;
		}
		free(fnames[i]);
	}
	free(fnames);

	/* Set initial values */
	f->no_part = UINT64_C(0);

	return f;
} /* io_mcubep3m_open */

extern void
io_mcubep3m_close(io_logging_t  log,
                  io_mcubep3m_t *f)
{
	int32_t i;

	/* Catch NULLs */
	if ((f == NULL) || (*f == NULL))
		return;

	/* Close */
	for (i = 0; i < (*f)->numfiles; i++)
		io_cubep3m_close(log, &(((*f)->files)[i]));

	/* Free */
	free((*f)->files);
	free((*f)->stem);
	free((*f)->path);
	free(*f);
	*f = NULL;

	return;
}

extern void
io_mcubep3m_init(io_logging_t  log,
                 io_mcubep3m_t f)
{
	int32_t i;

	if (f == NULL)
		return;

	if (f->files[0]->mode != IO_FILE_READ) {
		io_logging_warn(log, INT32_C(1),
		                "%s (first file of %" PRIi32
		                ") is not opened for reading. Will do nothing.",
		                f->files[0]->fname, f->numfiles);
		return;
	}

	/* Sum up the particle count */
	f->no_part = UINT64_C(0);
	for (i = 0; i < f->numfiles; i++) {
		io_cubep3m_init(log, (f->files)[i]);
		f->no_part += f->files[i]->no_part;
	}

	return;
}

extern uint64_t
io_mcubep3m_readpart(io_logging_t          log,
                     io_mcubep3m_t         f,
                     uint64_t              pskip,
                     uint64_t              pread,
                     io_file_strg_struct_t strg)
{
	uint64_t part_read, tmp;

	/* First read the particles unscaled. */
	part_read = io_mcubep3m_readpart_raw(log, f, pskip, pread, strg);
	if (part_read != pread) {
		return UINT64_C(0);
	}

	/* And do the scaling */
	tmp = io_cubep3m_scale_particles(log,
	                                 f->files[0]->header->boxsize,
	                                 f->files[0]->header->a,
	                                 f->files[0]->header->lunit,
	                                 f->files[0]->header->vunit,
	                                 part_read, strg);
	if (tmp != part_read) {
		return tmp;
	}

	/* Wow, we are done! */
	return part_read;
}

extern uint64_t
io_mcubep3m_readpart_raw(io_logging_t          log,
                         io_mcubep3m_t         f,
                         uint64_t              pskip,
                         uint64_t              pread,
                         io_file_strg_struct_t strg)
{
	long     tmp;
	uint64_t partread, pread_file, pread_done, partinfile;
	uint64_t pskip_file, pskip_done;
	int32_t  i;

	/* See if there is anything to do */
	if ((f == NULL) || (f->files == NULL))
		return UINT64_C(0);

	/* Set extreme position detectors */
//	f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e40;
//	f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e40;

	/* Initialize accounting of skipping and reading */
	pskip_done = pread_done = UINT64_C(0);

	/* Read the particles from the different files */
	for (i = 0; i < f->numfiles; i++) {
		/*
		 * First figure out how many particles are in the file (use the
		 * temporary long variable for that and copy it over to the 64
		 * bit integer afterwards.
		 */
		if (io_cubep3m_get(log, f->files[i],
		                   IO_FILE_GET_NOPART_IN_FILE, &tmp)
		    != true) {
			io_logging_fatal(log,
			                 "Could not get number of particles from %s",
			                 f->files[i]->fname);
			exit(EXIT_FAILURE);
		}
		partinfile = (uint64_t)tmp;

		/*
		 * Then do some arithmetic and set the skipping and reading
		 * numbers correctly for the file
		 */
		pskip_file = (pskip_done < pskip) ? pskip - pskip_done : UINT64_C(0);
		pread_file = (pread_done < pread) ? pread - pread_done : UINT64_C(0);
		if (pskip_file > partinfile) {
			pskip_file = partinfile;
			pread_file = UINT64_C(0);
		} else {
			if (pread_file > partinfile - pskip_file)
				pread_file = partinfile - pskip_file;
		}

		/* Now read the particles */
		partread = io_cubep3m_readpart_raw(log, f->files[i],
		                                   pskip_file, pread_file,
		                                   strg);
		if (pread_file != partread) {
			io_logging_fatal(log,
			                 "Something went wrong. Wanted to read %"
			                 PRIu64 " particles, but got %" PRIu64
			                 ". Aborting.", pread_file, partread);
			return pread_done + partread;
		}

		/* Update the skipping arithmetic */
		pskip_done += pskip_file;
		pread_done += pread_file;

		/* Move the particle pointers */
		strg.posx.val = (void *)(((char *)strg.posx.val)
		                         + strg.posx.stride * partread);
		strg.posy.val = (void *)(((char *)strg.posy.val)
		                         + strg.posy.stride * partread);
		strg.posz.val = (void *)(((char *)strg.posz.val)
		                         + strg.posz.stride * partread);
		strg.momx.val = (void *)(((char *)strg.momx.val)
		                         + strg.momx.stride * partread);
		strg.momy.val = (void *)(((char *)strg.momy.val)
		                         + strg.momy.stride * partread);
		strg.momz.val = (void *)(((char *)strg.momz.val)
		                         + strg.momz.stride * partread);
		if (strg.weight.val != NULL)
			strg.weight.val = (void *)(((char *)strg.weight.val)
			                           + strg.weight.stride * partread);
		if (strg.id.val != NULL)
			strg.id.val = (void *)(((char *)strg.id.val)
			                       + strg.id.stride * partread);
	} /* End of particle reading loop*/

	return pread_done;
} /* io_mcubep3m_readpart_raw */

extern bool
io_mcubep3m_get(io_logging_t  log,
                io_mcubep3m_t f,
                io_file_get_t what,
                void          *res)
{
	if ((f == NULL) || (f->files[0]->header == NULL))
		return false;

	switch (what) {
	case IO_FILE_GET_NOPART_IN_FILE:
	case IO_FILE_GET_NOPART:
		*((long *)res) = (long)(f->no_part);
		break;
	case IO_FILE_GET_NOVPART:
		*((double *)res) = (double)f->no_part;
		break;
	case IO_FILE_GET_NOSPECIES:
		*((int *)res) = 1;
		break;
	case IO_FILE_GET_BOXSIZE:
	case IO_FILE_GET_PMASS:
	case IO_FILE_GET_ZINITIAL:
	case IO_FILE_GET_Z:
	case IO_FILE_GET_AINITIAL:
	case IO_FILE_GET_A:
	case IO_FILE_GET_OMEGA0:
	case IO_FILE_GET_OMEGAL:
	case IO_FILE_GET_H:
	case IO_FILE_GET_DOUBLE:
	case IO_FILE_GET_MMASS:
	case IO_FILE_GET_NOTSTEP:
	case IO_FILE_GET_TSTEP:
	case IO_FILE_GET_HEADERSTR:
		io_cubep3m_get(log, f->files[0], what, res);
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
} /* io_mcubep3m_get */

extern void
io_mcubep3m_log(io_logging_t log, io_mcubep3m_t f)
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
	for (i = 0; i < f->numfiles; i++) {
		io_logging_msg(log, INT32_C(5),
		               "  ---> File %" PRIi32 ":",
		               i);
		io_cubep3m_log(log, (f->files)[i]);
	}

	return;
} /* io_mcubep3m_log */

/**********************************************************************
 *    Implementation of local functions                               *
 **********************************************************************/
