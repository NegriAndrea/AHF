/**
 * \file io_cubep3m_header.c
 *
 * Provides functions for reading and writing the header of CUBEP3M
 * files.
 */


/**********************************************************************
 *    Includes                                                        *
 **********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#include "io_cubep3m_header.h"
#include "io_cubep3m_util.h"
#include "io_util.h"


/**********************************************************************
 *    Local defines, structure definitions and typedefs               *
 **********************************************************************/


/**********************************************************************
 *    Prototypes of local functions                                   *
 **********************************************************************/

/**
 * \brief  Resets the values in the header to safe values.
 *
 * \param[in,out]  header
 *                    The header to null.
 *
 * \return  Returns nothing.
 */
inline static void
local_null_header(io_cubep3m_header_t header);


/**
 * \brief  Tries to get the file number from the file name.
 *
 * This will identify the last number in the file-name and interpret
 * this as the file number, i.e for 8.397xv24.dat it would return 24.
 *
 * \param  *fname  The filename.
 *
 * \return  Returns hopefully the file number...
 */
inline static int
local_extract_file_number_from_name(const char *fname);


/**
 * \brief  Calculates the offset of the file.
 *
 * \param[in]   file_number
 *                 The number of the file to calculate the offset for.  The
 *                 file number must be a positive or zero integer smaller
 *                 than the total number of files, which is given as
 *                 \c nodes_dim cubed.
 * \param[in]   nodes_dim
 *                 The number of files per dimension.  Must be a positive
 *                 integer.
 * \param[in]   ngrid
 *                 The number of grid cells in one dimension.  This must be
 *                 larger than \c nodes_dim.
 * \param[out]  offset
 *                 This will receive the offset of the file.  This must be
 *                 an array of at least 3 elements.
 *
 * \return  Returns nothing.
 */
inline static void
local_calculate_file_offset(int      file_number,
                            int      nodes_dim,
                            uint64_t ngrid,
                            double   *offset);


/**
 * \brief  Tries to open the CubeP3M info.
 *
 * \param[in]      info_file_name
 *                    The file name of the info file to open.  Must be a
 *                    valid string.
 * \param[in,out]  log
 *                    The logging module to use for error messages.  Must
 *                    not be a valid logging module.
 *
 * \return  Returns the file handle of the just opened info file.
 */
inline static FILE *
local_open_info_file(const char *info_file_name, io_logging_t log);


/**********************************************************************
 *    Implementation of global functions                              *
 **********************************************************************/
extern io_cubep3m_header_t
io_cubep3m_header_get(io_logging_t log, io_cubep3m_t f)
{
	io_cubep3m_header_t dummy;

	if ((f == NULL) || (f->file == NULL))
		return NULL;

	if (f->header != NULL)
		return f->header;

	dummy = (io_cubep3m_header_t)malloc(sizeof(io_cubep3m_header_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "CUBEP3M header structure");
		return NULL;
	}
	local_null_header(dummy);

	if (f->isChunked) {
		fseek(f->file, CUBEP3M_HEADER_MAGIC_SIZE, SEEK_SET);
		io_cubep3m_header_read_chunk_info(log, f->file, f->swapped, dummy);
	} else {
		rewind(f->file);
	}
	io_cubep3m_header_read_basics(log, f->file, f->swapped, dummy);
	io_cubep3m_header_read_extras(log, "cubep3m.info", dummy);

	dummy->lunit = io_cubep3m_util_lunit(dummy->boxsize, dummy->ngrid);
	dummy->munit = io_cubep3m_util_munit(dummy->boxsize, dummy->omega0,
	                                     dummy->nptotal, dummy->mass_p);
	dummy->vunit = io_cubep3m_util_vunit(dummy->boxsize, dummy->omega0,
	                                     dummy->a, dummy->ngrid);

	dummy->file_number = local_extract_file_number_from_name(f->fname);
	if (f->isChunked) {
		dummy->offset[0] = 0.0;
		dummy->offset[1] = 0.0;
		dummy->offset[2] = 0.0;
	} else {
		local_calculate_file_offset(dummy->file_number, dummy->nodes_dim,
		                            dummy->ngrid, dummy->offset);
	}

	f->header = dummy;
	return dummy;
} /* io_cubep3m_header_get */

extern void
io_cubep3m_header_read_basics(io_logging_t        log,
                              FILE                *f,
                              io_file_swap_t      swapped,
                              io_cubep3m_header_t header)
{
	assert(log != NULL);
	assert(f != NULL);
	assert(header != NULL);

	io_util_readint32(f, &(header->np_local), swapped);
	io_util_readfloat(f, &(header->a), swapped);
	io_util_readfloat(f, &(header->t), swapped);
	io_util_readfloat(f, &(header->tau), swapped);
	io_util_readint32(f, &(header->nts), swapped);
	io_util_readfloat(f, &(header->dt_f_acc), swapped);
	io_util_readfloat(f, &(header->dt_pp_acc), swapped);
	io_util_readfloat(f, &(header->dt_c_acc), swapped);
	io_util_readint32(f, &(header->cur_checkpoint), swapped);
	io_util_readint32(f, &(header->cur_projection), swapped);
	io_util_readint32(f, &(header->cur_halofind), swapped);
	io_util_readfloat(f, &(header->mass_p), swapped);
}

extern void
io_cubep3m_header_read_chunk_info(io_logging_t        log,
                                  FILE                *f,
                                  io_file_swap_t      swapped,
                                  io_cubep3m_header_t header)
{
	assert(log != NULL);
	assert(f != NULL);
	assert(header != NULL);

	io_util_readint64(f, &(header->tot_np_in_chunk), swapped);
	io_util_readint64(f, &(header->np_in_chunk_file), swapped);
	io_util_readint(f, &(header->num_chunk_files), swapped);
	io_util_readdouble(f, header->chunk_offset, swapped);
	io_util_readdouble(f, header->chunk_offset + 1, swapped);
	io_util_readdouble(f, header->chunk_offset + 2, swapped);
	io_util_readdouble(f, header->chunk_start_full_data, swapped);
	io_util_readdouble(f, header->chunk_start_full_data + 1, swapped);
	io_util_readdouble(f, header->chunk_start_full_data + 2, swapped);
	io_util_readdouble(f, header->chunk_end_full_data, swapped);
	io_util_readdouble(f, header->chunk_end_full_data + 1, swapped);
	io_util_readdouble(f, header->chunk_end_full_data + 2, swapped);
	io_util_readdouble(f, header->chunk_start_real_data, swapped);
	io_util_readdouble(f, header->chunk_start_real_data + 1, swapped);
	io_util_readdouble(f, header->chunk_start_real_data + 2, swapped);
	io_util_readdouble(f, header->chunk_end_real_data, swapped);
	io_util_readdouble(f, header->chunk_end_real_data + 1, swapped);
	io_util_readdouble(f, header->chunk_end_real_data + 2, swapped);
}

extern void
io_cubep3m_header_read_extras(io_logging_t        log,
                              const char          *extras_file_name,
                              io_cubep3m_header_t header)
{
	FILE *fInfo;
	char cubep3mline[2048];

	fInfo = local_open_info_file(extras_file_name, log);

	fgets(cubep3mline, 2048, fInfo);
	sscanf(cubep3mline, "%lf", &(header->omega0));

	fgets(cubep3mline, 2048, fInfo);
	sscanf(cubep3mline, "%lf", &(header->lambda0));

	fgets(cubep3mline, 2048, fInfo);
	sscanf(cubep3mline, "%lf", &(header->boxsize));

	fgets(cubep3mline, 2048, fInfo);
	sscanf(cubep3mline, "%" SCNu64 "", &(header->ngrid));
	header->nptotal = (header->ngrid / 2) * (header->ngrid / 2)
	                  * (header->ngrid / 2);

	fgets(cubep3mline, 2048, fInfo);
	sscanf(cubep3mline, "%i", &(header->nodes_dim));

	fclose(fInfo);
}

extern void
io_cubep3m_header_del(io_logging_t log, io_cubep3m_header_t *header)
{
	if ((header == NULL) || (*header == NULL))
		return;

	free(*header);
	*header = NULL;

	return;
}

extern void
io_cubep3m_header_write(io_logging_t        log,
                        io_cubep3m_header_t header,
                        io_cubep3m_t        f)
{
	if ((header == NULL) || (f == NULL))
		return;

	if (f->mode != IO_FILE_WRITE)
		return;

	if (f->file == NULL)
		return;

	if (header != f->header) {
		io_logging_msg(log, INT32_C(1),
		               "Writing a different header than stored in "
		               "the file object to the file.");
	}

	/* TODO: Write the header */
	return;
}

extern void
io_cubep3m_header_log(io_logging_t log, io_cubep3m_header_t header)
{
	io_logging_msg(log, INT32_C(5),
	               "   np_local:                     %" PRIi32,
	               header->np_local);
	io_logging_msg(log, INT32_C(5),
	               "   a:                            %e",
	               header->a);
	io_logging_msg(log, INT32_C(5),
	               "   t:                            %e",
	               header->t);
	io_logging_msg(log, INT32_C(5),
	               "   tau:                          %e",
	               header->tau);
	io_logging_msg(log, INT32_C(5),
	               "   nts:                          %" PRIi32,
	               header->nts);
	io_logging_msg(log, INT32_C(5),
	               "   dt_f_acc:                     %e",
	               header->dt_f_acc);
	io_logging_msg(log, INT32_C(5),
	               "   dt_pp_acc:                    %e",
	               header->dt_pp_acc);
	io_logging_msg(log, INT32_C(5),
	               "   dt_c_acc:                     %e",
	               header->dt_c_acc);
	io_logging_msg(log, INT32_C(5),
	               "   cur_checkpoint:               %" PRIi32,
	               header->cur_checkpoint);
	io_logging_msg(log, INT32_C(5),
	               "   cur_projection:               %" PRIi32,
	               header->cur_projection);
	io_logging_msg(log, INT32_C(5),
	               "   cur_halofind:                 %" PRIi32,
	               header->cur_halofind);
	io_logging_msg(log, INT32_C(5),
	               "   mass_p:                       %e",
	               header->mass_p);

	io_logging_msg(log, INT32_C(5),
	               "  Chunk information");
	io_logging_msg(log, INT32_C(5),
	               "   tot_np_in_chunk:              %" PRIi64,
	               header->tot_np_in_chunk);
	io_logging_msg(log, INT32_C(5),
	               "   np_in_chunk_file:             %" PRIi64,
	               header->np_in_chunk_file);
	io_logging_msg(log, INT32_C(5),
	               "   num_chunk_files:              %i",
	               header->num_chunk_files);
	io_logging_msg(log, INT32_C(5),
	               "   chunk_offset:                 %e %e %e",
	               header->chunk_offset[0],
	               header->chunk_offset[1],
	               header->chunk_offset[2]);
	io_logging_msg(log, INT32_C(5),
	               "   chunk_start_full_data:        %e %e %e",
	               header->chunk_start_full_data[0],
	               header->chunk_start_full_data[1],
	               header->chunk_start_full_data[2]);
	io_logging_msg(log, INT32_C(5),
	               "   chunk_end_full_data:          %e %e %e",
	               header->chunk_end_full_data[0],
	               header->chunk_end_full_data[1],
	               header->chunk_end_full_data[2]);
	io_logging_msg(log, INT32_C(5),
	               "   chunk_start_real_data:        %e %e %e",
	               header->chunk_start_real_data[0],
	               header->chunk_start_real_data[1],
	               header->chunk_start_real_data[2]);
	io_logging_msg(log, INT32_C(5),
	               "   chunk_end_real_data:          %e %e %e",
	               header->chunk_end_real_data[0],
	               header->chunk_end_real_data[1],
	               header->chunk_end_real_data[2]);

	io_logging_msg(log, INT32_C(5),
	               "  Derived information (and taken from the info file)");
	io_logging_msg(log, INT32_C(5),
	               "   omega0:                       %e",
	               header->omega0);
	io_logging_msg(log, INT32_C(5),
	               "   lambda0:                      %e",
	               header->lambda0);
	io_logging_msg(log, INT32_C(5),
	               "   boxsize:                      %e",
	               header->boxsize);
	io_logging_msg(log, INT32_C(5),
	               "   lunit:                        %e",
	               header->lunit);
	io_logging_msg(log, INT32_C(5),
	               "   vunit:                        %e",
	               header->vunit);
	io_logging_msg(log, INT32_C(5),
	               "   munit:                        %e",
	               header->munit);
	io_logging_msg(log, INT32_C(5),
	               "   nptotal:                      %" PRIu64,
	               header->nptotal);
	io_logging_msg(log, INT32_C(5),
	               "   ngrid:                        %" PRIu64,
	               header->ngrid);
	io_logging_msg(log, INT32_C(5),
	               "   nodes_dim:                    %i",
	               header->mass_p);
	io_logging_msg(log, INT32_C(5),
	               "  Derived information");
	io_logging_msg(log, INT32_C(5),
	               "   file_number:                  %i",
	               header->file_number);
	io_logging_msg(log, INT32_C(5),
	               "   offset:                       %e  %e  %e",
	               header->offset[0], header->offset[1],
	               header->offset[2]);

	return;
} /* io_cubep3m_header_log */

/**********************************************************************
 *    Implementation of local functions                               *
 **********************************************************************/

inline static void
local_null_header(io_cubep3m_header_t header)
{
	assert(header != NULL);

	header->np_local         = 0;
	header->a                = 0.0;
	header->t                = 0.0;
	header->tau              = 0.0;
	header->nts              = 0;
	header->dt_f_acc         = 0.0;
	header->dt_pp_acc        = 0.0;
	header->dt_c_acc         = 0.0;
	header->cur_checkpoint   = 0;
	header->cur_projection   = 0.0;
	header->cur_halofind     = 0.0;
	header->mass_p           = 0.0;

	header->tot_np_in_chunk  = 0;
	header->np_in_chunk_file = 0;
	header->num_chunk_files  = 0;
	for (int i = 0; i < 3; i++) {
		header->chunk_offset[i]          = 0.0;
		header->chunk_start_full_data[i] = 0.0;
		header->chunk_end_full_data[i]   = 0.0;
		header->chunk_start_real_data[i] = 0.0;
		header->chunk_end_real_data[i]   = 0.0;
	}

	header->omega0      = 0.0;
	header->lambda0     = 0.0;
	header->boxsize     = 0.0;
	header->lunit       = 0.0;
	header->vunit       = 0.0;
	header->munit       = 0.0;
	header->nptotal     = 0;
	header->ngrid       = 0;
	header->nodes_dim   = 0;
	header->file_number = 0;
	header->offset[0]   = 0.0;
	header->offset[1]   = 0.0;
	header->offset[2]   = 0.0;
} /* local_null_header */

inline static int
local_extract_file_number_from_name(const char *fname)
{
	int file_number = 0;
	int pos         = strlen(fname) - 1;

	/* Skip the ending of the file */
	while (pos > 0 && !(isdigit(fname[pos])))
		pos--;
	/* Now find the beginning of the *last* number in the filename */
	while (pos > 0 && isdigit(fname[pos - 1]))
		pos--;

	/* Parse that number as the file number */
	sscanf(fname + pos, "%i", &file_number);

	return file_number;
}

inline static void
local_calculate_file_offset(int      file_number,
                            int      nodes_dim,
                            uint64_t ngrid,
                            double   *offset)
{
	assert(nodes_dim > 0);
	assert(file_number >= 0);
	assert(file_number < nodes_dim * nodes_dim * nodes_dim);
	assert(ngrid > nodes_dim);
	assert(offset != NULL);

	int    pos[3];
	double delta;

	pos[2]       = file_number / (nodes_dim * nodes_dim);
	file_number -= pos[2] * (nodes_dim * nodes_dim);
	pos[1]       = file_number / nodes_dim;
	file_number -= pos[1] * nodes_dim;
	pos[0]       = file_number;
	delta        = ngrid / ((double)nodes_dim);

	offset[0]    = pos[0] * delta;
	offset[1]    = pos[1] * delta;
	offset[2]    = pos[2] * delta;
}

inline static FILE *
local_open_info_file(const char *info_file_name, io_logging_t log)
{
	FILE *f = NULL;

	f = fopen(info_file_name, "r");

	if (f == NULL) {
		io_logging_fatal(log,
		                 "Could not open cubep3m.info containing the "
		                 "following information:");
		io_logging_fatal(log,
		                 "0.3               CUBEP3M_OMEGA0");
		io_logging_fatal(log,
		                 "0.7               CUBEP3M_LAMBDA0");
		io_logging_fatal(log,
		                 "20.0              CUBEP3M_BOXSIZE [Mpc/h]");
		io_logging_fatal(log,
		                 "512               CUBEP3M_NGRID");
		io_logging_fatal(log,
		                 "10                CUBEP3M_NODES_DIM");
		io_logging_fatal(log, "Please generate this file. Aborting now!");
		exit(EXIT_FAILURE);
	}

	return f;
}
