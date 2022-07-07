/* $Id: io_ares_header.c,v 1.5 2008/07/22 12:26:40 knolli Exp $ */

/**
 * \file io_ares_header.c
 *
 * Provides functions for reading and writing the header of ARES
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdint.h>

#include "io_ares_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_ares_header_t
io_ares_header_get(io_logging_t log, io_ares_t f)
{
	io_ares_header_t dummy;
	int32_t trash;

	/* Some sanity checks */
	if ((f == NULL) || (f->file == NULL))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

	/* Create the header structure */
	dummy = (io_ares_header_t)malloc(sizeof(io_ares_header_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "ARES header structure");
		return NULL;
	}

	/* Go to the header */
	rewind(f->file);

	/* Now start reading the header */
	if (f->swapped == IO_FILE_UNKOWN_SWAPPING) {
		uint32_t len_swapped, len_unswapped;

		f->swapped = IO_FILE_ISNOT_SWAPPED;
		io_util_readuint32(f->file, &(len_unswapped), f->swapped);
		rewind(f->file);
		f->swapped = IO_FILE_IS_SWAPPED;
		io_util_readuint32(f->file, &(len_swapped), f->swapped);
		if (len_unswapped < len_swapped) {
			 f->swapped = IO_FILE_ISNOT_SWAPPED;
			dummy->len_id = len_unswapped;
		} else {
			 f->swapped = IO_FILE_IS_SWAPPED;
			dummy->len_id = len_swapped;
		}
	} else {
		io_util_readuint32(f->file, &(dummy->len_id), f->swapped);
	}
	if (dummy->len_id != ARES_LENGTH_IDENTIFIER) {
		io_logging_warn(log, INT32_C(1),
		                "The length of the file id is %" PRIu32
		                " but it should be %" PRIu32 ". Wrong file?");
	}
	fread(dummy->id, sizeof(char), ARES_LENGTH_IDENTIFIER, f->file);
	dummy->id[ARES_LENGTH_IDENTIFIER] = '\0';
	io_util_readuint64(f->file, &(dummy->no_part), f->swapped);
	io_util_readuint64(f->file, &(dummy->no_part_in_file), f->swapped);
	io_util_readuint32(f->file, &(dummy->bytes_float), f->swapped);
	io_util_readuint32(f->file, &(dummy->bytes_int), f->swapped);
	io_util_readuint64(f->file, &(dummy->no_species), f->swapped);
	io_util_readuint64(f->file, &(dummy->has_mass), f->swapped);
	io_util_readuint64(f->file, &(dummy->has_gas), f->swapped);
	io_util_readdouble(f->file, &(dummy->no_vpart), f->swapped);
	io_util_readdouble(f->file, &(dummy->boxsize), f->swapped);
	io_util_readdouble(f->file, &(dummy->omega0), f->swapped);
	io_util_readdouble(f->file, &(dummy->lambda0), f->swapped);
	io_util_readdouble(f->file, &(dummy->pmass), f->swapped);
	io_util_readdouble(f->file, &(dummy->minweight), f->swapped);
	io_util_readdouble(f->file, &(dummy->maxweight), f->swapped);
	io_util_readdouble(f->file, &(dummy->a_initial), f->swapped);
	io_util_readdouble(f->file, &(dummy->a_current), f->swapped);
	io_util_readdouble(f->file, &(dummy->timestep), f->swapped);
	io_util_readuint64(f->file, &(dummy->minkey), f->swapped);
	io_util_readuint64(f->file, &(dummy->maxkey), f->swapped);
	io_util_readint(f->file, &(dummy->lb_level), f->swapped);
	io_util_readint(f->file, &(dummy->rank), f->swapped);
	io_util_readint(f->file, &(dummy->size), f->swapped);

	/* And set the header */
	f->header = dummy;

	return dummy;
}

extern io_ares_header_t
io_ares_header_new(io_logging_t log)
{
	io_ares_header_t dummy;

	dummy = (io_ares_header_t)calloc(sizeof(io_ares_header_struct_t),
	                                 (size_t)1);
    if (dummy == NULL) {
		io_logging_memfatal(log, "ARES header structure");
		return NULL;
	}
	dummy->len_id = ARES_LENGTH_IDENTIFIER;
	strncpy(dummy->id, ARES_IDENTIFIER_STRING, ARES_LENGTH_IDENTIFIER);
	dummy->id[ARES_LENGTH_IDENTIFIER] = '\0';

	return dummy;
}

extern void
io_ares_header_del(io_logging_t log, io_ares_header_t *header)
{
	if ( (header == NULL) || (*header == NULL) )
		return;

	free(*header);

	*header = NULL;

	return;
}

extern void
io_ares_header_write(io_logging_t log,
                     io_ares_header_t header,
                     io_ares_t f)
{
	int32_t tmp, i;
	char pad[ARES_HEADER_FILLHEADER];

	/* Sanity checks */
	if ( (header == NULL) || (f == NULL) )
		return;

	if (f->mode != IO_FILE_WRITE)
		return;

	if (f->file == NULL)
		return;

	if (header != f->header) {
		io_logging_msg(log, INT32_C(1),
		               "Writing a different header than stored in "
		               "the file object to the file");
	}

	/* Go to the right spot */
	rewind(f->file);

	/* Start writing */
	fwrite(&header->len_id, sizeof(uint32_t), 1, f->file);
	fwrite(header->id, sizeof(char), ARES_LENGTH_IDENTIFIER, f->file);
	fwrite(&header->no_part, sizeof(uint64_t), 1, f->file);
	fwrite(&header->no_part_in_file, sizeof(uint64_t), 1, f->file);
	fwrite(&header->bytes_float, sizeof(uint32_t), 1, f->file);
	fwrite(&header->bytes_int, sizeof(uint32_t), 1, f->file);
	fwrite(&header->no_species, sizeof(uint64_t), 1, f->file);
	fwrite(&header->has_mass, sizeof(uint64_t), 1, f->file);
	fwrite(&header->has_gas, sizeof(uint64_t), 1, f->file);
	fwrite(&header->no_vpart, sizeof(double), 1, f->file);
	fwrite(&header->boxsize, sizeof(double), 1, f->file);
	fwrite(&header->omega0, sizeof(double), 1, f->file);
	fwrite(&header->lambda0, sizeof(double), 1, f->file);
	fwrite(&header->pmass, sizeof(double), 1, f->file);
	fwrite(&header->minweight, sizeof(double), 1, f->file);
	fwrite(&header->maxweight, sizeof(double), 1, f->file);
	fwrite(&header->a_initial, sizeof(double), 1, f->file);
	fwrite(&header->a_current, sizeof(double), 1, f->file);
	fwrite(&header->timestep, sizeof(double), 1, f->file);
	fwrite(&header->minkey, sizeof(uint64_t), 1, f->file);
	fwrite(&header->maxkey, sizeof(uint64_t), 1, f->file);
	fwrite(&header->lb_level, sizeof(int), 1, f->file);
	fwrite(&header->rank, sizeof(int), 1, f->file);
	fwrite(&header->size, sizeof(int), 1, f->file);

	/* Padding to the right size now */
	for (i=0 ; i<ARES_HEADER_FILLHEADER; i++)
		pad[i] = 0xff;
	fwrite(pad, sizeof(char), ARES_HEADER_FILLHEADER, f->file);

	/* Done writing the header */

	return;
}

extern void
io_ares_header_log(io_logging_t log, io_ares_header_t header)
{
	io_logging_msg(log, INT32_C(5),
	               "Headerobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Length of identifier       : %" PRIu32,
	               header->len_id);
	io_logging_msg(log, INT32_C(5),
	               "  Identifier                 : %s",
	               header->id);
	io_logging_msg(log, INT32_C(5),
	               "  Number of particles        : %" PRIu64,
	               header->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  Number of particles in file: %" PRIu64,
	               header->no_part_in_file);
	io_logging_msg(log, INT32_C(5),
	               "  Bytes per float in file    : %" PRIu32,
	               header->bytes_float);
	io_logging_msg(log, INT32_C(5),
	               "  Bytes per integer in file  : %" PRIu32,
	               header->bytes_int);
	io_logging_msg(log, INT32_C(5),
	               "  Number of species          : %" PRIu64,
	               header->no_species);
	io_logging_msg(log, INT32_C(5),
	               "  Has mass                   : %" PRIu64,
	               header->has_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Has gas                    : %" PRIu64,
	               header->has_gas);
	io_logging_msg(log, INT32_C(5),
	               "  Number of virtual particles: %e",
	               header->no_vpart);
	io_logging_msg(log, INT32_C(5),
	               "  Boxsize                    : %e",
	               header->boxsize);
	io_logging_msg(log, INT32_C(5),
	               "  Omega0                     : %e",
	               header->omega0);
	io_logging_msg(log, INT32_C(5),
	               "  lambda0                    : %e",
	               header->lambda0);
	io_logging_msg(log, INT32_C(5),
	               "  pmass                      : %e",
	               header->pmass);
	io_logging_msg(log, INT32_C(5),
	               "  Minimal weight             : %e",
	               header->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  Maximal weight             : %e",
	               header->maxweight);
	io_logging_msg(log, INT32_C(5),
	               "  a_initial                  : %e",
	               header->a_initial);
	io_logging_msg(log, INT32_C(5),
	               "  a_current                  : %e",
	               header->a_current);
	io_logging_msg(log, INT32_C(5),
	               "  Timestep                   : %e",
	               header->timestep);
	io_logging_msg(log, INT32_C(5),
	               "  Minimal SFC key            : %" PRIu64,
	               header->minkey);
	io_logging_msg(log, INT32_C(5),
	               "  Maximal SFC key            : %" PRIu64,
	               header->maxkey);
	io_logging_msg(log, INT32_C(5),
	               "  Loadbalancing level        : %i",
	               header->lb_level);
	io_logging_msg(log, INT32_C(5),
	               "  Rank of writing process    : %i",
	               header->rank);
	io_logging_msg(log, INT32_C(5),
	               "  Size of MPI domain         : %i",
	               header->size);

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
