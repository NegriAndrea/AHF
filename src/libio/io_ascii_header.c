/**
 * \file io_ascii_header.c
 *
 * Provides functions for reading and writing the header of ASCII
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

#include "io_ascii_header.h"
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
extern io_ascii_header_t
io_ascii_header_get(io_logging_t log, io_ascii_t f)
{
	io_ascii_header_t dummy;

	/* Some sanity checks */
	if ((f == NULL) || (f->file == NULL))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

	/* Create the header structure */
	dummy = malloc(sizeof(io_ascii_header_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "ASCII header structure");
		return NULL;
	}

	/* Go to the header */
	rewind(f->file);

	/* Now start reading the header */
	fscanf(f->file, "# Header = %s \n", dummy->header);
	fscanf(f->file, "# MultiMass = %"SCNi32" \n", &(dummy->multi_mass));
	fscanf(f->file, "# Boxsize = %lf \n", &(dummy->boxsize));
	fscanf(f->file, "# NumberOfPart = %li \n", &(dummy->no_part));
	fscanf(f->file, "# NumberOfSpecies = %li \n", &(dummy->no_species));
	fscanf(f->file, "# TotalMass = %lf \n", &(dummy->total_mass));
	fscanf(f->file, "# NumberOfTimestep = %"SCNi32" \n", &(dummy->no_timestep));
	fscanf(f->file, "# Omega0 = %lf \n", &(dummy->omega0));
	fscanf(f->file, "# Lambda0 = %lf \n", &(dummy->lambda0));
	fscanf(f->file, "# PMass = %lf \n", &(dummy->pmass));
	fscanf(f->file, "# ZInitial = %lf \n", &(dummy->a_initial));
	fscanf(f->file, "# ZCurrent = %lf \n", &(dummy->a_current));

  dummy->a_initial = 1./(1.+dummy->a_initial);
  dummy->a_current = 1./(1.+dummy->a_current);
  
//  fprintf(stderr,"%s\n",dummy->header);
//  fprintf(stderr,"%"PRIi32"\n",dummy->multi_mass);
//  fprintf(stderr,"%lf\n",dummy->boxsize);

	/* Set the header */
	f->header = dummy;

	/* And return it */
	return dummy;
}

extern io_ascii_header_t
io_ascii_header_new(io_logging_t log)
{
	io_ascii_header_t dummy;

	dummy = malloc(sizeof(io_ascii_header_struct_t));
    if (dummy == NULL) {
		io_logging_memfatal(log, "ASCII header structure");
		return NULL;
	}

	return dummy;
}

extern void
io_ascii_header_del(io_logging_t log, io_ascii_header_t *header)
{
	if ( (header == NULL) || (*header == NULL) )
		return;

	free(*header);

	*header = NULL;

	return;
}

extern void
io_ascii_header_write(io_logging_t log,
                      io_ascii_header_t header,
                      io_ascii_t f)
{
	/* NOT IMPLEMENTED */

	return;
}

extern void
io_ascii_header_log(io_logging_t log, io_ascii_header_t header)
{
	io_logging_msg(log, INT32_C(5),
	               "Headerobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Header:                      %s",
	               header->header);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass:                   %" PRIi32,
	               header->multi_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Number of particles:         %li",
	               header->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  Number of mass species:      %li",
	               header->no_species);
	io_logging_msg(log, INT32_C(5),
	               "  Total mass:                  %e",
	               header->total_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Number of timestep         : %" PRIi32,
	               header->no_timestep);
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
	               "  a_initial                  : %e",
	               header->a_initial);
	io_logging_msg(log, INT32_C(5),
	               "  a_current                  : %e",
	               header->a_current);

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
