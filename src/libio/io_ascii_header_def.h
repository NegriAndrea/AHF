#ifndef IO_ASCII_HEADER_DEF_H 
#define IO_ASCII_HEADER_DEF_H

/**
 * \file io_ascii_header_def.h
 *
 * Provides the structure definition for the ASCII header
 * structure. Including useful typedefs.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/** Defines the length for the identifying header string */
#define ASCII_HEADER_HEADERSTRING 256

/**
 * The header structure itself
 */
struct io_ascii_header_struct {
	char          header[ASCII_HEADER_HEADERSTRING];
	int32_t       multi_mass;
	long          no_part;
	long          no_species;
	double        total_mass;
	int32_t       no_timestep;
	double        boxsize;
	double        omega0;
	double        lambda0;
	double        pmass;
	double        a_initial;
	double        a_current;
};

/** Convenient typedef */
typedef struct io_ascii_header_struct io_ascii_header_struct_t;

/** Convenient typedef */
typedef io_ascii_header_struct_t *io_ascii_header_t;


#endif /* IO_ASCII_HEADER_DEF_H */
