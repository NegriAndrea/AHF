#ifndef IO_AMIGA_DEF_H
#define IO_AMIGA_DEF_H

/* $Id: io_amiga_def.h,v 1.5 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file io_amiga_def.h
 *
 * Provides the structure definition for the AMIGA file structure.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdint.h>

#include "io_amiga_header_def.h"
#include "io_file.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/**
 * The file structure itself
 */
struct io_amiga_struct {
	/** Holds the filetype ID, *must* be first in the structure */
	io_file_type_t ftype;
#ifdef WITH_MPI
	/** The global rank of the process */
	int rank;
	/** The size of the global communicator */
	int size;
	/** Stores the communicator used for intra libio communication */
	MPI_Comm mycomm;
	/** The size of the intra-library communicator */
	int size_mycomm;
	/** The rank of the local process */
	int rank_mycomm;
#endif
	/** Holds the file stream */
	FILE *file;
	/** Holds the filename */
	char *fname;
	/** In which mode should this be opened */
	io_file_mode_t mode;
	/** Is this file byteswapped? */
	io_file_swap_t swapped;
	/** How many bytes are long in the file? */
	int32_t file_sizeof_long;
	/** Holds the header information */
	io_amiga_header_t header;
	/** Stores the minimal weight of the particles */
	double minweight;
	/** Stores the maximal weight of the particles */
	double maxweight;
};

/** Convenient typedef */
typedef struct io_amiga_struct io_amiga_struct_t;

/** Convenient typedef */
typedef io_amiga_struct_t *io_amiga_t;


#endif /* IO_AMIGA_DEF_H */
