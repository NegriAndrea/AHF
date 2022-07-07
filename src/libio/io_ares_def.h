#ifndef IO_ARES_DEF_H
#define IO_ARES_DEF_H

/* $Id: io_ares_def.h,v 1.1 2007/11/30 01:21:06 knolli Exp $ */

/**
 * \file io_ares_def.h
 *
 * Provides the structure definition for the ARES file structure.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdint.h>

#include "io_ares_header_def.h"
#include "io_file.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/**
 * The file structure itself
 */
struct io_ares_struct {
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
	/** Holds the header information */
	io_ares_header_t header;
};

/** Convenient typedef */
typedef struct io_ares_struct io_ares_struct_t;

/** Convenient typedef */
typedef io_ares_struct_t *io_ares_t;


#endif /* IO_ARES_DEF_H */
