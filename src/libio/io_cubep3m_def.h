#ifndef IO_CUBEP3M_DEF_H
#define IO_CUBEP3M_DEF_H

/**
 * \file io_cubep3m_def.h
 *
 * Provides the structure definition for the CUBEP3M file structure.
 */


/**********************************************************************\
 *    Includes                                                        *
 \**********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "io_cubep3m_header_def.h"
#include "io_file.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              *
 \**********************************************************************/

/**
 * The file structure itself
 */
struct io_cubep3m_struct {
	/** Holds the filetype ID, *must* be first in the structure */
	io_file_type_t ftype;
#ifdef WITH_MPI
	/** The global rank of the process */
	int      rank;
	/** The size of the global communicator */
	int      size;
	/** Stores the communicator used for intra libio communication */
	MPI_Comm mycomm;
	/** The size of the intra-library communicator */
	int      size_mycomm;
	/** The rank of the local process */
	int      rank_mycomm;
#endif
	/** Holds the file stream of the opened file */
	FILE                *file;
	/** Holds the filename */
	char                *fname;
	/** In which mode should this be opened */
	io_file_mode_t      mode;
	/** Is this file byteswapped? */
	io_file_swap_t      swapped;
	/** Toggles between normal and chunked CubeP3M files */
	bool                isChunked;
	/** Stores the total number of particles in the file */
	uint64_t            no_part;
	/** Stores the minimal position values in file units */
	double              minpos[3];
	/** Stores the maximal position values in file units */
	double              maxpos[3];
	/** Holds the header information for this file */
	io_cubep3m_header_t header;
};

/** Convenient typedef */
typedef struct io_cubep3m_struct io_cubep3m_struct_t;

/** Convenient typedef */
typedef io_cubep3m_struct_t      *io_cubep3m_t;


#endif /* IO_CUBEP3M_DEF_H */
