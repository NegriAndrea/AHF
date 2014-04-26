#ifndef IO_MCUBEP3M_DEF_H
#define IO_MCUBEP3M_DEF_H


/**
 * \file io_mcubep3m_def.h
 *
 * Provides the structure definition for the Multiple CubeP3M file
 * structure.
 */


/**********************************************************************
 *    Includes                                                        *
 **********************************************************************/
#include <stdio.h>
#include <stdint.h>

#include "io_cubep3m_def.h"
#include "io_file.h"


/**********************************************************************
 *    Global defines, structure definitions and typedefs              *
 **********************************************************************/

/**
 * The file structure itself
 */
struct io_mcubep3m_struct {
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
	/** The Path */
	char         *path;
	/** The filename stem */
	char         *stem;
	/** Holds the number of files in this set */
	int32_t      numfiles;
	/** Holds the array of CubeP3M files */
	io_cubep3m_t *files;
	/** Stores the total number of particles in the file */
	uint64_t     no_part;
};

/** Convenient typedef */
typedef struct io_mcubep3m_struct io_mcubep3m_struct_t;

/** Convenient typedef */
typedef io_mcubep3m_struct_t      *io_mcubep3m_t;


#endif /* IO_MCUBEP3M_DEF_H */
