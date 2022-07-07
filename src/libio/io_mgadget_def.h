#ifndef IO_MGADGET_DEF_H
#define IO_MGADGET_DEF_H

/* $Id: io_mgadget_def.h,v 1.4 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file io_mgadget_def.h
 *
 * Provides the structure definition for the Multiple Gadget file
 * structure.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdint.h>

#include "io_gadget_def.h"
#include "io_file.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/**
 * The file structure itself
 */
struct io_mgadget_struct {
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
	/** The Path */
	char *path;
	/** The filename stem */
	char *stem;
	/** Holds the number of files in this set */
	int32_t numfiles;
	/** Holds the array of Gadget files (incl. the actual FILE pointers) */
	io_gadget_t *files;
	/** Stores the total number of particles in the file */
	uint64_t no_part;
	/** Flag if this is a multimass file */
	int8_t multimass;
	/** Stores the minimal weight of halo particles in file units */
	double mmass;
	/** Stores the minimal particle weight in file units */
	double minweight;
	/** Stores the maximal particle weight in file units */
	double maxweight;
	/** Stores the sum over all particle weights in file units */
	double sumweight;
	/** Stores the number of weight types */
	int32_t no_species;
	/** Stores the minimal position values in file units */
	double minpos[3];
	/** Stores the maximal position values in file units */
	double maxpos[3];
	/** Scaling value to convert from internal units to Mpc */
	double posscale;
	/** Scaling value to convert from internal mass to Msun */
	double weightscale;
};

/** Convenient typedef */
typedef struct io_mgadget_struct io_mgadget_struct_t;

/** Convenient typedef */
typedef io_mgadget_struct_t *io_mgadget_t;


#endif /* IO_MGADGET_DEF_H */
