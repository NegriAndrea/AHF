// Copyright (C) 2011, Steffen Knollmann
// Released under the terms of the GNU General Public License version 3.
// This file is part of `AHF'.


#ifndef IO_ARTNEW_DEF_H
#define IO_ARTNEW_DEF_H

/*--- Doxygen file description ------------------------------------------*/

/**
 * @file  io_artnew_def.h
 *
 * Provides the main structure of the ARTnew reader.
 */


/*--- Includes ----------------------------------------------------------*/
struct io_artnew_struct {
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
	int rank_mycomm;
#endif
	/** The ART file object. */
	art_t       art;
	/** Keeps a reference to the header (taken from io_artnew_struct::art) */
	artHeader_t header;
	/** The particles mass of particles with weight 1. */
	double      pmass;
};


#endif
