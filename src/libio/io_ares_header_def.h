#ifndef IO_ARES_HEADER_DEF_H 
#define IO_ARES_HEADER_DEF_H

/* $Id: io_ares_header_def.h,v 1.4 2008/07/22 12:26:40 knolli Exp $ */

/**
 * \file io_ares_header_def.h
 *
 * Provides the structure definition for the ARES header
 * structure. Including useful typedefs.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/** 
 * The size (in bytes) reserved at the beginning of a file for
 * the header 
 */
#define ARES_HEADER_SIZE 512

/** The length of the file identifier */
#define ARES_LENGTH_IDENTIFIER 16

/** The used bytes of the header. */
#define ARES_HEADER_HEADERSIZE (   ARES_LENGTH_IDENTIFIER*sizeof(char) \
                                 +  7*sizeof(uint64_t)\
                                 +  3*sizeof(uint32_t)\
                                 + 10*sizeof(double)\
                                 +  3*sizeof(int))

/** The unused bytes of the header */
#define ARES_HEADER_FILLHEADER (ARES_HEADER_SIZE - ARES_HEADER_HEADERSIZE)

/** The file identifier string */
#define ARES_IDENTIFIER_STRING "ARES v/02.00.01"


/**
 * The header structure itself
 */
struct io_ares_header_struct {
	/** Used to figure out the swapping, contains the length of the id */
	uint32_t len_id;
	/** 
	 * Identifier string, in memory we store the final \0, but this is
	 * not written to the file, hence the +1 here.
	 */
	char id[ARES_LENGTH_IDENTIFIER+1];
	/** The total number of particles in the simulation */
	uint64_t no_part;
	/** The total number of particles in this file */
	uint64_t no_part_in_file;
	/** The number of bytes used to store floats in this file */
	uint32_t bytes_float;
	/** The number of bytes used to store integers in this file */
	uint32_t bytes_int;
	/** The number of particle species in the simulation */
	uint64_t no_species;
	/**
	 * Set to anything larger than 0 to indicate that the particles have
	 * a mass.
	 */
	uint64_t has_mass;
	/**
	 * Set to anything larger than 0 to indicate that gas particle
	 * energies are in the file
	 */
	uint64_t has_gas;
	/** The number of virtual particles (total mass in the simulation) */
	double no_vpart;
	/** The boxsize (in Mpc) of the simulation */
	double boxsize;
	double omega0;
	double lambda0;
	/** The mass of one particle of weight 1 */
	double pmass;
	/** The minimal weight of a particle occuring in the simulation */
	double minweight;
	/** The maximal weight of a particle occuring in the simulation */
	double maxweight;
	/** Initial expansion factor */
	double a_initial;
	/** Current expansion factor */
	double a_current;
	/** Current timestep */
	double timestep;
	/** 
	 * Stores the minimal key that belonged to the process writing this
	 * file
	 */
	uint64_t minkey;
	/** 
	 * Stores the maximal key that belonged to the process writing this
	 * file
	 */
	uint64_t maxkey;
	/**
	 * Stores the level use for the keys; total number of keys is then
	 * (2^level)^3.
	 */
	int lb_level;
	/** Stores the rank of the mpi process that wrote this file */
	int rank;
	/** Store the size of the MPI domain used to generate the file(s) */
	int size;
};

/** Convenient typedef */
typedef struct io_ares_header_struct io_ares_header_struct_t;

/** Convenient typedef */
typedef io_ares_header_struct_t *io_ares_header_t;


#endif /* IO_ARES_HEADER_DEF_H */
