#ifndef IO_GADGET_HEADER_DEF_H 
#define IO_GADGET_HEADER_DEF_H

/* $Id: io_gadget_header_def.h,v 1.2 2006/11/13 15:22:07 knolli Exp $ */

/**
 * \file io_gadget_header_def.h
 *
 * Provides the structure definition for the Gadget header
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
#define GADGET_HEADER_SIZE 256

/** Defines the total used size (in bytes) of the header */
#define GADGET_HEADER_HEADERSIZE ( 13*sizeof(int32_t)\
                                  +12*sizeof(uint32_t)\
                                  +12*sizeof(double))

/** Defines the number of unused bytes in the header */
#define GADGET_HEADER_FILLHEADER (  GADGET_HEADER_SIZE \
                                  - GADGET_HEADER_HEADERSIZE)

/**
 * The header structure itself
 */
struct io_gadget_header_struct {
	/** The number of particles in the file */
	int32_t np[6];
	/** The mass array */
	double massarr[6];
	double expansion;
	double redshift;
	int32_t flagsfr;
	int32_t flagfeedback;
	/** Total number of particles in the simulation */
	uint32_t nall[6];
	int32_t flagcooling;
	/** The number of files */
	int32_t numfiles;
	/** The boxsize */
	double boxsize;
	double omega0;
	double omegalambda;
	double hubbleparameter;
	int32_t flagstellarage;
	int32_t flagmetals;
	uint32_t nallhighw[6];
	int32_t flagentropyu;
	char unused[GADGET_HEADER_FILLHEADER+1];
};

/** Convenient typedef */
typedef struct io_gadget_header_struct io_gadget_header_struct_t;

/** Convenient typedef */
typedef io_gadget_header_struct_t *io_gadget_header_t;


#endif /* IO_GADGET_HEADER_DEF_H */
