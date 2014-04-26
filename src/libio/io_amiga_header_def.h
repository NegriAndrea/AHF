#ifndef IO_AMIGA_HEADER_DEF_H 
#define IO_AMIGA_HEADER_DEF_H

/* $Id: io_amiga_header_def.h,v 1.5 2007/02/07 14:31:03 knolli Exp $ */

/**
 * \file io_amiga_header_def.h
 *
 * Provides the structure definition for the AMIGA header
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
#define AMIGA_HEADER_SIZE 2048

/** Defines the length for the identifying header string */
#define AMIGA_HEADER_HEADERSTRING 256

/** Defines the total used size (in bytes) of the header */
#define AMIGA_HEADER_HEADERSIZE ( AMIGA_HEADER_HEADERSTRING*sizeof(char)\
                                 +2*sizeof(long)+4*sizeof(int32_t)\
                                 +41*sizeof(double))

/** Defines the number of unused bytes in the header */
#define AMIGA_HEADER_FILLHEADER (  AMIGA_HEADER_SIZE \
                                 - AMIGA_HEADER_HEADERSIZE)

/**
 * The header structure itself
 */
struct io_amiga_header_struct {
	char          header[AMIGA_HEADER_HEADERSTRING];

	int32_t       multi_mass;
	int32_t       double_precision;

	long          no_part;
	long          no_species;
	double        no_vpart;

	double        timestep;
	int32_t       no_timestep;

	double        boxsize;
	double        omega0;
	double        lambda0;
	double        pmass;

	double        cur_reflevel;
	double        cur_frcres;

	double        a_initial;
	double        a_current;

	double        K_initial;
	double        K_current;
	double        U_initial;
	double        U_current;
	double        Eintegral;
	double        Econst;

	double        paramNSTEPS;
	double        paramNGRID_DOM;
	double        paramNth_dom;
	double        paramNth_ref;
	double        paramE_UPDATE;
	double        paramCELLFRAC_MAX;
	double        paramCELLFRAC_MIN;
	double        paramCA_CRIT;
	double        paramMAX_L1DIM;
	double        paramDOMSWEEPS;
	double        paramREFSWEEPS;

	double        paramAHF_MINPART;
	double        paramAHF_VTUNE;
	double        paramAHF_RISE;
	double        paramAHF_SLOPE;
	double        paramAHF_MAXNRISE;

   double        min_weight;
   double        max_weight;
   double        t_unit;
   double        B_init;                  // relevant for MHD solver
   double        param_dummy5;
   double        param_dummy6;
   double        param_dummy7;
   double        param_dummy8;            // empty space that once was used (downwards compatibility!)

   double        version;
	 int32_t       build;

	char          dummy[AMIGA_HEADER_FILLHEADER];
};

/** Convenient typedef */
typedef struct io_amiga_header_struct io_amiga_header_struct_t;

/** Convenient typedef */
typedef io_amiga_header_struct_t *io_amiga_header_t;


#endif /* IO_AMIGA_HEADER_DEF_H */
