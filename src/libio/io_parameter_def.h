#ifndef IO_PARAMETER_DEF_H
#define IO_PARAMETER_DEF_H

/**
 * \file io_parameter_def.h
 *
 * Provides functions for reading in AMIGA parameter files.
 */


/*--- Includes ----------------------------------------------------------*/
#include "io_file.h"
#include <stdint.h>


/*--- Defines -----------------------------------------------------------*/
#define IO_PARAMETER_MAXSTRING 1024


/*--- Main structure definition -----------------------------------------*/

/**
 * The header structure itself
 */
struct io_parameter_struct {
	/* Required parameter */
	char           *icfile_name;
	io_file_type_t ic_filetype;
	uint32_t       reader;
	char           *outfile_prefix;
	int            NGRID_DOM;
  int            NGRID_MAX;
	double         Nth_dom;
	double         Nth_ref;
	int            UseRhoBack;
	double         UserDvir;
	double         MaxGatherRad;
	int            lb_level;
  double         AHF_VTUNE;
  int            AHF_MINPART;
    
  double         GADGET_m2Msunh;
  double         GADGET_l2Mpch;  // the GADGET units

#ifdef AHF_LRSI
	double         lrsi_beta;
	double         lrsi_r_s;
#endif
  
#if (defined AHFmixHaloIDandSnapID || defined SUSSING2013)
  uint64_t       isnap;
#endif

#ifdef DARK_ENERGY
	char          *defile_name;
#endif
};

/** Convenient typedef */
typedef struct io_parameter_struct io_parameter_struct_t;

/** Convenient typedef */
typedef io_parameter_struct_t      *io_parameter_t;


#endif /* IO_PARAMETER_DEF_H */
