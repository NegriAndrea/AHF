#ifndef COMMON_H
#define COMMON_H

/* $Id: common.h,v 1.22 2008/04/29 12:47:22 aknebe Exp $ */

/**
 * \file common.h
 *
 * Defines all global variables. This header should be included by all
 * files that whish to manipulate the global variables.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include "tdef.h"
#include "param.h"
#include "libio/io.h"

#ifdef WITH_MPI
#	include <mpi.h>
#	include "libutility/loadbalance.h"
#endif
#	include "libsfc/sfc.h"


/***********************************************************************\
 *    Global typedefs                                                  * 
\***********************************************************************/

#ifdef DOUBLE
/** Sets the type of stored floating point values */
typedef double fpv_t;
#else
/** Sets the type of stored floating point values */
typedef float fpv_t;
#endif

/** Defines a vector of floating point values of the stored type */
typedef fpv_t vec_t[NDIM];

/** Specifically defines a vector of floats */
typedef float fvec_t[NDIM];

/** Specifically defines a vector of doubles */
typedef double dvec_t[NDIM];

typedef struct global_io_struct global_io_t;
typedef struct global_info_struct global_info_t;

#ifdef WITH_MPI
/** Convenient typedef */
typedef struct global_mpi_struct global_mpi_t;
#endif


/***********************************************************************\
 *    Global structure definitions                                     * 
\***********************************************************************/

/**
 * \brief Describes the IO related global information
 */
struct global_io_struct {
	/** Holds the log module */
	io_logging_t log;
	/** Holds the information read from a parameter file */
	io_parameter_t params;
	/** Holds the file module */
	io_file_t file;
};

/**
 * \brief Describes general global information
 */
struct global_info_struct {
	/** Gives access to all particles */
	partptr fst_part;
	/** Remembers how many particles are attached to fst_part */
	uint64_t no_part;
#ifdef WITH_MPI
	/** Stores the current loadbalancing scheme */
	loadbalance_t loadbal;
#endif
	sfc_key_t minkey;
	sfc_key_t maxkey;
	sfc_curve_t ctype;
	int level; 
	int rank;
	int size;
};

#ifdef WITH_MPI
/**
 * \brief Global structure holding the MPI informations
 */
struct global_mpi_struct {
	/** The number of this process */
	int rank;
	/** The total number of processes */
	int size;
	/** A summation operation that works with long long datatypes */
	MPI_Op op_sum;
	/** Will hold the MPI realisation of an int32_t */
	MPI_Datatype dt_int32;
	/** Will hold the MPI realisation of an uint32_t */
	MPI_Datatype dt_uint32;
	/** Will hold the MPI realisation of an int64_t */
	MPI_Datatype dt_int64;
	/** Will hold the MPI realisation of an uint64_t */
	MPI_Datatype dt_uint64;
#	ifdef MPI_TIMING
	/** Keeps track of the timing */
	double start;
	/** Keeps track of the timing */
	double stop;
#	endif
};
#endif


/***********************************************************************\
 *    Global variables                                                 * 
\***********************************************************************/
/* information about timings                    */
extern info_timing   timing;
/* important information about the domain grid  */
extern info_global   global;
/* all information regarding output files       */
extern info_io       io;
/* information about AHF                        */
extern info_ahf      ahf;
/* Layzer-Irvine energy check data              */ 
extern info_energy   energy;
/* on-the-fly power spectrum P(k)               */
extern info_Pk       PkSpectrum;
/* information regarding GADGET files           */
extern info_gadget   gadget;
/* contains simulation parameter */
extern param_simu    simu;
/* unit matrix */
extern int           unit_matrix[NDIM][NDIM];

extern global_io_t global_io;
extern global_info_t global_info;

#ifdef WITH_MPI
/* stores the MPI related quantities            */ 
extern global_mpi_t  global_mpi;
#endif


/***********************************************************************\
 *    Global functions                                                 * 
\***********************************************************************/
#ifdef WITH_MPI
/**
 * \brief Fires up the MPI environment and sets the global_mpi things.
 *
 * \param int *argc     Pointer to the argc parameter.
 * \param char ***argv  Pointer to the argv parameter.
 *
 * \return Nothing.
 */
extern void
common_initmpi(int *argc, char ***argv);
#endif

/**
 * \brief Function to terminate the program at any point. Basically a
 *        wrapper around exit() but taking into account that in MPI mode
 *        some other things need to be done too.
 *
 * If the supplied error code is something else then 0, an abnormal
 * program termination is suspected and hence in MPI mode an abort is
 * called. Otherwise the MPI environment is finalized.
 *
 * This function should mainly be used to abnormally terminate the
 * program, as in MPI mode it is not working well for as a debugging
 * method.
 *
 * \param errcode An error code.
 *
 * \return Never returns.
 */
extern void
common_terminate(int errcode);


#endif /* COMMON_H */
