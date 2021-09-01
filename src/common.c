/* $Id: common.c,v 1.16 2008/04/29 12:47:22 aknebe Exp $ */

/**
 * \file common.c
 *
 * Implements all global variables.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <stdlib.h>

#ifndef COMMON_INCLUDED
#define COMMON_INCLUDED
#include "common.h"


/***********************************************************************\
 *    Implementation of global variables                               * 
\***********************************************************************/
info_timing   timing;
info_global   global;
info_io       io;            
info_ahf      ahf;
info_energy   energy;        
info_Pk       PkSpectrum;        
info_gadget   gadget;        
param_simu    simu;
int           unit_matrix[NDIM][NDIM];

global_io_t global_io;
global_info_t global_info;

#ifdef WITH_MPI
global_mpi_t global_mpi;


/***********************************************************************\
 *    Prototypes of local functions                                    * 
\***********************************************************************/
static void
local_mpi_sum(void *invec,
              void *inoutvec,
              int *len,
              MPI_Datatype *datatype);
#endif


/***********************************************************************\
 *    Implementation of global functions                               * 
\***********************************************************************/

#ifdef WITH_MPI
extern void
common_initmpi(int *argc, char ***argv)
{
	/* Do the standard MPI_Init */
	MPI_Init(argc, argv);

	/* And set the global information */
	MPI_Comm_size(MPI_COMM_WORLD, &global_mpi.size);
	MPI_Comm_rank(MPI_COMM_WORLD, &global_mpi.rank);

	/* Generate MPI type for the extended integers */
	MPI_Type_contiguous(4, MPI_CHAR, &(global_mpi.dt_int32));
	MPI_Type_contiguous(4, MPI_CHAR, &(global_mpi.dt_uint32));
	MPI_Type_contiguous(8, MPI_CHAR, &(global_mpi.dt_int64));
	MPI_Type_contiguous(8, MPI_CHAR, &(global_mpi.dt_uint64));

	/* Activate them */
	MPI_Type_commit(&(global_mpi.dt_int32));
	MPI_Type_commit(&(global_mpi.dt_uint32));
	MPI_Type_commit(&(global_mpi.dt_int64));
	MPI_Type_commit(&(global_mpi.dt_uint64));

	/* Create needed operation for the new datatypes */
	MPI_Op_create(local_mpi_sum, 1, &(global_mpi.op_sum));

	return;
}
#endif /* WITH_MPI */

extern void
common_terminate(int errcode)
{
#	ifdef WITH_MPI
	MPI_Op_free(&(global_mpi.op_sum));

	MPI_Type_free(&(global_mpi.dt_uint64));
	MPI_Type_free(&(global_mpi.dt_int64));
	MPI_Type_free(&(global_mpi.dt_uint32));
	MPI_Type_free(&(global_mpi.dt_int32));

	if (errcode != 0)
		MPI_Abort(MPI_COMM_WORLD, errcode);
	else
		MPI_Finalize();
#	endif

	exit(errcode);
}


#ifdef WITH_MPI
/***********************************************************************\
 *    Implementation of local functions                                * 
\***********************************************************************/
static void
local_mpi_sum(void *invec,
              void *inoutvec,
              int *len,
              MPI_Datatype *datatype)
{
	int i;

	if  (*datatype == global_mpi.dt_int32) {
#		ifdef MPI_DEBUG
		io_logging_msg(global_io.log, INT32_C(10),
		               "Summing up %i int32_t things.",
		               *len);
#		endif
		for (i=0; i<*len; i++)
			((int32_t *)inoutvec)[i] += ((int32_t *)invec)[i];
		return;
	}
	if  (*datatype == global_mpi.dt_uint32) {
#		ifdef MPI_DEBUG
		io_logging_msg(global_io.log, INT32_C(10),
		               "Summing up %i uint32_t things.",
		               *len);
#		endif
		for (i=0; i<*len; i++)
			((uint32_t *)inoutvec)[i] += ((uint32_t *)invec)[i];
		return;
	}
	if  (*datatype == global_mpi.dt_int64) {
#		ifdef MPI_DEBUG
		io_logging_msg(global_io.log, INT32_C(10),
		               "Summing up %i int64_t things.",
		               *len);
#		endif
		for (i=0; i<*len; i++)
			((int64_t *)inoutvec)[i] += ((int64_t *)invec)[i];
		return;
	}
	if  (*datatype == global_mpi.dt_uint64) {
#		ifdef MPI_DEBUG
		io_logging_msg(global_io.log, INT32_C(10),
		               "Summing up %i uint64_t things.",
		               *len);
#		endif
		for (i=0; i<*len; i++)
			((uint64_t *)inoutvec)[i] += ((uint64_t *)invec)[i];
		return;
	}

	io_logging_fatal(global_io.log,
	                 "%s not defined for this dataype",
	                 __func__);
	common_terminate(EXIT_FAILURE);

	return;
}
#endif /* WITH_MPI */

#endif /* COMMON_INCLUDED */
