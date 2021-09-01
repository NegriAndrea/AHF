// Copyright (C) 2010, Steffen Knollmann
// Released under the terms of the GNU General Public License version 3.
// This file is part of `ginnungagap'.


/*--- Includes ----------------------------------------------------------*/
#include "timer.h"
#include <stdio.h>
#include <assert.h>
#if (defined WITH_MPI)
#  include <mpi.h>
#elif (defined _OPENMP)
#  include <omp.h>
#else
#  include <time.h>
#endif


/*--- Local variables ---------------------------------------------------*/
#if (!defined WITH_MPI && !defined _OPENMP)
static double CPS_INV = 1. / ((double)CLOCKS_PER_SEC);
#endif


/*--- Local defines -----------------------------------------------------*/


/*--- Prototypes of local functions -------------------------------------*/


/*--- Implementations of exported functios ------------------------------*/
extern double
timer_getTime(void)
{
	double timing;

#if (defined WITH_MPI)
	timing = MPI_Wtime();
#elif (defined _OPENMP)
	timing = omp_get_wtime();
#else
	timing = clock() * CPS_INV;
#endif

	return timing;
}

/*--- Implementations of local functions --------------------------------*/
