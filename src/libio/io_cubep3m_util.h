#ifndef IO_CUBEP3M_UTIL__H
#define IO_CUBEP3M_UTIL_H

/**
 * \file io_cubep3m_util.h
 *
 * Provides utility functions for CubeP3M related things.
 */


/***********************************************************************
 *    Includes                                                         *
 ***********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "io_logging.h"


/***********************************************************************
 *    Required defines                                                 *
 ***********************************************************************/
#ifndef RHO_CRIT_0
#  define RHO_CRIT_0 2.7755397e11   // [h^2*Msun]/[Mpc^3]
#endif

#ifndef H0
#  define H0 100.                   // [h*km]/[sec*Mpc]
#endif


/***********************************************************************
 *    Prototypes of global functions                                   *
 ***********************************************************************/


/***********************************************************************
 *    Inline functions                                                 *
 ***********************************************************************/

/**
 * \brief  Calculates the conversion factor for positions from internal
 *         CubeP3M units to Mpc/h.
 *
 * \param[in]  boxsize
 *                The boxsize in Mpc/h.
 * \param[in]  ngrid
 *                The size of the grid.
 *
 * \return  Returns the conversion factor.
 */
inline static double
io_cubep3m_util_lunit(double boxsize, uint64_t ngrid)
{
	return boxsize / ((double)ngrid);
}

/**
 * \brief  Calculates the conversion factor for masses from internal CubeP3M
 *         units to M_sun/h.
 *
 * \param[in]  boxsize
 *                The boxsize in Mpc/h.
 * \param[in]  omega0
 *                The matter density in units of the critical density.
 * \param[in]  nptotal
 *                The total number of particles used to cover the full
 *                simulation box.
 * \param[in]  weight
 *                The weight of one particle, this should normally be 8. as
 *                one particle is formed of 8 cells.
 *
 * \return  Returns the converison factor,
 */
inline static double
io_cubep3m_util_munit(double   boxsize,
                      double   omega0,
                      uint64_t nptotal,
                      double   weight)
{
	return boxsize * boxsize * boxsize * omega0 * RHO_CRIT_0
	       / ((double)(nptotal) * weight);
}

/**
 * \brief  Calculates the conversion factor for velocities from internal
 *         CubeP3M units to km/s.
 *
 * \param[in]  boxsize
 *                The boxsize in Mpc/h.
 * \param[in]  omega0
 *                The matter density in units of the critical density.
 * \param[in]  a
 *                The expansion factor.
 * \param[in]  ngrid
 *                The size of the grid.
 *
 * \return  Returns the converison factor,
 */
inline static double
io_cubep3m_util_vunit(double   boxsize,
                      double   omega0,
                      double   a,
                      uint64_t ngrid)
{
	return boxsize * 1.5 * sqrt(omega0) * H0 / (ngrid * a);
}

#endif
