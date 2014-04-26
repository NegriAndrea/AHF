/**
 * \file comm.h
 *
 * Provides functionality for communication between processes.
 */

/* Only needed in MPI mode */
#ifdef WITH_MPI

#ifndef COMM_H
#define COMM_H

/**********************************************************************\
 *    Includes                                                        *
\**********************************************************************/
#include <stdint.h>
#include "tdef.h"
#include "libamr_serial/get_nnodes.h"
#include "libutility/loadbalance.h"
#include "libio/io_logging.h"
#include "libgravity/gravity.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of global functions                                  * 
\**********************************************************************/

/**
 * \brief Will distribute the particles according to a given
 *        loadbalancing
 *
 * \param log        A logging module.
 * \param *fst_part  Pointer to the current particle storage.
 * \param *no_part   Pointer to the variable holding the local number of
 *                   particles (the length of the array, fst_part is
 *                   pointing to).
 * \param lb         The loadbalancing scheme after according to which
 *                   the distribution is to be done.
 *
 * \return Nothing, but the particle storage (and number) might be
 *         changed.
 */
extern void
comm_dist_part(io_logging_t log,
               partptr *fst_part,
               uint64_t *no_part,
               loadbalance_t lb);

/**
 * \brief Will distribute the particles in the boundary, hence
 *        will duplicate particles. This is needed for AHF.
 *
 * \param log        A logging module.
 * \param *fst_part  Pointer to the current particle storage.
 * \param *no_part   Pointer to the variable holding the local number of
 *                   particles (the length of the array, fst_part is
 *                   pointing to).
 * \param lb         The loadbalancing scheme after according to which
 *                   the distribution is to be done.
 *
 * \return Returns the number of the local duplicated particles,
 *         meaning those which do not actually belong to this
 *         CPU.
 */
extern uint64_t
comm_dist_part_ahf(io_logging_t log,
                   partptr *fst_part,
                   uint64_t *no_part,
                   loadbalance_t lb);

/**
 * \brief Will update the boundary values.
 *
 * \param log       A logging module.
 * \param lb        The current loadbalancing.
 * \param *curgrid  A pointer to the grid structure.
 * \param offset    Number of bytes to add to the pointer to the node
 *                  structure to reach the storage for the element that
 *                  should be updated.
 *
 * \return Nothing.
 */
extern void
comm_update_bound(io_logging_t log,
                  loadbalance_t lb,
                  gridls *curgrid,
                  int offset);

# ifdef WITH_FFTW
/**
 * \brief Will fill the fft slab correctly.
 *
 * \param log       A logging module.
 * \param *curgrid  A pointer to the grid structure.
 * \param fft       The FFT structure.
 * \param lb        The current loadbalancing.
 *
 * \return Nothing.
 */
extern void 
comm_fill_fft_slab(io_logging_t log,
                   gridls *curgrid,
                   fft_t fft,
                   loadbalance_t lb);

/**
 * \brief Will copy the potential from the FFT slab back to the grid
 *        structure.
 *
 * \param log       A logging module.
 * \param *curgrid  A pointer to the grid structure.
 * \param fft       The FFT structure.
 * \param lb        The current loadbalancing.
 *
 * \return Nothing.
 */
extern void 
comm_read_fft_slab(io_logging_t log,
                   gridls *curgrid,
                   fft_t fft,
                   loadbalance_t lb);
# endif /* WITH_FFTW */

#endif /* COMM_H */

#endif /* WITH_MPI */
