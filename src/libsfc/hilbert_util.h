#ifndef HILBERT_UTIL_H
#define HILBERT_UTIL_H

/* $Id: hilbert_util.h,v 1.4 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file hilbert_util.h
 *
 * Useful Hilbert Index functions.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include "hilbert.h"
#include <stdint.h>


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/

/** Define the datatype for the Hilbert Key */
typedef bitmask_t hikey_t;


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Contract a given Hilbert index of a given depth to a lower
 *        depth.
 *
 * Done via fast bitshift operations.
 *
 * \param trgt_level  The requested level to contract to.
 * \param src_level   The current level of the input Hilbert Key.
 * \param hikey       The key.
 *
 * \return Returns the contracted Hilbert Key
 */
extern hikey_t
hilbert_util_contract(unsigned int trgt_level,
                      unsigned int src_level,
                      hikey_t hikey);

/**
 * \brief Prolongs a given Hilber index of a given depth to a higher
 *        depth and returns the smallest key that would contracted back
 *        to the original depths yields the initial key.
 *
 * \param trgt_level  The requested level to prolong to.
 * \param src_level   The current level of the input Key.
 * \param key         The key.
 *
 * \return Returns the contracted SFC Key.
 */
extern hikey_t
hilbert_util_prolongMin(unsigned int trgt_level,
                        unsigned int src_level,
                        hikey_t key);

/**
 * \brief Prolongs a given Hilbert index of a given depth to a higher
 *        depth and returns the largest key that would contracted back
 *        to the original depths yields the initial key.
 *
 * \param trgt_level  The requested level to prolong to.
 * \param src_level   The current level of the input Key.
 * \param ctype       The type of the space filling curve.
 * \param key         The key.
 *
 * \return Returns the contracted SFC Key.
 */
extern hikey_t
hilbert_util_prolongMax(unsigned int trgt_level,
                        unsigned int src_level,
                        hikey_t key);


/**
 * \brief Takes the position of a particle and projects that to the
 *        according Hilbert Index at a given level.
 *
 * \param x     The x coordinate, this must be in [0,1[
 * \param y     The y coordinate, this must be in [0,1[
 * \param z     The z coordinate, this must be in [0,1[
 * \param bits  The number of bits used per dimension for the Hilbert
 *              Key; this depends on the exact definition of hikey_t and
 *              hikey_t must be sufficiently large to accomodate 3*bits.
 *              Hence for a 64bit hikey_t a maximum of 21 bits can be
 *              used per dimension.
 *
 * \return Returns the according Hilbert Key.
 */
extern hikey_t
hilbert_util_calcHikey(double x,
                       double y,
                       double z,
                       unsigned bits);

/**
 * \brief Takes a grid position (x,y,z) and converts that to the
 *        according Hilbert key index at a given grid level.
 *
 * \param x     The x coordinate, must be in [0,1<<bits[
 * \param y     The y coordinate, must be in [0,1<<bits[
 * \param z     The z coordinate, must be in [0,1<<bits[
 * \param bits  The number of bits used per dimension for the Hilbert
 *              Key, this gives effectively the grid level.
 */
extern hikey_t
hilbert_util_calcHikey_grid(uint32_t x,
                            uint32_t y,
                            uint32_t z,
                            unsigned bits);

/**
 * \brief Converts a given key to positiono integers.
 *
 * \param key   The Key to convert.
 * \param bits  Number of bits used per dimension.
 * \param *pos  Pointer to the array to hold the position integers.
 *
 * \return Nothing, but the array pos is pointing to will be filled with
 *         the position integers.
 */
extern void
hilbert_util_calcPos(hikey_t key, unsigned bits, unsigned *pos);

/**
 * \brief Takes a Hilbert Key and gives back the Hilbert Keys of all
 *        neighbouring cells.
 *
 * \param base    The Hilbert Key of the center cell
 * \param *shell  A pointer to a 27 element long array to hold all
 *                Hilbert Keys of the shell cells including the center
 *                shell. The ordering is as follows: It starts at the
 *                lower left front corner (x0-1, x1-1, x2-1) and then
 *                proceeds first in the x2 direction, then the x1 and
 *                finally the x0 direction. So the first few elemets
 *                are [(x0-1,x0-1,x2-1), (x0-1,x1-1,x2),
 *                (x0-1,x1-1,x2+1), (x0-1,x1,x2-1), (x0-1,x1,x2),...].
 * \param bits    The number of bits used per dimension for the Hilbert
 *                Key.
 *
 * \return Nothing.
 */
extern void
hilbert_util_getShell(hikey_t base, hikey_t *shell, unsigned bits);


#endif /* HILBERT_UTIL_H */
