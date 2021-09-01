#ifndef SFC_CURVE_H
#define SFC_CURVE_H

/* $Id: sfc_curve.h,v 1.3 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file sfc_curve.h
 *
 * Definitions of the space filling curve functions and types.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdint.h>


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/** Descriptive string of the Hilbert curve type */
#define SFC_CURVE_HILBERT_STR "Hilbert curve"

/** Descriptive string for an unkown curve type */
#define SFC_CURVE_UNKOWN_STR "Unkown curve type"

/** This defines the type of the space filling curve */
typedef enum {
	/** For Hilbert Curves */
	SFC_CURVE_HILBERT = 0,
	SFC_CURVE_PLAIN = 1
} sfc_curve_t;

/** Used to print a SFC key */
#define SFC_PRIkey PRIu64

/** The SFC key type */
typedef uint64_t sfc_key_t;


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/
/**
 * \brief Returns a string describing the curve type.
 *
 * \param ctype  The curve type to describe.
 *
 * \return A static string describing the curve type. The calling
 *         function must not try to change the string.
 */
extern const char*
sfc_curve_typestr(sfc_curve_t ctype);

/**
 * \brief Contract a given SFC index of a given depth to a lower
 *        depth.
 *
 * \param trgt_level  The requested level to contract to.
 * \param src_level   The current level of the input Key.
 * \param ctype       The type of the space filling curve.
 * \param key         The key.
 *
 * \return Returns the contracted SFC Key.
 */
extern sfc_key_t
sfc_curve_contract(uint32_t trgt_level,
                   uint32_t src_level,
                   sfc_curve_t ctype,
                   sfc_key_t key);

/**
 * \brief Prolongs a given SFC index of a given depth to a higher
 *        depth and returns the smallest key that would contracted back
 *        to the original depths yields the initial key.
 *
 * \param trgt_level  The requested level to prolong to.
 * \param src_level   The current level of the input Key.
 * \param ctype       The type of the space filling curve.
 * \param key         The key.
 *
 * \return Returns the contracted SFC Key.
 */
extern sfc_key_t
sfc_curve_prolongMin(uint32_t trgt_level,
                     uint32_t src_level,
                     sfc_curve_t ctype,
                     sfc_key_t key);

/**
 * \brief Prolongs a given SFC index of a given depth to a higher
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
extern sfc_key_t
sfc_curve_prolongMax(uint32_t trgt_level,
                     uint32_t src_level,
                     sfc_curve_t ctype,
                     sfc_key_t key);

/**
 * \brief Takes the position of a particle and projects that to the
 *        according SFC Index at a given level.
 *
 * \param ctype  The type of the space filling curve.
 * \param x      The x coordinate, this must be in [0,1[
 * \param y      The y coordinate, this must be in [0,1[
 * \param z      The z coordinate, this must be in [0,1[
 * \param bits   The number of bits used per dimension for the SFC
 *               Key; this depends on the exact definition of sfc_key_t:
 *               sfc_key_t must be sufficiently large to accomodate
 *               3*bits. Hence for a 64bit sfc_key_t a maximum of
 *               21 bits can be used per dimension.
 *
 * \return Returns the according SFC Key.
 */
extern sfc_key_t
sfc_curve_calcKey(sfc_curve_t ctype,
                  double x,
                  double y,
                  double z,
                  uint32_t bits);
/**
 * \brief Takes the grid position and converts that to the
 *        according SFC Index at a given level.
 *
 * \param ctype  The type of the space filling curve.
 * \param x      The x coordinate, must be in [0,1<<bits[
 * \param y      The y coordinate, must be in [0,1<<bits[
 * \param z      The z coordinate, must be in [0,1<<bits[
 * \param bits   The number of bits used per dimension for the SFC Key.
 *               This is effectively the grid level.
 *
 * \return Returns the according SFC Key.
 */
extern sfc_key_t
sfc_curve_calcKey_grid(sfc_curve_t ctype,
                       uint32_t x,
                       uint32_t y,
                       uint32_t z,
                       uint32_t bits);

/**
 * \brief Takes a key and converts that to position integers.
 *
 * \param ctype  The type of the space filling curve.
 * \param key    The SFC key which is to be converted to position
 *               integers.
 * \param bits   The number of bits used per dimension for the SFC
 *               Key; this depends on the exact definition of sfc_key_t:
 *               sfc_key_t must be sufficiently large to accomodate
 *               3*bits. Hence for a 64bit sfc_key_t a maximum of
 *               21 bits can be used per dimension.
 * \param *pos   Pointer to the array supposed to hold the resultin
 *               position integers. Must be at least 3 elements long.
 *
 * \return Nothing, but the array pos is pointing to will be filled with
 *         the position integers.
 */
extern void
sfc_curve_calcPos(sfc_curve_t ctype,
                  sfc_key_t key,
                  uint32_t bits,
                  uint32_t *pos);

/**
 * \brief Takes a SFC Key and gives back the Keys of all neighbouring
 *        cells.
 *
 * \param base    The SFC Key of the center cell
 * \param *shell  A pointer to a 27 element long array to hold all
 *                Keys of the shell cells including the center shell.
 *                The ordering is as follows: It starts at the lower
 *                left front corner (x0-1, x1-1, x2-1) and then
 *                proceeds first in the x2 direction, then the x1 and
 *                finally the x0 direction. So the first few elemets
 *                are [(x0-1,x0-1,x2-1), (x0-1,x1-1,x2),
 *                (x0-1,x1-1,x2+1), (x0-1,x1,x2-1), (x0-1,x1,x2),...].
 * \param bits    The number of bits used per dimension for the SFC
 *                Key.
 *
 * \return Nothing.
 */
extern void
sfc_curve_getShell(sfc_curve_t ctype,
                   sfc_key_t base,
                   sfc_key_t *shell,
                   uint32_t bits);

/**
 * \brief A helper function to compare two keys.
 *
 * \param *k1  The first key.
 * \param *k2  The second key.
 *
 * \return If k1<k2, -1 is returned, if k1 and k2 match, 0 is returned,
 *         1 otherwise.
 */
extern int
sfc_curve_comp_key(const void *k1, const void *k2);


#endif /* SFC_CURVE_H */
