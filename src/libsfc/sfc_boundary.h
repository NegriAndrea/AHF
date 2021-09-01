#ifndef SFC_BOUNDARY_H
#define SFC_BOUNDARY_H

/* $Id: sfc_boundary.h,v 1.4 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file sfc_boundary.h
 *
 * Provides a way to figure out all bounding cells of a given space
 * filling curve segment.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <stdint.h>
#include <stdbool.h>
#include "../libio/io_logging.h"
#include "sfc_curve.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/

/** Descriptive string of the inner boundary type */
#define SFC_BOUNDARY_TYPE_INNER_STR "Inner Boundary"

/** Descriptive string of the inner boundary type */
#define SFC_BOUNDARY_TYPE_OUTER_STR "Outer Boundary"

/** Descriptive string for an unknown boundary type */
#define SFC_BOUNDARY_TYPE_UNKNOWN_STR "Unknown Boundary Type"

/** Defines the type of the boundary */
typedef enum {
	/**
	 * This type contains all cells belonging to a given segment and
	 * situated on the surface of the volume filled by the segment. So
	 * all thees boundary cells have at least one neighbour which does
	 * not belong to the segment.
	 */
	SFC_BOUNDARY_TYPE_INNER = 0,
	/**
	 * This type contains all cells not belonging to a given segment
	 * but situated on the surface of the volume filled by the
	 * segment. So each cell of these boundary cells have at least one
	 * neighbour which belongs to the segment.
	 */
	SFC_BOUNDARY_TYPE_OUTER = 1
} sfc_boundary_type_t;

/** The boundary structure */
struct sfc_boundary_struct {
	/** The type of the boundary */
	sfc_boundary_type_t btype;
	/** Number of boundary elements */
	uint64_t num;
	/** Length of boundary element list, that might be more than num */
	uint64_t len;
	/** The increment by which to increase the boundary array */
	uint32_t inc;
	/** The keys of all boundary cells sorted from small to large */
	sfc_key_t *bound;
	/** The smallest key that still belongs to the segment */
	sfc_key_t minkey;
	/** The largest key that still belongs to the segment. */
	sfc_key_t maxkey;
	/** Remember for what sfc type this is the boundary */
	sfc_curve_t ctype;
	/** Stores the number of bits used per dimension for the SFC Key */
	uint32_t bits;
};

/** Convenient typedef */
typedef struct sfc_boundary_struct sfc_boundary_struct_t;

/** Convenient typedef */
typedef struct sfc_boundary_struct *sfc_boundary_t;

/** The 2nd boundary structure that describes everything needed */
struct sfc_boundary_2_struct {
	/** The outer boundary */
	sfc_boundary_t outer;
	/** The number of inner boundaries (one for each CPU) */
	uint32_t num_inner;
	/** The inner boundaries */
	sfc_boundary_t inner;
};

/** Convenient typedef */
typedef struct sfc_boundary_2_struct sfc_boundary_2_struct_t;

/** Convenient typedef */
typedef sfc_boundary_2_struct_t *sfc_boundary_2_t;


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Returns a string describing the boundary type.
 *
 * \param btype  The boundary type to describe.
 *
 * \return A static string describing the boundary type. The calling
 *         function must not try to change the string.
 */
extern const char*
sfc_boundary_typestr(sfc_boundary_type_t btype);

/**
 * \brief Finds all boundary cells of a given space filling
 *        curve segment.
 *
 * \param log     A logging modul.
 * \param btype   The type of boundary to create.
 * \param minkey  The smallest key that still belongs to the
 *                segment.
 * \param maxkey  The largest key that still belongs to the segment.
 * \param bits    The number of bits used per dimension for the SFC key.
 * \param ctype   The type of the SFC curve.
 *
 * \return Returns the freshly created boundary of the given space
 *         filling curve segment.
 */
extern sfc_boundary_t
sfc_boundary_get(io_logging_t log,
                 sfc_boundary_type_t btype,
                 sfc_key_t minkey,
                 sfc_key_t maxkey,
                 uint32_t bits,
                 sfc_curve_t ctype);

/**
 * \brief Delete the boundary structure.
 *
 * \param *bound  Pointer to the variable holding the boundary. Will be
 *                set to NULL;
 */
extern void
sfc_boundary_del(sfc_boundary_t *bound);

/**
 * \brief Finds all boundary cells of a given space filling segment,
 *        taking into account that a inner boundary cell can be an outer
 *        boundary cell on more than one CPU.
 *
 * \param log      A logging module.
 * \param rank     The CPU number for which the boundary information
 *                 needs to be found. Selects the SFC segment out of the
 *                 lstkey and fstkey arrays.
 * \param ncpu     The total number of CPU involved in the game, hence
 *                 the length of the following arrays.
 * \param *fstkey  An array with all the first keys according to the
 *                 domain decomposition sorted by CPU number.
 * \param *lstkey  An array with all the last keys according to the
 *                 domain decomposition sorted by CPU number.
 * \param bits     The number of bits used per dimension for the SFC
 *                 key.
 * \param ctype    The type of the SFC curve.
 *
 * \return Returns the freshly created boundary of the given space
 *         filling curve segment.
 */
extern sfc_boundary_2_t
sfc_boundary_2_get(io_logging_t log,
                   uint32_t rank,
                   uint32_t ncpu,
                   sfc_key_t *fstkey,
                   sfc_key_t *lstkey,
                   uint32_t bits,
                   sfc_curve_t ctype);
/**
 * \brief Delete the boundary structure 2nd order.
 *
 * \param *bound Pointer to the variable holding the boundary. Will be
 *                set to NULL.
 */
extern void
sfc_boundary_2_del(sfc_boundary_2_t *bound);

/**
 * \brief Delete the boundary structure of 2nd order.
 *
 * \param *bound  Pointer to the variable holding the boundary. Will be
 *                set to NULL.
 */

/**
 * \brief Adds a key to the boundary
 *
 * This is done by searching through the whole key list and inserting
 * the key at the right place. If the key is already in there, nothing
 * happens.
 *
 * \param log    A logging module.
 * \param bound  The boundary to which to add the key.
 * \param key    The key to be added.
 *
 * \return  Returns the updated boundary, or NULL if the resizing of the
 *          key array did not succeed, in which case the original
 *          boundary structure remains intact.
 */
extern sfc_boundary_t
sfc_boundary_addKey(io_logging_t log,
                    sfc_boundary_t bound,
                    sfc_key_t key);

/**
 * \brief Finds out where a given key is in the boundary, or where it
 *        should be entered.
 *
 * \param log      A logging module.
 * \param bound    The boundary to look in.
 * \param key      The key to look for.
 * \param *retpos  A pointer to the variable that will recieve the
 *                 position the key was found, or, if it was not found,
 *                 the position where it should be entered. This can be
 *                 NULL.
 * 
 * \return Returns 'true' if the key is in the boundary, 'false'
 *         otherwise. Depending on that, *retpos will have a diferent
 *         meaning.
 */
extern bool
sfc_boundary_findKeyPos(io_logging_t log,
                        sfc_boundary_t bound,
                        sfc_key_t key,
                        uint64_t *retpos);

/**
 * \brief Returns a box which includes all points in the boundars
 *
 * The box is given by the position of two points, consisting of the
 * minimal and maximal, repsectively, coordinate.
 *
 * \param bound  The boundary to be boxed, can be more than one but must
 *               be than consecutive in memory (bound, bound+1, bound+2,
 *               ..., bound+(num-1)).
 * \param num    The number of boundaries.
 * \param *p1    First point (smallest coordinates), must be a pointer
 *               to an array large enough to hold the coordinates.
 * \param *p2    Second point (largest coordinates), size as *p1.
 *
 * \return Nothing, but the external variables p1 and p2 are pointing to
 *         will be changed.
 */
extern void
sfc_boundary_getBox(sfc_boundary_t bound,
                    uint32_t num,
                    uint32_t *p1,
                    uint32_t *p2);

/**
 * \brief Write the boundary to the logfile.
 *
 * \param log    The logfile.
 * \param bound  The boundary.
 */
extern void
sfc_boundary_log(io_logging_t log,
                 sfc_boundary_t bound);


#endif /* SFC_BOUNDARY_H */
