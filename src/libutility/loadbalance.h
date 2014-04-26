#ifndef LOADBALANCE_H
#define LOADBALANCE_H

/**
 * \file loadbalance.h
 *
 * Provides functionality for loadbalancing.
 */

/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdint.h>
#include "../tdef.h"
#include "../libsfc/sfc.h"
#include "../libio/io.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/** Descriptive string of the scheme 0 */
#define LOADBALANCE_EQUALVOL_STR \
	"equalvol: Same number of cells per cpu"

/** Descriptive string of the scheme 1 */
#define LOADBALANCE_EQUALPART_STR \
	"equalpart: Same number of particles per cpu"

/** This defines the loadbalancing scheme employed */
typedef enum {
	/** Selects equal volume scheme */
	LOADBALANCE_EQUALVOL = 0,
	/** Selects equal particle scheme */
	LOADBALANCE_EQUALPART = 1
} loadbalance_scheme_t;

/**
 * \brief Holds detailed information about the loadbalancing
 */
struct loadbalance_struct {
	/** Which scheme are we using */
	loadbalance_scheme_t scheme;
	/** What kind of space filling curve is currently used */
	sfc_curve_t ctype;
	/** Level on which the loadbalancing happen */
	int level;
	/** First SFC key occuring (meant to allow for segemented loadbal) */
	sfc_key_t startkey;
	/** Number of total keys of the SFC and the length of the next list */
	uint64_t totkeys;
	/** This holds the global information how many particles are where */
	uint32_t *bf;
	/** This holds the local information how many particles are where */
	uint32_t *loc_bf;
	/** Number of CPUs (and hence the length of the following lists) */
	int ncpu;
	/** Stores for every CPU the first key it is supposed to hold */
	sfc_key_t *fstkey;
	/** Stores for every CPU the last key it is supposed to hold */
	sfc_key_t *lstkey;
	/** Stores for every CPU the number of particles it is supposed to *\
	\*  hold                                                           */
	uint64_t *no_parts;
	/** Stores for every CPU the number of particles we hold locally */
	uint64_t *no_parts_loc;
	/** Stores the boundary */
	sfc_boundary_2_t bound;
};

/** Convenient short name */
typedef struct loadbalance_struct * loadbalance_t;


/**********************************************************************\
 *    Prototypes of global functions                                  * 
\**********************************************************************/

/**
 * \brief Returns a string describing the scheme.
 *
 * \param scheme  The scheme to describe.
 *
 * \return A static string describing the scheme. The calling function
 *         must not try to change the string.
 */
extern const char *
loadbalance_schemestr(loadbalance_scheme_t scheme);

/**
 * \brief Function to create a loadbalancing structure, on which the
 *        different schemes will then work.
 *
 * \param log     The logging module.
 * \param scheme  The scheme to be employed.
 * \param ctype   The type of the SFC.
 * \param level   The loadbalance level.
 * \param ncpu    The number of CPUs to do the loadbalance for.
 *
 * \return A freshly allocated loadbalance structure.
 */
extern loadbalance_t
loadbalance_new(io_logging_t log,
                loadbalance_scheme_t scheme,
                sfc_curve_t ctype,
                int level,
                int ncpu);

/**
 * \brief Function for disposing a loadbalance structure cleanly.
 *
 * A pointer to the loadbalance structure must be given, as the function
 * takes care, that the reference gets set to NULL to prevent usage of
 * freed memory location.
 *
 * \param *loadbal  Pointer to the variable holding the loadbalance
 *                  structure. This will be set to NULL. If a NULL
 *                  pointer is given here, nothing will be done, also
 *                  nothing will be done, if the loadbalance structure
 *                  is already empty (*loadbal == NULL).
 *
 * \return Nothing.
 */
extern void
loadbalance_del(loadbalance_t *loadbal);

/**
 * \brief Updates the keys of the particles and counts them according to
 *        the loadbalancing level. Then a new distribution of the keys
 *        is calculated.
 *
 * \param log       A logging module.
 * \param loadbal   The loadbalance structure to work with.
 * \param part      The first particle.
 * \param no_parts  The total number of particles (length of the array).
 *
 * \return Returns nothing, but will update the loadbalance structure
 *         (the counting array specifically) and the SFC keys of the
 *         particles.
 */
extern void
loadbalance_update(io_logging_t log,
                   loadbalance_t loadbal,
                   partptr part,
                   uint64_t no_parts);

/**
 * \brief Updates the keys of the particles and counts them according to
 *        the loadbalancing level. The distribution of keys is NOT
 *        changed.
 *
 * \param log       A logging module.
 * \param loadbal   The loadbalance structure to work with.
 * \param part      The first particle.
 * \param no_parts  The total number of particles (length of the array).
 *
 * \return Returns nothing, but will update the loadbalance structure
 *         (the counting array specifically) and the SFC keys of the
 *         particles.
 */
extern void
loadbalance_update_plain(io_logging_t log,
                         loadbalance_t loadbal,
                         partptr part,
                         uint64_t no_parts);

/**
 * \brief Puts the current loadbalance information to the logfile.
 *
 * \param log  The logging module.
 * \param lb   The loadbalance structure to log.
 *
 * \return Nothing.
 */
extern void
loadbalance_log(io_logging_t log, loadbalance_t lb);

/**
 * \brief Will free the memory of the counting arrays to save memory.
 *
 * \param log  A logging module.
 * \param lb   The loadbalance structure to minimize.
 *
 * \return Returns nothing.
 */
extern void
loadbalance_minimalMemory(io_logging_t log, loadbalance_t lb);


#endif /* LOADBALANCE_H */
