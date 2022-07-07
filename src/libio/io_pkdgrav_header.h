#ifndef IO_PKDGRAV_HEADER_H
#define IO_PKDGRAV_HEADER_H

/* $Id: io_pkdgrav_header.h,v 1.2 2006/11/13 15:22:07 knolli Exp $ */

/**
 * \file io_pkdgrav_header.h
 *
 * Provides functions for reading and writing the header of pkdgrav
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include "io_pkdgrav_header_def.h"
#include "io_pkdgrav_def.h"
#include "io_logging.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of global functions                                  * 
\**********************************************************************/

/**
 * \brief Reads a pkdgrav header.
 *
 * The function will first rewind the file pointer to the start
 * of the header and then read in everything. Then it will
 * position the file pointer to the start of the particle date.
 *
 * \param log  A logging object.
 * \param f    A pkdgrav file object.
 *
 * \return A freshly filled header, or NULL, in case of memory
 *         problems.
 */
extern io_pkdgrav_header_t
io_pkdgrav_header_get(io_logging_t log, io_pkdgrav_t f);

/**
 * \brief This will delete a pkdgrav.header object.
 *
 * \param log      The logging object.
 * \param *header  A pointer to the variable holding the header object.
 *                 This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_pkdgrav_header_del(io_logging_t log, io_pkdgrav_header_t *header);

/**
 * \brief Writes a header to the file.
 *
 * \param log     The logging object.
 * \param header  The header to write.
 * \param f       The file the header will be written to.
 *
 * \return Nothing.
 */
extern void
io_pkdgrav_header_write(io_logging_t log,
                       io_pkdgrav_header_t header,
                       io_pkdgrav_t f);

/**
 * \briefs Writes the header information to the logfile.
 *
 * \param log     The logging object.
 * \param header  The header object.
 *
 * \return Nothing.
 */
extern void
io_pkdgrav_header_log(io_logging_t log, io_pkdgrav_header_t header);


#endif /* IO_PKDGRAV_HEADER_H */
