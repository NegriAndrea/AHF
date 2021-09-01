#ifndef IO_AMIGA_HEADER_H
#define IO_AMIGA_HEADER_H

/* $Id: io_amiga_header.h,v 1.5 2006/11/10 14:54:15 knolli Exp $ */

/**
 * \file io_amiga_header.h
 *
 * Provides functions for reading and writing the header of AMIGA
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include "io_amiga_header_def.h"
#include "io_amiga_def.h"
#include "io_logging.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of global functions                                  * 
\**********************************************************************/

/**
 * \brief Reads a AMIGA header.
 *
 * The function will first rewind the file pointer to the start
 * of the header and then read in everything. Then it will
 * position the file pointer to the start of the particle date.
 *
 * \param log  A logging object.
 * \param f    An AMIGA file object.
 *
 * \return A freshly filled header, or NULL, in case of memory
 *         problems.
 */
extern io_amiga_header_t
io_amiga_header_get(io_logging_t log, io_amiga_t f);

/**
 * \brief Generates an empty header object.
 *
 * \param log  A logging object.
 *
 * \return A freshly allocated header, or NULL, in case of memory
 *         problems.
 */
extern io_amiga_header_t
io_amiga_header_new(io_logging_t log);

/**
 * \brief This will delete an amiga_header object.
 *
 * \param log      The logging object.
 * \param *header  A pointer to the variable holding the header object.
 *                 This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_amiga_header_del(io_logging_t log, io_amiga_header_t *header);

/**
 * \brief Writes the header to the file
 *
 * \param log     A logging object.
 * \param header  The header to write.
 * \param f       The file the header will be written to.
 *
 * \return Nothing.
 */
extern void
io_amiga_header_write(io_logging_t log,
                      io_amiga_header_t header,
                      io_amiga_t f);

/**
 * \briefs Writes the header information to the logfile.
 *
 * \param log     The logging object.
 * \param header  The header object.
 *
 * \return Nothing.
 */
extern void
io_amiga_header_log(io_logging_t log, io_amiga_header_t header);


#endif /* IO_AMIGA_HEADER_H */
