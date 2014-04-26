#ifndef IO_TIPSY_HEADER_H
#define IO_TIPSY_HEADER_H

/**
 * \file io_tipsy_header.h
 *
 * Provides functions for reading and writing the header of TIPSY
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include "io_tipsy_header_def.h"
#include "io_tipsy_def.h"
#include "io_logging.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of global functions                                  * 
\**********************************************************************/

/**
 * \brief Reads a TIPSY header.
 *
 * The function will first rewind the file pointer to the start
 * of the header and then read in everything. Then it will
 * position the file pointer to the start of the particle date.
 *
 * \param log  A logging object.
 * \param f    A TIPSY file object.
 *
 * \return A freshly filled header, or NULL, in case of memory
 *         problems.
 */
extern io_tipsy_header_t
io_tipsy_header_get(io_logging_t log, io_tipsy_t f);

/**
 * \brief This will delete a tipsy.header object.
 *
 * \param log      The logging object.
 * \param *header  A pointer to the variable holding the header object.
 *                 This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_tipsy_header_del(io_logging_t log, io_tipsy_header_t *header);

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
io_tipsy_header_write(io_logging_t log,
                       io_tipsy_header_t header,
                       io_tipsy_t f);

/**
 * \briefs Writes the header information to the logfile.
 *
 * \param log     The logging object.
 * \param header  The header object.
 *
 * \return Nothing.
 */
extern void
io_tipsy_header_log(io_logging_t log, io_tipsy_header_t header);


#endif /* IO_TIPSY_HEADER_H */
