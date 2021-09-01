#ifndef IO_CUBEP3M_HEADER_H
#define IO_CUBEP3M_HEADER_H

/**
 * \file io_cubep3m_header.h
 *
 * Provides functions for reading and writing the header of CUBEP3M
 * files.
 */


/**********************************************************************\
 *    Includes                                                        *
 \**********************************************************************/
#include "io_cubep3m_header_def.h"
#include "io_cubep3m_def.h"
#include "io_logging.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              *
 \**********************************************************************/


/**********************************************************************\
 *    Prototypes of global functions                                  *
 \**********************************************************************/

/**
 * \brief Reads a CUBEP3M header.
 *
 * The function will first rewind the file pointer to the start
 * of the header and then read in everything. Then it will
 * position the file pointer to the start of the particle date.
 *
 * \param log  A logging object.
 * \param f    A CUBEP3M file object.
 *
 * \return A freshly filled header, or NULL, in case of memory
 *         problems.
 */
extern io_cubep3m_header_t
io_cubep3m_header_get(io_logging_t log, io_cubep3m_t f);


/**
 * \brief This will delete a cubep3m header object.
 *
 * \param log      The logging object.
 * \param *header  A pointer to the variable holding the header object.
 *                 This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_cubep3m_header_del(io_logging_t log, io_cubep3m_header_t *header);


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
io_cubep3m_header_write(io_logging_t        log,
                        io_cubep3m_header_t header,
                        io_cubep3m_t        f);


/**
 * \briefs Writes the header information to the logfile.
 *
 * \param log     The logging object.
 * \param header  The header object.
 *
 * \return Nothing.
 */
extern void
io_cubep3m_header_log(io_logging_t log, io_cubep3m_header_t header);


/**********************************************************************\
 *    Prototypes of protected functions                               *
 \**********************************************************************/

/**
 * \brief  Reads and fills the actual header values in the file for a
 *         CubeP3M file.
 *
 * This function should normally not be used directly, it is called from
 * io_cubep3m_header_get().
 *
 * \param[in,out]  log
 *                    The logging module.
 * \param[in,out]  *f
 *                    The file pointer to use.  The file must be opened for
 *                    binary reading and be positioned at the beginning of
 *                    the CubeP3M header block.
 * \param[in]      swapped
 *                    Toggles between byteswapping and no byteswapping for
 *                    endian adjustment.
 * \param[out]     header
 *                    This is the (already existing) CubeP3M header that
 *                    should be filled with the header values from the file.
 *
 * \return  Returns nothing.
 */
extern void
io_cubep3m_header_read_basics(io_logging_t        log,
                              FILE                *f,
                              io_file_swap_t      swapped,
                              io_cubep3m_header_t header);

extern void
io_cubep3m_header_read_chunk_info(io_logging_t        log,
                                  FILE                *f,
                                  io_file_swap_t      swapped,
                                  io_cubep3m_header_t header);


/**
 * \brief  Reads and fills the extra information for the CubeP3M header from
 *         an info file.
 *
 * This function should not be used directly, it is called from
 * io_cubep3m_header_get().
 *
 * \param[in,out]  log
 *                    The logging module that should be used.
 * \param[in]      extras_file_name
 *                    The name of the file from which to read the extra
 *                    information.
 * \param[out]     header
 *                    The (already) existing header that should be filled
 *                    with the extra values.
 *
 * \return  Returns nothing.
 */
extern void
io_cubep3m_header_read_extras(io_logging_t        log,
                              const char          *extras_file_name,
                              io_cubep3m_header_t header);

#endif /* IO_CUBEP3M_HEADER_H */
