#ifndef IO_DEVA_HEADER_H
#define IO_DEVA_HEADER_H

/**
 * \file io_deva_header.h
 *
 * Provides functions for reading and writing the header of DEVA
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include "io_deva_header_def.h"
#include "io_deva_def.h"
#include "io_logging.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/** Holds the size of integers in the file */
extern size_t deva_intsize;
/** Holds the size of reals in the file */
extern size_t deva_realsize;

/**********************************************************************\
 *    Prototypes of global functions                                  * 
\**********************************************************************/

/**
 * \brief Reads a DEVA header.
 *
 * The function will first rewind the file pointer to the start
 * of the header and then read in everything. Then it will
 * position the file pointer to the start of the particle date.
 *
 * \param log  A logging object.
 * \param f    A DEVA file object.
 *
 * \return A freshly filled header, or NULL, in case of memory
 *         problems.
 */
extern io_deva_header_t
io_deva_header_get(io_logging_t log, io_deva_t f);

/**
 * \brief Reads a DEVA header.
 *
 * The function will first rewind the file pointer to the start
 * of the header and then read in everything. Then it will
 * position the file pointer to the start of the particle date.
 *
 * \param log  A logging object.
 * \param f    A DEVA file object.
 *
 * \return A freshly filled header, or NULL, in case of memory
 *         problems.
 */
extern io_deva_header_t
io_deva_header_get_native(io_logging_t log, io_deva_t f);

/**
 * \brief This will delete a deva.header object.
 *
 * \param log      The logging object.
 * \param *header  A pointer to the variable holding the header object.
 *                 This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_deva_header_del(io_logging_t log, io_deva_header_t *header);

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
io_deva_header_write(io_logging_t log,
                       io_deva_header_t header,
                       io_deva_t f);

/**
 * \briefs Writes the header information to the logfile.
 *
 * \param log     The logging object.
 * \param header  The header object.
 *
 * \return Nothing.
 */
extern void
io_deva_header_log(io_logging_t log, io_deva_header_t header);


/**
 * \brief Subroutine that reads an (optionally) segmented block of data written by a Fortran program
 *
 * Please note that the file must be opened in read mode and enough memory be allocated for buf.
 * This implies you have to know beforehand the size of your data (the same happens in Fortran).
 *
 * \param f	The file object.
 * \param buf0	Pointer to where the data will be written.
 * \param ts	Size of the data type to be read (this is needed for byte-swapping).
 * \param len	Number of elements that will be read.
 * \param skip	Number of elements to skip since the beginning.
 * \param swap	Indicates if swapping should be performed.
 *
 * \return Number of elements actually read.
 *
 * See http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/bldaps_for/common/bldaps_rectypes.htm for a descriotion of Variable-Length Records
 *
 * Francisco Martinez-Serrano, 2008, 2010
 */

extern int
io_util_readfortran(FILE * const f, void * const buf0, size_t const ts, size_t const len, size_t const skip, int const swap);

#endif /* IO_DEVA_HEADER_H */
