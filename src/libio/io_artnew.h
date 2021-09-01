// Copyright (C) 2011, Steffen Knollmann
// Released under the terms of the GNU General Public License version 3.
// This file is part of `AHF'.

#ifndef IO_ARTNEW_H
#define IO_ARTNEW_H


/*--- Doxygen file description ------------------------------------------*/

/**
 * @file  io_artnew.h
 *
 * Provides functionality to read from ART files using as the backend
 * the ART file object.
 */

/*--- Includes ----------------------------------------------------------*/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#ifdef WITH_MPI
#  include <mpi.h>
#endif
#include "io_file.h"
#include "io_logging.h"


/*--- ADT handle --------------------------------------------------------*/

/**
 * @brief  Defines the handle for an ART object.
 */
typedef struct io_artnew_struct *io_artnew_t;


/*--- Prototypes of exported functions ----------------------------------*/

/**
 * @brief  Creates a new ART file object for use with AHF.
 *
 * This only works for reading ART files (in single or multiple files).
 * No file operation on the actual data file(s) or the header file is
 * done.  This function only creates the according io_file object and
 * sets all required values for use with io_artnew_init().
 *
 * @param[in,out]  log
 *                    The logging object to use
 * @param[in]      *fname
 *                    The filename of the extra file that provides the
 *                    construction information for the ART reader.
 * @param[in]      swapped
 *                    This is ignored, when reading the endianess of the
 *                    file will always be autodetected and every read from
 *                    the file will convert to the system endianess.
 * @param[in]      mode
 *                    This is ignored, the module only works for reading
 *                    at the moment.
 * @param[in]      reader
 *                    Number of processes reading. Only important if in
 *                    MPI mode, otherwise it will be forced to 1.
 *
 * @return Returns a partially initialized file object, or NULL if the
 *         file could not be opened.
 */
extern io_artnew_t
io_artnew_open(io_logging_t   log,
               char           *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t       reader);

/**
 * @brief  Closes the io_file object and frees all associated memory.
 *
 * @param[in,out]  log
 *                    The logging object to use.
 * @param[in,out]  *f
 *                    Pointer to the external variable holding the ART
 *                    file object that should be closed.  After deletion
 *                    and freeing of all memory, the external variable
 *                    will be set to @c NULL.
 *
 * @return  Returns nothing.
 */
extern void
io_artnew_close(io_logging_t log,
                io_artnew_t  *f);

/**
 *
 */
extern void
io_artnew_init(io_logging_t log,
               io_artnew_t  f);

extern uint64_t
io_artnew_readpart(io_logging_t          log,
                   io_artnew_t           f,
                   uint64_t              pskip,
                   uint64_t              pread,
                   io_file_strg_struct_t strg);

extern bool
io_artnew_get(io_logging_t  log,
              io_artnew_t   f,
              io_file_get_t what,
              void          *res);

extern void
io_artnew_log(io_logging_t log, io_artnew_t f);


#endif
