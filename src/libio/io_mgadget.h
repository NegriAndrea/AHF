#ifndef IO_MGADGET_H
#define IO_MGADGET_H

/* $Id: io_mgadget.h,v 1.9 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file io_mgadget.h
 *
 * Provides functions for reading and writing multiple Gadget files.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "io_mgadget_def.h"
#include "io_file.h"
#include "io_logging.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Just tries to open multiple Gadget file.
 *
 * This can be used either to read from an existing file, or to open a
 * new file for reading. See the mode flag.
 *
 * \param log       The logging object.
 * \param *fname    The filename stem of the multiple Gadget file. This
 *                  string will be divided into PATH and STEM and then a
 *                  search in PATH is performed to find everything
 *                  matching STEM.
 * \param swapped   If this is IO_FILE_ISNOT_SWAPPED or
 *                  IO_FILE_IS_SWAPPED an according Gadget file is
 *                  assumed and no checks will be done. Hence
 *                  IO_FILE_UNKOWN_SWAPPING should be supplied when
 *                  opening a Gadget file as the interal mechanism will
 *                  normally detect the right state.
 * \param mode      Tells if the file should be opened for reading or for
 *                  writing. If opened for writing, the value for swapped
 *                  will be ignored.
 * \param reader    Number of processes reading. Only important if in
 *                  MPI mode, otherwise it will be forced to 1.
 *
 * \return Returns a partially initialized file object, or NULL if the
 *         file could not be opened.
 */
extern io_mgadget_t
io_mgadget_open(io_logging_t log,
                char *fname,
                io_file_swap_t swapped,
                io_file_mode_t mode,
                uint32_t reader);

/**
 * \brief This will close and finalize multiple Gadget files.
 *
 * \param log  The logging object.
 * \param *f   Pointer to the variable holding the Multiple Gadget file
 *             object.
 *             This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_mgadget_close(io_logging_t log,
                 io_mgadget_t *f);

/**
 * \brief Initializes opened for reading Gadget files.
 *
 * \param log  The logging object.
 * \param f    The file object to be initialized.
 *
 * \return Nothing.
 */
extern void
io_mgadget_init(io_logging_t log,
                io_mgadget_t f);

/**                          
 * \brief Reads from an opened Multiple Gadget file all particle
 *        information.
 *
 * This functions requires the file object to be opened by
 * io_gadget_open and inititalized by io_amiga_init. It also requires
 * pointer to beginning of the particle array, which must be large
 * enough to accomodate all particles (can be check by evaluating the
 * number of particles given in the file header).
 * 
 * The particle structure can be arbitrarily arranged, the function only
 * needs to know where within in the structure the components of
 * position and momentum and the weight are stored and how large the
 * particle structure is.  This is described in the strg parameter, see
 * io_file_aux.h for the definition.
 *
 * \param log    The logging object.
 * \param f      The initialized file object.
 * \param pskip  Number of particles to skip.
 * \param pread  Number of particles to read.
 * \param strg   The abstract description of the external storage.
 *                                 
 * \return Returns the number of particles read from the file. If this
 *         is not the number of particles given as the pread parameter,
 *         something went wrong. The calling function hence should check
 *         the return value.
 */
extern uint64_t
io_mgadget_readpart(io_logging_t log,
                    io_mgadget_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg);

/**                          
 * \brief Reads from an opened Multiple Gadget file all particle
 *        information without scaling them to AMIGA units.
 *
 * See io_mgadget_readpart() for more details.
 *
 * \param log    The logging object.
 * \param f      The initialized file object.
 * \param pskip  Number of particles to skip.
 * \param pread  Number of particles to read.
 * \param strg   The abstract description of the external storage.
 *                                 
 * \return Returns the number of particles read from the file. If this
 *         is not the number of particles given as the pread parameter,
 *         something went wrong. The calling function hence should check
 *         the return value.
 */
extern uint64_t
io_mgadget_readpart_raw(io_logging_t log,
                        io_mgadget_t f,
                        uint64_t pskip,
                        uint64_t pread,
                        io_file_strg_struct_t strg);

/**
 * \brief Generic get-function to retrieve things from the file header.
 *
 * \param log   The logging module.
 * \param f     The file.
 * \param what  What should be returned.
 * \param *res  A pointer to the place where the result will be stored.
 *
 * \return True if the parameter could be read, false if not.
 */
extern bool
io_mgadget_get(io_logging_t log,
               io_mgadget_t f,
               io_file_get_t what,
               void *res);

/**
 * \brief Generic get-function to set things in the file header.
 *
 * \param log   The logging module.
 * \param f     The file.
 * \param what  What should be set
 * \param *res  Pointer to the value that will be stored.
 *
 * \return True if the parameter could be set, false if not.
 */
extern bool
io_mgadget_set(io_logging_t log,
               io_mgadget_t f,
               io_file_get_t what,
               void *res);

/**
 * \brief Writes the file information to the logfile.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Nothing.
 */
extern void
io_mgadget_log(io_logging_t log, io_mgadget_t f);

/**
 * \brief Resets the position and weight scales to given values
 *
 * \param log          The logging object.
 * \param f            The file object.
 * \param posscale     The new scale to translate from Gadget file
 *                     units to Mpc.
 * \param weightscale  The new scale to translate from Gadget file
 *                     to Msun.
 */
extern void
io_mgadget_resetscale(io_logging_t log,
                      io_mgadget_t f,
                      double posscale,
                      double weightscale);


#endif /* IO_MGADGET_H */
