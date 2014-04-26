#ifndef IO_MCUBEP3M_H
#define IO_MCUBEP3M_H

/**
 * \file io_mcubep3m.h
 *
 * Provides functions for reading and writing multiple CubeP3M files.
 */


/***********************************************************************
 *    Includes                                                         *
 ***********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "io_mcubep3m_def.h"
#include "io_file.h"
#include "io_logging.h"


/***********************************************************************
 *    Global defines, structure definitions and typedefs               *
 ***********************************************************************/


/***********************************************************************
 *    Prototypes of global functions                                   *
 ***********************************************************************/

/**
 * \brief Just tries to open multiple CubeP3M file.
 *
 * This can be used either to read from an existing file, or to open a
 * new file for reading. See the mode flag.
 *
 * \param log       The logging object.
 * \param *fname    The filename stem of the multiple CubeP3M file. This
 *                  string will be divided into PATH and STEM and then a
 *                  search in PATH is performed to find everything
 *                  matching STEM.
 * \param swapped   If this is IO_FILE_ISNOT_SWAPPED or
 *                  IO_FILE_IS_SWAPPED an according CubeP3m file is
 *                  assumed and no checks will be done. Hence
 *                  IO_FILE_UNKOWN_SWAPPING should be supplied when
 *                  opening a CubeP3M file as the interal mechanism will
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
extern io_mcubep3m_t
io_mcubep3m_open(io_logging_t   log,
                 char           *fname,
                 io_file_swap_t swapped,
                 io_file_mode_t mode,
                 uint32_t       reader);


/**
 * \brief This will close and finalize multiple CubeP3M files.
 *
 * \param log  The logging object.
 * \param *f   Pointer to the variable holding the Multiple CubeP3M file
 *             object.
 *             This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_mcubep3m_close(io_logging_t  log,
                  io_mcubep3m_t *f);


/**
 * \brief Initializes opened for reading CubeP3M files.
 *
 * \param log  The logging object.
 * \param f    The file object to be initialized.
 *
 * \return Nothing.
 */
extern void
io_mcubep3m_init(io_logging_t  log,
                 io_mcubep3m_t f);


/**
 * \brief Reads from an opened Multiple CubeP3M file all particle
 *        information.
 *
 * This functions requires the file object to be opened by
 * io_mcubep3m_open and inititalized by io_mcubep3m_init. It also requires
 * a pointer to beginning of the particle array, which must be large
 * enough to accomodate all particles (can be checked by evaluating the
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
io_mcubep3m_readpart(io_logging_t          log,
                     io_mcubep3m_t         f,
                     uint64_t              pskip,
                     uint64_t              pread,
                     io_file_strg_struct_t strg);


/**
 * \brief Reads from an opened Multiple CubeP3M file all particle
 *        information without scaling them to AMIGA units.
 *
 * See io_mcubep3m_readpart() for more details.
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
io_mcubep3m_readpart_raw(io_logging_t          log,
                         io_mcubep3m_t         f,
                         uint64_t              pskip,
                         uint64_t              pread,
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
io_mcubep3m_get(io_logging_t  log,
                io_mcubep3m_t f,
                io_file_get_t what,
                void          *res);


/**
 * \brief Writes the file information to the logfile.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Nothing.
 */
extern void
io_mcubep3m_log(io_logging_t log, io_mcubep3m_t f);


#endif /* IO_MCUBEP3M_H */
