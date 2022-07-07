#ifndef IO_CUBEP3M_H
#define IO_CUBEP3M_H

/**
 * \file io_cubep3m.h
 *
 * Provides functions for reading and writing CUBEP3M files.
 */


/***********************************************************************
 *    Includes                                                         *
 ***********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#ifdef WITH_MPI
#  include <mpi.h>
#endif

#include "io_cubep3m_header_def.h"
#include "io_cubep3m_def.h"
#include "io_file.h"
#include "io_logging.h"


/***********************************************************************
 *    Global defines, structure definitions and typedefs               *
 ***********************************************************************/


/***********************************************************************
 *    Prototypes of global functions                                   *
 ***********************************************************************/

/**
 * \brief Just tries to open a CUBEP3M file.
 *
 * This can be used either to read from an existing file, or to open a
 * new file for reading. See the mode flag.
 *
 * \param log       The logging object.
 * \param *fname    The filename of the CUBEP3M file.
 * \param swapped   If this is IO_FILE_ISNOT_SWAPPED or
 *                  IO_FILE_IS_SWAPPED an according CUBEP3M file is
 *                  assumed and no checks will be done. Hence
 *                  IO_FILE_UNKOWN_SWAPPING should be supplied when
 *                  opening a CUBEP3M file as the interal mechanism will
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
extern io_cubep3m_t
io_cubep3m_open(io_logging_t   log,
                char           *fname,
                io_file_swap_t swapped,
                io_file_mode_t mode,
                uint32_t       reader);


/**
 * \brief This will close and finalize an CUBEP3M file.
 *
 * \param log  The logging object.
 * \param *f   Pointer to the variable holding the CUBEP3M file object.
 *             This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_cubep3m_close(io_logging_t log,
                 io_cubep3m_t *f);


/**
 * \brief Initializes an opened for reading CUBEP3M file.
 *
 * \param log  The logging object.
 * \param f    The file object to be initialized.
 *
 * \return Nothing.
 */
extern void
io_cubep3m_init(io_logging_t log,
                io_cubep3m_t f);


/**
 * \brief Sets the number of files the cubep3m set consists of.
 *
 * \param log        The logging object.
 * \param f          The file object.
 * \param num_files  The number of files in this cubep3m set.
 *
 * \return  Nothing.
 */
extern void
io_cubep3m_init_num_files(io_logging_t log,
                          io_cubep3m_t f,
                          int          num_files);


/**
 * \brief Reads from an opened CUBEP3M file all particle information and
 *        converts them to AMIGA units.
 *
 * This functions requires the file object to be opened by
 * io_cubep3m_open and inititalized by io_amiga_init. It also requires
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
io_cubep3m_readpart(io_logging_t          log,
                    io_cubep3m_t          f,
                    uint64_t              pskip,
                    uint64_t              pread,
                    io_file_strg_struct_t strg);


/**
 * \brief Reads from an opened CUBEP3M file all particle information
 *        without converting to AMIGA units.
 *
 * Otherwise the io_cubep3m_readpart();
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
io_cubep3m_readpart_raw(io_logging_t          log,
                        io_cubep3m_t          f,
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
io_cubep3m_get(io_logging_t  log,
               io_cubep3m_t  f,
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
io_cubep3m_log(io_logging_t log, io_cubep3m_t f);


/**
 * \brief Does the scaling of particles.
 *
 * \param log             A logging module.
 * \param boxsize         Boxsize (in Mpc/h).
 * \param expansion       Expansion factor.
 * \param lunit           Length Unit (factor to go from file units to
 *                        Mpc/h).
 * \param vunit           Velocity Unit (factor to file units to km/s).
 * \param particles_read  The number of particles to scale.
 * \param strg            The storage structure.
 *
 * \return Returns the number of scaled particles, which should be
 *         exactly particles_read.
 */
extern uint64_t
io_cubep3m_scale_particles(io_logging_t          log,
                           double                boxsize,
                           double                expansion,
                           double                lunit,
                           double                vunit,
                           uint64_t              particles_read,
                           io_file_strg_struct_t strg);

#endif /* IO_CUBEP3M_H */
