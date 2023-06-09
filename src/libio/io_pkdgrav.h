#ifndef IO_PKDGRAV_H
#define IO_PKDGRAV_H

/**
 * \file io_pkdgrav.h
 *
 * Provides functions for reading and writing pkdgrav files.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif

#include "io_pkdgrav_header_def.h"
#include "io_pkdgrav_def.h"
#include "io_file.h"
#include "io_logging.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Just tries to open a pkdgrav file.
 *
 * This can be used either to read from an existing file, or to open a
 * new file for reading. See the mode flag.
 *
 * \param log       The logging object.
 * \param *fname    The filename of the pkdgrav file.
 * \param swapped   If this is IO_FILE_ISNOT_SWAPPED or
 *                  IO_FILE_IS_SWAPPED an according pkdgrav file is
 *                  assumed and no checks will be done. Hence
 *                  IO_FILE_UNKOWN_SWAPPING should be supplied when
 *                  opening a pkdgrav file as the interal mechanism will
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
extern io_pkdgrav_t
io_pkdgrav_open(io_logging_t log,
               char *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t reader);

/**
 * \brief This will close and finalize an pkdgrav file.
 *
 * \param log  The logging object.
 * \param *f   Pointer to the variable holding the pkdgrav file object.
 *             This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_pkdgrav_close(io_logging_t log,
                io_pkdgrav_t *f);

/**
 * \brief Initializes an opened for reading pkdgrav file.
 *
 * \param log  The logging object.
 * \param f    The file object to be initialized.
 *
 * \return Nothing.
 */
extern void
io_pkdgrav_init(io_logging_t log,
               io_pkdgrav_t f);

/**                          
 * \brief Reads from an opened pkdgrav file all particle information and
 *        converts them to AMIGA units.
 *
 * This functions requires the file object to be opened by
 * io_pkdgrav_open and inititalized by io_amiga_init. It also requires
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
io_pkdgrav_readpart(io_logging_t log,
                   io_pkdgrav_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg);

/**                          
 * \brief Reads from an opened pkdgrav file all particle information
 *        without converting to AMIGA units.
 *
 * Otherwise the io_pkdgrav_readpart();
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
io_pkdgrav_readpart_raw(io_logging_t log,
                       io_pkdgrav_t f,
                       uint64_t pskip,
                       uint64_t pread,
                       io_file_strg_struct_t strg);

/**
 * \brief Writes the particles to a pkdgrav binary file
 *
 * The file object given to the function needs to be opened for
 * writing.
 *
 * \param log     The logging object.
 * \param f       The initialized file object.
 * \param pskip   Number of particles in the file to skip.
 * \param pwrite  Number of particles to write.
 * \param strg    The particle storage.
 *
 * \return Returns the number of particles written to the file. This
 *         should correspond to the number of particles given in the
 *         header.
 */
extern uint64_t
io_pkdgrav_writepart(io_logging_t log,
                    io_pkdgrav_t f,
                    uint64_t pskip,
                    uint64_t pwrite,
                    io_file_strg_struct_t strg);

/**
 * \brief Writes the particles to a pkdgrav binary file in an ordered
 *        way.
 *
 * The file object given to the function needs to be opened for
 * writing.
 *
 * \param log        The logging object.
 * \param f          The initialized file object.
 * \param pskip      Number of particles in the file to skip.
 * \param pwrite     Number of particles to write.
 * \param *nxt_part  Pointer to the storage that holds the the pointer
 *                   to the next particle.
 * \param strg       The particle storage.
 *
 * \return Returns the number of particles written to the file. This
 *         should correspond to the number of particles given in the
 *         header.
 */

extern uint64_t
io_pkdgrav_writepart_ord(io_logging_t log,
                        io_pkdgrav_t f,
                        uint64_t pskip,
                        uint64_t pwrite,
                        void *nxt_part,
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
io_pkdgrav_get(io_logging_t log,
              io_pkdgrav_t f,
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
io_pkdgrav_set(io_logging_t log,
              io_pkdgrav_t f,
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
io_pkdgrav_log(io_logging_t log, io_pkdgrav_t f);


/**
 * \brief Does the scaling of particles.
 *
 * \param log
 * \param maxpos[]
 * \param minpos[]
 * \param *boxsize
 * \param expansion
 * \param posscale
 * \param mmass
 * \param particles_read
 * \param strg
 *
 * \return Returns the number of scaled particles, which should be
 *         exactly particles_read.
 */
extern uint64_t
io_pkdgrav_scale_particles(io_logging_t log,
                         double maxpos[],
                         double minpos[],
                         double boxsize,
                         double expansion,
                         double KpcUnit,
                         double MsolUnit,
                         double ErgPerGmUnit,
                         double KmPerSecUnit,
                         uint64_t particles_read,
                         io_file_strg_struct_t strg);

#ifdef WITH_MPI
/**
 * \brief Establishes the global min and max values needed for scaling.
 *        Only available when in MPI mode.
 *
 * \param log      A logging module.
 * \param comm     The communicator.
 * \param *maxpos  Array of maximal positions, will be updated to the
 *                 global values.
 * \param *minpos  Array of minimal positions, will be updated to the
 *                 global values.
 * \param *mmass   Minimal mass, will be updated to the global values.
 *
 * \return Nothing.
 */
extern void
io_pkdgrav_scale_global(io_logging_t log,
                       MPI_Comm comm,
                       double *maxpos,
                       double *minpos,
                       double *mmass);
#endif


#endif /* IO_PKDGRAV_H */
