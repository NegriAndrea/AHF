#ifndef IO_DEVA_H
#define IO_DEVA_H

/**
 * \file io_deva.h
 *
 * Provides functions for reading and writing DEVA files.
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

#include "io_deva_header_def.h"
#include "io_deva_def.h"
#include "io_file.h"
#include "io_logging.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Just tries to open a DEVA file.
 *
 * This can be used either to read from an existing file, or to open a
 * new file for reading. See the mode flag.
 *
 * \param log       The logging object.
 * \param *fname    The filename of the DEVA file.
 * \param swapped   If this is IO_FILE_ISNOT_SWAPPED or
 *                  IO_FILE_IS_SWAPPED an according DEVA file is
 *                  assumed and no checks will be done. Hence
 *                  IO_FILE_UNKOWN_SWAPPING should be supplied when
 *                  opening a DEVA file as the interal mechanism will
 *                  normally detect the right state.
 * \param mode      Tells if the file should be opened for reading or for
 *                  writing. If opened for writing, the value for swapped
 *                  will be ignored.
 * \param reader    Number of processes reading. Only important if in
 *                  MPI mode, otherwise it will be forced to 1.
 * \param native    Native or converted format.
 *
 * \return Returns a partially initialized file object, or NULL if the
 *         file could not be opened.
 */
extern io_deva_t
io_deva_open(io_logging_t log,
               char *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t reader,
               int native);

/**
 * \brief This will close and finalize an DEVA file.
 *
 * \param log  The logging object.
 * \param *f   Pointer to the variable holding the DEVA file object.
 *             This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_deva_close(io_logging_t log,
                io_deva_t *f);

/**
 * \brief Initializes an opened for reading DEVA file.
 *
 * \param log     The logging object.
 * \param f       The file object to be initialized.
 * \param native  Native or converted format
 *
 * \return Nothing.
 */
extern void
io_deva_init(io_logging_t log,
               io_deva_t f,
               int native);

/**                          
 * \brief Reads from an opened DEVA file all particle information and
 *        converts them to AMIGA units.
 *
 * This functions requires the file object to be opened by
 * io_deva_open and inititalized by io_amiga_init. It also requires
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
 * \param native Native or converted format
 *                                 
 * \return Returns the number of particles read from the file. If this
 *         is not the number of particles given as the pread parameter,
 *         something went wrong. The calling function hence should check
 *         the return value.
 */
extern uint64_t
io_deva_readpart(io_logging_t log,
                   io_deva_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg,
                   int native);

/**                          
 * \brief Reads from an opened DEVA file all particle information
 *        without converting to AMIGA units.
 *
 * Otherwise the io_deva_readpart();
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
io_deva_readpart_raw(io_logging_t log,
                       io_deva_t f,
                       uint64_t pskip,
                       uint64_t pread,
                       io_file_strg_struct_t strg);

/**                          
 * \brief Reads from an opened DEVA file all particle information
 *        without converting to AMIGA units.
 *
 * Otherwise the io_deva_readpart();
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
io_deva_readpart_raw_native(io_logging_t log,
                       io_deva_t f,
                       uint64_t pskip,
                       uint64_t pread,
                       io_file_strg_struct_t strg);

/**
 * \brief Writes the particles to a DEVA binary file
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
io_deva_writepart(io_logging_t log,
                    io_deva_t f,
                    uint64_t pskip,
                    uint64_t pwrite,
                    io_file_strg_struct_t strg);

/**
 * \brief Writes the particles to a DEVA binary file in an ordered
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
io_deva_writepart_ord(io_logging_t log,
                        io_deva_t f,
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
io_deva_get(io_logging_t log,
              io_deva_t f,
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
io_deva_set(io_logging_t log,
              io_deva_t f,
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
io_deva_log(io_logging_t log, io_deva_t f);

/**
 * \brief Resets the position and weight scales to given values
 *
 * \param log          The logging object.
 * \param f            The file object.
 * \param posscale     The new scale to translate from DEVA file
 *                     units to Mpc.
 * \param weightscale  The new scale to translate from DEVA file units
 *                     to Msun.
 */
extern void
io_deva_resetscale(io_logging_t log,
                     io_deva_t f,
                     double posscale,
                     double weightscale);

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
io_deva_scale_particles(io_logging_t log,
                          double maxpos[],
                          double minpos[],
                          float *boxsize,
                          float expansion,
                          float h0t0,
                          double posscale,
                          double mmass,
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
io_deva_scale_global(io_logging_t log,
                     MPI_Comm comm,
                     double *maxpos,
                     double *minpos,
                     double *mmass,
                     double *minweight,
                     double *maxweight,
                     double *sumweight);
#endif


#endif /* IO_DEVA_H */
