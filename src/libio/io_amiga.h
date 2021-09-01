#ifndef IO_AMIGA_H
#define IO_AMIGA_H

/**
 * \file io_amiga.h
 *
 * Provides functions for reading and writing AMIGA files.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "io_amiga_header_def.h"
#include "io_amiga_def.h"
#include "io_file.h"
#include "io_file_aux.h"
#include "io_logging.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Just tries to open an AMIGA file.
 *
 * This can be used either to read from an existing file, or to open a
 * new file for reading. See the mode flag.
 *
 * \param log       The logging object.
 * \param *fname    The filename of the AMIGA file.
 * \param swapped   If this is IO_AMIGA_ISNOT_SWAPPED or
 *                  IO_AMIGA_IS_SWAPPED an according AMIGA file is
 *                  assumed and no checks will be done. Hence
 *                  IO_AMIGA_UNKOWN_SWAPPING should be supplied when
 *                  opening an AIMGA file as the interal mechanism will
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
extern io_amiga_t
io_amiga_open(io_logging_t log,
              char *fname,
              io_file_swap_t swapped,
              io_file_mode_t mode,
              uint32_t reader);

/**
 * \brief This will close and finalize an AMIGA file.
 *
 * \param log  The logging object.
 * \param *f   Pointer to the variable holding the AMIGA file object.
 *             This variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_amiga_close(io_logging_t log,
               io_amiga_t *f);

/**
 * \brief Initializes an opened for reading amiga file.
 *
 * \param log  The logging object.
 * \param f    The file object to be initialized.
 *
 * \return Nothing.
 */
extern void
io_amiga_init(io_logging_t log,
              io_amiga_t f);

/**
 * \brief Reads from an opend amiga file all particle information.
 *
 * This functions requires the file object to be opened by io_amiga_open
 * and inititalized by io_amiga_init. It also requires an abstract
 * description of the external particle storage, cf. io_file_aux.h for
 * this.
 *
 * If the setting of the weight is requested on single mass files, a
 * warning will be written to the log file and the weight will be set to
 * 1.
 *
 * \param log     The logging object.
 * \param f       The initialized file object.
 * \param pskip   Number of particles to skip.
 * \param pread   Number of particles to read.
 * \param strg    Abstract description of the particle storage.
 *
 * \return Returns the number of particles read from the file. If this
 *         is not the number of particles given in the header, something
 *         went wrong. The calling function hence should check the
 *         return value.
 */
extern uint64_t
io_amiga_readpart(io_logging_t log,
                  io_amiga_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg);

/**
 * \brief Writes the particles to an AMIGA binary file
 *
 * The file object given to the function needs to be opened for
 * writing.
 *
 * \param log      The logging object.
 * \param f        The initialized file object.
 * \param pskip    Number of particles in the file to skip.
 * \param pwrite   Number of particles to write.
 * \param strg     Abstract description of the particle storage.
 *
 * \return Returns the number of particles written to the file. This
 *         should correspond to the number of particles given in the
 *         header.
 */
extern uint64_t
io_amiga_writepart(io_logging_t log,
                   io_amiga_t f,
                   uint64_t pskip,
                   uint64_t pwrite,
                   io_file_strg_struct_t strg);

/**
 * \brief Writes the particles to an AMIGA binary file in an ordered
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
 * \param strg      Abstract description of the particle storage.
 *
 * \return Returns the number of particles written to the file. This
 *         should correspond to the number of particles given in the
 *         header.
 */
extern uint64_t
io_amiga_writepart_ord(io_logging_t log,
                       io_amiga_t f,
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
io_amiga_get(io_logging_t log,
             io_amiga_t f,
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
io_amiga_set(io_logging_t log,
             io_amiga_t f,
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
io_amiga_log(io_logging_t log, io_amiga_t f);


#endif /* IO_AMIGA_H */
