#ifndef IO_FILE_H
#define IO_FILE_H

/* $Id: io_file.h,v 1.17 2007/12/10 13:13:08 knolli Exp $ */

/**
 * \file io_file.h
 *
 * Provides general definitions and function for reading and writing
 * files.
 */

/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif
#include "io_file_aux.h"
#include "io_logging.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/

/** The mode string for fopen for reading files */
#define IO_FILE_MODE_READ "rb"

/** The mode string for fopen for writing files */
#define IO_FILE_MODE_WRITE "wb"

/** Descriptive string of the AMIGA file type */
#define IO_FILE_AMIGA_STR "AMIGA binary"

/** Descriptive string of the PAMIGA file type */
#define IO_FILE_PAMIGA_STR "PAMIGA binary"

/** Descriptive string of the ARE file type */
#define IO_FILE_ARES_STR "ARES binary"

/** Descriptive string of the MLAPM file type */
#define IO_FILE_MLAPM_STR "MLAPM binary"

/** Descriptive string of the ASCII file type */
#define IO_FILE_ASCII_STR "ASCII"

/** Descriptive string of the CubeP3M file type */
#define IO_FILE_CUBEP3M_STR "CubeP3M binary"

/** Descriptive string of the Multiple CubeP3M file type */
#define IO_FILE_MCUBEP3M_STR "Multiple CubeP3M binary"

/** Descriptive string of the Gadget file type */
#define IO_FILE_GADGET_STR "Gadget binary"

/** Descriptive string of the Multiple Gadget file type */
#define IO_FILE_MGADGET_STR "Multiple Gadget binary"

/** Descriptive string of the DEVA file type */
#define IO_FILE_DEVA_STR "DEVA binary fixed-length records"

/** Descriptive string of the DEVA file type */
#define IO_FILE_DEVANATIVE_STR "DEVA binary native"

/** Descriptive string of the TIPSY file type */
#define IO_FILE_TIPSY_STR "TIPSY binary"

/** Descriptive string of the ART file type */
#define IO_FILE_ART_STR "ART binary"

/** Descriptive string of the Unkown file type */
#define IO_FILE_UNKOWN_STR "Unkown (trying to autoguess)"

/** Descriptive string of the Empty file type */
#define IO_FILE_EMPTY_STR "Empty file, dummy thing for parallel reads"

/** This defines the byte swapping of the file */
typedef enum {
	/* If the file is not swapped */
	IO_FILE_ISNOT_SWAPPED = 0,
	/* If the file is swapped */
	IO_FILE_IS_SWAPPED = 1,
	/* Convenient state if the swapping is unkown */
	IO_FILE_UNKOWN_SWAPPING = 4
} io_file_swap_t;

/** This defines the modes in which to open the files */
typedef enum {
	/* Only for reading */
	IO_FILE_READ = 1,
	/* For writing */
	IO_FILE_WRITE = 2
} io_file_mode_t;

/** This defines the type of the file(s) */
typedef enum {
	/** For AMIGA format */
	IO_FILE_AMIGA = 0,
	/** For PAMIGA format */
	IO_FILE_PAMIGA = 1,
	/** For MLAPM format */
	IO_FILE_MLAPM = 2,
	/** For ARES format */
	IO_FILE_ARES = 5,
	/** For ASCII format */
	IO_FILE_ASCII = 10,
	/** For CUBEP3M format */
	IO_FILE_CUBEP3M = 20,
	/** For MCubeP3M format */
	IO_FILE_MCUBEP3M = 21,
	/** For Gadget format */
	IO_FILE_GADGET = 60,
	/** For Gadget format in multiple files */
	IO_FILE_MGADGET = 61,
	/** For ART format */
	IO_FILE_ART = 70,
	/** For ART format, with new reader. */
	IO_FILE_ARTNEW = 71,
	/** For DEVA format */
	IO_FILE_DEVA = 80,
	/** For DEVA format */
	IO_FILE_DEVANATIVE = 81,
	/** For TIPSY format */
	IO_FILE_TIPSY = 90,
	/** Unkown format, try to guess (not implemented..) */
	IO_FILE_UNKOWN = 100,
	/** For an empty file handler */
	IO_FILE_EMPTY = 9999
} io_file_type_t;

/** Defines different things that can be requested with io_file_get() */
typedef enum {
	/** Requests the number of particles in the simulation, requires a 
	 *  long storage pointer */
	IO_FILE_GET_NOPART,
	/** 
	 * Requests the number of particle in the file (as opposed to 'in
	 * the simulation'. Requires a long storage pointer.
	 */
	IO_FILE_GET_NOPART_IN_FILE,
	/** Requests the number of virtual particles in the file, requires
	 *  a double storage pointer */
	IO_FILE_GET_NOVPART,
	/** Requests the number of mass species in the file, requires 
	 *  a int storage pointer */
	IO_FILE_GET_NOSPECIES,
	/** Requests the boxsize, requires a double storage pointer */
	IO_FILE_GET_BOXSIZE,
	/** Requests the particle mass, requires a double storage pointer */
	IO_FILE_GET_PMASS,
	/** Requests the initial redshift, requires a double storage pointer */
	IO_FILE_GET_ZINITIAL,
	/** Requests the redshift, requires a double storage pointer */
	IO_FILE_GET_Z,
	/** Requests the expansion factor at the beginning, requires a
	 *  double storage pointer */
	IO_FILE_GET_AINITIAL,
	/** Requests the expansion, requires a double storage pointer */
	IO_FILE_GET_A,
	/** Requests Omega0, requires a double storage pointer */
	IO_FILE_GET_OMEGA0,
	/** Requests OmegaLambda, requires a double storage pointer */
	IO_FILE_GET_OMEGAL,
	/** Requests the Hubble parameter, requires a double storage pointer */
	IO_FILE_GET_H,
	/** Asks if the floating point values are stored as doubles,
	 *  requires an int storage pointer */
	IO_FILE_GET_DOUBLE,
	/** Asks if the file is a multimass file, requires an int storage
	 *  pointer */
	IO_FILE_GET_MMASS,
	/** Asks for the number of timesteps, requires an int32_t storage
	 *  pointer */
	IO_FILE_GET_NOTSTEP,
	/** Asks for the timestep, requires a double storage pointer */
	IO_FILE_GET_TSTEP,
	/** Request the header string, requires an char pointer storage
	 *  pointer, not that a reference is returned here, not a copy */
	IO_FILE_GET_HEADERSTR,
	/** Requests the minimal weight in the simulation, requires a double
	 *  pointer */
	IO_FILE_GET_MINWEIGHT,
	/** Requests the maximal weight in the simulation, requires a double
	 *  pointer */
	IO_FILE_GET_MAXWEIGHT
} io_file_get_t;

/** The generic file object */
struct io_file_struct {
	/** This will store the file type */
	io_file_type_t ftype;
#ifdef WITH_MPI
	/** The global rank of the process */
	int rank;
	/** The size of the global communicator */
	int size;
	/** Stores the communicator used for intra libio communication */
	MPI_Comm mycomm;
	/** The size of the intra-library communicator */
	int size_mycomm;
	/** The rank of the local process */
	int rank_mycomm;
#endif
};

/** Convenient typedef */
typedef struct io_file_struct io_file_struct_t;

/** Convenient typedef */
typedef io_file_struct_t *io_file_t;


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Returns a string describing the file type.
 *
 * \param type  The filetype to describe.
 *
 * \return A static string describing the file type. The calling
 *         function must not try to change this string.
 */
extern const char*
io_file_typestr(io_file_type_t type);

/**
 * \brief Tries to open a file of given type.
 *
 * This can be used to read from an existing file (read mode), or to
 * open a new file for writing (write mode). When opening an existing
 * file in write mode, the file will be overwritten.
 *
 * Basically this function is a wrapper around a fopen()-call. It will
 * create the file object but will not do any reading of data, hence
 * header information an such are not set. To do this, a subsequent call
 * to io_file_init() is required.
 *
 * \param log       The logging object.
 * \param *fname    The filename, or for multiple files the stem.
 * \param type      The type of file to open.
 * \param swap      The byte-swapping of the file to be opened. The
 *                  interpretation of this depends on the filetype.
 * \param mode      Tells if a file should be opened for reading or
 *                  writing.
 * \param parallel  Signals if the IO should be done in parallel.
 * \param reader    Number of processes reading. Only important if
 *                  parallel is true, otherwise reader is forced to 1.
 *
 * \return A partially initialized file object will be return. If the
 *         file could not be opened, NULL will be returned.
 */
extern io_file_t
io_file_open(io_logging_t log,
             char *fname,
             io_file_type_t type,
             io_file_swap_t swapped,
             io_file_mode_t mode,
             uint32_t reader);

/**
 * \brief Closes and finalizes a opened file.
 *
 * \param log  The logging object.
 * \param *f   Pointer to the variable holding the file object. This
 *             variable will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_file_close(io_logging_t log,
              io_file_t *f);

/**
 * \brief Initializes an opened for reading file. This will do nothing
 *        on files opened for writing.
 *
 * This is required to setup the file object correctly and it will read
 * in the header information. Nothing will be done, if the provided file
 * object is NULL.
 *
 * \param log  The logging object.
 * \param f    The file object to be initialized.
 *
 * \return Nothing.
 */
extern void
io_file_init(io_logging_t log,
             io_file_t f);

/**
 * \brief Reads dark matter particles from an opened file.
 *
 * This functions requires the file to be opened for reading and being
 * initialized. If theses criteria are not met, the function will exit,
 * returning 0 particles read.
 *
 * The external particle structure can be arranged as wished, the
 * function requires the abstract description encoded in the
 * strg-structure to learn how to advance in the storage and how to
 * store values.
 *
 * Nothing will be done, if the provided file object is NULL.
 *
 * \param log    The logging object.
 * \param f      The initialized file object.
 * \param pskip  The number of particles to skip before starting to
 *               read.
 * \param pread  The number of particles to read. Giving 0 here means
 *               to read all particles.
 * \param strg   The abstract description of the external storage.
 *
 * \return Returns the number of particles read from the file. If this
 *         is not the number of particles given to read, something went
 *         wrong; this can be either an improper initialized file
 *         object, a reading error or simply the file did not contain as
 *         many particles as were requested to read. The calling
 *         function should check the return value.
 */
extern uint64_t
io_file_readpart(io_logging_t log,
                 io_file_t f,
                 uint64_t pskip,
                 uint64_t pread,
                 io_file_strg_struct_t strg);

/**
 * \brief Write particles to a file.
 *
 * Writing particles to a file always updates the header of the file
 * with the number of particles in the file. For the requirements on the
 * external structure holding the particle information, see description
 * of io_file_readpart().
 *
 * Nothing will be done, if the provided file object is NULL.
 *
 * \param log      The logging object.
 * \param f        The initialized file object.
 * \param pskip    The number of particles to skip in the file before
 *                 starting to write. If there are not enough particles
 *                 in the file to skip, the function will return with 0.
 *                 If more particles are already in the file, they will
 *                 be overwritten.
 * \param pwrite   The number of particles to write. 
 * \param strg     The external particle storage.
 *
 * \return Returns the number of particles written to the file. If this
 *         is not the number of particles given to write, something went
 *         wrong.
 */
extern uint64_t
io_file_writepart(io_logging_t log,
                  io_file_t f,
                  uint64_t pskip,
                  uint64_t pwrite,
                  io_file_strg_struct_t strg);

/**
 * \brief Write particles to a file in an ordered way.
 *
 * Writing particles to a file always updates the header of the file
 * with the number of particles in the file. 
 *
 * Nothing will be done, if the provided file object is NULL.
 *
 * \param log        The logging object.
 * \param f          The initialized file object.
 * \param pskip      The number of particles to skip in the file before
 *                   starting to write. If there are not enough
 *                   particles in the file to skip, the function will
 *                   return with 0. If more particles are already in the
 *                   file, they will be overwritten.
 * \param pwrite     The number of particles to write.
 * \param *nxt_part  Pointer to the next particle.
 * \param strg       The external particle storage. Note that the stride
 *                   information in the storage description will not be
 *                   used.
 *
 * \return Returns the number of particles written to the file. If this
 *         is not the number of particles given to write, something went
 *         wrong.
 */
extern uint64_t
io_file_writepart_ord(io_logging_t log,
                      io_file_t f,
                      uint64_t pskip,
                      uint64_t pwrite,
                      void *nxt_part,
                      io_file_strg_struct_t strg);

/**
 * \brief Returns the number of particles which will be read with the
 *        given choice of skipping and reading.
 *
 * This function should be used before using the particle reading
 * routines, as here the skipping and reading are adjusted and the
 * amount of particles which will be returned by the reading routines
 * are returned.
 *
 * This can also be used to get the total number of particles in the
 * file, if pskip is 0 and pread is very large (use UINT64_MAX). After
 * the call, pread will contain the actual number of particles in the
 * file.
 *
 * The return value of this function is pread, iff the file gets read by
 * only one process. Otherwise it will be the according number the
 * calling process will read.
 *
 * If the provided file object is NULL, the function will return 0.
 *
 * \param log    The logging module.
 * \param f      The file object.
 * \param *pskip The number of particles to skip, will be adjusted to
 *               the number of particles in the file if too big.
 * \param *pread The number of particles to read, will be adjusted to 
 *               the number of particles left in file after skipping, if
 *               initially too big.
 *
 * \return The number of particles.
 */
extern uint64_t
io_file_get_numpart(io_logging_t log,
                    io_file_t f,
                    uint64_t *pskip,
                    uint64_t *pread);

/**
 * \brief Generic get-function to retrieve things from the file header.
 *
 * Nothing will be done, if the provided file object is NULL.
 *
 * \param log   The logging module.
 * \param f     The file.
 * \param what  What should be returned.
 * \param *res  A pointer to the place where the result will be stored.
 *
 * \return True if the parameter could be read, false if not.
 */
extern bool
io_file_get(io_logging_t log,
            io_file_t f,
            io_file_get_t what,
            void *res);

/**
 * \brief Generic get-function to set things in the file header.
 *
 * Nothing will be done, if the provided file object is NULL.
 *
 * \param log   The logging module.
 * \param f     The file.
 * \param what  What should be set
 * \param *res  Pointer to the value that will be stored.
 *
 * \return True if the parameter could be set, false if not.
 */
extern bool
io_file_set(io_logging_t log,
            io_file_t f,
            io_file_get_t what,
            void *res);

/**
 * \brief Writes the file information to the logfile.
 *
 * Nothing will be done, if the provided file object is NULL.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Nothing.
 */
extern void
io_file_log(io_logging_t log,
            io_file_t f);


#endif /* IO_FILE_H */
