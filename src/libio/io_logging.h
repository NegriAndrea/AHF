#ifndef IO_LOGGING_H
#define IO_LOGGING_H

/* $Id: io_logging.h,v 1.10 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file io_logging.h
 *
 * Provides functions for logging.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/

/** 
 * In MPI mode don't write a logfile unless this process is the master
 * thread.
 */
#define IO_LOGGING_FLAGS_NO_OUTPUT_IF_NOT_MASTER 1

/**
 * Duplicate critical messages to the logfile.
 */
#define IO_LOGGING_FLAGS_DUPLICATE_CRITICAL 2

/**
 * Use stdout for critical messages as opposed to stderr normally used.
 */
#define IO_LOGGING_FLAGS_USE_STDOUT_FOR_CRITICAL 4

/**
 * Use this as base name if nothing is given as the stem.
 */
#define IO_LOGGING_STEM_STANDARD "standard-logfile"

/**
 * Use this as a suffix for the log-filename.
 */
#define IO_LOGGING_SUFFIX ".log"

/**
 * The logging object.
 */
struct io_logging_struct {
	/** Holds the logfile stream */
	FILE *logfile;
	/** Holds the filename of the logfile */
	char *fname;
	/** Holds the critical output stream */
	FILE *critical;
	/** Stores the verbosity */
	int32_t verbosity;
	/** Stores the flags */
	int32_t flags;
	/** Stores a counter for parts */
	int32_t cnt_part;
	/** Stores a counter for the sections */
	int32_t cnt_sec;
	/** Stores a counter for the subsections */
	int32_t cnt_subsec;
};

/** Convenient typedef */
typedef struct io_logging_struct io_logging_struct_t;

/** Convenient typedef */
typedef io_logging_struct_t *io_logging_t;


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Starts the logging.
 *
 * \param *stem      The base for the filename of the log-file.
 *                   Normally only .log will get appended but in MPI
 *                   mode, an additional integer will be added to
 *                   encode the process ID.
 * \param verbosity  Sets the basic level of verbosity. Only
 *                   messages with a verbosity value equal or lower
 *                   than this value will get printed.
 * \param flags      A combination of the IO_LOGGING_FLAGS-defines.
 *                   Normally 0 is fine here.
 *
 * \return Returns a logging object that must be given for all further
 *         calls to logging functions. This might be NULL if an error
 *         occured.
 */
extern io_logging_t
io_logging_start(char *stem,
                 int32_t verbosity,
                 int32_t flags);

/**
 * \brief Stops logging and disposes the logging structure.
 *
 * \param *log  A pointer to the external variable storing the logging
 *              module. The external variable will be set to NULL.
 *              If this is NULL, or the external variable points to
 *              NULL, the function does nothing.
 *
 * \return Nothing.
 */
extern void
io_logging_stop(io_logging_t *log);

/**
 * \brief Starts a new numbered part in the logfile.
 *
 * \param log    The log-module.
 * \param *frmt  The format of the output string, like in printf.
 * \param ...    List of arguments.
 *
 * \return Nothing.
 */
extern void
io_logging_part(io_logging_t log,
                const char *frmt,
                ...);

/**
 * \brief Starts a new numbered section in the logfile.
 *
 * \param log    The log-module.
 * \param *frmt  The format of the output string, like in printf.
 * \param ...    List of arguments.
 *
 * \return Nothing.
 */
extern void
io_logging_section(io_logging_t log,
                   const char *frmt,
                   ...);

/**
 * \brief Starts a new numbered subsetion in the logfile.
 *
 * \param log    The log-module.
 * \param *frmt  The format of the output string, like in printf.
 * \param ...    List of arguments.
 *
 * \return Nothing.
 */
extern void
io_logging_subsection(io_logging_t log,
                      const char *frmt,
                      ...);

/**
 * \brief Writes a greeting message to the log.
 *
 * \param log       The log-module to be used.
 * \param *version  The version string.
 *
 * \return Nothing.
 */
extern void
io_logging_hello(io_logging_t log, float version, int build);

/**
 * \brief Writes an identification message of the library to the
 *        logfile. Includes information of the processors.
 *
 * \param log   The log-module to be used.
 *
 * \return Nothing.
 */
extern void
io_logging_identify(io_logging_t log);

/**
 * \brief Writes a message to the logfile.
 *
 * Given that the message verbosity is smaller or equal to the
 * verbosity of the used log-module, the given string will be
 * written to the logfile. Some information will be prepended
 * and a newline added at the end.
 *
 * \param log        The log-module.
 * \param verbosity  The verbosity level.
 * \param *frmt      The format of the output string, like in printf.
 * \param ...        List of argmunets
 *
 * \return Nothing.
 */
extern void
io_logging_msg(io_logging_t log,
               int32_t verbosity,
               const char *frmt,
               ...);

/**
 * \brief Writes a message to the logfile without adding anything.
 *
 * \param log        The log-module.
 * \param verbosity  The verbosity level.
 * \param *frmt      The format of the output string, like in printf.
 * \param ...        List of argmunets
 *
 * \return Nothing.
 */
extern void
io_logging_msgplain(io_logging_t log,
                    int32_t verbosity,
                    const char *frmt,
                    ...);

/**
 * \brief Write a warning to the logfile.
 *
 * This does essentially the same as the msg-function, but adds
 * a 'WARNING' at the beginning.
 *
 * \param log        The log module.
 * \param verbosity  The verbosity level.
 * \param *frmt      The format of the output string, like in printf.
 * \param ...        List of argmunets
 *
 * \return Nothing.
 */
extern void
io_logging_warn(io_logging_t log,
                int32_t verbosity,
                const char *frmt,
                ...);

/**
 * \brief Writes a fatal message to the screen.
 *
 * If the DUPLICATE-flag is set, the message will also be
 * written to the logfile; otherwise the output is only to
 * stderr, or stdout. This message will be printed no matter
 * what verbosity value the log-module has. Also, 'FATAL' will
 * be prepended and newlines appended.
 *
 * \param log    The log module.
 * \param *frmt  The format of the output string, like in printf.
 * \param ...    List of argmunets
 *
 * \return Nothing.
 */
extern void
io_logging_fatal(io_logging_t log,
                 const char *frmt,
                 ...);

/**
 * \brief Writes a fatal message in the case of memory error.
 *
 * Does the same as the fatal-function, but 'FATAL: Could not
 * allocate memory for ' will be appended to the given str.
 * Hence this should only contain a description for what failed.
 * Again, the message will be terminated properly.
 *
 * \param log  The log module.
 * \param *frmt  The format of the output string, like in printf.
 * \param ...    List of argmunets
 *
 * \return Nothing.
 */
extern void
io_logging_memfatal(io_logging_t log,
                    const char *frmt,
                    ...);

/**
 * \brief Will flush the logfile stream.
 *
 * \param log  The log module.
 *
 * \return Nothing.
 */
extern void
io_logging_flush(io_logging_t log);


#endif /* IO_LOGGING_H */
