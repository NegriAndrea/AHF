/* $Id: io_logging.c,v 1.13 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file io_logging.c
 *
 * Provides functions for logging.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif

#include "io_logging.h"
#include "io_logging_defs.h"
#include "io_defs.h"

#define AHF2_overwrite_logfiles

/***********************************************************************\
 *    Local defines, structure definitions and typedefs                * 
\***********************************************************************/


/***********************************************************************\
 *    Prototypes of local functions                                    * 
\***********************************************************************/
/**
 * \brief Counts the digits of a integer. Only needed in the MPI
 *        version.
 *
 * \param i  The integer to be counted.
 *
 * \return The number of used digits. Will be at least 1.
 */
#ifdef WITH_MPI
static int32_t
local_countdigits(int64_t i);
#endif

/**
 * \brief Checks wether the given verbosity is to be printed by the
 *        given logging module.
 *
 * \param log        The logging module.
 * \param verbosity  The verbosity to check.
 *
 * \return Returns -1 if the message should not be printed, 0 otherwise.
 */
inline static int8_t
local_logcheck(io_logging_t log, int32_t verbosity);


/***********************************************************************\
 *    Implementation of global functions                               * 
\***********************************************************************/
extern io_logging_t
io_logging_start(char *stem,
                 int32_t verbosity,
                 int32_t flags)
{
	io_logging_t log;
	char *fname, *suffix;
	int8_t appending = 0;
#	ifdef WITH_MPI
	int32_t size, rank, digits;
#	endif

	/* Get memory for logging structure */
	log = (io_logging_t)malloc(sizeof(io_logging_struct_t));
	if (log == NULL) {
		fprintf(stderr,
		        "%s, %u: Could not allocate memory for logging "
		        "structure.\n",
		        __func__, (unsigned int)__LINE__);
		return NULL;
	}

	/* Initialise the counter settings */
	log->cnt_part = (uint32_t)('A' - 1);
	log->cnt_sec = INT32_C(0);
	log->cnt_subsec = INT32_C(0);

	/* First have a look if we have to go through all the hassle */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (    ( (flags & IO_LOGGING_FLAGS_NO_OUTPUT_IF_NOT_MASTER) != 0 )
	     && (rank > 0) ) {
		/* Cool. We are not supposed to write anything at all. */
		log->fname = NULL;
		log->logfile = NULL;
		log->verbosity = verbosity;
		log->flags = flags;
		if ((flags & IO_LOGGING_FLAGS_USE_STDOUT_FOR_CRITICAL) != 0) {
			log->critical = stdout;
		} else {
			log->critical = stderr;
		}

		/* That's it! */
		return log;
	}
#	endif

	/* Check parameter and make a local copy of the filename */
	if (stem == NULL) {
		fprintf(stderr,
		        "%s, %u: No base filename given. Using standard '%s'.\n",
		        __func__, (unsigned int)__LINE__,
		        IO_LOGGING_STEM_STANDARD);
		fname = (char *)malloc(sizeof(char) * (1 +
		                              strlen(IO_LOGGING_STEM_STANDARD)));
	} else {
		fname = (char *)malloc(sizeof(char) * (1 + strlen(stem)));
	}
	if (fname == NULL) {
		fprintf(stderr,
		        "%s, %u: Could not allocate memory for the stem.\n",
		        __func__, (unsigned int)__LINE__);
		free(log);
		return NULL;
	}
	if (stem == NULL) {
		strcpy(fname, IO_LOGGING_STEM_STANDARD);
	} else {
		strcpy(fname, stem);
	}

	/* Get memory for the suffix */
#	ifdef WITH_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* How many digits are needed to distinguish all processes? Make
	 * it one more */
	digits = local_countdigits((int64_t)size) + 1;

	/* Get memory for the suffix string, +2 because of the *\
	 * terminating NULL character and a dot separating the *
	\* process ID and the suffix                           */
	suffix = (char *)malloc(sizeof(char) * (  digits + 2
	                                        + strlen(IO_LOGGING_SUFFIX)));
#	else
	suffix = (char *)malloc(sizeof(char) * (  1 
	                                        + strlen(IO_LOGGING_SUFFIX)));
#	endif /* WITH_MPI */
	if (suffix == NULL) {
		fprintf(stderr,
	    	    "%s, %u: Could not allocate memory for suffix.\n",
		        __func__, (unsigned int)__LINE__);
		free(fname);
		free(log);
		return NULL;
	}

	/* Create the suffix */
#	ifdef WITH_MPI
	sprintf(suffix,
	        ".%0*" PRIi32 "%s",
	        (int)digits, rank, IO_LOGGING_SUFFIX);
#	else
	sprintf(suffix, "%s", IO_LOGGING_SUFFIX);
#	endif /* WITH_MPI */

	/* Concatenate the filename and the suffix, first create a   *\
	\* sufficiently large string                                 */
	log->fname = (char *)malloc(sizeof(char)*(  strlen(fname) + 1
	                                          + strlen(suffix)));
	if (log->fname == NULL) {
		fprintf(stderr,
	    	    "%s, %u: Could not allocate memory for the filename.\n",
		        __func__, (unsigned int)__LINE__);
		free(fname);
		free(suffix);
		free(log);
		return NULL;
	}
	strcpy(log->fname, fname);
	strcat(log->fname, suffix);

	log->verbosity = verbosity;
	log->flags = flags;

	if ((flags & IO_LOGGING_FLAGS_USE_STDOUT_FOR_CRITICAL) != 0) {
		log->critical = stdout;
	} else {
		log->critical = stderr;
	}

	/* Clean up */
	free(suffix);
	free(fname);

	/* Check if the logfile already exists */
	log->logfile = fopen(log->fname, "r");
	if (log->logfile == NULL) {
		/* Does not exist yet, create a new one */
		log->logfile = fopen(log->fname, "w+");
		appending = 0;
	} else {
		/* Does exist, we want to append then */
		fclose(log->logfile);
		log->logfile = fopen(log->fname, "a");
		appending = 1;
    
#ifdef AHF2_overwrite_logfiles
    // this opens a new logfile irrespective of whether it exists or not
		log->logfile = fopen(log->fname, "w+");
		appending = 0;
#endif
    
	}
	if (log->logfile == NULL) {
		fprintf(stderr,
		        "%s, %u: Could not open %s for reading/writing.\n",
		        __func__, __LINE__, log->fname);
		free(log->fname);
		free(log);
		return NULL;
	}

	/* Make it clear that we are appending */
	if (appending == 1)
		fprintf(log->logfile, IO_LOGGING_NEW_ENTRY);

	/* Done! */
	return log;
}

extern void
io_logging_stop(io_logging_t *log)
{
	if (log == NULL || *log == NULL)
		return;

	if ( (*log)->logfile != NULL ) {
		fprintf( (*log)->logfile, IO_LOGGING_STOP);
		fclose( (*log)->logfile );
	}
	if ( (*log)->fname != NULL )
		free( (*log)->fname );

	free(*log);
	*log = NULL;

	return;
}

extern void
io_logging_hello(io_logging_t log, float version, int build)
{

	if (log->logfile != NULL) {
		fprintf(log->logfile, IO_LOGGING_BAR "\n");
		fprintf(log->logfile, IO_LOGGING_LOGO);
		fprintf(log->logfile, "(%3.1f/%03d)\n", version,build);
#		ifdef WITH_MPI
		{
			int rank, size;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			fprintf(log->logfile,
			        "\t\tProcess %i. Total number of processes: %i\n",
			        rank, size);
#			ifdef DEBUG
			fflush(log->logfile);
#			endif
		}
#		endif
		fprintf(log->logfile, IO_LOGGING_BAR "\n\n");
	}

	return;
}

extern void
io_logging_identify(io_logging_t log)
{
	int rank = 0;
	int size = 1;
	if (log->logfile != NULL) {
#		ifdef WITH_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
#		endif
		fprintf(log->logfile,
		        "This is libio version %i.%i.%i "
		        "running on %i processes. I have rank %i.\n",
		        IO_VERSION_MAJOR, IO_VERSION_MINOR,
		        IO_VERSION_MICRO, size, rank);
	}

	return;
}

extern void
io_logging_part(io_logging_t log,
                const char *frmt,
                ...)
{
	va_list lst;

	va_start(lst,frmt);

	/* Update and reset the counter */
	log->cnt_part++;
	log->cnt_sec = INT32_C(0);
	log->cnt_subsec = INT32_C(0);

	/* Print the part blobb */
	fprintf(log->logfile, IO_LOGGING_PART_PRE);
	fprintf(log->logfile,
	        IO_LOGGING_PART_NUM_PRE "%c" IO_LOGGING_PART_NUM_POST,
	        ((int)log->cnt_part));
	vfprintf(log->logfile, frmt, lst);
	fprintf(log->logfile, IO_LOGGING_PART_POST);
	fflush(log->logfile);

	va_end(lst);

	return;
}

extern void
io_logging_section(io_logging_t log,
                   const char *frmt,
                   ...)
{
	va_list lst;

	va_start(lst,frmt);

	/* Update and reset the counter */
	log->cnt_sec++;
	log->cnt_subsec=INT32_C(0);

	/* Print the section blobb */
	fprintf(log->logfile, IO_LOGGING_SECTION_PRE);
	fprintf(log->logfile,
	        IO_LOGGING_SECTION_NUM_PRE
	        "%c-%" PRIi32
	        IO_LOGGING_SECTION_NUM_POST,
	        (int)(log->cnt_part), log->cnt_sec);
	vfprintf(log->logfile, frmt, lst);
	fprintf(log->logfile, IO_LOGGING_SECTION_POST);
	fflush(log->logfile);

	va_end(lst);

	return;
}

extern void
io_logging_subsection(io_logging_t log,
                      const char *frmt,
                      ...)
{
	va_list lst;

	va_start(lst,frmt);

	/* Update the counter */
	log->cnt_subsec++;

	/* Print the section blobb */
	fprintf(log->logfile, IO_LOGGING_SUBSECTION_PRE);
	fprintf(log->logfile,
	        IO_LOGGING_SUBSECTION_NUM_PRE
	        "%c-%" PRIi32 ".%" PRIi32
	        IO_LOGGING_SUBSECTION_NUM_POST,
	        (int)(log->cnt_part), log->cnt_sec, log->cnt_subsec);
	vfprintf(log->logfile, frmt, lst);
	fprintf(log->logfile, IO_LOGGING_SUBSECTION_POST);
	fflush(log->logfile);

	va_end(lst);

	return;
}

extern void
io_logging_msg(io_logging_t log,
               int32_t verbosity,
               const char *frmt,
               ...)
{
	va_list lst;

	if (local_logcheck(log, verbosity) != 0)
		return;

	va_start(lst, frmt);

#	ifdef PRINT_VERBOSE_LEVEL
	fprintf(log->logfile, "%02" PRIi32 ": ", verbosity);
#	endif
	vfprintf(log->logfile, frmt, lst);
	fprintf(log->logfile, "\n");
	fflush(log->logfile);

	va_end(lst);

	return;
}

extern void
io_logging_msgplain(io_logging_t log,
                    int32_t verbosity,
                    const char *frmt,
                    ...)
{
	va_list lst;

	if (local_logcheck(log, verbosity) != 0)
		return;

	va_start(lst, frmt);

	vfprintf(log->logfile, frmt, lst);
#	ifdef DEBUG
	fflush(log->logfile);
#	endif

	va_end(lst);

	return;
}

extern void
io_logging_warn(io_logging_t log,
                int32_t verbosity,
                const char *frmt,
                ...)
{
	va_list lst;

	if (local_logcheck(log, verbosity) != 0)
		return;

	va_start(lst, frmt);

#	ifdef PRINT_VERBOSE_LEVEL
	fprintf(log->logfile, "%02" PRIi32 ": WARNING: ", verbosity);
#	else
	fprintf(log->logfile, "WARNING: ");
#	endif
	vfprintf(log->logfile, frmt, lst);
	fprintf(log->logfile, "\n");
#	ifdef DEBUG
	fflush(log->logfile);
#	endif

	va_end(lst);

	return;
}

extern void
io_logging_fatal(io_logging_t log,
                 const char *frmt,
                 ...)
{
	va_list lst;

	if (local_logcheck(log, INT32_C(0)) != 0)
		return;


	/* Write to logfile if requested */
	if (log->flags & IO_LOGGING_FLAGS_DUPLICATE_CRITICAL) {
		va_start(lst, frmt);

#		ifdef PRINT_VERBOSE_LEVEL
		fprintf(log->logfile, "%02" PRIi32 ": FATAL: ", INT32_C(0));
#		else
		fprintf(log->logfile, "FATAL: ");
#		endif
		vfprintf(log->logfile, frmt, lst);
		fprintf(log->logfile, "\n");
		fflush(log->logfile);

		va_end(lst);
	}

	va_start(lst, frmt);

	/* Write to the critical stream */
#	ifdef PRINT_VERBOSE_LEVEL
	fprintf(log->critical, "%02" PRIi32 ": FATAL: ", INT32_C(0));
#	else
	fprintf(log->critical, "FATAL: ");
#	endif
	vfprintf(log->critical, frmt, lst);
	fprintf(log->critical, "\n");
	fflush(log->critical);

	va_end(lst);

	return;
}

extern void
io_logging_memfatal(io_logging_t log,
                    const char *frmt,
                    ...)
{
	va_list lst;

	if (local_logcheck(log, INT32_C(0)) != 0)
		return;

	va_start(lst, frmt);

	/* Write to logfile if requested */
	if ((log->flags & IO_LOGGING_FLAGS_DUPLICATE_CRITICAL)==1) {
#		ifdef PRINT_VERBOSE_LEVEL
		fprintf(log->logfile,
		        "%02" PRIi32 ": FATAL: Could not allocate memory for ",
		        INT32_C(0));
#		else
		fprintf(log->logfile, "FATAL: Could not allocate memory for ");
#		endif
		vfprintf(log->logfile, frmt, lst);
		fprintf(log->logfile, ".\n");
#		ifdef DEBUG
		fflush(log->logfile);
#		endif
	}

	/* Write to the critical stream */
#	ifdef PRINT_VERBOSE_LEVEL
	fprintf(log->critical,
	        "%02" PRIi32 ": FATAL: Could not allocate memory for ",
	        INT32_C(0));
#	else
	fprintf(log->critical, "FATAL: Could not allocate memory for ");
#	endif
	vfprintf(log->critical, frmt, lst);
	fprintf(log->critical, "\n");
	fflush(log->critical);

	va_end(lst);

	return;
}

extern void
io_logging_flush(io_logging_t log)
{
	fflush(log->logfile);
}


/***********************************************************************\
 *    Implementation of local functions                                * 
\***********************************************************************/
#ifdef WITH_MPI
static int32_t
local_countdigits(int64_t i)
{
	int32_t dig;
	uint64_t div;

	dig = 1;
	div = UINT64_C(10);

	while ( i / div != 0) {
		dig++;
		div *= UINT64_C(10);
	}

	return dig;
}
#endif /* WITH_MPI */

inline static int8_t
local_logcheck(io_logging_t log, int32_t verbosity)
{
	if (log == NULL)
		return -1;

	if (log->verbosity < verbosity)
		return -1;

	return 0;
}
