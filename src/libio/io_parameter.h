#ifndef IO_PARAMETER_H
#define IO_PARAMETER_H

/**
 * \file io_parameter.h
 *
 * Provides functions for reading in AMIGA parameter files.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include "io_parameter_def.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/

/** The standard name for parameter file */
#define IO_PARAMETER_FNAME "amiga.input"


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Gets the parameters.
 *
 * When in MPI mode this function will check if the parameters
 * are going to be read from a file and only then proceed. Keep
 * in  mind that even though you can thus call the function from
 * each MPI task, the parameter files needs to be readable from
 * every process. And you might produce a lot of traffic.
 *
 * \param *fname  The filename of the parameter file. May be
 *                NULL, in which case a standarized filename is
 *                used.
 *
 * \return Returns a structure holding the information present
 *         in the file. Might be NULL in the case of errors.
 */
extern io_parameter_t
io_parameter_get(char *fname);

/**
 * \brief Disposes a parameter object.
 *
 * \param *param  Pointer to the variable holding the parameter
 *                object. This will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_parameter_del(io_parameter_t *params);

#endif /* IO_PARAMETER_H */
