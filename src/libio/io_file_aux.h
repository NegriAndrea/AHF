#ifndef IO_FILE_AUX_H
#define IO_FILE_AUX_H

/* $Id: io_file_aux.h,v 1.2 2007/02/11 20:12:20 knolli Exp $ */

/**
 * \file io_file_aux.h
 *
 * Provides some useful auxiliary datatypes for reading files.
 */

/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <inttypes.h>
#include <stddef.h>
#include "io_logging.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/

/** A generic element of an array */
struct io_file_val_struct {
	/** This points to the values that is stored */
	void *val;
	/** By how many bytes should be advanced to reach the next value */
	ptrdiff_t stride;
};

/** Convenient typedef */
typedef struct io_file_val_struct io_file_val_struct_t;

/** A generic description of the particle storage */
struct io_file_strg_struct {
	io_file_val_struct_t posx;
	io_file_val_struct_t posy;
	io_file_val_struct_t posz;
	io_file_val_struct_t momx;
	io_file_val_struct_t momy;
	io_file_val_struct_t momz;
	io_file_val_struct_t weight;
	io_file_val_struct_t id;
	io_file_val_struct_t u;
#ifdef METALHACK
	io_file_val_struct_t z;
	io_file_val_struct_t age;
#endif
	uint32_t bytes_float;
	uint32_t bytes_int;
};

/** Convenient typedef */
typedef struct io_file_strg_struct io_file_strg_struct_t;


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Will print information about the value to the logfile.
 *
 * \param log  The logging module.
 * \param val  The value.
 *
 * \return Nothing.
 */
#define io_file_val_log(log, value) {\
	io_logging_msg(log, INT32_C(3), "  val    = %p", value.val); \
	io_logging_msg(log, INT32_C(3), "  stride = %ti", value.stride); \
}

/**
 * \brief Will print information about the storage to the logfile.
 *
 * \param log   The logging module.
 * \param strg  The storage description.
 *
 * \return Nothing.
 */
#define io_file_strg_log(log, strg) {\
	io_logging_msg(log, INT32_C(3), "posx:"); \
	io_file_val_log(log, strg.posx); \
	io_logging_msg(log, INT32_C(3), "posy:"); \
	io_file_val_log(log, strg.posy); \
	io_logging_msg(log, INT32_C(3), "posz:"); \
	io_file_val_log(log, strg.posz); \
	io_logging_msg(log, INT32_C(3), "momx:"); \
	io_file_val_log(log, strg.momx); \
	io_logging_msg(log, INT32_C(3), "momy:"); \
	io_file_val_log(log, strg.momy); \
	io_logging_msg(log, INT32_C(3), "momz:"); \
	io_file_val_log(log, strg.momz); \
	io_logging_msg(log, INT32_C(3), "weight:"); \
	io_file_val_log(log, strg.weight); \
	io_logging_msg(log, INT32_C(3), "id:"); \
	io_file_val_log(log, strg.id); \
	io_logging_msg(log, INT32_C(3), "u:"); \
	io_file_val_log(log, strg.u); \
	io_logging_msg(log, INT32_C(3), \
	               "bytes_float: %" PRIu32, strg.bytes_float); \
	io_logging_msg(log, INT32_C(3), \
	               "bytes_int:   %" PRIu32, strg.bytes_int); \
}


#endif /* IO_FILE_AUX_H */
