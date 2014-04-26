/* $Id: io_amiga_header.c,v 1.11 2007/02/11 19:42:42 knolli Exp $ */

/**
 * \file io_amiga_header.c
 *
 * Provides functions for reading and writing the header of AMIGA
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

#include "io_amiga_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/
/**
 * \brief Helper function to read a long variable, which might be 32bit
 *        or 64bit from a file.
 *
 * Only needed when reading the old header format. The new one only has
 * well defined integer variables.
 *
 * \param log    The logging module.
 * \param f      The file.
 * \param *trgt  The read variable.
 */
inline static void
local_getreadlong(io_logging_t log, io_amiga_t f, long *trgt);


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_amiga_header_t
io_amiga_header_get(io_logging_t log, io_amiga_t f)
{
	io_amiga_header_t dummy;
	int32_t trash;

	/* Some sanity checks */
	if ((f == NULL) || (f->file == NULL))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

	/* Create the header structure */
	dummy = (io_amiga_header_t)malloc(sizeof(io_amiga_header_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "AMIGA header structure");
		return NULL;
	}

	/* Go to the header */
	fseek(f->file, 4L, SEEK_SET);

	/* Now start reading the header */
	io_util_readstring(f->file,
	                   dummy->header,
	                   AMIGA_HEADER_HEADERSTRING);
	io_util_readint32(f->file, &(dummy->multi_mass), f->swapped);
	io_util_readint32(f->file, &(dummy->double_precision), f->swapped);

	local_getreadlong(log, f, (long *)&(dummy->no_part));
	local_getreadlong(log, f, (long *)&(dummy->no_species));

	io_util_readdouble(f->file, &(dummy->no_vpart), f->swapped);

	io_util_readdouble(f->file, &(dummy->timestep), f->swapped);
	io_util_readint32(f->file, &(dummy->no_timestep), f->swapped);
	if (f->file_sizeof_long == 8) {
		io_util_readint32(f->file, &trash, f->swapped);
	}

	io_util_readdouble(f->file, &(dummy->boxsize), f->swapped);
	io_util_readdouble(f->file, &(dummy->omega0), f->swapped);
	io_util_readdouble(f->file, &(dummy->lambda0), f->swapped);
	io_util_readdouble(f->file, &(dummy->pmass), f->swapped);

	io_util_readdouble(f->file, &(dummy->cur_reflevel), f->swapped);
	io_util_readdouble(f->file, &(dummy->cur_frcres), f->swapped);

	io_util_readdouble(f->file, &(dummy->a_initial), f->swapped);
	io_util_readdouble(f->file, &(dummy->a_current), f->swapped);

	io_util_readdouble(f->file, &(dummy->K_initial), f->swapped);
	io_util_readdouble(f->file, &(dummy->K_current), f->swapped);
	io_util_readdouble(f->file, &(dummy->U_initial), f->swapped);
	io_util_readdouble(f->file, &(dummy->U_current), f->swapped);
	io_util_readdouble(f->file, &(dummy->Eintegral), f->swapped);
	io_util_readdouble(f->file, &(dummy->Econst), f->swapped);

	io_util_readdouble(f->file, &(dummy->paramNSTEPS), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramNGRID_DOM), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramNth_dom), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramNth_ref), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramE_UPDATE), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramCELLFRAC_MAX), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramCELLFRAC_MIN), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramCA_CRIT), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramMAX_L1DIM), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramDOMSWEEPS), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramREFSWEEPS), f->swapped);

	io_util_readdouble(f->file, &(dummy->paramAHF_MINPART), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramAHF_VTUNE), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramAHF_RISE), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramAHF_SLOPE), f->swapped);
	io_util_readdouble(f->file, &(dummy->paramAHF_MAXNRISE), f->swapped);

	io_util_readdouble(f->file, &(dummy->min_weight), f->swapped);
	io_util_readdouble(f->file, &(dummy->max_weight), f->swapped);
	io_util_readdouble(f->file, &(dummy->t_unit), f->swapped);
	io_util_readdouble(f->file, &(dummy->B_init), f->swapped);
	io_util_readdouble(f->file, &(dummy->param_dummy5), f->swapped);
	io_util_readdouble(f->file, &(dummy->param_dummy6), f->swapped);
	io_util_readdouble(f->file, &(dummy->param_dummy7), f->swapped);
	io_util_readdouble(f->file, &(dummy->param_dummy8), f->swapped);
	io_util_readdouble(f->file, &(dummy->version), f->swapped);
	io_util_readint32(f->file, &(dummy->build), f->swapped);
	if (f->file_sizeof_long == 8) {
		io_util_readint32(f->file, &trash, f->swapped);
	}

	/* Read the fillheader correcting for a possible offset */
	io_util_readstring(f->file,
	                   dummy->dummy,
	                     AMIGA_HEADER_FILLHEADER
	                   - 2*(sizeof(uint64_t) - f->file_sizeof_long));

	f->header = dummy;

	return dummy;
}

extern io_amiga_header_t
io_amiga_header_new(io_logging_t log)
{
	io_amiga_header_t dummy;

	dummy = (io_amiga_header_t)calloc(sizeof(io_amiga_header_struct_t),
	                                  (size_t)1);
    if (dummy == NULL) {
		io_logging_memfatal(log, "AMIGA header structure");
		return NULL;
	}

	return dummy;
}

extern void
io_amiga_header_del(io_logging_t log, io_amiga_header_t *header)
{
	if ( (header == NULL) || (*header == NULL) )
		return;

	free(*header);

	*header = NULL;

	return;
}

extern void
io_amiga_header_write(io_logging_t log,
                      io_amiga_header_t header,
                      io_amiga_t f)
{
	int32_t tmp;

	if ( (header == NULL) || (f == NULL) )
		return;

	if (f->mode != IO_FILE_WRITE)
		return;

	if (f->file == NULL)
		return;

	if (header != f->header) {
		io_logging_msg(log, INT32_C(1),
		               "Writing a different header than stored in "
		               "the file object to the file");
	}

	fseek(f->file, 0L, SEEK_SET);

	tmp = sizeof(long);
	fwrite(&tmp, sizeof(int32_t), 1, f->file);
	fwrite(header, sizeof(io_amiga_header_struct_t), 1, f->file);

	return;
}

extern void
io_amiga_header_log(io_logging_t log, io_amiga_header_t header)
{
	io_logging_msg(log, INT32_C(5),
	               "Headerobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Header:                      %s",
	               header->header);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass:                   %" PRIi32,
	               header->multi_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Double precision:            %" PRIi32,
	               header->double_precision);
	io_logging_msg(log, INT32_C(5),
	               "  Number of particles:         %li",
	               header->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  Number of mass species:      %li",
	               header->no_species);
	io_logging_msg(log, INT32_C(5),
	               "  Number of virtual particles: %e",
	               header->no_vpart);
	io_logging_msg(log, INT32_C(5),
	               "  Timestep                   : %e",
	               header->timestep);
	io_logging_msg(log, INT32_C(5),
	               "  Number of timestep         : %" PRIi32,
	               header->no_timestep);
	io_logging_msg(log, INT32_C(5),
	               "  Boxsize                    : %e",
	               header->boxsize);
	io_logging_msg(log, INT32_C(5),
	               "  Omega0                     : %e",
	               header->omega0);
	io_logging_msg(log, INT32_C(5),
	               "  lambda0                    : %e",
	               header->lambda0);
	io_logging_msg(log, INT32_C(5),
	               "  pmass                      : %e",
	               header->pmass);
	io_logging_msg(log, INT32_C(5),
	               "  cur_reflevel               : %e",
	               header->cur_reflevel);
	io_logging_msg(log, INT32_C(5),
	               "  cur_frcres                 : %e",
	               header->cur_frcres);
	io_logging_msg(log, INT32_C(5),
	               "  a_initial                  : %e",
	               header->a_initial);
	io_logging_msg(log, INT32_C(5),
	               "  a_current                  : %e",
	               header->a_current);
	io_logging_msg(log, INT32_C(5),
	               "  K_initial                  : %e",
	               header->K_initial);
	io_logging_msg(log, INT32_C(5),
	               "  K_current                  : %e",
	               header->K_current);
	io_logging_msg(log, INT32_C(5),
	               "  U_initial                  : %e",
	               header->U_initial);
	io_logging_msg(log, INT32_C(5),
	               "  U_current                  : %e",
	               header->U_current);
	io_logging_msg(log, INT32_C(5),
	               "  Eintegral                  : %e",
	               header->Eintegral);
	io_logging_msg(log, INT32_C(5),
	               "  Econst                     : %e",
	               header->Econst);
	io_logging_msg(log, INT32_C(5),
	               "  paramNSTEPS                : %e",
	               header->paramNSTEPS);
	io_logging_msg(log, INT32_C(5),
	               "  paramNGRID_DOM             : %e",
	               header->paramNGRID_DOM);
	io_logging_msg(log, INT32_C(5),
	               "  paramNth_dom               : %e",
	               header->paramNth_dom);
	io_logging_msg(log, INT32_C(5),
	               "  paramNth_ref               : %e",
	               header->paramNth_ref);
	io_logging_msg(log, INT32_C(5),
	               "  paramE_UPDATE              : %e",
	               header->paramE_UPDATE);
	io_logging_msg(log, INT32_C(5),
	               "  paramCELLFRAC_MAX          : %e",
	               header->paramCELLFRAC_MAX);
	io_logging_msg(log, INT32_C(5),
	               "  paramCELLFRAC_MIN          : %e",
	               header->paramCELLFRAC_MIN);
	io_logging_msg(log, INT32_C(5),
	               "  paramCA_CRIT               : %e",
	               header->paramCA_CRIT);
	io_logging_msg(log, INT32_C(5),
	               "  paramMAX_L1DIM             : %e",
	               header->paramMAX_L1DIM);
	io_logging_msg(log, INT32_C(5),
	               "  paramDOMSWEEPS             : %e",
	               header->paramDOMSWEEPS);
	io_logging_msg(log, INT32_C(5),
	               "  paramREFSWEEPS             : %e",
	               header->paramREFSWEEPS);
	io_logging_msg(log, INT32_C(5),
	               "  paramAHF_MINPART           : %e",
	               header->paramAHF_MINPART);
	io_logging_msg(log, INT32_C(5),
	               "  paramAHF_VTUNE             : %e",
	               header->paramAHF_VTUNE);
	io_logging_msg(log, INT32_C(5),
	               "  paramAHF_RISE              : %e",
	               header->paramAHF_RISE);
	io_logging_msg(log, INT32_C(5),
	               "  paramAHF_SLOPE             : %e",
	               header->paramAHF_SLOPE);
	io_logging_msg(log, INT32_C(5),
	               "  paramAHF_MAXNRISE          : %e",
	               header->paramAHF_MAXNRISE);
	io_logging_msg(log, INT32_C(5),
	               "  min_weight                 : %e",
	               header->min_weight);
	io_logging_msg(log, INT32_C(5),
	               "  max_weight                 : %e",
	               header->max_weight);
	io_logging_msg(log, INT32_C(5),
	               "  t_unit                     : %e",
	               header->t_unit);
	io_logging_msg(log, INT32_C(5),
	               "  B_init                     : %e",
	               header->B_init);
	io_logging_msg(log, INT32_C(5),
	               "  param_dummy5               : %e",
	               header->param_dummy5);
	io_logging_msg(log, INT32_C(5),
	               "  param_dummy6               : %e",
	               header->param_dummy6);
	io_logging_msg(log, INT32_C(5),
	               "  param_dummy7               : %e",
	               header->param_dummy7);
	io_logging_msg(log, INT32_C(5),
	               "  param_dummy8               : %e",
	               header->param_dummy8);
	io_logging_msg(log, INT32_C(5),
	               "  version                    : %e",
	               header->version);
	io_logging_msg(log, INT32_C(5),
	               "  build                      : %" PRIi32,
	               header->build);

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
inline static void
local_getreadlong(io_logging_t log, io_amiga_t f, long *trgt)
{
	if (f->file_sizeof_long == sizeof(long)) {
		/* Easy here */
		if (sizeof(long) == 4)
			io_util_readint32(f->file, trgt, f->swapped);
		else
			io_util_readuint64(f->file, trgt, f->swapped);
	} else {
		if (f->file_sizeof_long <= (long)sizeof(long)) {
			/* Okay, not so bad, we can store it easily */
			if (f->file_sizeof_long == 4) {
				int32_t tmp;
				io_util_readint32(f->file, &tmp, f->swapped);
				*trgt = (long)tmp;
			} else {
				int64_t tmp;
				io_util_readuint64(f->file, &tmp, f->swapped);
				*trgt = (long)tmp;
			}
		} else {
			if (f->file_sizeof_long == sizeof(uint64_t)) {
				uint64_t tmp;
				io_logging_warn(log, INT32_C(0),
				                "Asked to read a long with "
				                "%i bytes but I can only store "
				                "%i bytes per long.",
				                f->file_sizeof_long, sizeof(long));
				io_util_readuint64(f->file, &tmp, f->swapped);
				*trgt = (long)tmp;
			} else {
				/* This should never happen. */
				io_logging_fatal(log,
				                 "Asked to read a long with "
				                 "%i bytes this is not supported. ",
				                f->file_sizeof_long);
			}
		}
	}

	return;
}
