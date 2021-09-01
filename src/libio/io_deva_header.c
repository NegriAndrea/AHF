
/**
 * \file io_deva_header.c
 *
 * Provides functions for reading and writing the header of DEVA
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

#include "io_deva_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/
#define SKIP(f) {fseek(f, 4L, SEEK_CUR);}

#define io_util_readintd(f,d,b) {if(deva_intsize == sizeof(int64_t)){int64_t dummyint0; io_util_readfortran(f, &dummyint0, deva_intsize, 1, 0, b);(d) = dummyint0;}\
				else if (deva_intsize == sizeof(int32_t)){int32_t dummyint0; io_util_readfortran(f, &dummyint0, deva_intsize, 1, 0, b); (d) = dummyint0;}}

#define io_util_readreald(f,d,b) {if(deva_realsize == sizeof(double)){double dummyreal0; io_util_readfortran(f, &dummyreal0, deva_realsize, 1, 0, b); (d) = dummyreal0;}\
				else if (deva_realsize == sizeof(float)){float dummyreal0; io_util_readfortran(f, &dummyreal0, deva_realsize, 1, 0, b); (d) = dummyreal0;}}
			


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/

/** Holds the size of integers in the file, global variable */
size_t deva_intsize;
/** Holds the size of reals in the file, global variable */
size_t deva_realsize;

extern int io_util_readfortran(FILE * const f, void * const buf0, size_t const ts, size_t const len, size_t const skip, int const swap);

/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_deva_header_t
io_deva_header_get_native(io_logging_t log, io_deva_t f)
{
	io_deva_header_t dummy;
	long skipsize;
  int iloop;

	/* Some sanity checks */
	if ((f == NULL) || (f->file == NULL))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

	/* Create the header structure array */
	dummy = (io_deva_header_t)malloc((size_t)sizeof(io_deva_header_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "DEVA header structure");
		return NULL;
	}

	/* just to be on the safe side: rewind the file to the start... */
	rewind(f->file);
	
	/* Now start reading the header */

	/* First find out integer size */
	int32_t deva_intsize0;
	f->swapped=IO_FILE_ISNOT_SWAPPED;
	io_util_readint32(f->file, &deva_intsize0, f->swapped);
	if(deva_intsize0>16){
		f->swapped=IO_FILE_IS_SWAPPED;
		io_util_sexchange(&deva_intsize0, sizeof(int32_t));
	}
	rewind(f->file);
	deva_intsize=deva_intsize0;



	io_util_readintd(f->file, dummy->itime, f->swapped);
	io_util_readintd(f->file, dummy->itstop, f->swapped);
	io_util_readintd(f->file, dummy->itdump, f->swapped);
	io_util_readintd(f->file, dummy->iout, f->swapped);
	io_util_readintd(f->file, dummy->nsformed, f->swapped);
	io_util_readintd(f->file, dummy->nsdead, f->swapped);
	io_util_readintd(f->file, dummy->irun, f->swapped);
	io_util_readintd(f->file, dummy->nobj, f->swapped);
	io_util_readintd(f->file, dummy->ngas, f->swapped);
	io_util_readintd(f->file, dummy->ndark, f->swapped);
	io_util_readintd(f->file, dummy->L, f->swapped);
	io_util_readintd(f->file, dummy->CHEMEVOL, f->swapped);
	io_util_readintd(f->file, dummy->ANOTHER, f->swapped);
	io_util_readintd(f->file, dummy->COOL, f->swapped);
	io_util_readintd(f->file, dummy->REFINEMENT, f->swapped);
	io_util_readintd(f->file, dummy->DEVA_HYDRO, f->swapped);
	io_util_readintd(f->file, dummy->GRAVITY, f->swapped);
	io_util_readintd(f->file, dummy->ISOLATED, f->swapped);
	io_util_readintd(f->file, dummy->EXPAND, f->swapped);
	io_util_readintd(f->file, dummy->COMOVING, f->swapped);
	io_util_readintd(f->file, dummy->STARFORM, f->swapped);
	io_util_readintd(f->file, dummy->GRADH, f->swapped);
	io_util_readintd(f->file, dummy->INITIALCOND, f->swapped);
	io_util_readintd(f->file, dummy->nstar, f->swapped);
	io_util_readintd(f->file, dummy->iseed1, f->swapped);
	io_util_readintd(f->file, dummy->ispec, f->swapped);
	io_util_readintd(f->file, dummy->indxsp, f->swapped);
	io_util_readintd(f->file, dummy->n_neigh, f->swapped);
	io_util_readintd(f->file, dummy->lastbar, f->swapped);
	for(iloop=0; iloop<100-29; iloop++)
		io_util_readintd(f->file, dummy->fill1[iloop], f->swapped);

	
	int32_t deva_realsize32;
	io_util_readint32(f->file, &deva_realsize32, f->swapped);
	deva_realsize = deva_realsize32;

	fseek(f->file,-sizeof(int32_t),SEEK_CUR);

	io_util_readreald(f->file, dummy->time, f->swapped);
	io_util_readreald(f->file, dummy->atime, f->swapped);
	io_util_readreald(f->file, dummy->htime, f->swapped);
	io_util_readreald(f->file, dummy->dtime, f->swapped);
	io_util_readreald(f->file, dummy->E_init, f->swapped);
	io_util_readreald(f->file, dummy->E_kin, f->swapped);
	io_util_readreald(f->file, dummy->E_ther, f->swapped);
	io_util_readreald(f->file, dummy->E_pot, f->swapped);
	io_util_readreald(f->file, dummy->Radiation, f->swapped);
	io_util_readreald(f->file, dummy->Esum, f->swapped);
	io_util_readreald(f->file, dummy->Rsum, f->swapped);
	io_util_readreald(f->file, dummy->cpu, f->swapped);
	io_util_readreald(f->file, dummy->time_end, f->swapped);
	io_util_readreald(f->file, dummy->tout, f->swapped);
	io_util_readreald(f->file, dummy->padding, f->swapped);
	io_util_readreald(f->file, dummy->Tlost, f->swapped);
	io_util_readreald(f->file, dummy->Qlost, f->swapped);
	io_util_readreald(f->file, dummy->Ulost, f->swapped);
	io_util_readreald(f->file, dummy->delta_min, f->swapped);
	io_util_readreald(f->file, dummy->delta_max, f->swapped);
	io_util_readreald(f->file, dummy->T_min, f->swapped);
	io_util_readreald(f->file, dummy->avisc, f->swapped);
	io_util_readreald(f->file, dummy->bvisc, f->swapped);
	io_util_readreald(f->file, dummy->eta2, f->swapped);
	io_util_readreald(f->file, dummy->rho_star, f->swapped);
	io_util_readreald(f->file, dummy->c_star, f->swapped);
	io_util_readreald(f->file, dummy->rmtot, f->swapped);
	io_util_readreald(f->file, dummy->rmsep, f->swapped);
	io_util_readreald(f->file, dummy->dnthres, f->swapped);
	io_util_readreald(f->file, dummy->sft0, f->swapped);
	io_util_readreald(f->file, dummy->sftmin, f->swapped);
	io_util_readreald(f->file, dummy->sftmax, f->swapped);
	io_util_readreald(f->file, dummy->h100, f->swapped);
	io_util_readreald(f->file, dummy->box100, f->swapped);
	io_util_readreald(f->file, dummy->rmgas, f->swapped);
	io_util_readreald(f->file, dummy->rmdark, f->swapped);
	io_util_readreald(f->file, dummy->omega0, f->swapped);
	io_util_readreald(f->file, dummy->xlambda0, f->swapped);
	io_util_readreald(f->file, dummy->h0t0, f->swapped);
	io_util_readreald(f->file, dummy->omegab0, f->swapped);
	io_util_readreald(f->file, dummy->sigma80, f->swapped);
	io_util_readreald(f->file, dummy->ztime0, f->swapped);
	io_util_readreald(f->file, dummy->e0, f->swapped);
	for(iloop=0; iloop<100-43; iloop++)
		io_util_readreald(f->file, dummy->fill2[iloop], f->swapped);

	f->header = dummy;

	return dummy;
}

extern io_deva_header_t
io_deva_header_get(io_logging_t log, io_deva_t f)
{
	io_deva_header_t dummy;
	long skipsize;
  int iloop;

	/* Some sanity checks */
	if ((f == NULL) || (f->file == NULL))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

	/* Create the header structure array */
	dummy = (io_deva_header_t)malloc((size_t)sizeof(io_deva_header_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "DEVA header structure");
		return NULL;
	}

	/* just to be on the safe side: rewind the file to the start... */
  rewind(f->file);
	
  /* Now start reading the header */
  io_util_readint32(f->file, &(dummy->itime), f->swapped);
  io_util_readint32(f->file, &(dummy->itstop), f->swapped);
  io_util_readint32(f->file, &(dummy->itdump), f->swapped);
  io_util_readint32(f->file, &(dummy->iout), f->swapped);
  io_util_readint32(f->file, &(dummy->nsformed), f->swapped);
  io_util_readint32(f->file, &(dummy->nsdead), f->swapped);
  io_util_readint32(f->file, &(dummy->irun), f->swapped);
  io_util_readint32(f->file, &(dummy->nobj), f->swapped);
  io_util_readint32(f->file, &(dummy->ngas), f->swapped);
  io_util_readint32(f->file, &(dummy->ndark), f->swapped);
  io_util_readint32(f->file, &(dummy->L), f->swapped);
  io_util_readfloat(f->file, &(dummy->CHEMEVOL), f->swapped);
  io_util_readfloat(f->file, &(dummy->ANOTHER), f->swapped);
  io_util_readfloat(f->file, &(dummy->COOL), f->swapped);
  io_util_readfloat(f->file, &(dummy->REFINEMENT), f->swapped);
  io_util_readfloat(f->file, &(dummy->DEVA_HYDRO), f->swapped);
  io_util_readfloat(f->file, &(dummy->GRAVITY), f->swapped);
  io_util_readint32(f->file, &(dummy->ISOLATED), f->swapped);
  io_util_readfloat(f->file, &(dummy->EXPAND), f->swapped);
  io_util_readfloat(f->file, &(dummy->COMOVING), f->swapped);
  io_util_readfloat(f->file, &(dummy->STARFORM), f->swapped);
  io_util_readfloat(f->file, &(dummy->GRADH), f->swapped);
  io_util_readint32(f->file, &(dummy->INITIALCOND), f->swapped);
  io_util_readint32(f->file, &(dummy->nstar), f->swapped);
  io_util_readint32(f->file, &(dummy->iseed1), f->swapped);
  io_util_readint32(f->file, &(dummy->ispec), f->swapped);
  io_util_readint32(f->file, &(dummy->indxsp), f->swapped);
  io_util_readint32(f->file, &(dummy->n_neigh), f->swapped);
  io_util_readint32(f->file, &(dummy->lastbar), f->swapped);
  for(iloop=0; iloop<100-29; iloop++)
    io_util_readint32(f->file, &(dummy->fill1[iloop]), f->swapped);

  io_util_readfloat(f->file, &(dummy->time), f->swapped);
  io_util_readfloat(f->file, &(dummy->atime), f->swapped);
  io_util_readfloat(f->file, &(dummy->htime), f->swapped);
  io_util_readfloat(f->file, &(dummy->dtime), f->swapped);
  io_util_readfloat(f->file, &(dummy->E_init), f->swapped);
  io_util_readfloat(f->file, &(dummy->E_kin), f->swapped);
  io_util_readfloat(f->file, &(dummy->E_ther), f->swapped);
  io_util_readfloat(f->file, &(dummy->E_pot), f->swapped);
  io_util_readfloat(f->file, &(dummy->Radiation), f->swapped);
  io_util_readfloat(f->file, &(dummy->Esum), f->swapped);
  io_util_readfloat(f->file, &(dummy->Rsum), f->swapped);
  io_util_readfloat(f->file, &(dummy->cpu), f->swapped);
  io_util_readfloat(f->file, &(dummy->time_end), f->swapped);
  io_util_readfloat(f->file, &(dummy->tout), f->swapped);
  io_util_readfloat(f->file, &(dummy->padding), f->swapped);
  io_util_readfloat(f->file, &(dummy->Tlost), f->swapped);
  io_util_readfloat(f->file, &(dummy->Qlost), f->swapped);
  io_util_readfloat(f->file, &(dummy->Ulost), f->swapped);
  io_util_readfloat(f->file, &(dummy->delta_min), f->swapped);
  io_util_readfloat(f->file, &(dummy->delta_max), f->swapped);
  io_util_readfloat(f->file, &(dummy->T_min), f->swapped);
  io_util_readfloat(f->file, &(dummy->avisc), f->swapped);
  io_util_readfloat(f->file, &(dummy->bvisc), f->swapped);
  io_util_readfloat(f->file, &(dummy->eta2), f->swapped);
  io_util_readfloat(f->file, &(dummy->rho_star), f->swapped);
  io_util_readfloat(f->file, &(dummy->c_star), f->swapped);
  io_util_readfloat(f->file, &(dummy->rmtot), f->swapped);
  io_util_readfloat(f->file, &(dummy->rmsep), f->swapped);
  io_util_readfloat(f->file, &(dummy->dnthres), f->swapped);
  io_util_readfloat(f->file, &(dummy->sft0), f->swapped);
  io_util_readfloat(f->file, &(dummy->sftmin), f->swapped);
  io_util_readfloat(f->file, &(dummy->sftmax), f->swapped);
  io_util_readfloat(f->file, &(dummy->h100), f->swapped);
  io_util_readfloat(f->file, &(dummy->box100), f->swapped);
  io_util_readfloat(f->file, &(dummy->rmgas), f->swapped);
  io_util_readfloat(f->file, &(dummy->rmdark), f->swapped);
  io_util_readfloat(f->file, &(dummy->omega0), f->swapped);
  io_util_readfloat(f->file, &(dummy->xlambda0), f->swapped);
  io_util_readfloat(f->file, &(dummy->h0t0), f->swapped);
  io_util_readfloat(f->file, &(dummy->omegab0), f->swapped);
  io_util_readfloat(f->file, &(dummy->sigma80), f->swapped);
  io_util_readfloat(f->file, &(dummy->ztime0), f->swapped);
  io_util_readfloat(f->file, &(dummy->e0), f->swapped);
  for(iloop=0; iloop<100-43; iloop++)
    io_util_readfloat(f->file, &(dummy->fill2[iloop]), f->swapped);

	f->header = dummy;

	return dummy;
}

extern void
io_deva_header_del(io_logging_t log, io_deva_header_t *header)
{
	if ( (header == NULL) || (*header == NULL) )
		return;

	free(*header);

	*header = NULL;

	return;
}

extern void
io_deva_header_write(io_logging_t log,
                       io_deva_header_t header,
                       io_deva_t f)
{
	if ( (header == NULL) || (f == NULL))
		return;

	if (f->mode != IO_FILE_WRITE)
		return;

	if (f->file == NULL)
		return;

	if (header != f->header) {
		io_logging_msg(log, INT32_C(1),
		               "Writing a different header than stored in "
		               "the file object to the file.");
	}
/*  // just to be on the safe side: rewind the file to the start... 
  rewind(f->file);
	
  // Now start reading the header 
  io_util_writeint32d(f->file, &(dummy->itime), f->swapped);
  io_util_writeint32d(f->file, &(dummy->itstop), f->swapped);
  io_util_writeint32d(f->file, &(dummy->itdump), f->swapped);
  io_util_writeint32d(f->file, &(dummy->iout), f->swapped);
  io_util_writeint32d(f->file, &(dummy->nsformed), f->swapped);
  io_util_writeint32d(f->file, &(dummy->nsdead), f->swapped);
  io_util_writeint32d(f->file, &(dummy->irun), f->swapped);
  io_util_writeint32d(f->file, &(dummy->nobj), f->swapped);
  io_util_writeint32d(f->file, &(dummy->ngas), f->swapped);
  io_util_writeint32d(f->file, &(dummy->ndark), f->swapped);
  io_util_writeint32d(f->file, &(dummy->L), f->swapped);
  io_util_writeint32d(f->file, &(dummy->CHEMEVOL), f->swapped);
  io_util_writeint32d(f->file, &(dummy->ANOTHER), f->swapped);
  io_util_writeint32d(f->file, &(dummy->COOL), f->swapped);
  io_util_writeint32d(f->file, &(dummy->REFINEMENT), f->swapped);
  io_util_writeint32d(f->file, &(dummy->DEVA_HYDRO), f->swapped);
  io_util_writeint32d(f->file, &(dummy->GRAVITY), f->swapped);
  io_util_writeint32d(f->file, &(dummy->ISOLATED), f->swapped);
  io_util_writeint32d(f->file, &(dummy->EXPAND), f->swapped);
  io_util_writeint32d(f->file, &(dummy->COMOVING), f->swapped);
  io_util_writeint32d(f->file, &(dummy->STARFORM), f->swapped);
  io_util_writeint32d(f->file, &(dummy->GRADH), f->swapped);
  io_util_writeint32d(f->file, &(dummy->INITIALCOND), f->swapped);
  io_util_writeint32d(f->file, &(dummy->nstar), f->swapped);
  io_util_writeint32d(f->file, &(dummy->iseed1), f->swapped);
  io_util_writeint32d(f->file, &(dummy->ispec), f->swapped);
  io_util_writeint32d(f->file, &(dummy->indxsp), f->swapped);
  io_util_writeint32d(f->file, &(dummy->n_neigh), f->swapped);
  io_util_writeint32d(f->file, &(dummy->lastbar), f->swapped);
  for(iloop=0; iloop<100-29; iloop++)
    io_util_writeint32d(f->file, &(dummy->fill1[iloop]), f->swapped);

  io_util_writereald(f->file, &(dummy->time), f->swapped);
  io_util_writereald(f->file, &(dummy->atime), f->swapped);
  io_util_writereald(f->file, &(dummy->htime), f->swapped);
  io_util_writereald(f->file, &(dummy->dtime), f->swapped);
  io_util_writereald(f->file, &(dummy->E_init), f->swapped);
  io_util_writereald(f->file, &(dummy->E_kin), f->swapped);
  io_util_writereald(f->file, &(dummy->E_ther), f->swapped);
  io_util_writereald(f->file, &(dummy->E_pot), f->swapped);
  io_util_writereald(f->file, &(dummy->Radiation), f->swapped);
  io_util_writereald(f->file, &(dummy->Esum), f->swapped);
  io_util_writereald(f->file, &(dummy->Rsum), f->swapped);
  io_util_writereald(f->file, &(dummy->cpu), f->swapped);
  io_util_writereald(f->file, &(dummy->time_end), f->swapped);
  io_util_writereald(f->file, &(dummy->tout), f->swapped);
  io_util_writereald(f->file, &(dummy->padding), f->swapped);
  io_util_writereald(f->file, &(dummy->Tlost), f->swapped);
  io_util_writereald(f->file, &(dummy->Qlost), f->swapped);
  io_util_writereald(f->file, &(dummy->Ulost), f->swapped);
  io_util_writereald(f->file, &(dummy->delta_min), f->swapped);
  io_util_writereald(f->file, &(dummy->delta_max), f->swapped);
  io_util_writereald(f->file, &(dummy->T_min), f->swapped);
  io_util_writereald(f->file, &(dummy->avisc), f->swapped);
  io_util_writereald(f->file, &(dummy->bvisc), f->swapped);
  io_util_writereald(f->file, &(dummy->eta2), f->swapped);
  io_util_writereald(f->file, &(dummy->rho_star), f->swapped);
  io_util_writereald(f->file, &(dummy->c_star), f->swapped);
  io_util_writereald(f->file, &(dummy->rmtot), f->swapped);
  io_util_writereald(f->file, &(dummy->rmsep), f->swapped);
  io_util_writereald(f->file, &(dummy->dnthres), f->swapped);
  io_util_writereald(f->file, &(dummy->sft0), f->swapped);
  io_util_writereald(f->file, &(dummy->sftmin), f->swapped);
  io_util_writereald(f->file, &(dummy->sftmax), f->swapped);
  io_util_writereald(f->file, &(dummy->h100), f->swapped);
  io_util_writereald(f->file, &(dummy->box100), f->swapped);
  io_util_writereald(f->file, &(dummy->rmgas), f->swapped);
  io_util_writereald(f->file, &(dummy->rmdark), f->swapped);
  io_util_writereald(f->file, &(dummy->omega0), f->swapped);
  io_util_writereald(f->file, &(dummy->xlambda0), f->swapped);
  io_util_writereald(f->file, &(dummy->h0t0), f->swapped);
  io_util_writereald(f->file, &(dummy->omegab0), f->swapped);
  io_util_writereald(f->file, &(dummy->sigma80), f->swapped);
  io_util_writereald(f->file, &(dummy->ztime0), f->swapped);
  io_util_writereald(f->file, &(dummy->e0), f->swapped);
  for(iloop=0; iloop<100-43; iloop++)
    io_util_writereald(f->file, &(dummy->fill2[iloop]), f->swapped);*/

	/* TODO: Write the header */

	return;
}


extern void
io_deva_header_log(io_logging_t log, io_deva_header_t header)
{
  io_logging_msg(log, INT32_C(5),
                 "  itime:                         %" PRIi32,
                 header->itime);
  io_logging_msg(log, INT32_C(5),
                 "  irun:                          %" PRIi32,
                 header->irun);
  io_logging_msg(log, INT32_C(5),
                 "  nobj:                          %" PRIi32,
                 header->nobj);
  io_logging_msg(log, INT32_C(5),
                 "  ndark:                         %" PRIi32,
                 header->ndark);
  io_logging_msg(log, INT32_C(5),
                 "  ngas:                          %" PRIi32,
                 header->ngas);
  io_logging_msg(log, INT32_C(5),
                 "  nstar:                         %" PRIi32,
                 header->nstar);
  io_logging_msg(log, INT32_C(5),
                 "  ztime:                         %e",
                 1./header->atime-1.);
  io_logging_msg(log, INT32_C(5),
                 "  atime:                         %e",
                 header->atime);
  io_logging_msg(log, INT32_C(5),
                 "  b100:                          %e",
                 header->box100);
  io_logging_msg(log, INT32_C(5),
                 "  omega0:                        %e",
                 header->omega0);
  io_logging_msg(log, INT32_C(5),
                 "  xlambda0:                      %e",
                 header->xlambda0);
  io_logging_msg(log, INT32_C(5),
                 "  omegab0:                       %e",
                 header->omegab0);
  io_logging_msg(log, INT32_C(5),
                 "  h100:                          %e",
                 header->h100);
  io_logging_msg(log, INT32_C(5),
                 "  sigma80:                       %e",
                 header->sigma80);
  io_logging_msg(log, INT32_C(5),
                 "  h0t0:                          %e",
                 header->h0t0);
  io_logging_msg(log, INT32_C(5),
                 "  rmtot:                         %e",
                 header->rmtot);
  io_logging_msg(log, INT32_C(5),
                 "  rmdark:                        %e",
                 header->rmdark);
  io_logging_msg(log, INT32_C(5),
                 "  rmgas:                         %e",
                 header->rmgas);
  
  return;
}


#define MIN(A,B)        ((A)<(B)?(A):(B))
#define ABS(A)        ((A)< 0?-(A):(A))

extern int
io_util_readfortran(FILE * const f, void * const buf0, size_t const ts, size_t const len, size_t const skip, int const swap)
{
	int32_t s1,s1p,s2,s2p; //Fortran heading and trailing record sizes
	size_t i,read=0,readt=0,remaining=len*ts,toRead,toSkip,remainingSkip=skip*ts; //counter
	unsigned char * buf= (unsigned char *) buf0;

	do
	{
		if (fread(&s1,sizeof(int32_t),1,f) != 1) return readt;
		//		printf("record size: %d ",s1);
		s1p=ABS(s1);
		if(swap) io_util_sexchange(&s1,sizeof(int32_t));
		toRead = s1p;
		toSkip = MIN(remainingSkip,toRead);
		//		printf("will skip: %ld ",toSkip);
		fseek(f,toSkip,SEEK_CUR);
		remainingSkip -= toSkip;
		toRead -= toSkip;
		toRead = MIN(toRead,remaining);
		//		printf("will read: %ld ",toRead);
		if(toRead > 0){
			if (fread(buf,1,toRead,f) != toRead) return readt;
			buf += toRead;
			read += toRead;
			readt=read/ts;
			remaining -= toRead;
		}
		//		printf("will skip %ld until the end. ",(s1p-(toSkip+toRead)));
		fseek(f,(s1p-(toSkip+toRead)),SEEK_CUR); //Skip to the end of the current slice
		if(fread(&s2,sizeof(int32_t),1,f)!=1) return readt;
		s2p=ABS(s2);
		//		printf("record size: %d ",s2);
		if(swap) io_util_sexchange(&s2,sizeof(int32_t));
		if(s1p != s2p){
			printf("Something went wrong: |%d| != |%d|",s1,s2);
			return readt; //Sanity check
		}

	} while (s1 < 0 /*&& toRead > 0*/);

	if(swap){
		buf=buf0;
		for(i=0;i<len;i++) 
		{
			io_util_sexchange(buf,ts);
			buf += ts;
		}
	}

	if(readt!=len)
		printf("read != len, %ld != %ld\n",readt,len);
	return readt;

}

/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
