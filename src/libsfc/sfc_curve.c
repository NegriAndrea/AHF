/* $Id: sfc_curve.c,v 1.2 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file sfc_curve.h
 *
 * Implementation of the space filling curve functions.
 */

/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include "sfc_curve.h"
#include "hilbert_util.h"


/***********************************************************************\
 *    Implemenation of global functions                                * 
\***********************************************************************/
extern const char*
sfc_curve_typestr(sfc_curve_t ctype)
{
	switch(ctype) {
		case SFC_CURVE_HILBERT:
			return SFC_CURVE_HILBERT_STR;
		default:
			return SFC_CURVE_UNKOWN_STR;
	}
}

extern sfc_key_t
sfc_curve_contract(uint32_t trgt_level,
                   uint32_t src_level,
                   sfc_curve_t ctype,
                   sfc_key_t key)
{
	switch(ctype) {
		case SFC_CURVE_HILBERT:
			return (sfc_key_t)hilbert_util_contract(trgt_level,
			                                        src_level,
			                                        (hikey_t)key);
		default:
			return 0;
	}
}

extern sfc_key_t
sfc_curve_prolongMin(uint32_t trgt_level,
                     uint32_t src_level,
                     sfc_curve_t ctype,
                     sfc_key_t key)
{
	switch(ctype) {
		case SFC_CURVE_HILBERT:
			return (sfc_key_t)hilbert_util_prolongMin(trgt_level,
			                                          src_level,
			                                          (hikey_t)key);
		default:
			return 0;
	}
}

extern sfc_key_t
sfc_curve_prolongMax(uint32_t trgt_level,
                     uint32_t src_level,
                     sfc_curve_t ctype,
                     sfc_key_t key)
{
	switch(ctype) {
		case SFC_CURVE_HILBERT:
			return (sfc_key_t)hilbert_util_prolongMax(trgt_level,
			                                          src_level,
			                                          (hikey_t)key);
		default:
			return 0;
	}
}

extern sfc_key_t
sfc_curve_calcKey(sfc_curve_t ctype,
                  double x,
                  double y,
                  double z,
                  uint32_t bits)
{
	switch(ctype) {
		case SFC_CURVE_HILBERT:
			return (sfc_key_t)hilbert_util_calcHikey(x, y, z, bits);
		default:
			return 0;
	}
}

extern sfc_key_t
sfc_curve_calcKey_grid(sfc_curve_t ctype,
                       uint32_t x,
                       uint32_t y,
                       uint32_t z,
                       uint32_t bits)
{
	switch(ctype) {
		case SFC_CURVE_HILBERT:
			return (sfc_key_t)hilbert_util_calcHikey_grid(x, y, z, bits);
		default:
			return 0;
	}
}

extern void
sfc_curve_calcPos(sfc_curve_t ctype,
                  sfc_key_t key,
                  uint32_t bits,
                  uint32_t *pos)
{
	switch(ctype) {
		case SFC_CURVE_HILBERT:
			hilbert_util_calcPos((hikey_t)key, (unsigned)bits,
			                     (unsigned *)pos);
			break;
		default:
			;
	}

	return;
}


extern void                                                                 
sfc_curve_getShell(sfc_curve_t ctype,
                   sfc_key_t base,
                   sfc_key_t *shell,
                   uint32_t bits)
{
	switch(ctype) {
		case SFC_CURVE_HILBERT:
			hilbert_util_getShell((hikey_t)base,
			                      (hikey_t *) shell, 
			                      bits);
			break;
		default:
			;
	}

	return;
}

extern int
sfc_curve_comp_key(const void *k1, const void *k2)
{
	if ( (*(const sfc_key_t *)k1) < (*(const sfc_key_t *)k2) )
		return -1;
	
	return ((*(const sfc_key_t *)k1) == (*(const sfc_key_t *)k2) ? 0 : 1);
}

