/**
 * \file hilbert_util.c
 *
 * Useful Hilbert Index functions.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include <math.h>
#include "hilbert_util.h"


/***********************************************************************\
 *    Local defines, structure definitions and typedefs                * 
\***********************************************************************/


/***********************************************************************\
 *    Prototypes of local functions                                    * 
\***********************************************************************/


/***********************************************************************\
 *    Implementation of global functions                               * 
\***********************************************************************/
extern hikey_t
hilbert_util_contract(unsigned int trgt_level,
                      unsigned int src_level,
                      hikey_t hikey)
{
	hikey_t rtrn;

	rtrn = hikey;
	rtrn >>= (3*(src_level - trgt_level));

	return rtrn;
}

extern hikey_t
hilbert_util_prolongMin(unsigned int trgt_level,
                        unsigned int src_level,
                        hikey_t hikey)
{
	hikey_t rtrn;

	rtrn = hikey;
	rtrn <<= (3*(trgt_level - src_level));

	return rtrn;
}

extern hikey_t
hilbert_util_prolongMax(unsigned int trgt_level,
                        unsigned int src_level,
                        hikey_t hikey)
{
	hikey_t rtrn;
   uint64_t one=UINT64_C(1);
  
	rtrn = hilbert_util_prolongMin(trgt_level, src_level, hikey);
	rtrn += (((one)<<(3*(trgt_level - src_level))) - 1);

	return rtrn;
}

extern hikey_t
hilbert_util_calcHikey(double x,
                       double y,
                       double z,
                       unsigned bits)
{
	hikey_t coord[3];
	hikey_t max = 1L<<bits;

	/* This will scale the real positions to the maximal dynamical range
	 * given by max (number of cell pers dimension). It will then be
	 * truncated to the next lowest integer and also a check will be
	 * performed whether the calculated grid position is valid, this is
	 * needed, as it can happen, that the position is too close to the
	 * upper boundary to so that (x*max < max) might be violated.
	 */
	coord[0] = (hikey_t)trunc(x * max);
	coord[0] = (coord[0] == max) ? max-1 : coord[0];
	coord[1] = (hikey_t)trunc(y * max);
	coord[1] = (coord[1] == max) ? max-1 : coord[1];
	coord[2] = (hikey_t)trunc(z * max);
	coord[2] = (coord[2] == max) ? max-1 : coord[2];
	
	return (hikey_t)hilbert_c2i(3, bits, coord);
}

extern hikey_t
hilbert_util_calcHikey_grid(uint32_t x,
                            uint32_t y,
                            uint32_t z,
                            unsigned bits)
{
	hikey_t coord[3];

	coord[0] = (hikey_t)x;
	coord[1] = (hikey_t)y;
	coord[2] = (hikey_t)z;

	return (hikey_t)hilbert_c2i(3, bits, coord);
}

extern void
hilbert_util_calcPos(hikey_t key, unsigned bits, unsigned *pos)
{
	bitmask_t coord[3];

	hilbert_i2c(3, bits, (bitmask_t)key, coord);

	pos[0] = (unsigned)coord[0];
	pos[1] = (unsigned)coord[1];
	pos[2] = (unsigned)coord[2];

	return;
}

extern void
hilbert_util_getShell(hikey_t base, hikey_t *shell, unsigned bits)
{
	hikey_t basecoord[3];
	hikey_t coord[3];
	hikey_t num1dim;
	int i, j, k;
	int count;

	/* Calculate the number of possible Key values per dimension */
	num1dim = (1L << bits);

	/* Convert the Hilbert Key to a position */
	hilbert_i2c(3, bits, base, basecoord);

	/*****************************************************************\
	 * Loop over all cells within a cube of sidelength 3 centered on *
	 * the base cell. Starting at the lower left front corner, first *
	 * going in the x2 direction, than x1 and finally x0             *
	\*****************************************************************/
	count = 0;
	for (i=num1dim-1; i<=num1dim+1; i++) {
		coord[0] = (basecoord[0] + i) % num1dim;
		for (j=num1dim-1; j<=num1dim+1; j++) {
			coord[1] = (basecoord[1] + j) % num1dim;
			for (k=num1dim-1; k<=num1dim+1; k++) {
				coord[2] = (basecoord[2] + k) % num1dim;
				shell[count] = hilbert_c2i(3, bits, coord);
				count++;
			}
		}
	}

	return;
}


/***********************************************************************\
 *    Implementation of local functions                                * 
\***********************************************************************/
