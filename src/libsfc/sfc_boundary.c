/* $Id: sfc_boundary.c,v 1.8 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file sfc_boundary.c
 *
 * Provides a way to figure out all bounding cells of a given space
 * filling curve segment.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include "sfc_boundary.h"
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>


/***********************************************************************\
 *    Local defines, structure definitions and typedefs                * 
\***********************************************************************/
#define SFC_BOUNDARY_MIN_SIZE 15
#define SFC_BOUNDARY_MIN_INC 5


/***********************************************************************\
 *    Prototypes of local functions                                    * 
\***********************************************************************/
/**
 * \brief Does a linear search on the key arrays to localize on which
 *        CPU a given key is located.
 *
 * \param key      The key to search for.
 * \param ncpu     The number of CPUs, hence the length of the following
 *                 arrays.
 * \param *fstkey  The array of all first keys sorted by CPU.
 * \param *lstkey  The array of all last keys sorted by CPU.
 *
 * \return Returns the CPU number on which the key is located.
 */
inline static uint32_t
local_localize_key(sfc_key_t key,
                   uint32_t ncpu,
                   sfc_key_t *fstkey,
                   sfc_key_t *lstkey);


/***********************************************************************\
 *    Implementation of global functions                               * 
\***********************************************************************/
extern const char*
sfc_boundary_typestr(sfc_boundary_type_t btype)
{
	switch(btype) {
		case SFC_BOUNDARY_TYPE_INNER:
			return SFC_BOUNDARY_TYPE_INNER_STR;
		case SFC_BOUNDARY_TYPE_OUTER:
			return SFC_BOUNDARY_TYPE_OUTER_STR;
	}
	return SFC_BOUNDARY_TYPE_UNKNOWN_STR;
}

extern sfc_boundary_t
sfc_boundary_new(io_logging_t log,
                 uint32_t num,
                 sfc_boundary_type_t btype,
                 sfc_key_t minkey,
                 sfc_key_t maxkey,
                 uint32_t bits,
                 sfc_curve_t ctype)
{
	sfc_boundary_t dummy;
	uint32_t i;

	/* Generate memory for the boundarie */
	dummy = (sfc_boundary_t)malloc(sizeof(sfc_boundary_struct_t)*num);
	if (dummy == NULL) {
		io_logging_memfatal(log, "boundary array");
		return NULL;
	}

	/* Initialize the boundaries */
	for (i=0; i<num; i++) {
		(dummy+i)->num = UINT64_C(0);
		(dummy+i)->len = UINT64_C(0);
		(dummy+i)->inc = SFC_BOUNDARY_MIN_INC;
		(dummy+i)->minkey = minkey;
		(dummy+i)->maxkey = maxkey;
		(dummy+i)->btype = btype;
		(dummy+i)->ctype = ctype;
		(dummy+i)->bits = bits;
		(dummy+i)->bound = NULL;
	}

	return dummy;
}

extern sfc_boundary_2_t
sfc_boundary_2_get(io_logging_t log,
                   uint32_t rank,
                   uint32_t ncpu,
                   sfc_key_t *fstkey,
                   sfc_key_t *lstkey,
                   uint32_t bits,
                   sfc_curve_t ctype)
{
	sfc_boundary_2_t dummy;
	sfc_key_t key, shell[27];
	int i, j;
	uint32_t is_on_cpu;

	/* Aquire the memory for the structure */
	dummy = (sfc_boundary_2_t)malloc(sizeof(sfc_boundary_2_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "boundary structure");
		return NULL;
	}

	/* Generate the real boundaries (key arrays) */
	dummy->outer = sfc_boundary_new(log, 1, SFC_BOUNDARY_TYPE_OUTER,
	                                fstkey[rank], lstkey[rank], bits,
	                                ctype);
	dummy->num_inner = ncpu;
	dummy->inner = sfc_boundary_new(log, ncpu, SFC_BOUNDARY_TYPE_INNER,
	                                fstkey[rank], lstkey[rank], bits,
	                                ctype);

	/* Now loop over all keys and put them in the right boundary array */
	for (key=fstkey[rank]; key<=lstkey[rank]; key++) {
		sfc_curve_getShell(ctype, key, shell, bits);
		for (i=0; i<27; i++) {
			is_on_cpu = local_localize_key(shell[i], ncpu,
			                               fstkey, lstkey);
			if (is_on_cpu != rank) {
				sfc_boundary_addKey(log, dummy->outer, shell[i]);
				sfc_boundary_addKey(log, dummy->inner+is_on_cpu, key);
			}
		}
	}

	return dummy;
}

extern sfc_boundary_t
sfc_boundary_get(io_logging_t log,
                 sfc_boundary_type_t btype,
                 sfc_key_t minkey,
                 sfc_key_t maxkey,
                 uint32_t bits,
                 sfc_curve_t ctype)
{
	sfc_boundary_t dummy;
	sfc_key_t i;
	sfc_key_t shell[27];
	int32_t j;

	/* Sanity checks */
	if (    (btype != SFC_BOUNDARY_TYPE_INNER)
	     && (btype != SFC_BOUNDARY_TYPE_OUTER) ) {
		io_logging_fatal(log, "Unsupported boundary type: %s",
		                 sfc_boundary_typestr(btype));
		return NULL;
	}

	/* Get the structure */
	dummy = (sfc_boundary_t)malloc(sizeof(sfc_boundary_struct_t));
	if (dummy == NULL) {
		io_logging_memfatal(log, "boundary structure");
		return NULL;
	}

	/* Set general stuff */
	dummy->num = 0;
	dummy->len = (uint64_t)((maxkey-minkey)/UINT64_C(10));
	dummy->inc = (uint32_t)(dummy->len / 5);
	dummy->minkey = minkey;
	dummy->maxkey = maxkey;
	dummy->btype = btype;
	dummy->ctype = ctype;
	dummy->bits = bits;
	if (dummy->len < SFC_BOUNDARY_MIN_SIZE)
		dummy->len = SFC_BOUNDARY_MIN_SIZE;
	if (dummy->inc < SFC_BOUNDARY_MIN_INC)
		dummy->inc = SFC_BOUNDARY_MIN_INC;


	/* Get a first boundary array */
	dummy->bound = (sfc_key_t *)malloc(sizeof(sfc_key_t)*dummy->len);
	if (dummy->bound == NULL) {
		free(dummy);
		io_logging_memfatal(log, "boundary key array");
		return NULL;
	}

	/* Get the boundary */
	switch (btype) {
		case SFC_BOUNDARY_TYPE_OUTER:
			/* Check all shell keys and add those not in the segment */
			for (i=minkey; i<=maxkey; i++) {
				sfc_curve_getShell(ctype, i, shell, bits);
				for (j=0; j<27; j++) {
					if ( (shell[j] < minkey) || (shell[j] > maxkey) )
						sfc_boundary_addKey(log, dummy, shell[j]);
				}
			}
			break;
		case SFC_BOUNDARY_TYPE_INNER:
			/* Check the shell of all keys, if one of the shell keys
			 * does not belong to the segment, add the shell center to
			 * the boundary */
			for (i=minkey; i<=maxkey; i++) {
				sfc_curve_getShell(ctype, i, shell, bits);
				for (j=0; j<27; j++) {
					if ( (shell[j] < minkey) || (shell[j] > maxkey) ) {
						sfc_boundary_addKey(log, dummy, i);
						break;
					}
				}
			}
			break;
	}

	return dummy;
}

extern void
sfc_boundary_del(sfc_boundary_t *bound)
{
	if (bound == NULL || *bound == NULL)
		return;

	free((*bound)->bound);
	free(*bound);
	*bound = NULL;

	return;
}

extern void
sfc_boundary_2_del(sfc_boundary_2_t *bound)
{
	uint32_t i;
	sfc_boundary_t bi;

	
	free((*bound)->outer->bound);
	for (i=0; i<(*bound)->num_inner; i++) {
		free((*bound)->inner[i].bound);
	}
	free((*bound)->outer);
	free((*bound)->inner);
	free(*bound);

	*bound = NULL;

	return;
}

extern sfc_boundary_t
sfc_boundary_addKey(io_logging_t log,
                    sfc_boundary_t bound, 
                    sfc_key_t key)
{
	uint64_t pos;
	sfc_key_t *tmp;

	/* See if the key is in the boundary and where to put it, if not */
	if (sfc_boundary_findKeyPos(log, bound, key, &pos))
		return bound;

	/* Make sure that the array is large enough */
	if (bound->num >= bound->len) {
#		ifdef DEBUG
		io_logging_msg(log, INT32_C(10),
		               "Resizing array, old: %" PRIu64, bound->len);
		io_logging_msg(log, INT32_C(10),
		               "Resizing array, inc: %" PRIu64, bound->inc);
#		endif
		bound->len += bound->inc;
#		ifdef DEBUG
		io_logging_msg(log, INT32_C(10),
		               "Resizing array, new size: %" PRIu64, bound->len);
#		endif
		tmp = (sfc_key_t *)realloc((void *)(bound->bound),
		                           bound->len * (sizeof(sfc_key_t)));
		if (tmp == NULL) {
			io_logging_memfatal(log, "resized boundary array");
			return NULL;
		}
		bound->bound = tmp;
	}

	/* Now move everything up if needed */
	if (bound->num > pos) {
		memmove(bound->bound+pos+UINT64_C(1),
		        bound->bound+pos,
		        (bound->num-pos)*sizeof(sfc_key_t));
	}
	bound->bound[pos] = key;
	bound->num++;

	/* Done */
	return bound;
}

extern bool
sfc_boundary_findKeyPos(io_logging_t log,
                        sfc_boundary_t bound,
                        sfc_key_t key,
                        uint64_t *retpos)
{
	uint64_t pos;
	uint64_t step;
	bool done = false;

	/* Setting first try position */
	pos = bound->num / UINT64_C(2);
	step = bound->num / UINT64_C(4);
	if (step == 0)
		step = 1;
#	ifdef DEBUG
	io_logging_msg(log, INT32_C(15),
	               "Trying to find key %" SFC_PRIkey " at %" PRIu64,
	               key, pos);
#	endif

	/* See if there is anything to look at */
	if (bound->num == UINT64_C(0)) {
		if (retpos != NULL)
			*retpos = UINT64_C(0);
		return false;
	}

	/* Do a binary search to find the right position to enter */
	while (!done) {
		if (bound->bound[pos] == key) {
#			ifdef DEBUG
			io_logging_msg(log, INT32_C(15),
		                   "Key already in the boundary");
#			endif
			 if (retpos != NULL)
				*retpos = pos;
			return true;
		}
		if (bound->bound[pos] < key) {
			if (    (pos == bound->num - UINT64_C(1))
			     || (bound->bound[pos+UINT64_C(1)] > key) ) {
				/* We have to fit the key at pos + 1 */
				pos += UINT64_C(1);
#				ifdef DEBUG
			    io_logging_msg(log, INT32_C(15),
				               "Found right position at %" PRIu64
				               " (approaching from below)",
				              pos);
#				endif
				done = true;
			} else {
				/* Take another step */
				pos += step;
				if (step > 1)
					step /= UINT64_C(2);
#				ifdef DEBUG
			    io_logging_msg(log, INT32_C(50),
				               "Stepping up, new pos: %" PRIu64,
				               pos);
#				endif
			}
		} else {
			if (    (pos == UINT64_C(0))
			     || (bound->bound[pos-UINT64_C(1)] < key) ) {
				/* We have to fit the key at pos */
#				ifdef DEBUG
			    io_logging_msg(log, INT32_C(15),
				               "Found right position at %" PRIu64
				               " (approaching from above)",
				               pos);
#				endif
				done = true;
			} else {
				/* Next step */
				pos -= step;
				if (step > 1)
					step /= UINT64_C(2);
#				ifdef DEBUG
			    io_logging_msg(log, INT32_C(50),
				               "Stepping down, new pos: %" PRIu64,
				               pos);
#				endif
			}
		}	
	} /* End of while loop */

	/* Set the position where the key should go, if requested */
	if (retpos != NULL)
		*retpos = pos;

	/* The key was not in the boundary */
	return false;
}

extern void
sfc_boundary_getBox(sfc_boundary_t bound,
                    uint32_t num,
                    uint32_t *p1,
                    uint32_t *p2)
{
	uint64_t i;
	uint32_t pos[3], j;
	bool boundary_exists;

	/* Initialize the coordinates */
	boundary_exists = false;
	for (j=0; j<num; j++)
		if ((bound+j)->num != UINT64_C(0))
			boundary_exists = true;
	if (!boundary_exists) {
		p1[0] = p1[1] = p1[2] = UINT32_C(0);
		p2[0] = p2[1] = p2[2] = (1<<bound->bits)-1;
		return;
	}
	p1[0] = p1[1] = p1[2] = UINT32_MAX;
	p2[0] = p2[1] = p2[2] = UINT32_C(0);

	/* Loop over all boundaries */
	for (j=0; j<num; j++) {
		/* Start looping over all boundary keys */
		for (i=UINT64_C(0); i<(bound+j)->num; i++) {
			/* Get the position...*/
			sfc_curve_calcPos((bound+j)->ctype, (bound+j)->bound[i],
			                  (bound+j)->bits, pos);
			/* ...and update the corner points if needed */
			p1[0] = (pos[0] < p1[0]) ? pos[0] : p1[0];
			p2[0] = (pos[0] > p2[0]) ? pos[0] : p2[0];
			p1[1] = (pos[1] < p1[1]) ? pos[1] : p1[1];
			p2[1] = (pos[1] > p2[1]) ? pos[1] : p2[1];
			p1[2] = (pos[2] < p1[2]) ? pos[2] : p1[2];
			p2[2] = (pos[2] > p2[2]) ? pos[2] : p2[2];
		}
	}

	return;
}

extern void
sfc_boundary_log(io_logging_t log,                                       
                 sfc_boundary_t bound)
{
	sfc_key_t i;

	io_logging_msg(log, INT32_C(5),
	               "Boundary Informantion:");
	io_logging_msg(log, INT32_C(5),
	               "Type:      : %s", sfc_boundary_typestr(bound->btype));
	io_logging_msg(log, INT32_C(5),
	               "num        : %" SFC_PRIkey, bound->num);
	io_logging_msg(log, INT32_C(5),
	               "len        : %" SFC_PRIkey, bound->len);
	io_logging_msg(log, INT32_C(5),
	               "inc        : %" SFC_PRIkey, bound->inc);
	io_logging_msg(log, INT32_C(5),
	               "minkey     : %" SFC_PRIkey, bound->minkey);
	io_logging_msg(log, INT32_C(5),
	               "maxkey     : %" SFC_PRIkey, bound->maxkey);
	io_logging_msg(log, INT32_C(5),
	               "Curve type : %s", sfc_curve_typestr(bound->ctype));
	if (bound->num == 0)
		return;
	io_logging_msg(log, INT32_C(5),
	               "List of all keys in the boundary follows");
	for (i=0; i<bound->num-1; i++) {
		io_logging_msgplain(log, INT32_C(5),
		                    "%" SFC_PRIkey ", ", bound->bound[i]);
	}
	io_logging_msgplain(log, INT32_C(5),
	                    "%" SFC_PRIkey "\n", bound->bound[i]);
	io_logging_msg(log, INT32_C(5),
	               "Quality    : %" SFC_PRIkey "/%" SFC_PRIkey " = %g",
	               bound->maxkey - bound->minkey + 1, 
	               bound->num,
	                 ((double)(bound->maxkey - bound->minkey + 1))
	               / ((double)(bound->num)));

	return;
}


/***********************************************************************\
 *    Implementation of local functions                                * 
\***********************************************************************/
inline static uint32_t
local_localize_key(sfc_key_t key,
                   uint32_t ncpu,
                   sfc_key_t *fstkey,
                   sfc_key_t *lstkey)
{
	uint32_t i;

	i=0;
	do {
		if (fstkey[i] <= key && lstkey[i] >= key)
			break;
		i++;
	} while (i<ncpu);

	return i;
}
