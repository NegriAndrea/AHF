/* See ahf_halos_sfc.h for a general description. */

#ifdef AHF2

/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
/** Okay, we need those to catch all the defines etc. */
#include "../define.h"
#include "../param.h"
#include "../common.h"

/** Standard includes */
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <math.h>
#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>
#include "ahf_halos_sfc.h"
#include "ahf_halos.h"

/** Needed for the extended bsearch and the required compare function */
#include "../libutility/util_bsearch.h"
#include "../libutility/specific.h"

/** This is required for the HALO structure */
#include "../tdef.h"

/** We need access to the SFC function here */
#include "../libsfc/sfc_curve.h"


/***********************************************************************\
 *    Local defines                                                    * 
\***********************************************************************/
/*
 * Multiply the actual gather radius by this factor to catch a bit more
 * particles before checking whether they are within the gathering
 * raduis.  This is meant to ensure that the initial particle selection
 * based on the SFC doesn't already, due to numerical effects, exclude
 * some particles from being considered to be halo particles.
 */ 
#define GATHERRAD_FAC 1.001

/*
 * The ipart array of the halo that holds its particles is increased by
 * this many elements every time the array end is reached.  This is
 * meant to minimize to calls to realloc().  Theoretically, this can be
 * chosen quite large as firstly (depending on the system) we can
 * actually overcommit memory and it will work as long as we don't use
 * it.  And secondly, at the end, the ipart array is resized to its
 * actual size anyways.  If you run into problems with memory
 * allocation, try to lower this number.  Note that BUFFERSIZE 1 will
 * call realloc() each time a particle is added, so that is the lowest
 * memory footprint that can be achieved here at the cost of *a lot* of
 * realloc() calls.
 */
#define BUFFERSIZE 100

/*
 * This is used to calculate the dx, dy and dz of two vectors taking
 * into account the periodicity of the box.  Note, a and b must be
 * arrays of at least 3 elements (and only the first three elements will
 * be used).
 */
#define LOCAL_CAL_DXYZ(a,b)\
	dx = fabs((a)[0] - (b)[0]); \
	dy = fabs((a)[1] - (b)[1]); \
	dz = fabs((a)[2] - (b)[2]); \
	if (isgreater(dx,0.5)) \
		dx = 1.0-dx; \
	if (isgreater(dy,0.5)) \
		dy = 1.0-dy; \
	if (isgreater(dz,0.5)) \
		dz = 1.0-dz;


/***********************************************************************\
 *    Global variables                                                 * 
\***********************************************************************/
extern double r_fac, x_fac, v_fac, m_fac, rho_fac, phi_fac, Hubble;



/***********************************************************************\
 *    Definitions of local functions                                   * 
\***********************************************************************/
static uint32_t
local_getSuitableBits(HALO *halo, double gatherRad);

static bool
local_checkCellDistance(sfc_key_t key,
                        double centre[],
                        double gatherRad,
                        sfc_curve_t ctype,
                        uint32_t bits);

static void
local_checkAndAddParticles(sfc_key_t key,
                           HALO *halo,
                           double centre[],
                           double gatherRad2,
                           sfc_curve_t ctype,
                           uint32_t bits);

static uint64_t
local_findMinOffset(partptr fstPart,
                    uint64_t numParts,
                    sfc_key_t minKey);


/***********************************************************************\
 *    Implemenation of exported functions                              * 
\***********************************************************************/
void
ahf_halos_sfc_constructHalo(HALO *halo)
{
  /* if the grid-tree did not collect any particles on spatialRef[][] it makes no sense to construct a halo */
  if(halo->npart == 0)
    return;
  
	/* Get all particles inside the gather radius */
	ahf_halos_sfc_gatherParts(halo);
	
	/* Sort the particles */
	sort_halo_particles(halo);

	/* Remove particle outside of Rvir (re-estimating Rvir), 1st call */
	rem_outsideRvir(halo, 0);

#ifdef AHFnoremunbound
   //rem_nothing(halo);
#else
	/* Remove the unbound particles */
	rem_unbound(halo);
#endif

	/* And remove particles again, 2nd call */
	rem_outsideRvir(halo, 1);

	/* Calculate the profiles */
	HaloProfiles(halo);

#	ifdef AHFphspdens
	/* Do the phase space stuff */
	HaloProfilesPhaseSpace(halo);
#	endif
  
#ifdef AHFdisks
	HaloProfilesDisk(halo);
#endif

	/* Done */
	return;
}

void
ahf_halos_sfc_gatherParts(HALO *halo)
{
	double centre[3];
	sfc_key_t centreKey;
	sfc_key_t shell[27];
	uint32_t bits;
	double gatherRad, gatherRad2;
	sfc_curve_t ctype;
	int i;

#ifdef VERBOSE2
	fprintf(io.logfile,"\n    ahf_halos_sfc_constructHalo:   x=%g y=%g z=%g r=%g (npart=%d)  ",
          halo->pos.x*x_fac,halo->pos.y*x_fac,halo->pos.z*x_fac,halo->gatherRad*x_fac,halo->npart);
	fflush(io.logfile);
#endif
  
	/* Init the halo structure */
	halo->nll   = 0;
	halo->ll    = NULL;
	halo->npart = 0;
	halo->ipart = NULL;
  
	/*--- INIT --------------------------------------------------------*/
	/* Do that once: Set the SFC type */
#	if (!defined WITH_MPI)
	ctype = global_info.ctype;
#	else
	ctype = global_info.loadbal->ctype;
#	endif
	/* Calculate/Copy the bounding box/centre positions */
	centre[0] = halo->pos.x;
	centre[1] = halo->pos.y;
	centre[2] = halo->pos.z;
   
	/* Compare to the squared gather radius so that no sqrt is needed */
	gatherRad  = halo->gatherRad;
	gatherRad2 = gatherRad*gatherRad;

	/*--- SETUP -------------------------------------------------------*/
	/* Now calculate the most suitable coarseness of the search grid */
	bits = local_getSuitableBits(halo, halo->gatherRad);
	/* Get the key of the cell containing the center */
	centreKey = sfc_curve_calcKey(ctype ,centre[0], centre[1],
	                              centre[2], bits);
	/* Get the shell around this cell */
	sfc_curve_getShell(ctype, centreKey, shell, bits);

	/*--- GATHERING ---------------------------------------------------*/
	/* Now loop over all possible cells and check which of their
	 * particles ought to be added to the halo. */
  
  if(bits == 1)
   {
    /* this loops over all cells inside the box */
    for (i=0; i<8; i++) {
      //if (local_checkCellDistance(i, centre, gatherRad, ctype,bits))
        local_checkAndAddParticles(i, halo, centre, gatherRad2, ctype, bits);
    }
   }
  else
   {
    /* this loops over all direct neighbouring cells */
    for (i=0; i<27; i++) {
      if (local_checkCellDistance(shell[i], centre, gatherRad, ctype,bits))
        local_checkAndAddParticles(shell[i], halo, centre, gatherRad2, ctype, bits);
    }
   }
    
    
	/*--- FINISHING ---------------------------------------------------*/
	/* Done we are */
#ifdef VERBOSE2
	fprintf(io.logfile,"  -> collected %ld particles\n", halo->npart);
	fflush(io.logfile);
#endif
  
	fflush(io.logfile);
	return;
}


/***********************************************************************\
 *    Implementation of local functions                                * 
\***********************************************************************/
static uint32_t
local_getSuitableBits(HALO *halo, double gatherRad)
{
	uint32_t bits = 1;

	while (    (GATHERRAD_FAC * gatherRad < 1./(1<<(bits+1))) 
	        && ((bits+1)<=BITS_PER_DIMENSION) )
		bits++;
//	if (bits == 1) {
//		fprintf(stderr,
//		        "This is not good, the search grid is the whole box:-(\n");
//	}
	return bits;
}

static bool
local_checkCellDistance(sfc_key_t key,
                        double centre[],
                        double gatherRad,
                        sfc_curve_t ctype,
                        uint32_t bits)
{
	uint32_t cellGridPos[3];
	double cellRealPos[3];
	double gridToRealFac;
	double cellLargestRadius;
	double dx, dy, dz, dist2;
   uint64_t one=UINT64_C(1);

	/*
	 * Calculate the grid scale and the largest radius in the cell,
	 * which is given as the distance between the center of the cell and
	 * one of the corner points.
	 */
	gridToRealFac = 1./(one<<bits);
	cellLargestRadius = 0.5*sqrt(3.)*gridToRealFac;

	/* Get the grid coordinate of the cell and convert it to real-space */
	sfc_curve_calcPos(ctype, key, bits, cellGridPos);
	cellRealPos[0] = (gridToRealFac*cellGridPos[0])+0.5*gridToRealFac;
	cellRealPos[1] = (gridToRealFac*cellGridPos[1])+0.5*gridToRealFac;
	cellRealPos[2] = (gridToRealFac*cellGridPos[2])+0.5*gridToRealFac;

	/* The squared distance between the cell center and the halo center */
	LOCAL_CAL_DXYZ(cellRealPos, centre);
	dist2 = dx*dx + dy*dy + dz*dz;

	/*
	 * If the halo center is further away than the gather radius
	 * enlarged by largest cell radius, then the cell is too far away to
	 * hold any halo particles, hence ignore the cell 
	 */
	if (isgreater(sqrt(dist2), cellLargestRadius+GATHERRAD_FAC*gatherRad))
		return false;

	/* Well... seems we need to deal with this cell after all */
	return true;
}

static void
local_checkAndAddParticles(sfc_key_t key,
                           HALO *halo,
                           double centre[],
                           double gatherRad2,
                           sfc_curve_t ctype,
                           uint32_t bits)
{
	sfc_key_t minKey, maxKey;
	uint64_t minOffset;
	double dist2, dx, dy, dz;
	uint64_t partIn, partOut, partTot, offset;
	unsigned long *tmp;

	/* Find out where the particles start and end */
	minKey = sfc_curve_prolongMin(BITS_PER_DIMENSION, bits, ctype, key);
	maxKey = sfc_curve_prolongMax(BITS_PER_DIMENSION, bits, ctype, key);
	minOffset = local_findMinOffset(global_info.fst_part,
	                                global_info.no_part,
	                                minKey);

	/* Now do check */
	partTot = partIn = partOut = 0;
	offset = minOffset;
   
	while (    (offset<(global_info.no_part)
	        && ((global_info.fst_part+offset)->sfckey <= maxKey)) ) {
		LOCAL_CAL_DXYZ((global_info.fst_part+offset)->pos, centre);
		dist2 = dx*dx + dy*dy + dz*dz;
      
		if (islessequal(dist2, gatherRad2)) {
			/* Add particle*/
      //fprintf(stderr,"\ndist2=%lf gatherRad2=%lf partIn=%ld",dist2,gatherRad2,partIn);
			if (partIn%BUFFERSIZE == 0) {
				tmp = realloc(halo->ipart,
				              (halo->npart+partIn+BUFFERSIZE)
				              *sizeof(unsigned long));
				if (tmp == NULL) {
					fprintf(stderr,
					        "Not enough memory to allocate ipart array "
					        "in %s.  Aborting.\n\n", __func__);
					common_terminate(EXIT_FAILURE);
				} else {
					halo->ipart = tmp;
				}
			}
			halo->ipart[halo->npart+partIn] = (unsigned long)offset;
			partIn++;
		} else {
			/* Do not add particle */
			partOut++;
		}
		partTot++;
		offset++;
	}

	/* Do maintenance work on the ipart array */
	if (partIn > 0) {
		halo->npart += partIn;
		halo->ipart = realloc(halo->ipart,
		                      (halo->npart*sizeof(unsigned long)));
	}

	/* Done */
	return;
}

static uint64_t
local_findMinOffset(partptr fstPart,
                    uint64_t numParts,
                    sfc_key_t minKey)
{
	uint64_t minOffset = UINT64_C(0);
	part keyPart; 
	partptr result;
	sfc_key_t keyFound;

	/* Find the a particle with this SFC key */
	keyPart.sfckey = minKey;
	(void)util_bsearch((void *)(&keyPart),
	                   (void *)fstPart,
	                   (size_t)numParts,
	                   sizeof(part),
	                   &cmp_sfckey_part,
	                   (void *)(&result));
	minOffset = (uint64_t)(result-fstPart);
	keyFound  = (global_info.fst_part+minOffset)->sfckey;

	/* This catches the situation that a key is searched which is larger
	 * than any particle key.  The result particle returned above would
	 * then be the last particle, but we need to return an offset that
	 * is larger than this. */
	if (minOffset == numParts - 1 && keyFound < minKey)
		return numParts;

	/* Verify the assumptions (only done for DEBUG) */
	assert(keyFound >= minKey);
	assert(minOffset >= 0 && minOffset < numParts);

	/* Make sure that we return the offset to the first particle with
	 * this key */
	while (    (minOffset>0)
	        && ((global_info.fst_part+minOffset-1)->sfckey == keyFound) )
		minOffset--;

	/* Verify that we return the proper thing (only done for DEBUG) */
	assert((global_info.fst_part+minOffset)->sfckey >= minKey);
	assert(minOffset >= 0 && minOffset < numParts);
	
	return minOffset;
}

#endif // AHF2

