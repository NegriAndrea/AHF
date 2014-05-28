#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"
#include "../libsfc/sfc.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libutility/specific.h"

#include "cubekey.h"
#include "patch.h"
#include "utilities.h"

#ifdef NEWAMR

//*********************************************************************************************************
// DEFINES
//*********************************************************************************************************


//*********************************************************************************************************
// PROTOTYPES
//*********************************************************************************************************
int32_t get_patch_trunk         (patch_t);
int     qcomparePatchNpart      (const void *, const void *);
int     check_merger            (patch_t *, patch_t *);

//*********************************************************************************************************
// set_patch_trunk():
//    assign a trunk to each patch from the list of daughters
//*********************************************************************************************************
void set_patch_trunk(patch_t **patch_tree, int64_t *n_patches)
{
  int32_t initial_level, final_level, ilevel;
  int64_t ipatch;
 
  // determine initial and final level
  get_patch_level_range(patch_tree, n_patches, &initial_level, &final_level);
  
  // it makes more sense to parallelise over the patches on each level than over the levels
  for(ilevel=initial_level; ilevel<final_level; ilevel++) {
    
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) private(ipatch) shared(patch_tree,n_patches,ilevel,stderr) default(none)
#endif
    for(ipatch=0; ipatch<n_patches[ilevel]; ipatch++) {
      
      patch_tree[ilevel][ipatch].trunk = get_patch_trunk(patch_tree[ilevel][ipatch]);
      
    } // ipatch
  } // ilevel
}

//*********************************************************************************************************
// get_patch_trunk():
//    assign a trunk to each patch from the list of daughters
//    (currently only PARDAU_PARTS is supported as this proved to be the choice for AHF1)
//*********************************************************************************************************
int32_t get_patch_trunk(patch_t patch)
{
  int       idaughter, merger;
  int32_t   trunk;
  uint64_t  Npart_daughter, max_Npart_daughter;
  patch_t  *daughter_patch;
  
  max_Npart_daughter = 0;
  trunk              = END_OF_TRUNK;
  
  // mulitple daughters => find trunk
  //----------------------------------
  if(patch.n_daughters > 1) {
    
    // not considered merger yet
    merger = FALSE;

    // sort daughters by Npart
    qsort((void *)patch.daughters, patch.n_daughters, sizeof(patch_t *), qcomparePatchNpart);
   
    // check the two largest daughters for a possible merger!
    merger = check_merger(patch.daughters[0],patch.daughters[1]);
    
    if(merger) {
      trunk = NO_SINGLE_TRUNK;
    }
    else {
#ifdef PARDAU_PARTS

      // we just sorted the daughters by Npart and hence the first daughter is the trunk
      trunk = 0;

#else // PARDAU_PARTS
      
      // loop over all daughters checking for whatever criterion
      for(idaughter=0; idaughter<patch.n_daughters; idaughter++) {
        
        daughter_patch = patch.daughters[idaughter];
        Npart_daughter = daughter_patch->Npart;
        
#ifdef PARDAU_DISTANCE   // base trunk on closest in distance
      // to be implemented...
#endif
        
#ifdef PARDAU_NODES     // base trunk on maximum number of cubes
      // to be implemented...
#endif
      } // for(idaughter)
      
#endif // PARDAU_PARTS
    } // if-else(merger)
    
  } // if(n_daughters > 1)
  
  // one daughter => that's the trunk then
  //---------------------------------------
  else if (patch.n_daughters == 1) {
    trunk = 0;
  }
  
  // no daughters => end of trunk
  //------------------------------
  else {
    //trunk = END_OF_TRUNK;  // already assigned to this above
  }
  
  return(trunk);
}


//*********************************************************************************************************
// qcomparePatchNpart():
//    compare patch Npart (used with qsort)
//*********************************************************************************************************
int qcomparePatchNpart(const void *patch1, const void *patch2)
{
	uint64_t n1, n2;
  
	n1 = ((patch_t *)patch1)->Npart;
	n2 = ((patch_t *)patch2)->Npart;
  
	return n1 < n2 ? -1 : (n1 > n2 ? 1 : 0);
}

//*********************************************************************************************************
// check_merger():
//    compare two patches whether they should be considered as an ongoing merger
//*********************************************************************************************************
int check_merger(patch_t *p1, patch_t *p2)
{
  if(fabs(p1->Npart      - p2->Npart)     /(double)p1->Npart      > (1-MERGER_NPART_FRAC) &&
     fabs(p1->n_subcubes - p2->n_subcubes)/(double)p1->n_subcubes > (1-MERGER_VOL_FRAC)      ) {
    return(TRUE);
  }
  else {
    return(FALSE);
  }
}

#endif // NEWAMR
