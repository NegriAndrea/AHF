#include <stddef.h>
#include <stdio.h>
#include <math.h>
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
int32_t get_patch_trunk(patch_t);

//*********************************************************************************************************
// set_patch_trunk():
//    assign a trunk to each patch from the list of daughters
//*********************************************************************************************************
void set_patch_trunk(patch_t **patch_tree, int *n_patches, int initial_depth)
{
  int ilevel, ipatch;
 
  // it makes more sense to parallelise over the patches on each level than over the levels
  for(ilevel=initial_depth; ilevel<MAX_DEPTH; ilevel++) {
    
#pragma omp parallel for schedule(dynamic) private(ipatch) shared(patch_tree,n_patches,ilevel) default(none)
    for(ipatch=0; ipatch<n_patches[ilevel]; ipatch++) {
      
      patch_tree[ilevel][ipatch].trunk = get_patch_trunk(patch_tree[ilevel][ipatch]);
      
    } // ipatch
  } // ilevel
}

//*********************************************************************************************************
// set_patch_trunk():
//    assign a trunk to each patch from the list of daughters
//    (currently only PARDAU_PARTS is supported as this proved to be the choice for AHF1)
//*********************************************************************************************************
int32_t get_patch_trunk(patch_t patch)
{
  int       idaughter;
  int32_t   trunk;
  uint64_t  Npart_daughter, max_Npart_daughter;
  patch_t  *daughter_patch;
  
  max_Npart_daughter = 0;
  trunk              = END_OF_TRUNK;
  
  // mulitple daughters => find trunk
  //----------------------------------
  if(patch.n_daughters > 1) {
    
    // loop over all daughters
    for(idaughter=0; idaughter<patch.n_daughters; idaughter++) {
     
      daughter_patch = patch.daughters[idaughter];
      Npart_daughter = daughter_patch->Npart;
      
#ifdef PARDAU_PARTS     // base trunk on maximum number of particles
      if(Npart_daughter > max_Npart_daughter) {
        trunk              = idaughter;
        max_Npart_daughter = Npart_daughter;
      }
#endif
      
#ifdef PARDAU_DISTANCE   // base trunk on closest in distance
      // to be implemented...
#endif
      
#ifdef PARDAU_NODES     // base trunk on maximum number of cubes
      // to be implemented...
#endif
      
      
    }
  }
  
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

#endif // NEWAMR