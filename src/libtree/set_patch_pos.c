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



//*********************************************************************************************************
// set_patch_pos():
//    decide which of the various centre definitions to eventually use as patch_tree[][].pos[]
//*********************************************************************************************************
void set_patch_pos(patch_t **patch_tree, int64_t *n_patches)
{
  int initial_level, final_level, ilevel, ipatch;
  
  // determine initial and final level
  get_patch_level_range(patch_tree, n_patches, &initial_level, &final_level);
  
  // loop over all patches on all levels
  for(ilevel=initial_level; ilevel<final_level; ilevel++) {
    
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) private(ipatch) shared(patch_tree,n_patches,ilevel,stderr) default(none)
#endif
    for(ipatch=0; ipatch<n_patches[ilevel]; ipatch++) {
    
#ifdef AHFcomcentre
      patch_tree[ilevel][ipatch].pos[0] = patch_tree[ilevel][ipatch].centre.part.com[0];
      patch_tree[ilevel][ipatch].pos[1] = patch_tree[ilevel][ipatch].centre.part.com[1];
      patch_tree[ilevel][ipatch].pos[2] = patch_tree[ilevel][ipatch].centre.part.com[2];
#endif
      
      // TODO: implement possible alternative centre definitions for patch
    
    } // ipatch
  } // ilevel  
}



#endif // NEWAMR