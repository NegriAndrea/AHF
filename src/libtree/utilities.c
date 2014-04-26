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

#ifdef NEWAMR 

//*********************************************************************************************************
// patch_property():
//    return some patch property that might be used to classify patches
//*********************************************************************************************************
double patch_property(patch_t patch)
{
  // mass comparisoin
  return((double)patch.Npart);
  
  // mean density comparison
  //return((double)patch.Npart/(double)patch.n_subcubes);
  
  // volume comparisoin
  return((double)patch.n_subcubes);
  
}

#endif