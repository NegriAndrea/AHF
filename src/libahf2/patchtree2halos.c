#ifdef AHF2

/***************************************************************************
 *   Includes                                                              *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <inttypes.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

#include "../libutility/utility.h"

#include "../libtree/tree.h"
#include "../libtree/utilities.h"


//*********************************************************************************************************
// PROTOTYPES
//*********************************************************************************************************
void generate_daughter_halos(HALO *halos, patch_t parent_patch, uint64_t ihost, uint64_t *ihalo);


//*********************************************************************************************************
// patchtree2halos()
//*********************************************************************************************************
halo_s* patchtree2halos (patch_t **patch_tree, uint64_t *n_patches)
{
  int      initial_level, final_level, ilevel;
  uint64_t ipatch, iparent;
  HALO    *halos;
  uint64_t ihalo, nhalos, ihost;
  halo_s  *return_values;
  patch_t  parent_patch;
  double   x, y, z;
  
  // determine initial and final level
  get_patch_level_range(patch_tree, n_patches, &initial_level, &final_level);
  
  //======================================================================================
  // first guess for number of haloes = all the end-leaves
  //======================================================================================
  nhalos = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for default(none) private(ilevel,ipatch) shared(patch_tree,n_patches,initial_level,final_level) schedule(dynamic) reduction( + : nhalos )
#endif
  for(ilevel=initial_level; ilevel<final_level; ilevel++) {
    for(ipatch=0; ipatch<n_patches[ilevel]; ipatch++) {
      if(patch_tree[ilevel][ipatch].n_daughters == 0) {
        nhalos++;
      } // if()
    } // ipatch
  } // ilevel

#ifdef VERBOSE
  fprintf(stderr,"  patchtre2halos: number of end-leaves in patch_tree[][] -> nhalos=%"PRIu64"\n",nhalos);
#endif
  fprintf(io.logfile,"  patchtre2halos: number of end-leaves in patch_tree[][] -> nhalos=%"PRIu64"\n",nhalos);
  fflush(io.logfile);
  
  //======================================================================================
  // allocate memory to store all the prospective halos
  //======================================================================================
  halos = (HALO *) calloc(nhalos, sizeof(HALO));
  if(halos == NULL) {
#ifdef VERBOSE
    fprintf(stderr,"  patchtre2halos: cannot allocate memory for halos[] array\n");
#endif
    fprintf(io.logfile,"  patchtre2halos: cannot allocate memory for halos[] array\n");
    fflush(io.logfile);
    exit(0);
  }
  

  //======================================================================================
  // pick each parent patch in initial_level and treat it as a halo
  //======================================================================================
  // use initial_level as starting point
  ilevel = initial_level;

  // accumulate haloes
  ihalo = 0;

  // this loop is intrinsically serial due to the use of 'ihalo' as the counter for the created halos
  for(iparent=0; iparent<n_patches[ilevel]; iparent++) {
    
    // set parent patch for easier access
    parent_patch = (patch_tree[ilevel][iparent]);
    
    // TODO: deal with possible mergers on the highest level
    if(parent_patch.trunk == NO_SINGLE_TRUNK) {
    }
  
    // this walks down the trunk to get the most refined centre
    get_patch_centre(parent_patch, &x, &y, &z);
    
#ifdef DEBUG_AHF2
    fprintf(stderr,"creating parent halo %"PRIu64" from patch-%"PRIi32"-%"PRIi64"\n",ihalo,parent_patch.level,parent_patch.id);
#endif
    
    // copy properties over to halos[] array
    halos[ihalo].npart         = parent_patch.Npart;
    halos[ihalo].numNodes      = parent_patch.n_subcubes;
    halos[ihalo].refLev        = parent_patch.level;
    halos[ihalo].spaRes        = 1./(pow(2.0,parent_patch.level));
    halos[ihalo].pos.x         = x;
    halos[ihalo].pos.y         = y;
    halos[ihalo].pos.z         = z;
    halos[ihalo].R_vir         = MIN(parent_patch.radius, simu.MaxGatherRad/simu.boxsize);
    halos[ihalo].gatherRad     = MIN(parent_patch.radius, simu.MaxGatherRad/simu.boxsize);
    
    halos[ihalo].hostHalo      = -1;
    halos[ihalo].hostHaloLevel = parent_patch.level;
    halos[ihalo].numSubStruct  = 0;    // will be calculated as we go along
    halos[ihalo].subStruct     = NULL;
    
    // keep track of the host halo
    ihost = ihalo;
    
    // generate daughter halos[]
    generate_daughter_halos(halos, parent_patch, ihost, &ihalo);
    
    // next parent halo
    ihalo++;
  }
  
#ifdef VERBOSE
  fprintf(stderr,"  patchtree2halos: filled halos[] with %"PRIu64" halos (nhalos=%"PRIu64")\n",ihalo,nhalos);
#endif
  fprintf(io.logfile,"  patchtree2halos: filled halos[] with %"PRIu64" halos (nhalos=%"PRIu64")\n",ihalo,nhalos);
  fflush(io.logfile);
  
  //======================================================================================
  // put everything into a halo_s container to be returned
  //======================================================================================
  return_values         = (halo_s *) calloc(1, sizeof(halo_s));
  return_values->halos  = halos;
  return_values->nhalos = nhalos;
  return(return_values);
}


//*********************************************************************************************************
// generate_daughter_halos()
//*********************************************************************************************************
void generate_daughter_halos(HALO *halos, patch_t parent_patch, uint64_t ihost, uint64_t *ihalo)
{
  uint64_t inewhost;
  int32_t  idaughter;
  double   x, y, z;
  patch_t  daughter_patch;
  
  // loop over all daughters
  for(idaughter=0; idaughter<parent_patch.n_daughters; idaughter++) {
    
    // direct access to the daughter patch
    daughter_patch = *(parent_patch.daughters[idaughter]);
    
    if(parent_patch.trunk == NO_SINGLE_TRUNK) {
#ifdef DEBUG_AHF2
      fprintf(stderr,"generate_daughter_halos(): untreated merger for ihalo=%"PRIu64" ihost=%"PRIu64" trunk=%d\n",*ihalo,ihost,parent_patch.trunk);
#endif
      parent_patch.trunk = 0;  // TODO: sub-mergers not treated correctly yet
    }
    
    // do not create the parent halo again
    if(idaughter != parent_patch.trunk) {
      
      // we are inserting halos[] even though they are "just" the daughters
      (*ihalo)++;
      
#ifdef DEBUG_AHF2
      fprintf(stderr,"creating daughter halo %"PRIu64" from patch-%"PRIi32"-%"PRIi64" (whose parent is patch-%"PRIi32"-%"PRIi64")\n",
              *ihalo,daughter_patch.level,daughter_patch.id,parent_patch.level,parent_patch.id);
#endif
      
      // add a pointer to the host
      halos[ihost].numSubStruct++;
      halos[ihost].subStruct = (int *) realloc(halos[ihost].subStruct, halos[ihost].numSubStruct*sizeof(int));
      if(halos[ihost].subStruct == NULL) {
        fprintf(stderr,"generate_daughter_halos: could not realloc memory for halos[].subStruct\n");
        exit(-1);
      }
      halos[ihost].subStruct[halos[ihost].numSubStruct-1] = (int) (*ihalo);
      
#ifdef DEBUG_AHF2
      fprintf(stderr,"   halo %"PRIu64" is now substructure of ihost=%"PRIu64": halos[%"PRIu64"].numSubStruct=%d \n",
              *ihalo,ihost,ihost,halos[ihost].numSubStruct);
#endif
      
      // this walks down the trunk to get the most refined centre
      get_patch_centre(daughter_patch, &x, &y, &z);
      
      // copy properties over to halos[] array
      halos[*ihalo].npart          = daughter_patch.Npart;
      halos[*ihalo].numNodes       = daughter_patch.n_subcubes;
      halos[*ihalo].refLev         = daughter_patch.level;
      halos[*ihalo].spaRes         = 1./(pow(2.0,daughter_patch.level));
      halos[*ihalo].pos.x          = x;
      halos[*ihalo].pos.y          = y;
      halos[*ihalo].pos.z          = z;
      halos[*ihalo].R_vir          = MIN(daughter_patch.radius, simu.MaxGatherRad/simu.boxsize);
      halos[*ihalo].gatherRad      = MIN(daughter_patch.radius, simu.MaxGatherRad/simu.boxsize);
      
      halos[*ihalo].hostHalo       = (int)ihost;
      halos[*ihalo].hostHaloLevel  = parent_patch.level;
      halos[*ihalo].numSubStruct   = 0;    // will be calculated as we go along
      halos[*ihalo].subStruct      = NULL;
      
    } // if(trunk)
    
    // simply follow the trunk, i.e. we do not need to do anything special here
    else {
#ifdef DEBUG_AHF2
      fprintf(stderr,"following the trunk of halo %"PRIu64" from patch-%"PRIi32"-%"PRIi64" (whose parent is patch-%"PRIi32"-%"PRIi64")\n",
              *ihalo,daughter_patch.level,daughter_patch.id,parent_patch.level,parent_patch.id);
#endif
    }

    // treat daughter as the host to new subhaloes
    inewhost = *ihalo;
    
    // generate daughter halos[] of this daughter
    generate_daughter_halos(halos, daughter_patch, inewhost, ihalo);
    
  } // for(idaughter)
}

#endif // AHF2

