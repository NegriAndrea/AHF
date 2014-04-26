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
#define DAUGHTER_DISTANCE_FRACTION    0.5

//*********************************************************************************************************
// PROTOTYPES
//*********************************************************************************************************
double parent_patch_radius      (patch_t patch);
double sub_patch_radius         (int, int, patch_t **, int *);
double patch_property           (patch_t );
void   get_patch_centre         (patch_t , double *, double *, double *);
void   set_daughter_patch_radii (patch_t );


//*********************************************************************************************************
// set_patch_radii():
//    for each patch on each level calculate the radius of the pathc to be later used
//    as an initial guess for the radius of the halo
//*********************************************************************************************************
void set_patch_radii(patch_t **patch_tree, int *n_patches, int initial_depth)
{
  int     ilevel, iparent;
  patch_t parent_patch;
  
  ilevel = initial_depth;
  
#pragma omp parallel for schedule(dynamic) private(iparent,parent_patch) shared(patch_tree,n_patches,ilevel) default(none)
  for(iparent=0; iparent<n_patches[ilevel]; iparent++) {
    
    // set parent patch for easier access
    parent_patch = patch_tree[ilevel][iparent];
    
    //---------------------------------------------------------------------------------------------------
    // parent patch => simply set radius to half of the diagonal of the patch-rectangular
    //---------------------------------------------------------------------------------------------------
    parent_patch.radius = parent_patch_radius(parent_patch);
  
    //---------------------------------------------------------------------------------------------------
    // daughter patches => recursively set radii to (fraction of) distance to nearest companion daughter
    //---------------------------------------------------------------------------------------------------
    set_daughter_patch_radii(parent_patch);
    
  } // iparent
  
}


//*********************************************************************************************************
// set_daughter_patch_radii():
//    for a given parent, loop over all daughters setting their radii
//*********************************************************************************************************
void set_daughter_patch_radii(patch_t parent_patch)
{
  patch_t *idaughter_patch, *jdaughter_patch;
  int      idaughter, jdaughter;
  double   xi, yi, zi, xj, yj, zj;
  double   dx, dy, dz, radius2, min_radius2;
  
  for(idaughter=0; idaughter<parent_patch.n_daughters; idaughter++) {
    
    // direct access to idaughter patch
    idaughter_patch = parent_patch.daughters[idaughter];
    
    // ignore the trunk (which is amongst the list of daughters, too)
    if(idaughter != parent_patch.trunk) {
      
      // get most credible centre for this patch (following its trunk)
      get_patch_centre(*idaughter_patch, &xi, &yi, &zi);
      
      // find minimum distance to any other daughter
      min_radius2 = pow2(2.0);
      
      // loop over all other daughters
      for(jdaughter=0; jdaughter<parent_patch.n_daughters; jdaughter++) {
        
        // exclude idaughter and trunk from this loop
        if(jdaughter != idaughter) {
          
          // maybe also exclude the trunk?
          //if(jdaughter != parent_patch.trunk)
         {
          // direct access to idaughter patch
          jdaughter_patch = parent_patch.daughters[jdaughter];
          
          // get most credible centre for this patch (following its trunk)
          get_patch_centre(*jdaughter_patch, &xj, &yj, &zj);
          
          // distance (periodic boundaries)
          dx = fabs(xi-xj);
          dy = fabs(yi-yj);
          dz = fabs(zi-zj);
          if(dx > 0.5) dx=1.0-dx;
          if(dy > 0.5) dy=1.0-dy;
          if(dz > 0.5) dz=1.0-dz;
          
          // record minimum distance as radius
          radius2 = pow2(dx)+pow2(dy)+pow2(dz);
          
          if(radius2 < min_radius2) min_radius2 = radius2;
          
         } // if(jdaughter != trunk)
          
        } // if(jdaughter != idaughter)
        
      } // jdaughter
      
      // set radius to DAUGHTER_DISTANCE_FRACTION times the distance to the nearest other daughter
      idaughter_patch->radius = DAUGHTER_DISTANCE_FRACTION * sqrt(min_radius2);
    }
    // if(idaughter == trunk) => inherit radius from parent
    else {
      idaughter_patch->radius = parent_patch.radius;
    }
    
    // treat every daughter now as the parent and assign radii to all its daughters
    set_daughter_patch_radii(*idaughter_patch);
    
  } // idaughter
}

//*********************************************************************************************************
// parent_patch_radius():
//    for an individual patch set the radius to maximum distance of the centre to one of the corners
//*********************************************************************************************************
double parent_patch_radius(patch_t patch)
{
  double  x,  y,  z;
  double dx, dy, dz, dist2, radius2;
  
  // patch centre
  get_patch_centre(patch, &x, &y, &z);
  
  // keep track of the maximum distance
  radius2 = -1.0;
  
  // distance to all the six corners of the rectangular encompassing the patch
  dx = fabs(x-patch.Xmin);
  dy = fabs(y-patch.Ymin);
  dz = fabs(z-patch.Zmin);
  if(dx > 0.5) dx=1.0-dx;
  if(pow2(dx)+pow2(dy)+pow2(dz) > radius2) radius2 = dist2;
  
  dx = fabs(x-patch.Xmax);
  dy = fabs(y-patch.Ymin);
  dz = fabs(z-patch.Zmin);
  if(dx > 0.5) dx=1.0-dx;
  if(pow2(dx)+pow2(dy)+pow2(dz) > radius2) radius2 = dist2;
  
  dx = fabs(x-patch.Xmin);
  dy = fabs(y-patch.Ymax);
  dz = fabs(z-patch.Zmin);
  if(dy > 0.5) dy=1.0-dy;
  if(pow2(dx)+pow2(dy)+pow2(dz) > radius2) radius2 = dist2;
  
  dx = fabs(x-patch.Xmin);
  dy = fabs(y-patch.Ymin);
  dz = fabs(z-patch.Zmax);
  if(dy > 0.5) dy=1.0-dy;
  if(pow2(dx)+pow2(dy)+pow2(dz) > radius2) radius2 = dist2;
  
  dx = fabs(x-patch.Xmin);
  dy = fabs(y-patch.Ymax);
  dz = fabs(z-patch.Zmax);
  if(dz > 0.5) dz=1.0-dz;
  if(pow2(dx)+pow2(dy)+pow2(dz) > radius2) radius2 = dist2;
  
  dx = fabs(x-patch.Xmax);
  dy = fabs(y-patch.Ymax);
  dz = fabs(z-patch.Zmax);
  if(dz > 0.5) dz=1.0-dz;
  if(pow2(dx)+pow2(dy)+pow2(dz) > radius2) radius2 = dist2;
  
  return(sqrt(radius2));
}

//*********************************************************************************************************
// get_patch_centre():
//    walk down the patch-trunk to locate most credible patch centre at the end of the trunk branch
//*********************************************************************************************************
void get_patch_centre(patch_t patch, double *x, double *y, double *z)
{
  int32_t  trunk;
  patch_t *cur_patch;
  
  // set current patch and its trunk
  trunk     = patch.trunk;
  cur_patch = &patch;
  
  // walk down the trunk branch of the patch
  while(trunk >= 0) {
    cur_patch = cur_patch->daughters[trunk];
    trunk     = cur_patch->trunk;
  }
  
  // set position
  *x = cur_patch->pos[0];
  *y = cur_patch->pos[1];
  *z = cur_patch->pos[2];
  
  // we should also deal with the situation NO_SINGLE_TRUNK
  //  if(trunk == NO_SINGLE_TRUNK && cur_patch->n_daughters > 0) {
  //    // we leave this to later when transferring patches to halos
  //  }
}


//*********************************************************************************************************
// patch_property():
//    we only take the distance to those patches into account in the radius calculation
//    that meet some criterion, e.g. that are more massive (Npart) or have higher mean density (Npart/Ncubes)
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

#endif // NEWAMR