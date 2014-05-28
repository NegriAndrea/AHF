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
#define DAUGHTER_DISTANCE_FRACTION    0.5

//*********************************************************************************************************
// PROTOTYPES
//*********************************************************************************************************
double parent_patch_radius      (patch_t patch);
double sub_patch_radius         (int, int, patch_t **, int *);
void   set_daughter_patch_radii (patch_t );
void   write_debuginfo          (patch_t, double, double, double, double );


//*********************************************************************************************************
// set_patch_radii():
//    for each patch on each level calculate the radius of the pathc to be later used
//    as an initial guess for the radius of the halo
//*********************************************************************************************************
void set_patch_radii(patch_t **patch_tree, int64_t *n_patches)
{
  int      initial_level, final_level, ilevel, iparent;
  patch_t *parent_patch;
  
  // determine initial and final level
  get_patch_level_range(patch_tree, n_patches, &initial_level, &final_level);

  // use initial_level as starting point
  ilevel = initial_level;
  
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) private(iparent,parent_patch) shared(patch_tree,n_patches,ilevel,stderr) default(none)
#endif
  for(iparent=0; iparent<n_patches[ilevel]; iparent++) {
    
    // set parent patch for easier access
    parent_patch = &(patch_tree[ilevel][iparent]);
    
    //---------------------------------------------------------------------------------------------------
    // parent patch => simply set radius to half of the diagonal of the patch-rectangular
    //---------------------------------------------------------------------------------------------------
    parent_patch->radius = parent_patch_radius(*parent_patch);

    //---------------------------------------------------------------------------------------------------
    // daughter patches => recursively set radii to (fraction of) distance to nearest companion daughter
    //---------------------------------------------------------------------------------------------------
    set_daughter_patch_radii(*parent_patch);

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
  
#ifdef DEBUG_AHF2libtree
 {
  FILE *fp;
  char parentname[2048], daughtername[2048];
  
  fp = fopen("radii.geom","a");
  
  if(parent_patch.n_daughters > 1) {
    sprintf(parentname,"parent-%05"PRIi32"-%05"PRIi64".geom",parent_patch.level,parent_patch.id);
    write_patch_geom(parentname,parent_patch);
    fprintf(fp,"s %lf %lf %lf %lf 0 0 0\n",parent_patch.pos[0],parent_patch.pos[1],parent_patch.pos[2],parent_patch.radius);
    for(idaughter=0; idaughter<parent_patch.n_daughters; idaughter++) {
      idaughter_patch = parent_patch.daughters[idaughter];
      sprintf(daughtername,"parent-%05"PRIi32"-%05"PRIi64"_daughter-%05"PRIi32"-%05"PRIi64".geom",parent_patch.level,parent_patch.id,idaughter_patch->level,idaughter_patch->id);
      write_patch_geom(daughtername, *idaughter_patch);
      fprintf(fp,"s %lf %lf %lf %lf 0 1 0\n",idaughter_patch->pos[0],idaughter_patch->pos[1],idaughter_patch->pos[2],idaughter_patch->radius);
    }
  }
  
  fclose(fp);
 }
#endif
  

}

//*********************************************************************************************************
// parent_patch_radius():
//    for an individual patch set the radius to maximum distance of the centre to one of the corners
//*********************************************************************************************************
#define CHECK_DISTANCE {\
if(dx > 0.5) dx=1.0-dx;\
if(dy > 0.5) dy=1.0-dy;\
if(dz > 0.5) dz=1.0-dz;\
dist2 = pow2(dx)+pow2(dy)+pow2(dz); \
if(dist2 > radius2) radius2 = dist2; \
}
double parent_patch_radius(patch_t patch)
{
  double  x,  y,  z;
  double dx, dy, dz, dist2, radius2;
  
  // patch centre
  get_patch_centre(patch, &x, &y, &z);
  
  // keep track of the maximum distance
  radius2 = -1.0;
  
  // distance to all the 8 corners of the rectangular encompassing the patch
  dx = fabs(x-patch.Xmin);
  dy = fabs(y-patch.Ymin);
  dz = fabs(z-patch.Zmin);
  CHECK_DISTANCE;
  
  dx = fabs(x-patch.Xmax);
  dy = fabs(y-patch.Ymin);
  dz = fabs(z-patch.Zmin);
  CHECK_DISTANCE;
  
  dx = fabs(x-patch.Xmin);
  dy = fabs(y-patch.Ymax);
  dz = fabs(z-patch.Zmin);
  CHECK_DISTANCE;
  
  dx = fabs(x-patch.Xmin);
  dy = fabs(y-patch.Ymin);
  dz = fabs(z-patch.Zmax);
  CHECK_DISTANCE;
  
  dx = fabs(x-patch.Xmin);
  dy = fabs(y-patch.Ymax);
  dz = fabs(z-patch.Zmax);
  CHECK_DISTANCE;
  
  dx = fabs(x-patch.Xmax);
  dy = fabs(y-patch.Ymin);
  dz = fabs(z-patch.Zmax);
  CHECK_DISTANCE;
  
  dx = fabs(x-patch.Xmax);
  dy = fabs(y-patch.Ymax);
  dz = fabs(z-patch.Zmin);
  CHECK_DISTANCE;
  
  dx = fabs(x-patch.Xmax);
  dy = fabs(y-patch.Ymax);
  dz = fabs(z-patch.Zmax);
  CHECK_DISTANCE;
  
  // diagonal of patch rectangle
#ifdef AHF2_hostradius_is_patchdiagonal
  dx = fabs(patch.Xmax-patch.Xmin);
  dy = fabs(patch.Ymax-patch.Ymin);
  dz = fabs(patch.Zmax-patch.Zmin);
  if(dx > 0.5) dx=1.0-dx;
  if(dy > 0.5) dy=1.0-dy;
  if(dz > 0.5) dz=1.0-dz;
  radius2 = (pow2(dx)+pow2(dy)+pow2(dz))/4.0;
#endif
  
#ifdef DEBUG_AHF2libtree
  if(patch.periodic[0] == 0) {
    if(x < patch.Xmin || x > patch.Xmax) {
      write_debuginfo(patch,x,y,z,radius2);
    }
  }
  
  if(patch.periodic[0] == 1) {
    if(x > patch.Xmin && x < patch.Xmax) {
      write_debuginfo(patch,x,y,z,radius2);
    }
  }

  if(patch.periodic[1] == 0) {
    if(y < patch.Ymin || y > patch.Ymax) {
      write_debuginfo(patch,x,y,z,radius2);
    }
  }
  
  if(patch.periodic[1] == 1) {
    if(y > patch.Ymin && y < patch.Ymax) {
      write_debuginfo(patch,x,y,z,radius2);
    }
  }

  if(patch.periodic[2] == 0) {
    if(z < patch.Zmin || z > patch.Zmax) {
      write_debuginfo(patch,x,y,z,radius2);
    }
  }
  
  if(patch.periodic[2] == 1) {
    if(z > patch.Zmin && z < patch.Zmax) {
      write_debuginfo(patch,x,y,z,radius2);
    }
  }
#endif

#ifdef DEBUG_AHF2libtree_x
  if(sqrt(radius2) > 0.0075) {
    write_debuginfo(patch,x,y,z,radius2);
  }
#endif
  
  return(sqrt(radius2));
}
#undef CHECK_DISTANCE





#ifdef DEBUG_AHF2libtree
//*********************************************************************************************************
// write_debuginfo()
//*********************************************************************************************************

void write_debuginfo(patch_t patch, double x, double y, double z, double radius2)
{
  // write some statistics about the centres
  fprintf(stderr,"periodic[0]=%"PRIi8" Xmin=%lf x=%lf Xmax=%lf radius=%lf Npart=%"PRIu64" Npart_stored=%"PRIu64"\n",
          patch.periodic[0],
          patch.Xmin,x,patch.Xmax,sqrt(radius2),
          patch.Npart,patch.Npart_stored);
  fprintf(stderr,"periodic[1]=%"PRIi8" Ymin=%lf y=%lf Ymax=%lf\n",
          patch.periodic[1],
          patch.Ymin,y,patch.Ymax);
  fprintf(stderr,"periodic[2]=%"PRIi8" Zmin=%lf z=%lf Zmax=%lf\n",
          patch.periodic[2],
          patch.Zmin,z,patch.Zmax);
  fprintf(stderr," patch.pos[0]=%lf, centre.cube.geom[0]=%lf centre.cube.wgeom[0]=%lf centre.part.com[0]=%lf centre.part.Nmax[0]=%lf\n",
          patch.pos[0],patch.centre.cube.geom[0],patch.centre.cube.wgeom[0],patch.centre.part.com[0],patch.centre.part.Nmax[0]);
  fprintf(stderr," patch.pos[1]=%lf, centre.cube.geom[1]=%lf centre.cube.wgeom[1]=%lf centre.part.com[1]=%lf centre.part.Nmax[1]=%lf\n",
          patch.pos[1],patch.centre.cube.geom[1],patch.centre.cube.wgeom[1],patch.centre.part.com[1],patch.centre.part.Nmax[1]);
  fprintf(stderr," patch.pos[2]=%lf, centre.cube.geom[2]=%lf centre.cube.wgeom[2]=%lf centre.part.com[2]=%lf centre.part.Nmax[2]=%lf\n",
          patch.pos[2],patch.centre.cube.geom[2],patch.centre.cube.wgeom[2],patch.centre.part.com[2],patch.centre.part.Nmax[2]);
  
  // and write the whole patch causing the trouble (cf. utilities.c)
  write_patch_geom("patch.geom",patch);
  
  exit(0);
}
#endif

#endif // NEWAMR