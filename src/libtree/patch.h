#ifndef _PATCH_H_
#define _PATCH_H_

//#include <math.h>
#include <stdint.h>
//#include "../common.h"
#include "subcube.h"
//#include "uthash.h"
#include "utarray.h"

typedef struct patch_s{
  // this is stuff needed by FC
  UT_array  *psubcubes;                //pointers to subcubes contained in the patch.
  uint64_t   n_subcubes;
  psubcube_t most_particles_subcube;
  int64_t    id;                 // allow for largest possible number (reserving highest bit for flagging things)
  int32_t    level;              // maximum = MAX_DEPTH


  // stuff needed by AK
  //====================

  // this is the patch tree
  struct     patch_s*  parent;     // there should only be one parent (otherwise the tree/simulation is screwed)
  uint32_t   n_daughters;          // 2^32 should be sufficient for daughters (any more casts doubt on the simulation)
  struct     patch_s** daughters;  // there are n_daughters

//#ifdef AHF2_read_spatialRef
  uint64_t  *idaughter;
//#endif
  
  int32_t    trunk;                // needed to follow the trunk of each patch
  
  int8_t     periodic[3];          // records whether the patch wraps around the periodic boundaries
  int8_t     periodic_shift[3];    // periodic shift to be applied (either +1 or -1, coordinates are in [0,1]!)
  
  //uint64_t   ncubes;      // this is n_subcubes defined above amongst the FC stuff
  uint64_t   Npart;         // total number of particles on patch
  uint64_t   Npart_stored;  // the actual number of particles we have access to (if there is a sub-patch, there will be a "particle hole")
  uint32_t   Nmax;          // number of particles in subcube with the most particles
  double     Xmin,Xmax,Ymin,Ymax,Zmin,Zmax; // dimensions of rectangle fully containing patch
  
  struct {
    
    struct {
      double  geom[3]; //          geometrical centre
      double wgeom[3]; // weighted geometrical centre
    } cube;
    
    struct {
      double com[3];  // centre-of-mass (ignoring particle masses though!)
      double Nmax[3]; // position of cube with most particles
    } part;
    
  } centre; // we calculate various possible centres, but only one will eventually be used (cf. pos[3] below)
  
  double     pos[3];  // this will be the centre of the patch to be used throughout the code
  double     radius;  // this will/might be used as the initial guess for the virial radius of the halo
  //====================
  // stuff needed by AK

  
  
} patch_t;

//ahf2_patches_t is the data structure storing the patches tree and the number of patches at each level
typedef struct ahf2_patches {
  patch_t    *tree[NLEVELS];
  int64_t    n_patches[NLEVELS];
  int64_t    n_rejected_patches[NLEVELS];
} ahf2_patches_t;

#define SIZEOF_PATCH (sizeof(patch_t))
#define SIZEOF_AHF2_PATCHES (sizeof(ahf2_patches_t))


/*
 * visitable_subcubes_t is the data structure that will store all the subcubes that should be
 * visited in the next iteration of a patch creation
 */
typedef struct visitable_subcubes {
  cubekey_t *list;
  int64_t next;
  int64_t size;
} visitable_subcubes_t;

//Block size for reallocating list of cubekeys. We want to avoid a call to realloc each time we add a new cubekey.
//Larger value will probably waste some memory. Smaller value will produce more realloc calls.
#define VISITABLE_SUBCUBES_BLOCK 10

void visitable_subcubes_create(visitable_subcubes_t**);
void visitable_subcubes_clear(visitable_subcubes_t*);
void visitable_subcubes_add(visitable_subcubes_t*, cubekey_t*);
int64_t visitable_subcubes_count(visitable_subcubes_t*);
cubekey_t* visitable_subcubes_get(visitable_subcubes_t*, int64_t);
void visitable_subcubes_free(visitable_subcubes_t**);

//Patch manipulation
//pointers to subcubes array manipulation
void patch_add_psubcube(patch_t*, psubcube_t*);

//psubcube_t* patch_pop_back(patch_t* patch);
psubcube_t* patch_next_psubcube(patch_t*, psubcube_t*);

void patch_create(patch_t**, int64_t, int32_t);
void patch_init(patch_t*, int64_t, int32_t);
void patch_clear(patch_t*);
void patch_free(patch_t**);
void patch_free_content(patch_t*);
uint64_t  patch_get_num_psubcubes(patch_t*);
void ahf2_patches_free(ahf2_patches_t**);
void patch_free_psubcubes_array(patch_t*);
void patch_connect_tree(patch_t**,subcube_t**, int64_t*);
int  patch_connect_my_parent(patch_t**, subcube_t**, patch_t*);
void patch_connection_review(ahf2_patches_t*);
void patch_formation_review(ahf2_patches_t*);
void patches_generation(ahf2_patches_t**, subcube_t**, int32_t, int32_t, uint64_t);
int64_t patches_generation_per_level(ahf2_patches_t*, subcube_t*, int32_t, uint64_t);

// routines added by AK
void init_patch_physics(patch_t *);
void add_patch_physics(patch_t *, psubcube_t );
void finish_patch_physics(patch_t *);

int patch_include_adjacent_subcubes(subcube_t*, subcube_t*, patch_t*);
#endif

