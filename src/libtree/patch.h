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
  int64_t    patch_id;                 // allow for largest possible number (reserving highest bit for flagging things)
  int8_t     patch_level;              // maximum = MAX_DEPTH


  // this is the patch tree
  struct     patch_s*  parent;     // there should only be one parent (otherwise the tree/simulation is screwed)
  uint32_t   n_daughters;          // 2^32 should be sufficient for daughters (any more casts doubt on the simulation)
  struct     patch_s** daughters;  // there are n_daughters

  int32_t    trunk;                // needed to follow the trunk of each patch
  
  
  // stuff needed by AK
  //====================
  //uint64_t   ncubes; // this is n_subcubes defined above amongst the FC stuff
  uint64_t   Npart;    // total number of particles on patch
  uint32_t   Nmax;     // number of particles in subcube with the most particles
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

#define SIZEOF_PATCH (sizeof(patch_t))

//Patch manipulation
//pointers to subcubes array manipulation
void patch_add_psubcube(patch_t*, psubcube_t*);

//psubcube_t* patch_pop_back(patch_t* patch);
psubcube_t* patch_next_psubcube(patch_t*, psubcube_t*);

void patch_create(patch_t**, int, int);
void patch_init(patch_t*, int, int);
void patch_clear(patch_t*);
int  patch_get_num_psubcubes(patch_t*);
void patch_free_psubcubes_array(patch_t*);

// routines added by AK
void init_patch_physics(patch_t *);
void add_patch_physics(patch_t *, psubcube_t );
void finish_patch_physics(patch_t *);

int patch_include_adjacent_subcubes(subcube_t*, subcube_t*, patch_t*);
#endif

