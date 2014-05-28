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
  
  // TODO: we should also deal with the situation NO_SINGLE_TRUNK
  //  if(trunk == NO_SINGLE_TRUNK && cur_patch->n_daughters > 0) {
  //    // we leave this to later when transferring patches to halos
  //  }
}


//*********************************************************************************************************
// get_patch_level_range():
//    return some patch property that might be used to classify patches
//*********************************************************************************************************
void get_patch_level_range(patch_t **patch_tree, int64_t *n_patches, int *initial_level, int *final_level)
{
  // determine initial level
  *initial_level = 0;
  while(patch_tree[*initial_level] == NULL && *initial_level < NLEVELS) {
    (*initial_level)++;
  }
  
  // determine final level
  *final_level = *initial_level;
  while(patch_tree[*final_level] != NULL && *final_level < NLEVELS) {
    (*final_level)++;
  }
}

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

//*********************************************************************************************************
// fprintf_cube_geom():
//*********************************************************************************************************
void fprintf_cube_geom(FILE *fp, flouble x, flouble y, flouble z, flouble L)
{
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x-L/2.,y-L/2.,z-L/2., x+L/2.,y-L/2.,z-L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x-L/2.,y-L/2.,z-L/2., x-L/2.,y+L/2.,z-L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x+L/2.,y-L/2.,z-L/2., x+L/2.,y+L/2.,z-L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x-L/2.,y+L/2.,z-L/2., x+L/2.,y+L/2.,z-L/2.);
  
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x-L/2.,y-L/2.,z+L/2., x+L/2.,y-L/2.,z+L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x-L/2.,y-L/2.,z+L/2., x-L/2.,y+L/2.,z+L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x+L/2.,y-L/2.,z+L/2., x+L/2.,y+L/2.,z+L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x-L/2.,y+L/2.,z+L/2., x+L/2.,y+L/2.,z+L/2.);
  
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x-L/2.,y-L/2.,z-L/2., x-L/2.,y-L/2.,z+L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x+L/2.,y-L/2.,z-L/2., x+L/2.,y-L/2.,z+L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x-L/2.,y+L/2.,z-L/2., x-L/2.,y+L/2.,z+L/2.);
  fprintf(fp,"l   %lf %lf %lf   %lf %lf %lf      0 0 1\n", x+L/2.,y+L/2.,z-L/2., x+L/2.,y+L/2.,z+L/2.);
  
  fprintf(fp,"P   %lf %lf %lf 1 0 0 5\n",x,y,z);
}

//*********************************************************************************************************
// write_patch_geom():
//*********************************************************************************************************
void write_patch_geom(char *fname, patch_t patch)
{
  FILE       *fp;
  psubcube_t *ppsubcube_aux;
  psubcube_t  psubcube_aux;
  flouble     x, y, z;     // centre coordinates
  flouble     xc, yc, zc;  // corner coordinates
  flouble     L;
  double      Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;

  Xmin =  2.0;
  Xmax = -1.0;
  Ymin =  2.0;
  Ymax = -1.0;
  Zmin =  2.0;
  Zmax = -1.0;
  
  fp = fopen(fname,"w");
  ppsubcube_aux = NULL;
  while((ppsubcube_aux=patch_next_psubcube(&patch, ppsubcube_aux)) != NULL) {
    
    psubcube_aux = *ppsubcube_aux;
    ck2coor(&xc, &yc, &zc, &L, psubcube_aux->cubekey);    // returns corner coordinates *and* edge size

    // centre of the cube
    x = xc+L/2.;
    y = yc+L/2.;
    z = zc+L/2.;
    
    if(patch.periodic[0] && x < 0.5) x += 1.0;
    if(patch.periodic[1] && y < 0.5) y += 1.0;
    if(patch.periodic[2] && z < 0.5) z += 1.0;
    fprintf_cube_geom(fp,x,y,z,L);
    
    if(x > Xmax) Xmax = x;
    if(y > Ymax) Ymax = y;
    if(z > Zmax) Zmax = z;
    if(x < Xmin) Xmin = x;
    if(y < Ymin) Ymin = y;
    if(z < Zmin) Zmin = z;
    
  }
  fclose(fp);

  fprintf(stderr,"Xmin = %lf Xmax = %lf\n",fmod(Xmin+1.0,1.0),fmod(Xmax+1.0,1.0));
  fprintf(stderr,"Ymin = %lf Ymax = %lf\n",fmod(Ymin+1.0,1.0),fmod(Ymax+1.0,1.0));
  fprintf(stderr,"Zmin = %lf Zmax = %lf\n",fmod(Zmin+1.0,1.0),fmod(Zmax+1.0,1.0));
}

#endif