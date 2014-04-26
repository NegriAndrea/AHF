/*
 * patch.c
 * 
 * Copyright 2014 F. Campos <fernandocdelpozo@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include "patch.h"
#include "subcube.h"

//Structure needed by UTARRAY to create a cubekey array
UT_icd patch_icd = {SIZEOF_PSUBCUBE, NULL, NULL, NULL };

void patch_create(patch_t** patch, int id, int level){
  if((*patch=(patch_t*)malloc(sizeof(patch_t)))==NULL){
    perror("malloc(patch_t): ");
    fprintf(stderr,"[patch_create] malloc for new patch structure failed\n");
    exit (-1);
  }
  //(*patch)->n_subcubes=0;
  (*patch)->patch_id=id;
  (*patch)->patch_level=level;
  utarray_new((*patch)->psubcubes, &patch_icd);


  // stuff needed by patch physics
  init_patch_physics(*patch);
}

void inline patch_init(patch_t* patch, int id, int level){
  //(*patch)->n_subcubes=0;
  patch_clear(patch);
  patch->patch_id=id;
  patch->patch_level=level;
  if(patch->psubcubes==NULL){
    utarray_new(patch->psubcubes, &patch_icd);
  }
  
  // stuff needed by patch physics
  //init_patch_physics(patch); // not needed as called inside of patch_clear() already
}
void inline patch_clear(patch_t* patch){
  patch->parent=NULL;
  patch->daughters=NULL;
  patch->patch_id=-1;
  patch->patch_level=-1;
  if(patch->psubcubes){
    utarray_clear(patch->psubcubes);
  }

  // stuff needed by patch physics
  init_patch_physics(patch);
}

void patch_add_psubcube(patch_t* patch, psubcube_t* psc){

  utarray_push_back(patch->psubcubes,psc);

  // stuff needed by patch physics
  add_patch_physics(patch, *psc);
}


//psubcube_t* patch_pop_back(patch_t* patch){
  //psubcube_t* psc_aux;
  //psc_aux=utarray_pop_back(patch->psubcubes);
  //return psc_aux;
//}

psubcube_t* patch_next_psubcube(patch_t* patch, psubcube_t* psubcube_it){
  return (psubcube_t*) utarray_next(patch->psubcubes,psubcube_it);
}

int patch_get_num_psubcubes(patch_t* patch){
  if(patch->psubcubes==NULL) return 0;
  return utarray_len(patch->psubcubes);
}
void patch_free_psubcubes_array(patch_t* patch){
  if(patch->psubcubes){
    utarray_free(patch->psubcubes);
    patch->psubcubes=NULL;
  }
}

int patch_include_adjacent_subcubes(subcube_t* sc_table, subcube_t* sc, patch_t* patch){
  //static long unsigned int calls, call_depth;
  //Look for adjacent subcubes
  ck_adjacents_t ck_adj;
  cubekey_t* pck_aux=NULL;
  subcube_t* sc_aux=NULL;
  int i,n_subcubes=0;
  subcube_t* valid_adjacents[NUM_ADJACENT]={NULL};
  
  //calls++;
  //call_depth++;
  //if(calls%100==0)
    //fprintf(stderr,"[patch_include_adjacent_subcubes] Recursive call number %lu, depth %lu\n",calls, call_depth);
    
  //Calculate 26 adjacent subcubes' cubekeys
  ck_get_adjacents(sc->cubekey,&ck_adj);
  
  //Iterate on the adjacent cubekeys with pck_aux (ckBBB is the first element in ck_adj. &ckBBB is equal to &ck_adj)
  //pck_aux=&ck_adj;
  pck_aux=&(ck_adj.ckBBB);
  
  for(i=0;i<NUM_ADJACENT;i++){
    if(*pck_aux!=0){
      //If the adjacent subcube exists, we include it in the current patch
      table_find_subcube(sc_table,pck_aux,&sc_aux);
      if(sc_aux!=NULL){
	if(sc_aux->patch_id==NO_PATCH){
	  valid_adjacents[i]=sc_aux;
	  n_subcubes++;
	  sc_aux->patch_id=patch->patch_id;
	  patch_add_psubcube(patch,&sc_aux);
	}
      }
    }
    pck_aux++;
  }
  for(i=0;i<NUM_ADJACENT;i++){
    //Recursive call for adjacent subcubes from the not previously visited subcubes
    if(valid_adjacents[i]!=NULL){
      n_subcubes+=patch_include_adjacent_subcubes(sc_table,valid_adjacents[i],patch);
    }
  }
  //call_depth--;
  return n_subcubes;
}

void init_patch_physics(patch_t *patch)
{
  patch->Npart                = 0;
  
  patch->centre.cube.geom[0]  = 0.0;
  patch->centre.cube.geom[1]  = 0.0;
  patch->centre.cube.geom[2]  = 0.0;
  patch->centre.cube.wgeom[0] = 0.0;
  patch->centre.cube.wgeom[1] = 0.0;
  patch->centre.cube.wgeom[2] = 0.0;
  patch->centre.part.com[0]   = 0.0;
  patch->centre.part.com[1]   = 0.0;
  patch->centre.part.com[2]   = 0.0;
  patch->centre.part.Nmax[0]  = 0.0;
  patch->centre.part.Nmax[1]  = 0.0;
  patch->centre.part.Nmax[2]  = 0.0;

  patch->Nmax                 = 0;
  patch->Xmin                 = +2.0;
  patch->Xmax                 = -1.0;
  patch->Ymin                 = +2.0;
  patch->Ymax                 = -1.0;
  patch->Zmin                 = +2.0;
  patch->Zmax                 = -1.0;
  
  patch->radius               = 0.0;
}

void add_patch_physics(patch_t *patch, psubcube_t psc)
{
  flouble  x, y, z, xpart, ypart, zpart;
  partptr *pp_part;          //subcube_next_particle() needs pointer to pointer to particle, since UTarray stores pointers to particle
  partptr  p_part;           //We will handle particles with a "direct" pointer to make syntax more clear (p_part=*pp_part)

  // obtain centre coordinates of subcube
  ck2coor_center(&x, &y, &z, psc->cubekey);
  
  // update patch physics
  patch->Npart                 += psc->nparticles;
  patch->n_subcubes            += 1;
  
  patch->centre.cube.geom[0]   += (double)x;
  patch->centre.cube.geom[1]   += (double)y;                    // to be divided by the total number of subcubes in the patch
  patch->centre.cube.geom[2]   += (double)z;
  patch->centre.cube.wgeom[0]  += (double)psc->nparticles*(double)x;
  patch->centre.cube.wgeom[1]  += (double)psc->nparticles*(double)y;   // to be divided by the total number of particles in the patch
  patch->centre.cube.wgeom[2]  += (double)psc->nparticles*(double)z;

  if(x < patch->Xmin) patch->Xmin = x;
  if(y < patch->Xmin) patch->Ymin = y;
  if(z < patch->Xmin) patch->Zmin = z;
  if(x > patch->Xmax) patch->Xmax = x;
  if(y > patch->Xmax) patch->Ymax = y;
  if(z > patch->Xmax) patch->Zmax = z;
  
  if(psc->nparticles > patch->Nmax) {
    patch->Nmax = psc->nparticles;
    
    patch->centre.part.Nmax[0] = (double)x;
    patch->centre.part.Nmax[1] = (double)y; // not to be divided by anything
    patch->centre.part.Nmax[2] = (double)z;
  }
  
	pp_part=NULL;
	while((pp_part=subcube_next_particle(psc, pp_part))){
    p_part=*pp_part;
    xpart = p_part->pos[0];
    ypart = p_part->pos[1];
    zpart = p_part->pos[2];
    
    patch->centre.part.com[0] += xpart;
    patch->centre.part.com[1] += ypart;  // to be divided by total number of particles in the patch
    patch->centre.part.com[2] += zpart;
  }

}

void finish_patch_physics(patch_t *patch)
{
  patch->centre.cube.geom[0] /= (double)patch->n_subcubes;
  patch->centre.cube.geom[1] /= (double)patch->n_subcubes;
  patch->centre.cube.geom[2] /= (double)patch->n_subcubes;

  patch->centre.cube.wgeom[0] /= (double)patch->Npart;
  patch->centre.cube.wgeom[1] /= (double)patch->Npart;
  patch->centre.cube.wgeom[2] /= (double)patch->Npart;
  
  patch->centre.part.com[0]   /= (double)patch->Npart;
  patch->centre.part.com[1]   /= (double)patch->Npart;
  patch->centre.part.com[2]   /= (double)patch->Npart;
 
}
