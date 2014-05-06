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
#include "../libutility/utility.h"

//Structure needed by UTARRAY to create a cubekey array
UT_icd patch_icd = {SIZEOF_PSUBCUBE, NULL, NULL, NULL };

void patch_create(patch_t** patch, int id, int level){
  if((*patch=(patch_t*)malloc(sizeof(patch_t)))==NULL){
    perror("malloc(patch_t): ");
    fprintf(stderr,"[patch_create] malloc for new patch structure failed\n");
    exit (-1);
  }
  //(*patch)->n_subcubes=0;
  (*patch)->id=id;
  (*patch)->level=level;
  (*patch)->parent=NULL;
  (*patch)->daughters=NULL;
  (*patch)->n_daughters=0;
  utarray_new((*patch)->psubcubes, &patch_icd);


  // stuff needed by patch physics
  init_patch_physics(*patch);
}

void inline patch_init(patch_t* patch, int id, int level){
  //(*patch)->n_subcubes=0;
  patch_clear(patch);
  patch->id=id;
  patch->level=level;
  if(patch->psubcubes==NULL){
    utarray_new(patch->psubcubes, &patch_icd);
  }
  
  // stuff needed by patch physics
  //init_patch_physics(patch); // not needed as called inside of patch_clear() already
}

void inline patch_clear(patch_t* patch){
  patch->parent=NULL;
  patch->daughters=NULL;
  patch->n_daughters=0;
  
  patch->id=-1;
  patch->level=-1;
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

uint64_t patch_get_num_psubcubes(patch_t* patch){
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
    
  //Calculate 26 adjacent subcubes' cubekeys (6 sides + 12 edges + 8 corners)
#ifdef AHF2_libtree_26neighbours
  ck_get_adjacents(sc->cubekey,&ck_adj);
#endif
  
  //Calculate 18 adjacent subcubes' cubekeys (6 sides + 12 edges)
#ifdef AHF2_libtree_18neighbours
  ck_get_adjacents_side_edge(sc->cubekey,&ck_adj);
#endif
  
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
          sc_aux->patch_id=patch->id;
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

void patch_connect_tree(patch_t** patch_tree, subcube_t** sc_table, uint64_t* n_patches){
  int initial_level=-1, final_level=-1;
  int patch_it, level_it;
  patch_t *p_patch_aux=NULL;
  
  for(level_it=0; level_it<NLEVELS; level_it++){
    if(patch_tree[level_it]!=NULL){
      //Check that n_patches is not 0 when there are patches at level_it
      if(n_patches[level_it]==0){
        fprintf(stderr, "ERROR: n_patches[%d] is 0, when patch_tree[%d] is not NULL. Aborting...\n", level_it,level_it);
        exit(-1);
      }
      //If initial_level equal to -1, we found the initial level
      if(initial_level==-1){
        initial_level=level_it;
        continue;
      }
    }//if(patch_tree[level_it]!=NULL
    else{
      //If patch_tree[level_it] is NULL, we passed by the deepest level
      if(initial_level!=-1){
        final_level=level_it-1;
        break;
      }
    }//else
  }//for

#ifdef  PATCH_DEBUG
  fprintf(stderr,"[patch_connect_tree] DEBUG: initial_level=%d, final_level=%d\n",initial_level, final_level);
#endif
  
  //We do not iterate for initial_level since the initial_level patches can not have any parents
  for(level_it=final_level;level_it>initial_level;level_it--){
    for(patch_it=0;patch_it<n_patches[level_it];patch_it++){
      p_patch_aux=&(patch_tree[level_it][patch_it]);
      patch_connect_my_parent(patch_tree, sc_table, p_patch_aux);
      #ifdef  PATCH_DEBUG
      if(p_patch_aux->parent==NULL){
        fprintf(stderr,"[patch_connect_tree] DEBUG: patch %d-%lu DIDN'T FOUND PARENT\n",p_patch_aux->level,p_patch_aux->id);
      }else{
        fprintf(stderr,"[patch_connect_tree] DEBUG: patch %d-%lu found his parent %d-%lu\n",p_patch_aux->level,p_patch_aux->id,p_patch_aux->parent->level,p_patch_aux->parent->id);
      }
      #endif
    }
  }
  
}

int patch_connect_my_parent(patch_t **patch_tree, subcube_t** sc_table, patch_t* patch){
#ifdef PATCH_DEBUG
  static FILE*   parent_and_daughtersF=NULL;
  static int try_fopen=1;
#endif
  cubekey_t      ck_patch=0, ck_parent=0;
  psubcube_t     psubcube_aux=NULL, psubcube_parent=NULL;
  psubcube_t*    ppsubcube_aux=NULL;
  
  //Extract first subcube from patch
  ppsubcube_aux=NULL;
  ppsubcube_aux=patch_next_psubcube(patch,ppsubcube_aux);
  psubcube_aux=*ppsubcube_aux;

  if(psubcube_aux==NULL){
    fprintf(stderr,"[patch_connect_my_parent] ERROR: patch_next_psubcube returned NULL. We cannot access to one subcube of the patch(%d-%"PRIu64")->n_subcubes=%"PRIu64"\n", patch->level, patch->id, patch->n_subcubes);
    return -1;
  }
  
  //Extract the cubekey of the first subcube extracted from the patch
  ck_patch=psubcube_aux->cubekey;
  
  //Calculate the subcube's cubekey of the parent of the first subcube extracted from the patch
  ck_get_parent(ck_patch, &ck_parent);
  if(ck_parent==0){
    fprintf(stderr,"[patch_connect_my_parent] ERROR: Parent cubekey equal to 0 when patch-cubekey(%d-%"PRIu64")/subcube->cubekey=%"PRIu64". We cannot search for the parent...\n", patch->level, patch->id, psubcube_aux->cubekey);
    return -1;
  }
  
  //Look for the parent subcube in the corresponding (level-1) subcube's hash table and point to it as your parent
  psubcube_parent=NULL;
  table_find_subcube(sc_table[patch->level-1], &ck_parent, &psubcube_parent);
  if(psubcube_parent==NULL){
    fprintf(stderr,"[patch_connect_my_parent] ERROR: patch(%d-%"PRIu64")/subcube->cubekey=%"PRIu64" didn't found his parent with the cubekey equal to %"PRIu64"...\n", patch->level, patch->id, psubcube_aux->cubekey, ck_parent);
    return -1;
  }
  
  //Check that psubcube_parent->patch_id is a valid id
  if(psubcube_parent->patch_id==-1){
    fprintf(stderr,"[patch_connect_my_parent] ERROR: the parent subcube of patch(%d-%"PRIu64") has NO_PATCH as patch_id. [n_part=%"PRIu64"]\n", patch->level, patch->id,patch_get_num_psubcubes(patch));
    return -1;
  }
  
  
  //========================================================
  // TO-DO: We should check psubcube_parent->patch_id < patches->n_patches[ck_get_depth(ck_parent)],
  //        but we do need to receive ahf2_patches_t instead of patch_tree...
  //========================================================

#ifdef PATCH_DEBUG
  if(parent_and_daughtersF==NULL && try_fopen==1){
    if((parent_and_daughtersF=fopen("parent_and_daughters.out", "w"))==NULL){
      perror("fopen(parent_and_daughters.out): ");
      fprintf(stderr, "Error opening parent_and_daughters.out\n");
      try_fopen=0;
    }
  }
  if(parent_and_daughtersF){
    fprintf(parent_and_daughtersF, "[patch_connect_my_parent] pointing patch(%d-%"PRIu64")->parent to patch(%d-%"PRIu64")\n",patch->level,patch->id,patch->level-1,psubcube_parent->patch_id);
    fflush(parent_and_daughtersF);
  }
#endif

  
  patch->parent=&(patch_tree[patch->level-1][psubcube_parent->patch_id]);
  
#ifdef PATCH_DEBUG
  if(parent_and_daughtersF){
    fprintf(parent_and_daughtersF, "[patch_connect_my_parent] realloc(patch(%d-%"PRIu64")->parent(%d-%"PRIu64")->daughters=%p),n_daughters=%"PRIu32"\n",patch->level,patch->id,patch->parent->level,patch->parent->id,patch->parent->daughters,patch->parent->n_daughters);
    fflush(parent_and_daughtersF);
  }
#endif
  
  //Register as a daughter patch in patch->parent->daughters
  //PROTECT daughters array access!!!
  patch->parent->daughters=realloc(patch->parent->daughters,(patch->parent->n_daughters+1)*sizeof(patch_t*));
  if(patch->parent->daughters==NULL){
    perror("realloc:");
    fprintf(stderr, "ERROR: realloc daughters array returned NULL. Aborting...\n");
    exit(-1);
  }
  patch->parent->daughters[patch->parent->n_daughters]=patch;
  patch->parent->n_daughters++;
  
  return 0;
}

void patch_connection_review(ahf2_patches_t* patches){
  FILE* patchtreeF=NULL;
  int it_level=0;
	uint64_t it_daughter=0, it_patch=0;
  patch_t *ppatch_aux=NULL, *ppatch_parent=NULL, *ppatch_daughter=NULL;
  
  //Review the patches->tree connections
	if((patchtreeF=fopen("patchtree.out", "w"))==NULL){
		perror("fopen(patchtree.out): ");
		fprintf(stderr, "Error opening patchtree.out\n");
		return;
	}
	fprintf(patchtreeF,"Let's review the patches->tree connections (Patch syntax is LEVEL-ID):\n");
	for(it_level=0;it_level<NLEVELS;it_level++){
		if(patches->tree[it_level]!=NULL){
			if(patches->n_patches[it_level]==0){
				fprintf(stderr,"[patch_connection_review] ERROR: patches->tree[%d] not NULL but patches->n_patches[%d] equal to 0!",it_level,it_level);
				continue;
			}
			fprintf(patchtreeF,"\nPatch_tree, at level %d, has %"PRIu64" patches (%"PRIu64" rejected)\n", it_level, patches->n_patches[it_level], patches->n_rejected_patches[it_level]);
			
			
			//Review all patches at level it_level
			for(it_patch=0;it_patch<patches->n_patches[it_level];it_patch++){
				ppatch_aux=&(patches->tree[it_level][it_patch]);
				fprintf(patchtreeF,"Patch %d-%"PRIu64" (Parent=%p, N_Subcubes=%"PRIu64", Ndaughters=%"PRIu32"):\n", ppatch_aux->level, ppatch_aux->id, ppatch_aux->parent, ppatch_aux->n_subcubes, ppatch_aux->n_daughters);
				//Check ppatch_aux PARENT
				ppatch_parent=ppatch_aux->parent;
				if(ppatch_parent==NULL){
					if(patches->tree[it_level-1]!=NULL){
						fprintf(patchtreeF,"\tERROR: patch %d-%"PRIu64" has no parent but patches->tree at level %d is not NULL!!\n", ppatch_aux->level, ppatch_aux->id, it_level-1);
					}
					//fprintf(patchtreeF,"\tPatch %d-%lu: NO parent (probably it's the initial level)\n", ppatch_aux->level, ppatch_aux->id);
				}//if
				else{
					fprintf(patchtreeF,"\tParent of %d-%"PRIu64" is patch %d-%lu\n", ppatch_aux->level, ppatch_aux->id, ppatch_parent->level, ppatch_parent->id);
				}
				
				//Check ppatch_aux DAUGHTERS
				if(ppatch_aux->daughters==NULL){
					if(it_level==MAX_DEPTH){
						fprintf(patchtreeF,"\tNo daughters because patch is at deepest posible level (%d)\n", MAX_DEPTH);
					}else{
						if(patches->tree[it_level+1]==NULL){
							fprintf(patchtreeF,"\tNo daughters because there are no patches at deeper level\n");
						}/*else{
							fprintf(patchtreeF,"\tWARNING: patch %d-%lu could have daughters, but patch->daughters pointer is NULL. (patch->n_daughters=%d)\n", ppatch_aux->level, ppatch_aux->id, ppatch_aux->n_daughters);
						}*/
					}
				}//if(ppatch_aux->daughters==NULL)
				else{
					//ppatch_daughter=ppatch_aux->daughters[0];
					for(it_daughter=0;it_daughter<ppatch_aux->n_daughters;it_daughter++){
						ppatch_daughter=ppatch_aux->daughters[it_daughter];
						if(ppatch_daughter==NULL){
							fprintf(patchtreeF,"\tERROR: pointer to daughter #%lu is NULL!!\n", it_daughter);
							continue; //Should we break??
						}
						fprintf(patchtreeF,"\tDaughter of %d-%lu is patch %d-%lu\n", ppatch_aux->level, ppatch_aux->id, ppatch_daughter->level, ppatch_daughter->id);
						fflush(patchtreeF);
					}
				}
			}//for(it_patch in patches->tree[it_level])
		}
	}
  if(patchtreeF){
    fclose(patchtreeF);
    patchtreeF=NULL;
  }
}

void patch_formation_review(ahf2_patches_t* patches){
  int it_level=0;
  FILE *patchF=NULL;
  
  if((patchF=fopen("patch.out", "w"))==NULL){
		perror("fopen(patch.out): ");
		fprintf(stderr, "Error opening patch.out\n");
		return;
	}

  /*
	 * Check and print the patches->tree formation
	 */
	fprintf(patchF,"LET'S CHECK PATCHES AT ALL LEVELS\n");
	for(it_level=0;it_level<NLEVELS;it_level++){
		if(patches->tree[it_level]==NULL){
			/*if(patches->n_rejected_patches[it_level]==0){
				break;
			}*/
			fprintf(patchF,"=> Level %2d, patches->tree[%d]==NULL . %lu invalid patches\n",it_level,it_level,patches->n_rejected_patches[it_level]);
			continue;
		}
		if(patches->n_patches[it_level]==0){
			fprintf(stderr,"[patch_tree_review]ERROR!: patches->tree[%d]!=NULL and patches->n_patches[%d]=%lu (Level %d)\n",it_level,it_level,patches->n_patches[it_level],it_level);
			continue;
		}
    
		fprintf(patchF,"====> Level %2d, %lu valid patches, %lu rejected patches\n",it_level,patches->n_patches[it_level],patches->n_rejected_patches[it_level]);
	}//for(depth) CHECKING PATCHES AT ALL LEVELS and printing to patchF

  if(patchF){
    fclose(patchF);
    patchF=NULL;
  }
}


void init_patch_physics(patch_t *patch)
{
  patch->Npart                = 0;
  patch->Npart_stored         = 0;
  patch->n_subcubes           = 0;
  
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
  
  patch->periodic[0]          = 0;
  patch->periodic[1]          = 0;
  patch->periodic[2]          = 0;
  patch->periodic_shift[0]    = 0;
  patch->periodic_shift[1]    = 0;
  patch->periodic_shift[2]    = 0;
}

void add_patch_physics(patch_t *patch, psubcube_t psc)
{
  flouble  xc, yc, zc, L;    // return values from ck2coor(): corner position and side length of cube
  partptr *pp_part;          // subcube_next_particle() needs pointer to pointer to particle, since UTarray stores pointers to particle
  partptr  p_part;           // we will handle particles with a "direct" pointer to make syntax more clear (p_part=*pp_part)
  uint64_t npart;
  int8_t   i, periodic[3];
  double   pos[3], posc[3], pospart[3];
  
  /*=============================================
   * Commented all FILE* patch_physicsF related code
   * because it's generating huge debug log file 4,1GB 
   * for FullBox_z0.0_0256
  #ifdef  PATCH_DEBUG
  static FILE* patch_physicsF=NULL;
  static int try_fopen=1;
  #endif
  ************************************************/
  
  // obtain lower-left corner coordinates of subcube
  ck2coor(&xc, &yc, &zc, &L, psc->cubekey);
  posc[0] = xc;
  posc[1] = yc;
  posc[2] = zc;
  
  // centre of subcube
  for(i=0; i<3; i++)
    pos[i] = posc[i]+L/2.;
  
  
  //===============================================================================================
  // check whether the patch itselfs starts to wrap periodically
  //===============================================================================================
  for(i=0; i<3; i++ ) {
    
    // at this moment we do not need periodic treatment of the present cube
    periodic[i] = 0;

    // the patch is not periodically wrapped in i-direction (yet)
    if(patch->periodic[i] == FALSE) {
      
      // does the present cube touch the right edge?
      if(posc[i]+1.5*L > 1) {
        patch->periodic_shift[i] = +1;    // we will move all "left"  cubes to the range [1.0,1.5]
        patch->periodic[i]       = TRUE;  // take note that the patch is periodically wrapped in i-direction
      }
      // does the present cube touch the left edge?
      else if(posc[i]-0.5*L < 0) {
        patch->periodic_shift[i] = -1;    // we will move all "right" cubes to the range [-0.5,0.0]
        patch->periodic[i]       = TRUE;  // take note that the patch is periodically wrapped in i-direction
      }
      else {
        // no cube has yet touched any boundary and hence we leave all flags as they are (i.e. at zero!)
        // patch->periodic_shift[i] = 0;
        // patch->periodic[i]       = 0;
        // periodic[i]              = 0;
      }
    } // if(patch->periodic[i])
    
    // the patch is periodically wrapped in x-direction
    else {
      // does the present cube require a shift?
      if((pos[i] > 0.5 && patch->periodic_shift[i] == -1) || (pos[i] <= 0.5 && patch->periodic_shift[i] == +1)) {
        periodic[i] = 1;  // signal to use patch->periodic_shift[0]
                          // note, perodic_shift[] encodes whether we first hit the right or left edge of the box!
      }
    } // else(patch->periodic[i])
    
    // periodic patch -> correct cube coordinates
    if(periodic[i])
      pos[i] += (double)patch->periodic_shift[i];
    
  } // for(i)
  
  /*=============================================
#ifdef  PATCH_DEBUG
  if(patch_physicsF==NULL && try_fopen==1){
    if((patch_physicsF=fopen("patch_physics.out", "w"))==NULL){
      perror("fopen(patch_physics.out): ");
      fprintf(stderr, "Error opening patch_physics.out\n");
      try_fopen=0;
    }
  }
  if(patch_physicsF){
    fprintf(patch_physicsF,"[add_patch_physics] Adding physics to patch(%d-%"PRIu64") from subcube (npart=%"PRIu64",npart_stored=%"PRIu64",cubekey=%"PRIu64")\n", patch->level, patch->id, subcube_get_num_particles(psc),subcube_get_num_stored_particles(psc), psc->cubekey);
  }
#endif
  ************************************************/
  
  // update patch physics
  patch->Npart                 += psc->nparticles;
  patch->n_subcubes            += 1;
  
  if(pos[0] < patch->Xmin) patch->Xmin = pos[0];
  if(pos[1] < patch->Ymin) patch->Ymin = pos[1];
  if(pos[2] < patch->Zmin) patch->Zmin = pos[2];
  
  if(pos[0] > patch->Xmax) patch->Xmax = pos[0];
  if(pos[1] > patch->Ymax) patch->Ymax = pos[1];
  if(pos[2] > patch->Zmax) patch->Zmax = pos[2];
  
  if(psc->nparticles > patch->Nmax) {
    patch->Nmax = psc->nparticles;
    
    patch->centre.part.Nmax[0] = (double)pos[0];
    patch->centre.part.Nmax[1] = (double)pos[1]; // not to be divided by anything
    patch->centre.part.Nmax[2] = (double)pos[2];
  }
  
	pp_part = NULL;
  npart   = 0;
	while((pp_part=subcube_next_particle(psc, pp_part))){
    npart++;
    p_part      = *pp_part;
    pospart[0]  =  p_part->pos[0];
    pospart[1]  =  p_part->pos[1];
    pospart[2]  =  p_part->pos[2];
    
    // periodic patch -> correct particle positions
    if(periodic[0]) pospart[0] += (double)patch->periodic_shift[0];
    if(periodic[1]) pospart[1] += (double)patch->periodic_shift[1];
    if(periodic[2]) pospart[2] += (double)patch->periodic_shift[2];
    
    patch->centre.part.com[0] += pospart[0];
    patch->centre.part.com[1] += pospart[1];  // to be divided by total number of stored(!!!) particles in the patch
    patch->centre.part.com[2] += pospart[2];
  }
  patch->Npart_stored += npart;

  /*=============================================
  #ifdef  PATCH_DEBUG
  if(patch_physicsF){
    if(npart != psc->nparticles) {
      fprintf(patch_physicsF,"[add_patch_physics] Adding REFINED subcube (psc->nparticles=%"PRIu64") to patch (%d-%"PRIu64")\n",psc->nparticles, patch->level,patch->id);
      //exit(0);
    }else{
      fprintf(patch_physicsF,"[add_patch_physics] Adding NOT REFINED subcube (psc->nparticles=%"PRIu64") to patch (%d-%"PRIu64")\n",psc->nparticles, patch->level,patch->id);
    }
  }
  #endif
  ============================================= */
}

void finish_patch_physics(patch_t *patch)
{
  double tmp;
  
  // normalize centres
  patch->centre.cube.geom[0] /= (double)patch->n_subcubes;
  patch->centre.cube.geom[1] /= (double)patch->n_subcubes;
  patch->centre.cube.geom[2] /= (double)patch->n_subcubes;

  patch->centre.cube.wgeom[0] /= (double)patch->Npart;
  patch->centre.cube.wgeom[1] /= (double)patch->Npart;
  patch->centre.cube.wgeom[2] /= (double)patch->Npart;
  
  patch->centre.part.com[0]   /= (double)patch->Npart_stored;
  patch->centre.part.com[1]   /= (double)patch->Npart_stored;
  patch->centre.part.com[2]   /= (double)patch->Npart_stored;
 
  // correct for periodicity
  patch->centre.cube.geom[0]  = fmod(patch->centre.cube.geom[0]  +1.0, 1.0);
  patch->centre.cube.geom[1]  = fmod(patch->centre.cube.geom[1]  +1.0, 1.0);
  patch->centre.cube.geom[2]  = fmod(patch->centre.cube.geom[2]  +1.0, 1.0);

  patch->centre.cube.wgeom[0] = fmod(patch->centre.cube.wgeom[0] +1.0, 1.0);
  patch->centre.cube.wgeom[1] = fmod(patch->centre.cube.wgeom[1] +1.0, 1.0);
  patch->centre.cube.wgeom[2] = fmod(patch->centre.cube.wgeom[2] +1.0, 1.0);

  patch->centre.part.com[0]   = fmod(patch->centre.part.com[0]   +1.0, 1.0);
  patch->centre.part.com[1]   = fmod(patch->centre.part.com[1]   +1.0, 1.0);
  patch->centre.part.com[2]   = fmod(patch->centre.part.com[2]   +1.0, 1.0);
  
  patch->centre.part.Nmax[0]  = fmod(patch->centre.part.Nmax[0]  +1.0, 1.0);
  patch->centre.part.Nmax[1]  = fmod(patch->centre.part.Nmax[1]  +1.0, 1.0);
  patch->centre.part.Nmax[2]  = fmod(patch->centre.part.Nmax[2]  +1.0, 1.0);
  
  
  patch->Xmin = fmod(patch->Xmin +1.0, 1.0);
  patch->Ymin = fmod(patch->Ymin +1.0, 1.0);
  patch->Zmin = fmod(patch->Zmin +1.0, 1.0);

  patch->Xmax = fmod(patch->Xmax +1.0, 1.0);
  patch->Ymax = fmod(patch->Ymax +1.0, 1.0);
  patch->Zmax = fmod(patch->Zmax +1.0, 1.0);
  
  if(patch->periodic[0]) SWAP(patch->Xmin,patch->Xmax,tmp);
  if(patch->periodic[1]) SWAP(patch->Ymin,patch->Ymax,tmp);
  if(patch->periodic[2]) SWAP(patch->Zmin,patch->Zmax,tmp);
}
