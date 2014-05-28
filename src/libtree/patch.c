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
#include <stdio.h>
#include <string.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "patch.h"
#include "subcube.h"
#include "../libutility/utility.h"

#ifdef PATCH_THREADS_LOG
FILE* patch_threadsF=NULL;
char patch_threads_name[256];
int try_patch_threadsF=1;
#endif

//Structure needed by UTARRAY to create a cubekey array
UT_icd patch_icd = {SIZEOF_PSUBCUBE, NULL, NULL, NULL };

void patch_create(patch_t** patch, int64_t id, int32_t level){
#ifdef WITH_OPENMP
  io_logging_msg(global_io.log, INT32_C(2),
      "[%s:%d] We should call patch_create(%2d-%ld) ONLY ONCE per thread building patches at every level. (thread %d/%d)", __FILE__, __LINE__, level, id, omp_get_thread_num(), omp_get_num_threads());
#else
  io_logging_msg(global_io.log, INT32_C(2),
      "[%s:%d] We should call patch_create(%2d-%ld) ONLY ONCE per thread building patches at every level.", __FILE__, __LINE__, level, id);
#endif

  if((*patch=(patch_t*)malloc(sizeof(patch_t)))==NULL){
    perror("malloc(patch_t): ");
    fprintf(stderr,"[patch_create] malloc for new patch structure failed\n");
    exit (EXIT_FAILURE);
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

void inline patch_init(patch_t* patch, int64_t id, int32_t level){
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
void inline patch_free(patch_t** patch){
  if(*patch==NULL) return;
  patch_free_content(*patch);
  free(*patch);
  *patch=NULL;
}

void inline patch_free_content(patch_t* patch){
  if(patch==NULL) return;
  if((patch)->psubcubes){
    patch_free_psubcubes_array(patch);
  }
  if(patch->daughters){
    free(patch->daughters);
    patch->daughters=NULL;
  }
}

void patch_add_psubcube(patch_t* patch, psubcube_t* psc){

  (*psc)->patch_id=patch->id;
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
  if(patch==NULL) return;

  if(patch->psubcubes){
    utarray_free(patch->psubcubes);
    patch->psubcubes=NULL;
  }
}

void ahf2_patches_free(ahf2_patches_t** patches){
  clock_t t_start, t_end;       //Timing variables for POSIX Clock
  double elapsed;
  int it_level, it_patch;
  if(*patches==NULL) return;
  t_start=clock();
  for(it_level=0; it_level<NLEVELS; it_level++){
    if((*patches)->tree[it_level]==NULL) continue;
    for(it_patch=0; it_patch<(*patches)->n_patches[it_level]; it_patch++){
      patch_free_content(&(*patches)->tree[it_level][it_patch]);
    }
    free((*patches)->tree[it_level]);
    (*patches)->tree[it_level]=NULL;
  }
  free(*patches);
  *patches=NULL;
  t_end=clock();
  elapsed=((double) (t_end - t_start)) / CLOCKS_PER_SEC;
  io_logging_msg(global_io.log, INT32_C(2),
      "[%s:%d] We freed ahf2_patches_structure in %f seconds.", __FILE__, __LINE__, elapsed);
}

int patch_include_adjacent_subcubes(subcube_t* sc_table, subcube_t* sc, patch_t* patch){
  //static long unsigned int calls, call_depth;
  //Look for adjacent subcubes
  ck_adjacents_t ck_adj;
  cubekey_t *pck_aux=NULL, *pck_aux_next=NULL;
  subcube_t* sc_aux=NULL;
  int64_t i,j;
  int64_t new_subcubes=0;
  //subcube_t* valid_adjacents[NUM_ADJACENT]={NULL};
  visitable_subcubes_t *v_s_current=NULL, *v_s_next=NULL, *v_s_SWAP=NULL;

  //for(i=0;i<NUM_ADJACENT;i++) valid_adjacents[i]=NULL;
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

  //Create visitable_subcubes (current and next) structures
  visitable_subcubes_create(&v_s_current);
  visitable_subcubes_create(&v_s_next);

  //Insert ck_adj in v_s_current
  pck_aux=&(ck_adj.ckBBB);
  new_subcubes=0;
  for(i=0;i<NUM_ADJACENT;i++){
      if(*pck_aux!=0){
        visitable_subcubes_add(v_s_current,pck_aux);
      }
      pck_aux++;
  }

//  pck_aux=&(ck_adj.ckBBB);
//  new_subcubes=0;

  /*
   * This while will keep searching adjacent subcubes until all posibilites have been checked
   */
  while(visitable_subcubes_count(v_s_current)>0){

    //Visit all visitable subcubes
    for(i=0;i<visitable_subcubes_count(v_s_current);i++){
      pck_aux=visitable_subcubes_get(v_s_current, i);

      //If the adjacent subcube exists, we include it in the current patch
      table_find_subcube(sc_table,pck_aux,&sc_aux);
      if(sc_aux!=NULL){
        //if the adjacent subcube has not been visited yet
        if(sc_aux->patch_id==NO_PATCH){
          //Include the subcube in the patch we are building
          patch_add_psubcube(patch,&sc_aux);
          new_subcubes++;

          //Include it's adjacents in the next visitable subcubes structure (v_s_next)
          //Calculate 26 adjacent subcubes' cubekeys (6 sides + 12 edges + 8 corners)
          #ifdef AHF2_libtree_26neighbours
          ck_get_adjacents(sc_aux->cubekey,&ck_adj);
          #endif
          //Calculate 18 adjacent subcubes' cubekeys (6 sides + 12 edges)
          #ifdef AHF2_libtree_18neighbours
          ck_get_adjacents_side_edge(sc_aux->cubekey,&ck_adj);
          #endif

          pck_aux_next=&(ck_adj.ckBBB);
          for(j=0;j<NUM_ADJACENT;j++){
            if(*pck_aux_next!=0){
              visitable_subcubes_add(v_s_next,pck_aux_next);
            }
            pck_aux_next++;
          }

        }//if(NO_PATCH)
        // Warning message if we visit a subcube WE have already visited previously
        /*else if(sc_aux->patch_id==patch->id){
          fprintf(stderr,"[patch_include_adjacent_subcubes] WARNING: we found a visited subcube (ck=%"PRIu64") with patch_id=%"PRIi64" when populating patch %d-%"PRIu64" from ck=%"PRIu64".\n", sc_aux->cubekey, sc_aux->patch_id, patch->level,patch->id, sc->cubekey);
        }*/
        // ERROR message if we visit a subcube with another patch_id. We should reach the current subcube sc when building the previous patch sc_aux->patch_id
        else if(sc_aux->patch_id!=patch->id){
          fprintf(stderr,"[%s:%d] ERROR in patch_include_adjacent_subcubes(): we found a visited subcube (ck=%"PRIck\
              ") with patch_id=%"PRId64" when populating patch %d-%"PRId64" from ck=%"PRIck".\n",
              __FILE__, __LINE__, sc_aux->cubekey, sc_aux->patch_id, patch->level,patch->id, sc->cubekey);
        }//else if ERROR
      }//if(sc_aux!=NULL)
    }//  for(i<visitable_subcubes_count(v_s_current))

    /*
     * In the next while iteration, we will need to use the next visitable_subcubes (v_s_next) as the current one.
     * And we will need another visitable_subcubes structure for the next visitable: let's swap and recycle!
     */
    v_s_SWAP=v_s_current;
    v_s_current=v_s_next;
    v_s_next=v_s_SWAP;

    //Clear the new v_s_next structure
    visitable_subcubes_clear(v_s_next);
  }//while(visitable_subcubes_count not EMPTY)

  //At this point, both current and next visitable_subcubes must be empty
  if(visitable_subcubes_count(v_s_current)!=0){
    fprintf(stderr,"[%s:%d] ERROR: In patch_include_adjacent_subcubes(), we exit the while loop and v_s_current is not empty!",
                __FILE__, __LINE__);
  }
  if(visitable_subcubes_count(v_s_next)!=0){
      fprintf(stderr,"[%s:%d] ERROR: In patch_include_adjacent_subcubes(), we exit the while loop and v_s_next is not empty!",
                  __FILE__, __LINE__);
  }
  //Free memory of both visitable_subcubes structures
  visitable_subcubes_free(&v_s_current);
  visitable_subcubes_free(&v_s_next);

  return new_subcubes;
}

void patch_connect_tree(patch_t** patch_tree, subcube_t** sc_table, int64_t* n_patches){
  int initial_level=-1, final_level=-1;
  int patch_it, level_it;
  patch_t *p_patch_aux=NULL;
  
  for(level_it=0; level_it<NLEVELS; level_it++){
    if(patch_tree[level_it]!=NULL){
      //Check that n_patches is not 0 when there are patches at level_it
      if(n_patches[level_it]==0){
        fprintf(stderr, "ERROR: n_patches[%d] is 0, when patch_tree[%d] is not NULL. Aborting...\n", level_it,level_it);
        exit(EXIT_FAILURE);
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
    fprintf(stderr,"[patch_connect_my_parent] ERROR: patch_next_psubcube returned NULL.\
         We cannot access to one subcube of the patch(%d-%"PRIu64")->n_subcubes=%"PRIu64"\n", patch->level, patch->id, patch->n_subcubes);
    return -1;
  }
  
  //Extract the cubekey of the first subcube extracted from the patch
  ck_patch=psubcube_aux->cubekey;
  
  //Calculate the subcube's cubekey of the parent of the first subcube extracted from the patch
  ck_get_parent(ck_patch, &ck_parent);
  if(ck_parent==0){
    fprintf(stderr,"[patch_connect_my_parent] ERROR: Parent cubekey equal to 0 when "\
        "patch-cubekey(%"PRId32"-%"PRId64")/subcube->cubekey=%"PRIck". We cannot search for the parent...\n",
        patch->level, patch->id, psubcube_aux->cubekey);
    return -1;
  }
  
  //Look for the parent subcube in the corresponding (level-1) subcube's hash table and point to it as your parent
  psubcube_parent=NULL;
  table_find_subcube(sc_table[patch->level-1], &ck_parent, &psubcube_parent);
  if(psubcube_parent==NULL){
    fprintf(stderr,"[patch_connect_my_parent] ERROR: patch(%"PRId32"-%"PRId64")/subcube->cubekey=%"PRIck" didn't found"\
        " his parent with the cubekey equal to %"PRIck"...\n",
        patch->level, patch->id, psubcube_aux->cubekey, ck_parent);
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
    exit(EXIT_FAILURE);
  }
  patch->parent->daughters[patch->parent->n_daughters]=patch;
  patch->parent->n_daughters++;
  
  return 0;
}

void patch_connection_review(ahf2_patches_t* patches){
  FILE* patchtreeF=NULL;
  char patchtree_name[256];
  int it_level=0;
  uint64_t it_daughter=0, it_patch=0;
  patch_t *ppatch_aux=NULL, *ppatch_parent=NULL, *ppatch_daughter=NULL;
  clock_t t_start, t_end;       //Timing variables for POSIX Clock
  double elapsed;

  t_start=clock();
  //Review the patches->tree connections
  strncpy(patchtree_name,global_io.params->outfile_prefix,256);
  strcat(patchtree_name, ".patchtree.log");

	if((patchtreeF=fopen(patchtree_name, "w"))==NULL){
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
        if(ppatch_aux->parent==NULL){
          fprintf(patchtreeF,"Patch %d-%"PRIu64" (Parent=NULL, N_Subcubes=%"PRIu64", Ndaughters=%"PRIu32"):\n", ppatch_aux->level, ppatch_aux->id, ppatch_aux->n_subcubes, ppatch_aux->n_daughters);
        }
        else{
          fprintf(patchtreeF,"Patch %d-%"PRIu64" (Parent=VALID, N_Subcubes=%"PRIu64", Ndaughters=%"PRIu32"):\n", ppatch_aux->level, ppatch_aux->id, ppatch_aux->n_subcubes, ppatch_aux->n_daughters);
        }
          
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
  t_end=clock();
  elapsed=((double) (t_end - t_start)) / CLOCKS_PER_SEC;
  io_logging_msg(global_io.log, INT32_C(2),
        "[%s:%d] patch_connection_review(). patchtree.out file generation took %f seconds.", __FILE__, __LINE__, elapsed);

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

void patches_generation(ahf2_patches_t** patches, subcube_t** sc_table, int32_t initial_depth, int32_t final_depth, uint64_t NminPerHalo){
  int i=0, level=0;
  uint64_t rc=0;
#ifdef WITH_OPENMP
  double start,end;
#endif

#ifdef PATCH_THREADS_LOG
  strncpy(patch_threads_name,global_io.params->outfile_prefix,256);
  strcat(patch_threads_name, ".patchthreads.log");
  if((patch_threadsF=fopen(patch_threads_name, "w"))==NULL) {
    perror("fopen(patch_thread.log): ");
    fprintf(stderr, "Error opening %s\n", patch_threads_name);
    exit(EXIT_FAILURE);
  }
#endif

  //Allocate memory for ahf2_patches_t structure and initialize it
  if ((*patches = (ahf2_patches_t*) malloc(SIZEOF_AHF2_PATCHES)) == NULL) {
    perror("malloc:");
    fprintf(stderr, "[%s:%d] ERROR: malloc for patches structure returned NULL. Aborting...\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  //Initialize patches->tree arrays and patches->n_patches counters
  for (i = 0; i < NLEVELS; i++) {
    (*patches)->tree[i] = NULL;
    (*patches)->n_patches[i] = 0;
    (*patches)->n_rejected_patches[i] = 0;
    //n_subcubes_in_patches[i]=0;
  }

  /*Generate patches at every level [PATCH-TREE]*/
  //We launch one thread per level to execute patches_generation_per_level()
#ifdef PATCH_THREADS_LOG
  fprintf(patch_threadsF, "[%s:%d] Launching OpenMP threads in for-loop calling patches_generation_per_level().\n", __FILE__, __LINE__); fflush(patch_threadsF);
#endif

#ifdef WITH_OPENMP
  io_logging_msg(global_io.log, INT32_C(2),
      "[%s:%d] Launching OpenMP threads in for-loop calling patches_generation_per_level().", __FILE__, __LINE__);
  start=omp_get_wtime();

  #ifdef PATCH_THREADS_LOG
    #pragma omp parallel for default(none) private(level,rc) \
    shared(initial_depth, final_depth, patches, sc_table, NminPerHalo, stderr, global_io, patch_threadsF) schedule(dynamic,1)

  #else
    #pragma omp parallel for default(none) private(level,rc) \
    shared(initial_depth, final_depth, patches, sc_table, NminPerHalo, stderr, global_io) schedule(dynamic,1)
  #endif //PATCH_THREADS_LOG
#endif//WITH_OPENMP
  for(level=initial_depth; level<=final_depth; level++){

    #ifdef WITH_OPENMP
    fprintf(stderr, "(thread %d/%d) patches_generation_per_level(level=%2d)\n", omp_get_thread_num(), omp_get_num_threads(), level);
    #endif

    #ifdef PATCH_THREADS_LOG
    fprintf(patch_threadsF, "\t(thread %d/%d) Calling patches_generation_per_level(level=%2d)\n", omp_get_thread_num(), omp_get_num_threads(), level); fflush(patch_threadsF);
    #endif

    rc=patches_generation_per_level(*patches, sc_table[level], level, NminPerHalo);
    //Sanity check
    if((*patches)->n_patches[level] != rc){
      fprintf(stderr, "[%s:%d] ERROR: returned value from patches_generation_per_level (%"PRId64") doesn't match patches->n_patches[level] (%"PRIu64")\n"
          , __FILE__, __LINE__, rc, (*patches)->n_patches[level]);
    }
  }//for


#ifdef PATCH_THREADS_LOG
  end=omp_get_wtime();
  fprintf(patch_threadsF, "(thread %d/%d) All (and ONLY) the calls to patches_generation_per_level() took %lf seconds\n", omp_get_thread_num(), omp_get_num_threads(), end-start);
  fclose(patch_threadsF);
  patch_threadsF=NULL;
#endif

}


int64_t patches_generation_per_level(ahf2_patches_t* patches, subcube_t* sc_table_level, int32_t level, uint64_t NminPerHalo){
  subcube_t *sc_aux = NULL, *sc_tmp = NULL;
//  uint64_t sc_count = 0, part_count = 0, sb_part_count = 0, initial_depth_nsubcubes = 0;
  int64_t current_patch_id=0;
  patch_t *ppatch_aux = NULL;
  int n_subcubes = 0;

#ifdef WITH_OPENMP
  double start,end;
  io_logging_msg(global_io.log, INT32_C(2), "(thread %d/%d) Building the patches at level %d", omp_get_thread_num(), omp_get_num_threads(), level);

  start=omp_get_wtime();
#else
  io_logging_msg(global_io.log, INT32_C(2), "Building the patches at level %d", level);
#endif
  if (sc_table_level != NULL) {
    table_iterate_subcube(sc_table_level, sc_aux, sc_tmp){
      //If subcube is not assigned to any patch, we include it in the new patch
      if (sc_aux->patch_id == NO_PATCH) {
        current_patch_id = patches->n_patches[level];
//        #ifdef WITH_OPENMP
//          fprintf(stderr,"[patches_generation_per_level] (thread %d/%d) Building =patch(%d-%ld)=\n",omp_get_thread_num(), omp_get_num_threads(), level,current_patch_id);
//        #else
//          fprintf(stderr,"[patches_generation_per_level] Building =patch(%d-%ld)=\n",level,current_patch_id);
//        #endif
        //If there is no patch allocated, we create a new one
        if (ppatch_aux == NULL) {
          patch_create(&ppatch_aux, current_patch_id, level);
        }
        //If there is one allocated patch, we recycle it with patch_init()
        else {
          patch_init(ppatch_aux, current_patch_id, level);
        }
        //sc_aux->patch_id=ppatch_aux->id; //->subcube.patch_id assigned inside patch_add_psubcube()
        patch_add_psubcube(ppatch_aux, &sc_aux);
        n_subcubes = 1;
        n_subcubes += patch_include_adjacent_subcubes(sc_table_level, sc_aux, ppatch_aux);
        //Sanity check
        if (n_subcubes != patch_get_num_psubcubes(ppatch_aux)) {
          fprintf(stderr, "[%s:%d] ERROR: patch(%2d-%2ld) n_subcubes %d != patch_get_num_psubcubes %lu\n", __FILE__, __LINE__, ppatch_aux->level, ppatch_aux->id, n_subcubes, patch_get_num_psubcubes(ppatch_aux));
        }


        //==============================================================================================
        // REJECTION OF USELESS PATCHES:
        //
        //   we want to avoid using patches stemming from Poisson noise in the particles distribution
        //
        //   problem:  we cannot reject patches based upon the number of subcubes as the patch tree
        //             is built starting from level 0 that only contains 1 subcube, etc.
        //
        //   solution: we know that each subcube contains at least Nth_dom particles;
        //             therefore, a patch of size MIN_NNODES should at least contain
        //             MIN_NNODES*Nth_dom particles
        //             if there are even more particles, it indicates either a large patch or
        //             a very high density patch which should be kept!
        //
        //==============================================================================================
        if (ppatch_aux->Npart >= NminPerHalo) {
          patches->n_patches[level]++;

          // finish all the cumulative physical properties of the just created patch
          finish_patch_physics(ppatch_aux);

          //Realloc for a new patch in patches->tree[depth] array
          patches->tree[level] = (patch_t*) realloc(patches->tree[level], patches->n_patches[level] * SIZEOF_PATCH);
          //Copy the just generated patch in the last position of patches->tree[deep] array
          memcpy(&(patches->tree[level][patches->n_patches[level] - 1]), ppatch_aux, SIZEOF_PATCH);

          //fprintf(generate_treeF,"Included patch(%d-%"PRIu64") in patches->tree\n", patches->tree[depth][patches->n_patches[depth]-1].level, patches->tree[depth][patches->n_patches[depth]-1].id);

          //ppatch_aux->psubcubes must be set to NULL to avoid corruption of valid patch copied to the patches->tree
          //The rest of the patch_t structure will be recycled in patch_init()
          ppatch_aux->psubcubes = NULL;
        }//if(NminPerHalo)
        else {
          patches->n_rejected_patches[level]++;

          /*
           * Flagging rejected subcubes is a waste of time, only interesting for debuging
           //Flag all subcubes of rejected patch with subcube.patch_id=PATCH_REJECTED
           while((ppsubcube_aux=patch_next_psubcube(ppatch_aux,ppsubcube_aux))!=NULL){
           psubcube_aux=*ppsubcube_aux;
           psubcube_aux->patch_id=PATCH_REJECTED;
           }
           */

        }//else
      }//if NO_PATCH
    }//table_iterate loop
  }//if(sc_table_level != NULL)

  //Free aux patch used to build the patches in patches->tree[level]
  if(ppatch_aux){
    patch_free(&ppatch_aux);
  }
  #ifdef WITH_OPENMP
  end=omp_get_wtime();
  io_logging_msg(global_io.log, INT32_C(2), "We built %lu patches (%lu rejected) at level %d in %lf seconds",
          patches->n_patches[level], patches->n_rejected_patches[level], level, end-start);

  #ifdef PATCH_THREADS_LOG
  fprintf(patch_threadsF, "\t(thread %d/%d) We built %lu patches (%lu rejected) at level %d in %lf seconds\n",
      omp_get_thread_num(), omp_get_num_threads(), patches->n_patches[level],
      patches->n_rejected_patches[level], level, end-start);
  fflush(patch_threadsF);
  #endif

  #else
  io_logging_msg(global_io.log, INT32_C(2), "We built %lu patches (%lu rejected) at level %d",
        patches->n_patches[level], patches->n_rejected_patches[level], level);
  #endif

  return patches->n_patches[level];
}

/*
 * Functions to manage the visitable_subcubes_t structure
 * We use memcpy and memcmp functions to avoid problems with 128bits
 * integer implementation and comparisons.
 */
//typedef struct visitable_subcubes {
//  cubekey_t *list;
//  int64_t next;
//  int64_t size;
//} visitable_subcubes_t;

//Allocate memory for visitable_subcubes_t structure and the internal cubekey_t list. Initialization of size and next.
void visitable_subcubes_create(visitable_subcubes_t** ppv_s){
  if((*ppv_s=(visitable_subcubes_t*)malloc(sizeof(visitable_subcubes_t)))==NULL){
    perror("malloc(visitable_subcubes_t): ");
    fprintf(stderr,"[visitable_subcubes_create] malloc for new visitable_subcubes_t structure failed\n");
    exit (EXIT_FAILURE);
  }
  (*ppv_s)->size=VISITABLE_SUBCUBES_BLOCK;
  (*ppv_s)->next=0;
  if(((*ppv_s)->list=(cubekey_t*)malloc(sizeof(cubekey_t)*((*ppv_s)->size)))==NULL){
    perror("malloc(visitable_subcubes_t->list): ");
    fprintf(stderr,"[visitable_subcubes_create] malloc for new visitable_subcubes_t->list array failed\n");
    exit (EXIT_FAILURE);
  }
}

inline void visitable_subcubes_clear(visitable_subcubes_t* pv_s){
  pv_s->next=0;
}

void visitable_subcubes_add(visitable_subcubes_t* pv_s, cubekey_t* pck){
  int64_t i;

  //Check pointer validity
  if(pv_s==NULL){
    fprintf(stderr,"[%s:%d] ERROR: In visitable_subcubes_add() pointer to visitable_subcubes_t is NULL!! Did you call visitable_subcubes_create()??\nAborting...",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  //We check if the cubekey we are about to add is already in the list
  for(i=0;i<pv_s->next;i++){
    //If the cubekey to be added is already stored in the list, we have nothing to do -> return
    if(memcmp(pck,&pv_s->list[i],SIZEOF_CUBEKEY)==0){
      return;
    }
  }

  //If the cubekey ck is not present in list, we add it
  //Check if the list is full (next==size). In that case, increase the list size reallocating memory
  if(pv_s->next==pv_s->size){
    pv_s->size+=VISITABLE_SUBCUBES_BLOCK;
    if((pv_s->list=(cubekey_t*)realloc(pv_s->list, SIZEOF_CUBEKEY*(pv_s->size)))==NULL){
      perror("realloc(visitable_subcubes_t->list): ");
      fprintf(stderr,"[%s:%d] ERROR: In visitable_subcubes_add() realloc of visitable_subcubes_t->list array failed\n. Aborting...", __FILE__, __LINE__);
      exit (EXIT_FAILURE);
    }
  }
  //Once we are sure we have memory available, we insert the cubekey in the list
  memcpy(&pv_s->list[pv_s->next],pck,SIZEOF_CUBEKEY);
  pv_s->next++;
}

//Return the number of valid cubekeys stored in the visitable_subcubes_t structure
inline int64_t visitable_subcubes_count(visitable_subcubes_t* pv_s){
  if(pv_s==NULL) return 0;
  else  return pv_s->next;
}

//Extract the cubekey of the visitable subcube at index position
inline cubekey_t* visitable_subcubes_get(visitable_subcubes_t* pv_s, int64_t index){
  //Check that index is not out of the limit
  if(index>=pv_s->next){
    return NULL;
  }
  return &(pv_s->list[index]);
}

inline void visitable_subcubes_free(visitable_subcubes_t** ppv_s){
  //Check pointer validity
  if(*ppv_s==NULL){
    fprintf(stderr,"[%s:%d] ERROR: In visitable_subcubes_free() pointer to visitable_subcubes_t is NULL!! We cannot free anything...",
            __FILE__, __LINE__);
    return;
  }
  //Free list memory
  if((*ppv_s)->list){
    free((*ppv_s)->list);
    (*ppv_s)->list=NULL;
  }
  //Free structure memory
  free(*ppv_s);
  *ppv_s=NULL;
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
  if(L<0){
    fprintf(stderr,"[%s:%d] ERROR: call to ck2coor() from add_patch_physics() failed for cubekey=%"PRIck"."\
        "\nPhysical results will be WRONG.\n",
        __FILE__, __LINE__, psc->cubekey);
  }
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
      //exit(EXIT_FAILURE);
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
