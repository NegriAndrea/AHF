#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
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
#include "subcube.h"
#include "patch.h"


inline int subcube_refine(subcube_t** sc_table, subcube_t* sc_src, partptr fst_part){
	unsigned int depth;
	int sc_generated=0;
	partptr *pp_part=NULL;     //We do need pointer to pointer to particle: 
                              // UTarray stores simple pointers 'partptr' type, we pass/receive "re-pointers" to/from routines
	cubekey_t cubekey;
	subcube_t *sc_aux=NULL;
  
	depth=ck_get_depth(sc_src->cubekey);
	if(depth==MAX_DEPTH){
		io_logging_msg(global_io.log, INT32_C(1), "We can not refine to a deeper level, MAX_DEPTH(%d) reached.", MAX_DEPTH);
		return 0;
	}
  
	pp_part=NULL;
	while((pp_part=subcube_next_particle(sc_src, pp_part))){
		/*
		 * TODO: implement a function to calculate the new cubekey for next refinement level from the original cubekey.
		 * Currently we recalculate the shifts at every refinement level with coor2ck.
		 * Proposed prototype: ck_refine(subcube_t* ck_new, subcube_t ck_original, partptr part)
		 * Proposed call line: ck_refine(&ck_new, sc_src->cubekey, fst_part+(*part_it));
		 */
		coor2ck(&cubekey, (*pp_part)->pos[0], (*pp_part)->pos[1], (*pp_part)->pos[2], depth+1);
		//Look for the new cubekey subcube
		table_find_subcube(sc_table[depth+1],&cubekey,&sc_aux);
    
		//If the corresponding subcube doesnt exist, we create it and add it to the hash-table
		if(sc_aux==NULL){
			subcube_create(&sc_aux,cubekey);
			table_add_subcube(&sc_table[depth+1], sc_aux);
			sc_generated++;
		}
    
		//Add the particle-ID (partptr; pointer to particle) to the particles array in the subcube structure
		subcube_add_particle(sc_aux,pp_part);
	}
	subcube_free_particles_array(sc_src);
	return sc_generated;
}


int initial_box_division(long unsigned npart, partptr fst_part, int initial_depth, subcube_t **sc_table){
	uint64_t   ipart=0;
	cubekey_t  cubekey;
	subcube_t *sc_aux=NULL;
	partptr    p_part=NULL;
	partptr*   pp_part=NULL; //We do need pointer to pointer to particle: 
                            // UTarray stores simple pointers 'partptr' type, we pass/receive "re-pointers" to/from routines
	
	if(initial_depth>MAX_DEPTH){
		fprintf(stderr,"ERROR: initial_depth (%d) cannot be larger than MAX_DEPTH (%d).\nAborting...\n",initial_depth,MAX_DEPTH);
		exit(-1);
	}
	
	for(ipart=0;ipart<npart; ipart++){
		p_part=fst_part+ipart;
		pp_part=&p_part;
				
		//Calculate cubekey for the particle positions and the given depth
		coor2ck(&cubekey, p_part->pos[0], p_part->pos[1], p_part->pos[2], initial_depth);
		//Look for the cubekey subcube
		table_find_subcube(sc_table[initial_depth],&cubekey,&sc_aux);
		//If the corresponding subcube doesnt exist, we create it and add it to the hash-table
		if(sc_aux==NULL){
			subcube_create(&sc_aux,cubekey);
			table_add_subcube(&sc_table[initial_depth], sc_aux);
		}
		//Add the particle-ID (position) to the particles array in the subcube structue
		subcube_add_particle(sc_aux,pp_part);
	}//for(ipart)
	return table_get_num_subcubes(sc_table[initial_depth]);
}


void	refine_box_division(int initial_depth, subcube_t **sc_table, int threshold, partptr fst_part){
	int depth;
	subcube_t *sc_aux, *sc_tmp;
	
	for(depth=initial_depth;depth<=MAX_DEPTH;depth++){
		sc_aux=NULL;
		sc_tmp=NULL;
    
		//If subcube-table is not empty, iterate all the subcubes at the given depth
		if(sc_table[depth]!=NULL){
			table_iterate_subcube(sc_table[depth], sc_aux, sc_tmp){
				//When # of particles in subcube > NtreeMin --> subcube_refine
				if(subcube_get_num_particles(sc_aux)>threshold){
					subcube_refine(sc_table,sc_aux,fst_part);
				}
			}//table_iterate_subcube
		}//if sc_table not NULL
		else{
			//If the subcube-table is empty, all the deeper levels will be empty as well
			if(depth>initial_depth)	break; //break the for depth loop
		}
	}
}

//generate_tree(global_info.no_part, global_info.fst_part, simu.Nth_dom);

void generate_tree(long unsigned npart, partptr fst_part, double Nth_dom){
	int32_t          NtreeMin, i=0, j=0;
	const int32_t   Nlevels=MAX_DEPTH+1; //We do need Nlevels=MAX_DEPTH+1 as the arrays size because there may exist level 0 and level MAX_DEPTH
	
	//sc_table[i] will be the subcubes hash table at refinement level i. Must be initialized to NULL
	subcube_t       *sc_table[Nlevels];
	//sc_aux is an auxiliar pointer to subcube. sc_tmp is used internally to iterate
	subcube_t       *sc_aux=NULL, *sc_tmp=NULL;
	uint64_t         sc_count=0, part_count=0, sb_part_count=0, initial_depth_nsubcubes=0;
	int32_t          depth=0, initial_depth=0, depth_aux=0;
	
	
	patch_t *patch_tree[Nlevels];
	patch_t *ppatch_aux=NULL;
	int n_patches[Nlevels];
	int n_invalid_patches[Nlevels];
	int n_valid_patches[Nlevels];       //This validity criteria will be applied later!! num_subcubes>=MIN_NNODES
	int n_subcubes_in_patches[Nlevels]; //sum of subcubes in all patches at each level
	int current_patch_id=0, n_subcubes=0;

	FILE* generate_treeF=NULL;
	FILE* patchF=NULL;

#ifdef DUMP_REFINEMENT
	char filename_aux[128]="";
	FILE* partCSV=NULL;
	FILE* subcubesCSV=NULL;
#endif
	
	NtreeMin = (int)(Nth_dom+0.5);
	//uthash pointers must be initialize to NULL before adding first structure
	for (i=0; i<Nlevels; i++)
		sc_table[i]=NULL;
  
	if((generate_treeF=fopen("generate_tree.out", "w"))==NULL){
		perror("fopen(generate_tree.out): ");
		fprintf(stderr, "Error opening generate_tree.out");
		return;
	}
	io_logging_msg(global_io.log, INT32_C(1), "Max particles per subcube: NtreeMin=%d",NtreeMin);
	
	
	initial_depth=ptree_min_level(Nth_dom);
	depth=initial_depth;
	io_logging_msg(global_io.log, INT32_C(1), "DEPTH returned by ptree_min_level = %d",initial_depth);
	fprintf(generate_treeF, "initial_depth returned by ptree_min_level = %d\n",depth);
	
	/*==================================================================================================================
	 * Generate particle-tree
	 *==================================================================================================================*/
	//Initial Box division at ptree_min_level depth
	initial_depth_nsubcubes=initial_box_division(npart,fst_part,initial_depth,sc_table);
	
	fprintf(generate_treeF, "initial_box_division at depth %u generated %lu subcubes\n", initial_depth, initial_depth_nsubcubes);
	io_logging_msg(global_io.log, INT32_C(1), "initial_box_division() at ptree_min_level depth (%u) generated %lu subcubes", initial_depth, initial_depth_nsubcubes);
	
	/*
	 * TO-DO: Move this printing to the AHF log file
	 */
	//Print the current refinement status
	fprintf(generate_treeF,"Refined at initial level %d\n",initial_depth);
  
	//Let's iterate the initial_depth level subcubes counting subcubes and particles
	sc_count=0;
	sc_aux=NULL;
	sc_tmp=NULL;
	sb_part_count=0;
	part_count=0;
	
	table_iterate_subcube(sc_table[initial_depth], sc_aux, sc_tmp){
		sc_count++;
		sb_part_count=subcube_get_num_particles(sc_aux);
		part_count+=sb_part_count;
	}
  
	fprintf(generate_treeF,"Depth level=%2d, sc_count=%lu / %lu, global_no_part=%lu, stored_particles=%7lu\n", depth_aux, sc_count, 1L<<3*initial_depth, global_info.no_part, part_count);
  
  
	//Once the initial refinement at initial_depth level is done, let's split subcubes when num_particles > threshold
	//Doing it in depth order could generate problems when code is parallelized (hopefully soon)
	refine_box_division(initial_depth,sc_table,NtreeMin,fst_part);
  
	//REFINEMENT IS OVER
	/*
	 * TO-DO: Move this printing to the AHF log file
	 */
	fprintf(generate_treeF,"\nREFINEMENT IS OVER\nLet's take a look at the subcubes formation:\n");
	io_logging_msg(global_io.log, INT32_C(1), "REFINEMENT IS OVER");
	for(depth=initial_depth;depth<Nlevels;depth++){
		sc_count=0;
		sc_aux=NULL;
		sc_tmp=NULL;
		sb_part_count=0;
		part_count=0;
		if(sc_table[depth]!=NULL){
			table_iterate_subcube(sc_table[depth], sc_aux, sc_tmp){
				sc_count++;
				sb_part_count=subcube_get_num_particles(sc_aux);
				part_count+=sb_part_count;
			}
		}
		else{
			if(depth>initial_depth)	break;
		}
		
		fprintf(generate_treeF,"Depth level=%2d (%5lu div per dim), sc_count=%7lu/%10lu (HASH_COUNT=%u), global_no_part=%lu, stored_particles=%7lu\n", depth, (uint64_t)1L<<depth, sc_count, (uint64_t)1L<<3*depth, HASH_COUNT(sc_table[depth]), global_info.no_part, part_count);
		io_logging_msg(global_io.log, INT32_C(1), "Depth level=%2d (%5lu div per dim), sc_count=%7lu/%10lu (HASH_COUNT=%u), global_no_part=%lu, stored_particles=%7lu", depth, (uint64_t)1L<<depth, sc_count, (uint64_t)1L<<3*depth, HASH_COUNT(sc_table[depth]), global_info.no_part, part_count);
	}
	
#ifdef DUMP_REFINEMENT
  
	fprintf(stderr,"RE-CHECK subcubes formation\n");
	//For each subcube, check particles position is inside subcube boundaries.
	//When particles_array is empty, check that the sum of particles_arrays' size in subcubes is equal to nparticles
	
	//particles.csv
	sprintf(filename_aux,"particles.csv");
	if((partCSV=fopen(filename_aux, "w"))==NULL){
		perror("fopen(particles.csv): ");
		fprintf(stderr, "Error opening/creating particles.csv");
		exit(-1);
	}
	fprintf(stderr, "Dumping %s file...\n",filename_aux);
	
	//subcubes.csv
	sprintf(filename_aux,"subcubes.csv");
	if((subcubesCSV=fopen(filename_aux, "w"))==NULL){
		perror("fopen(subcubes.csv): ");
		fprintf(stderr, "Error opening/creating subcubes.csv");
		exit(-1);
	}
	fprintf(stderr, "Dumping %s file...\n",filename_aux);
	
	//Print particles.csv header
	fprintf(partCSV,"#posX,posY,posZ,SFCkey,cubekey,level,partid\n");
	
	//Print subcubes.csv header
	fprintf(subcubesCSV,"#posX,posY,posZ,edge,level,cubekey,nparticles,ntotalsubparticles\n");
	
	for(depth=1;depth<Nlevels;depth++){
		sc_count=0;
		sc_aux=NULL;
		sc_tmp=NULL;
		sb_part_count=0;
		part_count=0;
		if(sc_table[depth]!=NULL){
			table_iterate_subcube(sc_table[depth], sc_aux, sc_tmp){
				sc_count++;
				sb_part_count=subcube_get_num_particles(sc_aux);
				part_count+=sb_part_count;
				pp_part_iter=NULL;
				ck2coor(&Xaux,&Yaux,&Zaux,&edge,sc_aux->cubekey);
				level=ck_get_depth(sc_aux->cubekey);
				fprintf(subcubesCSV,"%e,%e,%e,%e,%d,%lu,%lu,%lu\n",Xaux,Yaux,Zaux,edge,level,sc_aux->cubekey,sb_part_count,sc_aux->nparticles);
				while((pp_part_iter=subcube_next_particle(sc_aux, pp_part_iter))!=NULL){
					p_part_aux=(*pp_part_iter);
					fprintf(partCSV,"%e,%e,%e,%lu,%lu,%d,%lu\n", p_part_aux->pos[0], p_part_aux->pos[1], p_part_aux->pos[2], \
                  p_part_aux->sfckey, sc_aux->cubekey, level, p_part_aux->id);
				}
			}
		}
		else{
			if(depth>initial_depth)	break;
		}
	}//for(depth)
	fclose(subcubesCSV);
	subcubesCSV=NULL;
	fclose(partCSV);
	partCSV=NULL;
	
#endif //DUMP_REFINEMENT
	
	/*==================================================================================================================
	 * Generate patch-tree
	 *==================================================================================================================*/
	io_logging_msg(global_io.log, INT32_C(1), "Patch-tree building");
	
	if((patchF=fopen("patch.out", "w"))==NULL){
		perror("fopen(patch.out): ");
		fprintf(stderr, "Error opening patch.out");
		return;
	}
	//Initialize patch_tree arrays and n_patches counters
	for(i=0;i<Nlevels;i++){
		patch_tree[i]=NULL;
		n_patches[i]=0;
		n_invalid_patches[i]=0;
		n_subcubes_in_patches[i]=0;
	}
	
	/*
	 * Generate patches at every level [PATCH-TREE]
	 */
	for(depth=1;depth<Nlevels;depth++){
		sc_count=0;
		sc_aux=NULL;
		sc_tmp=NULL;
		sb_part_count=0;
		part_count=0;
		current_patch_id=0;
		
		if(sc_table[depth]!=NULL){
			table_iterate_subcube(sc_table[depth], sc_aux, sc_tmp){
				//If subcube is not assigned to any patch, we include it in the new patch
				if(sc_aux->patch_id==NO_PATCH){
					current_patch_id=n_patches[depth];
          
					//Realloc for a new patch
					//patch_tree[depth]=(patch_t*)realloc(patch_tree[depth],n_patches[depth]*sizeof(patch_t));
					if(ppatch_aux==NULL){
						patch_create(&ppatch_aux,current_patch_id,depth);
					}
					else{
						patch_init(ppatch_aux, current_patch_id, depth);
					}
          
					patch_add_psubcube(ppatch_aux,&sc_aux);
					n_subcubes=1;
					n_subcubes+=patch_include_adjacent_subcubes(sc_table[depth],sc_aux,ppatch_aux);
					//TO-DO: If, and only if, n_subcubes >= MIN_NNODES, add new patch to patch_tree
					if(n_subcubes!=patch_get_num_psubcubes(ppatch_aux)){
						fprintf(patchF,"WTF!! n_subcubes %d != patch_get_num_psubcubes %d\n",n_subcubes, patch_get_num_psubcubes(ppatch_aux));
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
					if(ppatch_aux->Npart >= MIN_NNODES*Nth_dom)
           {
            //Remove MIN_NNODES condition here. Re-think about patch_create/init/clear/whatever...
            fprintf(patchF, "patch at level %d, patch_id %d, n_subcubes %d\n",depth, current_patch_id,n_subcubes);
            n_patches[depth]++;
            
            // finish all the cumulative physical properties of the just created patch
            finish_patch_physics(ppatch_aux);
            
            //Realloc for a new patch
            patch_tree[depth]=(patch_t*)realloc(patch_tree[depth],n_patches[depth]*SIZEOF_PATCH);
            memcpy(&(patch_tree[depth][n_patches[depth]-1]),ppatch_aux,SIZEOF_PATCH);
            
            //ppatch_aux must be set to NULL to avoid corruption of valid patch copied to the patch_tree
            ppatch_aux=NULL;
           }
          else
           {
           }
				}//if NO_PATCH
			}//table_iterate
		}//if(sc_table[depth]!=NULL)
		else{
			fprintf(patchF,"Emtpy subcube table at depth %d (sc_table[%d])\n", depth, depth);
			if(depth>initial_depth)	break;
		}
	}//for depth
	
	/*
	 * Check and print the patch_tree formation
	 */
	fprintf(patchF,"LET'S ITERATE THE TREE\n");
	for(depth=initial_depth;depth<Nlevels;depth++){
		if(patch_tree[depth]==NULL){
			if(n_invalid_patches[depth]==0){
				break;
			}
			fprintf(patchF,"=> Level %2d, patch_tree[%d]==NULL . %d invalid patches\n",depth,depth,n_invalid_patches[depth]);
			continue;
		}
		if(n_patches[depth]==0){
			fprintf(patchF,"WTF?? patch_tree[%d]!=NULL and n_patches[%d]=%d (Level %d)\n",depth,depth,n_patches[depth],depth);
			continue;
		}
    
		n_valid_patches[depth]=0;
		for(j=0;j<n_patches[depth];j++){
				n_subcubes_in_patches[depth]+=patch_get_num_psubcubes(&(patch_tree[depth][j]));
				n_valid_patches[depth]++;
		}
		fprintf(patchF,"====> Level %2d, %d valid patches (%d stored patches, %d stored subcubes)\n",depth,n_patches[depth],n_valid_patches[depth],n_subcubes_in_patches[depth]);
	}
  
	
	
	/*
	 * TODO:
	 * Freeing hash-tables, particle-arrays and patches doesn't hurt. DO IT!!
	 */
  
	fclose(generate_treeF);
	generate_treeF=NULL;
  fclose(patchF);
	patchF=NULL;

	return;
}
