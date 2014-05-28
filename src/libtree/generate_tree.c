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
#include "set_patch_trunk.h"
#include "set_patch_radii.h"
#include "set_patch_pos.h"

//Internal functions declaration
int subcube_refine(subcube_t**, subcube_t*, partptr);
int initial_box_division(long unsigned, partptr, int, subcube_t**);
void refine_box_division(int, subcube_t**, int, partptr);
ahf2_patches_t* generate_tree(long unsigned, partptr, double, int);

inline int subcube_refine(subcube_t** sc_table, subcube_t* sc_src, partptr fst_part) {
  unsigned int depth;
  int sc_generated=0;
  partptr *pp_part=NULL;     //We do need pointer to pointer to particle:
                             // UTarray stores simple pointers 'partptr' type, we pass/receive "re-pointers" to/from routines
  cubekey_t cubekey;
  subcube_t *sc_aux=NULL;

  depth=ck_get_depth(sc_src->cubekey);
  if (depth==MAX_DEPTH){
    io_logging_msg(global_io.log, INT32_C(1), "We can not refine to a deeper level, MAX_DEPTH(%d) reached.", MAX_DEPTH);
    return 0;
  }

  pp_part=NULL;
  while ( (pp_part=subcube_next_particle(sc_src, pp_part))){
    /*
     * TODO: implement a function to calculate the new cubekey for next refinement level from the original cubekey.
     * Currently we recalculate the shifts at every refinement level with coor2ck.
     * Proposed prototype: ck_refine(subcube_t* ck_new, subcube_t ck_original, partptr part)
     * Proposed call line: ck_refine(&ck_new, sc_src->cubekey, fst_part+(*part_it));
     */
    coor2ck(&cubekey, (*pp_part)->pos[0], (*pp_part)->pos[1], (*pp_part)->pos[2], depth+1);
    //Look for the new cubekey subcube
    table_find_subcube(sc_table[depth+1], &cubekey, &sc_aux);

    //If the corresponding subcube doesnt exist, we create it and add it to the hash-table
    if (sc_aux==NULL){
      subcube_create(&sc_aux, cubekey);
      table_add_subcube(&sc_table[depth+1], sc_aux);
      sc_generated++;
    }

    //Add the particle-ID (partptr; pointer to particle) to the particles array in the subcube structure
    subcube_add_particle(sc_aux, pp_part);
  }
  subcube_free_particles_array(sc_src);
  return sc_generated;
}

int initial_box_division(long unsigned npart, partptr fst_part, int initial_depth, subcube_t **sc_table) {
  uint64_t ipart=0;
  cubekey_t cubekey;
  subcube_t *sc_aux=NULL;
  partptr p_part=NULL;
  partptr* pp_part=NULL; //We do need pointer to pointer to particle:
                         // UTarray stores simple pointers 'partptr' type, we pass/receive "re-pointers" to/from routines

  if(initial_depth>MAX_DEPTH){
    fprintf(stderr, "ERROR: initial_depth (%d) cannot be larger than MAX_DEPTH (%d).\nAborting...\n", initial_depth, MAX_DEPTH);
    exit(-1);
  }

  for (ipart=0; ipart<npart; ipart++){
    p_part=fst_part+ipart;
    pp_part=&p_part;

    //Calculate cubekey for the particle positions and the given depth
    coor2ck(&cubekey, p_part->pos[0], p_part->pos[1], p_part->pos[2], initial_depth);
    //Look for the cubekey subcube
    table_find_subcube(sc_table[initial_depth], &cubekey, &sc_aux);
    //If the corresponding subcube doesnt exist, we create it and add it to the hash-table
    if (sc_aux==NULL){
      subcube_create(&sc_aux, cubekey);
      table_add_subcube(&sc_table[initial_depth], sc_aux);
    }
    //Add the particle-ID (position) to the particles array in the subcube structue
    subcube_add_particle(sc_aux, pp_part);
  }                            //for(ipart)
  return table_get_num_subcubes(sc_table[initial_depth]);
}

void refine_box_division(int initial_depth, subcube_t **sc_table, int threshold, partptr fst_part) {
  int depth;
  subcube_t *sc_aux, *sc_tmp;

  for (depth=initial_depth; depth<NLEVELS; depth++){
    sc_aux=NULL;
    sc_tmp=NULL;

    //If subcube-table is not empty, iterate all the subcubes at the given depth
    if (sc_table[depth]!=NULL){
      table_iterate_subcube(sc_table[depth], sc_aux, sc_tmp){
        //When # of particles in subcube > NtreeMin --> subcube_refine
        if (subcube_get_num_particles(sc_aux)>threshold){
          subcube_refine(sc_table, sc_aux, fst_part);
        }
      }                            //table_iterate_subcube
    }                            //if sc_table not NULL
    else{
      //If the subcube-table is empty, all the deeper levels will be empty as well
      if (depth>initial_depth) break; //break the for depth loop
    }
  }
}

//generate_tree(global_info.no_part, global_info.fst_part, simu.Nth_dom);

ahf2_patches_t* generate_tree(long unsigned npart, partptr fst_part, double Nth_dom, int NminPerHalo) {
  int32_t NtreeMin=0, i=0, it_level=0, it_patch=0;

  //sc_table[i] will be the subcubes hash table at refinement level i. Must be initialized to NULL
  subcube_t *sc_table[NLEVELS];
  //sc_aux is an auxiliar pointer to subcube. sc_tmp is used internally to iterate
  subcube_t *sc_aux=NULL, *sc_tmp=NULL;
  uint64_t sc_count=0, part_count=0, initial_depth_nsubcubes=0, sb_part_count=0;
  int32_t depth=0, initial_depth=0, final_depth=0;

  ahf2_patches_t* patches=NULL;

  clock_t t_start, t_end;       //Timing variables for POSIX Clock
  double elapsed;

#ifdef GENERATE_TREE_LOG
  FILE* generate_treeF=NULL;
  char generate_tree_name[256];
#endif
  FILE* patchtreeF=NULL;

#ifdef DUMP_TREES
  char filename_aux[128]="";
  FILE* partCSV=NULL;
  FILE* subcubesCSV=NULL;
  FILE* patchesCSV=NULL;
#endif

  NtreeMin=(int) (Nth_dom+0.5);
  //UThash pointers must be initialize to NULL before adding first structure
  for (i=0; i<NLEVELS; i++)
    sc_table[i]=NULL;

#ifdef GENERATE_TREE_LOG
  strncpy(generate_tree_name,global_io.params->outfile_prefix,256);
  strcat(generate_tree_name, ".generate_tree.log");
  if((generate_treeF=fopen(generate_tree_name, "w"))==NULL){
    perror("fopen(generate_tree.out): ");
    fprintf(stderr, "Error opening generate_tree.out\n");
    return NULL;
  }
  fprintf(generate_treeF, "Memory size of basic structures:\n");
  fprintf(generate_treeF, "\tsize of particle ID: %lu bytes\n", sizeof(part_id_t));
  fprintf(generate_treeF, "\tsize of cubekey: %lu bytes\n", sizeof(cubekey_t));
  fprintf(generate_treeF, "\tsize of subcube structure: %lu bytes\n", sizeof(subcube_t));
  fprintf(generate_treeF, "\tsize of patch structure: %lu bytes\n\n", sizeof(patch_t));
#endif

  io_logging_msg(global_io.log, INT32_C(1), "Max particles per subcube: NtreeMin=%d", NtreeMin);

  initial_depth=ptree_min_level(Nth_dom);
  depth=initial_depth;
  io_logging_msg(global_io.log, INT32_C(1), "DEPTH returned by ptree_min_level = %d", initial_depth);
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "initial_depth returned by ptree_min_level = %d\n",depth);
#endif

  io_logging_msg(global_io.log, INT32_C(1), "Initial box division...");

#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "Initial box division... "); fflush(generate_treeF);
#endif

  /*==========================================================
   * We want to test the UThash capacity to store more than INT_MAX (limits.h) elements.
   * Let's try with 2^33 elements (8.589.934.592)...
   * Uncomment block of code below to run the test
   ==========================================================*/
  /***************************************************************
   fprintf(stderr, "Let's try to insert %lu elements in a UThash table...\n",(uint64_t)2L<<33);

   for(sc_count=0;sc_count<(uint64_t)2<<33;sc_count++){
   if(sc_count%10000000 == 0){
   fprintf(stderr, "%lu subcubes inserted\n",table_get_num_subcubes(sc_table[initial_depth]));
   }
   sc_aux=NULL;
   subcube_create(&sc_aux,sc_count);
   table_add_subcube(&sc_table[initial_depth], sc_aux);
   }

   for(sc_count=0;sc_count<(uint64_t)2<<33;sc_count+=(100000-1)){
   sc_aux=NULL;
   table_find_subcube(sc_table[initial_depth],&sc_count,&sc_aux);
   if(sc_aux==NULL){
   fprintf(stderr, "ERROR: Subcube not found for id %lu\n",sc_count);
   }
   else{
   fprintf(stderr, "FOUND subcube for id %lu\n",sc_count);
   }
   }
   exit(0);
   ***************************************************************/

  /*==================================================================================================================
   * Generate particle-tree
   *==================================================================================================================*/
  //Initial Box division at ptree_min_level depth
  t_start=clock();
  initial_depth_nsubcubes=initial_box_division(npart, fst_part, initial_depth, sc_table);
  t_end=clock();
  elapsed= ((double) (t_end-t_start))/CLOCKS_PER_SEC;
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "DONE\n");
  fflush(generate_treeF);
#endif

#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "\t=TIMING= initial_box_division() took %f seconds \n", elapsed);
  fprintf(generate_treeF, "initial_box_division at depth %u generated %lu subcubes\n", initial_depth, initial_depth_nsubcubes);
  fflush(generate_treeF);
#endif

  io_logging_msg(global_io.log, INT32_C(1), "initial_box_division() at ptree_min_level depth (%u) generated %lu subcubes", initial_depth, initial_depth_nsubcubes);
  io_logging_msg(global_io.log, INT32_C(1), "=timing= initial_box_division() took %f seconds", elapsed);

#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "Refine recursively... ");
  fflush(generate_treeF);
#endif

  io_logging_msg(global_io.log, INT32_C(1), "Refine recursively...");
  //Once the initial refinement at initial_depth level is done, let's split subcubes when num_particles > threshold
  //Doing it in depth order could generate problems when code is parallelized (hopefully soon)

  t_start=clock();
  refine_box_division(initial_depth, sc_table, NtreeMin, fst_part);
  t_end=clock();
  elapsed= ((double) (t_end-t_start))/CLOCKS_PER_SEC;

  //Calculate the final_depth level where there exist subcubes
  final_depth=-1;
  for (i=MAX_DEPTH; i>=initial_depth; i--){
    if (sc_table[i]!=NULL){
      final_depth=i;
      break;
    }
  }  //for
  if (final_depth==-1){
    fprintf(stderr, "[%s:%d] ERROR: Impossible to find the final_depth\n. Aborting...\n", __FILE__, __LINE__);
    exit(-1);
  }

  //REFINEMENT IS OVER
  io_logging_msg(global_io.log, INT32_C(1), "REFINEMENT IS OVER (initial_depth=%d, final_depth=%d).", initial_depth, final_depth);
  io_logging_msg(global_io.log, INT32_C(1), "=timing= refine_box_division() took %f seconds", elapsed);

#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "DONE\nREFINEMENT IS OVER (initial_depth=%d, final_depth=%d).\n", initial_depth, final_depth);
  fprintf(generate_treeF, "\t=TIMING= refine_box_division() took %f seconds \n", elapsed);
  fprintf(generate_treeF,"Subcubes formation:\n");
  fflush(generate_treeF);
#endif

  //Review Subcubes formation
  t_start=clock();
  for (depth=initial_depth; depth<NLEVELS; depth++){
    sc_count=0;
    sc_aux=NULL;
    sc_tmp=NULL;
    sb_part_count=0;
    part_count=0;
    if (sc_table[depth]!=NULL){
      table_iterate_subcube(sc_table[depth], sc_aux, sc_tmp)
      {
        sc_count++;
        //sb_part_count=subcube_get_num_particles(sc_aux);
        sb_part_count=subcube_get_num_stored_particles(sc_aux);
        part_count+=sb_part_count;
      }
    }  //if(sc_table[depth]!=NULL)
    else{
      if (depth>initial_depth){
        if (depth!=final_depth+1){
          fprintf(stderr, "[__FILE__:__LINE__] ERROR: We should exit from the loop when depth == final_depth+1, but depth=%d and final_depth=%d\n", depth, final_depth);
        }
        break;
      }
    }  //else

#ifdef GENERATE_TREE_LOG
    fprintf(generate_treeF,
        "\tDepth level=%2d (%10lu div per dim), sc_count=%10lu (HASH_COUNT=%10u), global_no_part=%10lu, stored_particles=%10lu\n",
        depth, (uint64_t) 1L << depth, sc_count, HASH_COUNT(sc_table[depth]),
        global_info.no_part, part_count);
#endif
    io_logging_msg(global_io.log, INT32_C(1), "Depth level=%2d (%10lu div per dim), sc_count=%10lu (HASH_COUNT=%10u), global_no_part=%10lu, stored_particles=%10lu", depth, (uint64_t)1L<<depth,
        sc_count, HASH_COUNT(sc_table[depth]), global_info.no_part, part_count);
  }  //for
  t_end=clock();
  elapsed= ((double) (t_end-t_start))/CLOCKS_PER_SEC;

  io_logging_msg(global_io.log, INT32_C(1), "=timing= Review Subcubes formation took %f seconds", elapsed);
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "\t=TIMING= Review Subcubes formation took %f seconds \n", elapsed);
  fflush(generate_treeF);
#endif

  /*==================================================================================================================
   * Generate patch-tree
   *==================================================================================================================*/
  io_logging_msg(global_io.log, INT32_C(1), "Building the patches...");
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "Building the patches... ");
  fflush(generate_treeF);
#endif

  t_start=clock();
  patches_generation(&patches, sc_table, initial_depth, final_depth, NminPerHalo);
  t_end=clock();
  elapsed= ((double) (t_end-t_start))/CLOCKS_PER_SEC;

  io_logging_msg(global_io.log, INT32_C(1), "=timing= patches_generation() took %f seconds", elapsed);
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "DONE\n");
  fprintf(generate_treeF, "\t=TIMING= patches_generation() took %f seconds \n", elapsed);
  fflush(generate_treeF);
#endif

  //patch.out file generation
  //fprintf(generate_treeF,"Let's review the patch_tree formation -not connected yet- with patch_tree_review() (generating patch.out log file...)\n"); fflush(generate_treeF);
  //patch_formation_review(patches);

#ifdef DUMP_TREES
  fprintf(stderr,"DUMP_TREES ENABLED: Generating CSV files for particles, subcubes and patches\n");
  //For each subcube, check particle position is inside subcube boundaries.
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

  //patches.csv
  sprintf(filename_aux,"patches.csv");
  if((patchesCSV=fopen(filename_aux, "w"))==NULL){
    perror("fopen(patches.csv): ");
    fprintf(stderr, "Error opening/creating patches.csv");
    exit(-1);
  }
  fprintf(stderr, "Dumping %s file...\n",filename_aux);

  //Print particles.csv header
  fprintf(partCSV,"#posX,posY,posZ,SFCkey,cubekey,level,partid\n");

  //Print subcubes.csv header
  fprintf(subcubesCSV,"#posX,posY,posZ,edge,level,cubekey,nparticles,ntotalsubparticles\n");

  for(depth=1;depth<NLEVELS;depth++){
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
          fprintf(partCSV,"%e,%e,%e,%lu,%lu,%d,%lu\n", p_part_aux->pos[0], p_part_aux->pos[1], p_part_aux->pos[2],
              p_part_aux->sfckey, sc_aux->cubekey, level, p_part_aux->id);
        }
      }
    }
    else{
      if(depth>initial_depth) break;
    }
  }	//for(depth)
  fclose(subcubesCSV);
  subcubesCSV=NULL;
  fclose(partCSV);
  partCSV=NULL;

#endif //DUMP_TREES

  //Connect parents and daughters...
  io_logging_msg(global_io.log, INT32_C(1), "Connecting patches between parents and daughters");
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "Connecting patches between parents and daughters... ");
  fflush(generate_treeF);
#endif

  t_start=clock();
  patch_connect_tree(patches->tree, sc_table, patches->n_patches);
  t_end=clock();
  elapsed= ((double) (t_end-t_start))/CLOCKS_PER_SEC;

  io_logging_msg(global_io.log, INT32_C(1), "Patches->tree connection has finished.");
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "DONE\n");
  fflush(generate_treeF);
#endif

  io_logging_msg(global_io.log, INT32_C(1), "=timing= patches_connect_tree() took %f seconds", elapsed);
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "\t=TIMING= patch_connect_tree() took %f seconds \n", elapsed);
  fflush(generate_treeF);
#endif

  //=======================================================================================
  // prepare things for patchtree2halos() of libahf2.a
  // (do not change the ordering, i.e. pos need to be set *before* trunk needs to be set *before* radii!)
  //=======================================================================================
  t_start=clock();
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "Calling to set_patch_pos()\n");
  fflush(generate_treeF);
#endif
  set_patch_pos(patches->tree, patches->n_patches);

#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "Calling to set_patch_trunk()\n");
  fflush(generate_treeF);
#endif
  set_patch_trunk(patches->tree, patches->n_patches);

#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "Calling to set_patch_radii()\n");
  fflush(generate_treeF);
#endif
  set_patch_radii(patches->tree, patches->n_patches);
  t_end=clock();
  elapsed= ((double) (t_end-t_start))/CLOCKS_PER_SEC;
  io_logging_msg(global_io.log, INT32_C(1), "=timing= set_patch_pos/trunk/radii took %f seconds", elapsed);
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "\t=TIMING= set_patch_pos/trunk/radii took %f seconds\n", elapsed);
  fflush(generate_treeF);
#endif

#ifdef DEBUG_AHF2libtree
  {
    FILE *fp;
    int initial_depth, final_depth, ilevel, ipatch;
    partptr cur_part;
    float r,g,b;

    fp = fopen("patchtree.geom","w");

    // determine initial level
    initial_depth = 0;
    while((patches->tree)[initial_depth] == NULL && initial_depth < MAX_DEPTH){
      initial_depth++;
    }
    // determine final level
    final_depth = initial_depth;
    while((patches->tree)[final_depth] != NULL && final_depth < MAX_DEPTH){
      final_depth++;
    }

    // add the simulation box to the geom file
    fprintf(fp,"l   0 0 0   1 0 0      0 0 1\n");
    fprintf(fp,"l   1 0 0   1 1 0      0 0 1\n");
    fprintf(fp,"l   1 1 0   0 1 0      0 0 1\n");
    fprintf(fp,"l   0 1 0   0 0 0      0 0 1\n");
    fprintf(fp,"l   0 0 1   1 0 1      0 0 1\n");
    fprintf(fp,"l   1 0 1   1 1 1      0 0 1\n");
    fprintf(fp,"l   1 1 1   0 1 1      0 0 1\n");
    fprintf(fp,"l   0 1 1   0 0 1      0 0 1\n");
    fprintf(fp,"l   0 0 0   0 0 1      0 0 1\n");
    fprintf(fp,"l   0 1 0   0 1 1      0 0 1\n");
    fprintf(fp,"l   1 0 0   1 0 1      0 0 1\n");
    fprintf(fp,"l   1 1 0   1 1 1      0 0 1\n");

    // add all patches
    for(ilevel=initial_depth; ilevel<final_depth; ilevel++){
      for(ipatch=0; ipatch<(patches->n_patches)[ilevel]; ipatch++){

        r = sqrt((float)(ilevel-initial_depth)/(float)(final_depth-initial_depth-1));
        g = 0.0;
        b = 1.0-(float)(ilevel-initial_depth)/(float)(final_depth-initial_depth-1);

        fprintf(fp,"s %f %f %f %f %f %f %f   %d %d  %lf %lf %lf %lf %lf %lf %"PRIu64"\n",
            (float)(patches->tree)[ilevel][ipatch].pos[0],
            (float)(patches->tree)[ilevel][ipatch].pos[1],
            (float)(patches->tree)[ilevel][ipatch].pos[2],
            //0.01,
            (float)(patches->tree)[ilevel][ipatch].radius,
            r,g,b,
            ilevel,ipatch,
            (patches->tree)[ilevel][ipatch].Xmin,(patches->tree)[ilevel][ipatch].Xmax,
            (patches->tree)[ilevel][ipatch].Ymin,(patches->tree)[ilevel][ipatch].Ymax,
            (patches->tree)[ilevel][ipatch].Zmin,(patches->tree)[ilevel][ipatch].Zmax,
            (patches->tree)[ilevel][ipatch].Npart
        );

      } // ipatch
    } // ilevel

    // add all particles
    for(cur_part=fst_part; cur_part<fst_part+npart; cur_part++){
      fprintf(fp,"p %f %f %f 0 1 0\n",
          (float)cur_part->pos[0],
          (float)cur_part->pos[1],
          (float)cur_part->pos[2]);
    }

    fclose(fp);
    exit(0);
  }
#endif

  /*===================================================================
   * FREE MEMORY OF ALL SUBCUBES, and the internal particles in UTarray's,
   * from the subcubes' hash tables. Hash table itself is sc_table[] structures,
   * locally declared inside this function.
   * We do also set patch->psubcubes pointers to NULL to avoid later access
   * to freed memory
   *===================================================================*/
#ifdef GENERATE_TREE_LOG
  fprintf(generate_treeF, "Freeing hash tables of subcubes, and arrays of pointers to subcubes from patches\n");
  fflush(generate_treeF);
#endif

  t_start=clock();
  for (it_level=0; it_level<NLEVELS; it_level++){
    if (sc_table[it_level]==NULL) continue;
    sc_aux=NULL;
    sc_tmp=NULL;
#ifdef GENERATE_TREE_LOG
    //Print count and overhead per table
    fprintf(generate_treeF, " -Table at level=%2d: COUNT=%10lu subcubes (%10lu bytes), OVERHEAD=%10lu bytes [TOTAL=%11lu]\n",
        it_level, table_get_num_subcubes(sc_table[it_level]),table_get_num_subcubes(sc_table[it_level])*sizeof(subcube_t),
        table_get_overhead(sc_table[it_level]),
        table_get_num_subcubes(sc_table[it_level])*sizeof(subcube_t)+table_get_overhead(sc_table[it_level]));
    fflush(generate_treeF);
#endif		
    table_iterate_subcube(sc_table[it_level],sc_aux,sc_tmp)
    {
      table_del_subcube(sc_table[it_level], sc_aux);
      subcube_free(&sc_aux);
    } //table_iterate_subcube
    if (patches->tree[it_level]!=NULL){
      for (it_patch=0; it_patch<patches->n_patches[it_level]; it_patch++){
        //Free and Clear pointers to subcubes (or UTarray's of subcubes) from inside the patches
        patch_free_psubcubes_array(&patches->tree[it_level][it_patch]);
        patches->tree[it_level][it_patch].most_particles_subcube=NULL;
      } //for(it_patch)
    } //if
  } //for(it_level)

  t_end=clock();
  elapsed= ((double) (t_end-t_start))/CLOCKS_PER_SEC;

  io_logging_msg(global_io.log, INT32_C(1), "=timing= Freeing memory of UThash tables of subcubes and internal UTarrays of particles took %f seconds", elapsed);
#ifdef GENERATE_TREE_LOG
  //Print count and overhead per table
  fprintf(generate_treeF, "\t=TIMING= Freeing memory of UThash tables of subcubes and internal UTarrays of particles took %f seconds\n", elapsed);
  fflush(generate_treeF);
#endif

  if (patchtreeF){
    fclose(patchtreeF);
    patchtreeF=NULL;
  }

#ifdef GENERATE_TREE_LOG
  if (generate_treeF){
    fclose(generate_treeF);
    generate_treeF = NULL;
  }
#endif

  return patches;
}
