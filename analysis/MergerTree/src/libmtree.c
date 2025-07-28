#include "include.h"
#include "common.h"
#include "libmerit.h"
#include "libsnapskipping.h"

/*==================================================================================================
 * particle_halo_mapping:
 *
 *  for each particle remember to which halo(s) it belongs
 *
 *==================================================================================================*/
int particle_halo_mapping(int isimu)
{
  int64_t  ihalo;         // the downwards for-loop is running until ihalo=-1
  uint64_t ipart, jpart, *position;
  clock_t  elapsed;
  
  elapsed = clock();
  
  PidMax_global = MAX(PidMax[0],PidMax_global);
  PidMax_global = MAX(PidMax[1],PidMax_global);
  fprintf(stderr,"  o creating particle<->halo mapping for file %d (PidMax=%"PRIu64", PidMax_global=%"PRIu64") ... ",isimu,PidMax[isimu],PidMax_global);
  
  parts[isimu] = (PARTptr) calloc((PidMax_global+1), sizeof(PARTS)); // +1 because we are accessing an array like [PidMax_global] !
  if(parts[isimu] == NULL) {
    fprintf(stderr,"\nCould not allocate memory for parts[] array of isimu = %d (size = %f GB)\nAborting\n",
            isimu,(float)((PidMax_global+1)*sizeof(PARTS)/1024./1024./1024.));
    exit(0);
  }
  
  /* recording every halo it belongs to: running from low to high mass objects to allow for unique assignment! */
#ifdef WITH_OPENMP_EXCLUSIVE // parallelisation not tested!
#pragma omp parallel for schedule(dynamic) private(ihalo,jpart,ipart) shared(nHalos,halos,parts,isimu)
#endif
  for(ihalo=nHalos[isimu]-1; ihalo>=0; ihalo--)
  {
    for(jpart=0; jpart<halos[isimu][ihalo].npart; jpart++)
    {
#ifdef USE_PIDMAP
      //  find halos[isimu][ihalo].Pid[jpart] in PidMap[isimu] to get the 'ipart' to be used with parts[isimu][ipart]
      position = (uint64_t *)bsearch(&(halos[isimu][ihalo].Pid[jpart]), PidMap[isimu], NPids[isimu], sizeof(uint64_t), pid_cmp);
      if(position != NULL) {
        ipart = (uint64_t) (position - PidMap[isimu]);
        
#ifdef EXCLUSIVE_PARTICLES
        if(parts[isimu][ipart].nhalos == 0)
#endif
        {
          parts[isimu][ipart].nhalos++;
          parts[isimu][ipart].Hid = (uint64_t *) realloc(parts[isimu][ipart].Hid, parts[isimu][ipart].nhalos*sizeof(uint64_t)); // valgrind says "16 bytes lost"
          if(parts[isimu][ipart].Hid == NULL) {
            fprintf(stderr,"particle_halo_mapping(): could not realloc() memory for Hid array\n");
            exit(0);
          }
          
          parts[isimu][ipart].Hid[parts[isimu][ipart].nhalos-1] = ihalo;
        }
      }
      else {
        // this particle must come from some previous simulation and does not exist here...
      }
#else // USE_PIDMAP
      // simply use the Pid as the 'ipart'
      ipart = halos[isimu][ihalo].Pid[jpart];

#ifdef EXCLUSIVE_PARTICLES
      if(parts[isimu][ipart].nhalos == 0)
#endif
      {
        parts[isimu][ipart].nhalos++;
        parts[isimu][ipart].Hid = (uint64_t *) realloc(parts[isimu][ipart].Hid, parts[isimu][ipart].nhalos*sizeof(uint64_t)); // valgrind says "16 bytes lost"
        if(parts[isimu][ipart].Hid == NULL) {
          fprintf(stderr,"particle_halo_mapping(): could not realloc() memory for Hid array\n");
          exit(0);
        }
        
        parts[isimu][ipart].Hid[parts[isimu][ipart].nhalos-1] = ihalo;
      } // if(nhalos == 0)
      
#endif // USE_PIDMAP
      
      
    } // for(jpart)
  } // for(ihalo)
  
  elapsed = clock()-elapsed;
  fprintf(stderr,"done in %4.2f sec.\n",(float)elapsed/CLOCKS_PER_SEC);
  
  return(1);
}

/*==================================================================================================
 * cross_correlation:
 *
 *  for each halo at isimu=0 figure out how many particles are in common with khalo at isimu=1
 *
 *==================================================================================================*/
int cross_correlation(char OutFile[MAXSTRING])
{
  uint64_t  ihalo;
  clock_t  elapsed;
  
  
  /*---------------------------------------------------------
   * backwards correlation
   *---------------------------------------------------------*/
  elapsed = clock();
  fprintf(stderr,"  o generating cross-correlation 0->1 for %"PRIu64" haloes ...",nHalos[0]);
#ifdef WITH_OPENMP
#  pragma omp parallel for schedule(dynamic) shared(nHalos, stderr) private(ihalo) default(none)
#endif
  /* cross-correlation simu0->simu1 */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    create_mtree(ihalo, 0, 1);
  }
  elapsed = clock()-elapsed;
  fprintf(stderr," done in %4.2f sec.\n", (double)elapsed/CLOCKS_PER_SEC);
  
  
  
  
#ifdef MTREE_BOTH_WAYS
  
  /*---------------------------------------------------------
   * forward correlation
   *---------------------------------------------------------*/
  elapsed = clock();
  fprintf(stderr,"  o generating cross-correlation 1->0 for %"PRIu64" haloes ...",nHalos[1]);
#ifdef WITH_OPENMP
#  pragma omp parallel for schedule(dynamic) shared(nHalos, stderr) private(ihalo) default(none)
#endif
  /* cross-correlation simu0<-simu1 */
  for(ihalo=0; ihalo<nHalos[1]; ihalo++) {
    create_mtree(ihalo, 1, 0);
  }
  elapsed = clock()-elapsed;
  fprintf(stderr," done in %4.2f sec.\n", (double)elapsed/CLOCKS_PER_SEC);
  
  
  
  elapsed = clock();
  fprintf(stderr,"  o removing network connections ...");
#ifdef WITH_OPENMP_CLEAN // not tested/checked properly, but this part is not taking up any time anyways...
#  pragma omp parallel for schedule (dynamic)	shared(nHalos) private(ihalo)
#endif
  /* clean connections simu0->simu1 */
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    clean_connection(ihalo, 0, 1);
  }
  elapsed = clock()-elapsed;
  fprintf(stderr," done in %4.2f sec.\n", (double)elapsed/CLOCKS_PER_SEC);
  
#endif // MTREE_BOTH_WAYS
  
  
#ifdef SNAPSKIPPING
  elapsed = clock();
  fprintf(stderr,"  o removing missing/uncredible connections ...");
  check_connections(); // this routine has OpenMP parallelisation in it...
  fprintf(stderr,"    -> removing missing/uncredible connections -> done in %4.2f sec.\n", (double)elapsed/CLOCKS_PER_SEC);
#endif
  
  
  write_mtree(OutFile);
  
  return(1);
}

/*==================================================================================================
 * clean_connection
 *==================================================================================================*/
void clean_connection(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t jhalo, icroco, ncroco_new;
  MTREE    *mtree;
#ifdef DEBUG
  uint64_t idesc;
#endif
  
  /* count number of new crocos */
  ncroco_new = 0;
  mtree      = NULL;
  
  /* loop over all cross-correlated haloes */
  for(icroco=0; icroco<halos[isimu0][ihalo].ncroco; icroco++) {
    jhalo = halos[isimu0][ihalo].mtree[icroco].id[1];
    
    /* check whether the present halo is the most likely descendant of this progenitor */
    if(max_merit(jhalo, isimu1) == ihalo) {
      // keep jhalo in mtree-list
      ncroco_new++;
      mtree = (MTREEptr) realloc(mtree, (ncroco_new+1)*sizeof(MTREE));
      
      mtree[ncroco_new-1].id[0]     = halos[isimu0][ihalo].mtree[icroco].id[0];
      mtree[ncroco_new-1].haloid[0] = halos[isimu0][ihalo].mtree[icroco].haloid[0];
      mtree[ncroco_new-1].npart[0]  = halos[isimu0][ihalo].mtree[icroco].npart[0];
      mtree[ncroco_new-1].common    = halos[isimu0][ihalo].mtree[icroco].common;
      mtree[ncroco_new-1].id[1]     = halos[isimu0][ihalo].mtree[icroco].id[1];
      mtree[ncroco_new-1].haloid[1] = halos[isimu0][ihalo].mtree[icroco].haloid[1];
      mtree[ncroco_new-1].npart[1]  = halos[isimu0][ihalo].mtree[icroco].npart[1];
      mtree[ncroco_new-1].merit     = halos[isimu0][ihalo].mtree[icroco].merit;

#ifdef DEBUG
      fprintf(stderr,"icroco=%ld (of %ld) for ihalo=%ld: jhalo=%ld is     a real progenitor of ihalo=%ld (jhalo has %ld descendants)\n",
              icroco,halos[isimu0][ihalo].ncroco,ihalo,
              halos[isimu1][jhalo].haloid,halos[isimu0][ihalo].haloid,
              halos[isimu1][jhalo].ncroco);
#endif
    }
    else {
      // remove jhalo from mtree-list and hence do not add it to the new mtree[] list
#ifdef DEBUG
      fprintf(stderr,"icroco=%ld (of %ld) for ihalo=%ld: jhalo=%ld is NOT a real progenitor of ihalo=%ld (jhalo has %ld descendants)\n",
              icroco,halos[isimu0][ihalo].ncroco,ihalo,
              halos[isimu1][jhalo].haloid,halos[isimu0][ihalo].haloid,
              halos[isimu1][jhalo].ncroco);
      for(idesc=0; idesc<halos[isimu1][jhalo].ncroco; idesc++) {
        fprintf(stderr,"    idesc=%ld haloidesc=%ld\n",idesc,halos[isimu1][jhalo].mtree[idesc].haloid[1]);
      }
#endif
    }
  } // for(icroco)
  
  /* replace halos[isimu0][ihalo].mtree[] with new structure array */
  if(halos[isimu0][ihalo].mtree != NULL) {
    free(halos[isimu0][ihalo].mtree);
    halos[isimu0][ihalo].mtree = NULL;
  }
  halos[isimu0][ihalo].ncroco = ncroco_new;
  halos[isimu0][ihalo].mtree  = mtree;
}

/*==================================================================================================
 * create_mtree
 *==================================================================================================*/
void create_mtree(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t  jhalo, khalo, ipart, jpart, ncroco, icroco;
  int64_t   jcroco;
  
  uint64_t      *common;
  MTREE         *mtree;
  double        *merit;
  long unsigned *idx;
#ifdef USE_PIDMAP
  uint64_t *position;
#endif
  
  // reset the actual mtree[] pointer
  if(halos[isimu0][ihalo].mtree != NULL)
    free(halos[isimu0][ihalo].mtree);
  halos[isimu0][ihalo].mtree = NULL;
  
  /* temporary array pointers */
  common = NULL;
  mtree  = NULL;
  merit  = NULL;
  idx    = NULL;
  
  /* common[] records how many particles ihalo(isimu0) has in common with khalo(isimu1) */
  common = (uint64_t *) calloc(nHalos[isimu1], sizeof(uint64_t));
  
  for(jpart=0; jpart<halos[isimu0][ihalo].npart; jpart++) {
#ifdef USE_PIDMAP
    //  find halos[isimu0][ihalo].Pid[jpart] in PidMap[isimu1] to get the actual ipart to be used with parts[isimu1][ipart]
    position = (uint64_t *) bsearch(&(halos[isimu0][ihalo].Pid[jpart]), PidMap[isimu1], NPids[isimu1], sizeof(uint64_t), pid_cmp);
    
    // when using the PidMap[] we need to make this additional check as this map for each isimu only contains the actual present particles
    if(position != NULL) {
      ipart = (uint64_t) (position - PidMap[isimu1]);
      
      /* ipart belongs to nhalos halos in isimu1 */
      for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
        khalo          = parts[isimu1][ipart].Hid[jhalo];
        common[khalo] += 1;
      }
    }
    else {
      // Pid[jpart] is not in [isimu1] and hence not shared with any halo from [isimu1]
    }
#else
    // simply use Pid as 'ipart'
    ipart = halos[isimu0][ihalo].Pid[jpart];
    
    // here we do not need that check as parts[][] contains PidMax_global entries, irrespective of isimu
    // (and .nhalos will be 0 for non-existent particles...)
    
    /* ipart belongs to nhalos halos in isimu1 */
    for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
      khalo          = parts[isimu1][ipart].Hid[jhalo];
      common[khalo] += 1;
    }
#endif
  }
  
  /* determine number of credible cross-correlations */
  ncroco = 0;
  for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
    if(common[khalo] > MINCOMMON)
      ncroco++;
  }
  halos[isimu0][ihalo].ncroco = ncroco;
  
  /* does not make sense to continue if there are no cross-correlations */
  if(ncroco > 0) {
    
    /* allocate memory for cross-correlations */
    halos[isimu0][ihalo].mtree  = (MTREEptr) calloc(ncroco, sizeof(MTREE));
    mtree  = (MTREEptr)        calloc(ncroco, sizeof(MTREE));
    idx    = (long unsigned *) calloc(ncroco, sizeof(long unsigned));
    merit  = (double *)        calloc(ncroco, sizeof(double));
    
    /* store cross-correlations temporarily for sorting by indexx() */
    icroco = 0;
    for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
      if(common[khalo] > MINCOMMON){
        mtree[icroco].id[0]     = ihalo;
        mtree[icroco].haloid[0] = halos[isimu0][ihalo].haloid;
        mtree[icroco].npart[0]  = halos[isimu0][ihalo].npart;
        mtree[icroco].common    = common[khalo];
        mtree[icroco].id[1]     = khalo;
        mtree[icroco].haloid[1] = halos[isimu1][khalo].haloid;
        mtree[icroco].npart[1]  = halos[isimu1][khalo].npart;
        
        // this is the applied merit function
        mtree[icroco].merit     = pow2((double)common[khalo])/((double)halos[isimu0][ihalo].npart*(double)halos[isimu1][khalo].npart);
        
        merit[icroco] = mtree[icroco].merit;
        
        icroco++;
      }
    }
    
    /* order by merit function */
    indexx((long unsigned)ncroco, merit-1, idx-1);
    
    /* store in halos[isimu0][*].mtree structure (descending order!) */
    for(jcroco=0; jcroco<ncroco; jcroco++) {
      icroco = idx[ncroco-1-jcroco]-1;
      
      /* store mtree[] inside halos[][] structure in correct order */
      halos[isimu0][ihalo].mtree[jcroco].id[0]    = mtree[icroco].id[0];
      halos[isimu0][ihalo].mtree[jcroco].haloid[0]= mtree[icroco].haloid[0];
      halos[isimu0][ihalo].mtree[jcroco].npart[0] = mtree[icroco].npart[0];
      halos[isimu0][ihalo].mtree[jcroco].common   = mtree[icroco].common;
      halos[isimu0][ihalo].mtree[jcroco].id[1]    = mtree[icroco].id[1];
      halos[isimu0][ihalo].mtree[jcroco].haloid[1]= mtree[icroco].haloid[1];
      halos[isimu0][ihalo].mtree[jcroco].npart[1] = mtree[icroco].npart[1];
      halos[isimu0][ihalo].mtree[jcroco].merit    = mtree[icroco].merit;
    }
    
    /* free temporary structures */
    if(mtree) {
      free(mtree);
      mtree = NULL;
    }
    if(idx) {
      free(idx);
      idx = NULL;
    }
    if(merit) {
      free(merit);
      merit = NULL;
    }
  } // if(ncroco)
  
  /* free temporary structures */
  if(common) {
    free(common);
    common = NULL;
  }
}

/*==================================================================================================
 * create_mtree_qsort
 *==================================================================================================*/
void create_mtree_qsort(uint64_t ihalo, int isimu0, int isimu1)
{
  uint64_t  jhalo, khalo, ipart, jpart, ncroco, icroco;
  int64_t   jcroco;
  uint64_t *common;
#ifdef USE_PIDMAP
  uint64_t *position;
#endif
  
  // reset the actual mtree[] pointer
  if(halos[isimu0][ihalo].mtree != NULL)
    free(halos[isimu0][ihalo].mtree);
  halos[isimu0][ihalo].mtree = NULL;
  
  /* temporary array pointers */
  common = NULL;
  
  /* common[] records how many particles ihalo(isimu0) has in common with khalo(isimu1) */
  common = (uint64_t *) calloc(nHalos[isimu1], sizeof(uint64_t));
  
  for(jpart=0; jpart<halos[isimu0][ihalo].npart; jpart++) {
#ifdef USE_PIDMAP
    //  find halos[isimu0][ihalo].Pid[jpart] in PidMap[isimu1] to get the actual ipart to be used with parts[isimu1][ipart]
    position = (uint64_t *) bsearch(&(halos[isimu0][ihalo].Pid[jpart]), PidMap[isimu1], NPids[isimu1], sizeof(uint64_t), pid_cmp);
    
    // when using the PidMap[] we need to make this additional check as this map for each isimu only contains the actual present particles
    if(position != NULL) {
      ipart = (uint64_t) (position - PidMap[isimu1]);
      
      /* ipart belongs to nhalos halos in isimu1 */
      for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
        khalo          = parts[isimu1][ipart].Hid[jhalo];
        common[khalo] += 1;
      }
    }
    else {
      // Pid[jpart] is not in [isimu1] and hence not shared with any halo from [isimu1]
    }
#else
    // simply use Pid as 'ipart'
    ipart = halos[isimu0][ihalo].Pid[jpart];
    
    // here we do not need that check as parts[][] contains PidMax_global entries, irrespective of isimu
    // (and .nhalos will be 0 for non-existent particles...)
    
    /* ipart belongs to nhalos halos in isimu1 */
    for(jhalo=0; jhalo<parts[isimu1][ipart].nhalos; jhalo++) {  // valgrind says "invalid read of size 4" here!?
      khalo          = parts[isimu1][ipart].Hid[jhalo];
      common[khalo] += 1;
    }
#endif
  }
  
  /* determine number of credible cross-correlations */
  ncroco = 0;
  for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
    if(common[khalo] > MINCOMMON)
      ncroco++;
  }
  halos[isimu0][ihalo].ncroco = ncroco;
  
  /* does not make sense to continue if there are no cross-correlations */
  if(ncroco > 0) {
    /* allocate memory for cross-correlations */
    halos[isimu0][ihalo].mtree  = (MTREEptr) calloc(ncroco, sizeof(MTREE));
    
    /* store cross-correlations temporarily for sorting by indexx() */
    icroco = 0;
    for(khalo=0; khalo<nHalos[isimu1]; khalo++) {
      if(common[khalo] > MINCOMMON){
        halos[isimu0][ihalo].mtree[icroco].id[0]     = ihalo;
        halos[isimu0][ihalo].mtree[icroco].haloid[0] = halos[isimu0][ihalo].haloid;
        halos[isimu0][ihalo].mtree[icroco].npart[0]  = halos[isimu0][ihalo].npart;
        halos[isimu0][ihalo].mtree[icroco].common    = common[khalo];
        halos[isimu0][ihalo].mtree[icroco].id[1]     = khalo;
        halos[isimu0][ihalo].mtree[icroco].haloid[1] = halos[isimu1][khalo].haloid;
        halos[isimu0][ihalo].mtree[icroco].npart[1]  = halos[isimu1][khalo].npart;
        halos[isimu0][ihalo].mtree[icroco].merit     = pow2((double)common[khalo])/((double)halos[isimu0][ihalo].npart*(double)halos[isimu1][khalo].npart);
        
        icroco++;
      }
    }
    
    /* qsort particles according to merit function */
    qsort((void *)halos[isimu0][ihalo].mtree, ncroco, sizeof(MTREE), merit_sort);
    
  } else { // if(ncroco)
    halos[isimu0][ihalo].mtree = NULL;
  }
  
  
  /* free temporary structures */
  if(common) {
    free(common);
    common = NULL;
  }
}


