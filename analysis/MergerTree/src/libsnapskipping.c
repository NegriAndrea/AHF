#include "include.h"
#include "common.h"

/*==================================================================================================
 * check for uncredible progenitors:
 *
 *   if a halo at [0] does not have a credible progenitor at [1] it will be added to the
 *   list of halos at [1] and considered when checking the connection [1]->[2]
 *==================================================================================================*/
int connectionrejection(HALOS halo, MTREE progmtree)
{
  double Mratio;
  
  //  Mratio = (double)halo.npart/(double)halo.mtree[0].npart[1];
  Mratio = (double)halo.npart/(double)progmtree.npart[1];
  if(Mratio<1) Mratio=1.0/Mratio;
  
#ifdef SNAPSKIPPING_UNCREDIBLEMASSRATIO
  if(Mratio > SNAPSKIPPING_UNCREDIBLEMASSRATIO)
    return(1);
  else
    return(0);
#else
  return(0);
#endif
}

void check_connections()
{
  MTREE    mtree_tmp;
  uint64_t ihalo, ipart, nHalos1, icroco, jcroco;
  int      rejected, missing;
  int     *remi_flag; // we generate an array as large as nHalos[0] that stores 1(remove) / 0(keep) for each ihalo
  
  PidMax_global = MAX(PidMax[0],PidMax_global);
  PidMax_global = MAX(PidMax[1],PidMax_global);
  nHalos1       = nHalos[1];
  
  fprintf(stderr,"(PidMax_global=%"PRIu64")",PidMax_global);
  
  remi_flag = (int *)calloc(nHalos[0],sizeof(int)); // name logic: REjectedMIssing
  
  // TODO: when using remi_flag[] we could parallize the following for(ihalo)
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) private(ihalo, rejected, missing, icroco, mtree_tmp) shared(halos, remi_flag, nHalos, stderr) default(none)
#endif
  // loop over all halos at [0]
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    
    // we assume that this halo is *not* yet missing/rejected
    rejected = 0;
    missing  = 0;
    
    // if a progenitor exists check for credibility...
    if(halos[0][ihalo].ncroco > 0) {
      
#ifdef SNAPSKIPPING_CONSIDERALLPROGENITORS
      // loop over all progenitors, checking the connectionrejection condition for halos[0][ihalo]
      for(icroco=0; icroco<halos[0][ihalo].ncroco; icroco++) {
        if(connectionrejection(halos[0][ihalo], halos[0][ihalo].mtree[icroco]) == FALSE) {
          // the first one not giving 'connectionrejection(halos[0][ihalo]) == TRUE' should be the new top-of-the-list
          break;
        }
      }
      
      // if all progenitors failed the test we cannot assign a credible progenitor
      if(icroco == halos[0][ihalo].ncroco) {
#ifdef DEBUG_SNAPSKIPPING
        fprintf(stderr,"   rejected connection:\n");
        fprintf(stderr,"      [0]: haloid=%12"PRIu64" npart=%12"PRIu64" ncroco=%12"PRIu64" mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64" mtree.common=%12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].npart, halos[0][ihalo].ncroco, halos[0][ihalo].mtree[0].haloid[0], halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].npart[0], halos[0][ihalo].mtree[0].common);
        fprintf(stderr,"      [1]:                                                               mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", halos[0][ihalo].mtree[0].haloid[1], halos[0][ihalo].mtree[0].id[1], halos[0][ihalo].mtree[0].npart[1]);
#endif
        rejected = 1;
        
        // flag halos[0][ihalo] to not be written (only halos with ncroco>0 will be written to file...)
        halos[0][ihalo].ncroco = 0;
      }
      // but the first one that complies with our criterion (i.e. icroco) will now become the top-of-the-list
      else {
#ifdef DEBUG_SNAPSKIPPING
        if(icroco>0) {
          fprintf(stderr,"   using minor branch connection:\n");
          fprintf(stderr,"      [0]: haloid=%12"PRIu64" npart=%12"PRIu64" ncroco=%12"PRIu64" mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].npart, halos[0][ihalo].ncroco, halos[0][ihalo].mtree[0].haloid[0], halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].npart[0]);
          fprintf(stderr,"      [1]: old =>                  merit=%lf                                mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", halos[0][ihalo].mtree[0].merit, halos[0][ihalo].mtree[0].haloid[1], halos[0][ihalo].mtree[0].id[1], halos[0][ihalo].mtree[0].npart[1]);
          fprintf(stderr,"      [1]: new => icroco=%"PRIu64" merit=%lf                       mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", icroco,halos[0][ihalo].mtree[icroco].merit,halos[0][ihalo].mtree[icroco].haloid[1], halos[0][ihalo].mtree[icroco].id[1], halos[0][ihalo].mtree[icroco].npart[1]);
        }
#endif
        
        // re-shuffle the mtree[] list up to this point
        mtree_tmp.haloid[0] = halos[0][ihalo].mtree[icroco].haloid[0];
        mtree_tmp.haloid[1] = halos[0][ihalo].mtree[icroco].haloid[1];
        mtree_tmp.id[0]     = halos[0][ihalo].mtree[icroco].id[0];
        mtree_tmp.id[1]     = halos[0][ihalo].mtree[icroco].id[1];
        mtree_tmp.npart[0]  = halos[0][ihalo].mtree[icroco].npart[0];
        mtree_tmp.npart[1]  = halos[0][ihalo].mtree[icroco].npart[1];
        mtree_tmp.common    = halos[0][ihalo].mtree[icroco].common;
        mtree_tmp.merit     = halos[0][ihalo].mtree[icroco].merit;
        for(jcroco=1; jcroco<=icroco; jcroco++) {
          halos[0][ihalo].mtree[jcroco].haloid[0] = halos[0][ihalo].mtree[jcroco-1].haloid[0];
          halos[0][ihalo].mtree[jcroco].haloid[1] = halos[0][ihalo].mtree[jcroco-1].haloid[1];
          halos[0][ihalo].mtree[jcroco].id[0]     = halos[0][ihalo].mtree[jcroco-1].id[0];
          halos[0][ihalo].mtree[jcroco].id[1]     = halos[0][ihalo].mtree[jcroco-1].id[1];
          halos[0][ihalo].mtree[jcroco].npart[0]  = halos[0][ihalo].mtree[jcroco-1].npart[0];
          halos[0][ihalo].mtree[jcroco].npart[1]  = halos[0][ihalo].mtree[jcroco-1].npart[1];
          halos[0][ihalo].mtree[jcroco].common    = halos[0][ihalo].mtree[jcroco-1].common;
          halos[0][ihalo].mtree[jcroco].merit     = halos[0][ihalo].mtree[jcroco-1].merit;
        }
        halos[0][ihalo].mtree[0].haloid[0] = mtree_tmp.haloid[0];
        halos[0][ihalo].mtree[0].haloid[1] = mtree_tmp.haloid[1];
        halos[0][ihalo].mtree[0].id[0]     = mtree_tmp.id[0];
        halos[0][ihalo].mtree[0].id[1]     = mtree_tmp.id[1];
        halos[0][ihalo].mtree[0].npart[0]  = mtree_tmp.npart[0];
        halos[0][ihalo].mtree[0].npart[1]  = mtree_tmp.npart[1];
        halos[0][ihalo].mtree[0].common    = mtree_tmp.common;
        halos[0][ihalo].mtree[0].merit     = mtree_tmp.merit;
      }
      
#else // SNAPSKIPPING_CONSIDERALLPROGENITORS
      
      if(connectionrejection(halos[0][ihalo], halos[0][ihalo].mtree[0]) == TRUE) {
#ifdef DEBUG_SNAPSKIPPING
        fprintf(stderr,"   rejected connection:\n");
        fprintf(stderr,"      [0]: haloid=%12"PRIu64" npart=%12"PRIu64" ncroco=%12"PRIu64" mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64" mtree.common=%12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].npart, halos[0][ihalo].ncroco, halos[0][ihalo].mtree[0].haloid[0], halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].npart[0], halos[0][ihalo].mtree[0].common);
        fprintf(stderr,"      [1]:                                                               mtree.haloid=%12"PRIu64" mtree.id=%12"PRIu64" mtree.npart=%12"PRIu64"\n", halos[0][ihalo].mtree[0].haloid[1], halos[0][ihalo].mtree[0].id[1], halos[0][ihalo].mtree[0].npart[1]);
#endif
        rejected = 1;
        
        // flag halos[0][ihalo] to not be written (only halos with ncroco>0 will be written to file...)
        halos[0][ihalo].ncroco = 0;
      }
#endif // SNAPSKIPPING_CONSIDERALLPROGENITORS
    }
    // if a progenitor has not been found, treat as rejected and keep searching...
    else {
#ifdef DEBUG_SNAPSKIPPING
      fprintf(stderr,"   missing connection:\n");
      fprintf(stderr,"      [0]: haloid=%12"PRIu64" npart=%12"PRIu64" ncroco=%12"PRIu64" mtree=%12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].npart, halos[0][ihalo].ncroco, halos[0][ihalo].mtree);
#endif
      missing = 1;
    }
    
    
    // now deal with the missing/rejected conncetion...
    if(missing || rejected) {
      
      // only try to follow the missing/rejected halo, if it has enough particles
      if(halos[0][ihalo].npart > MINCOMMON) {
        
        // REjectedMIssing = True
        remi_flag[ihalo] = 1;
        
        // in order to keep this loop parallel, we adjust the actual halos[1] array below in a separate, serial loop
        
      }// if(MINCOMMON)
    } // if (missing || rejected)
  } // for(ihalo)
    
  // use remi_flag[] to copy rejected/missing objects over to nHalos[1], using a new, serial for(ihalo) loop this time
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    if(remi_flag[ihalo] == 1) {
#ifdef DEBUG_SNAPSKIPPING
        fprintf(stderr,"           copying information:\n");
        fprintf(stderr,"             [0]:  haloid=%"PRIu64" npart=%"PRIu64" ncroco=%"PRIu64"\n",halos[0][ihalo].haloid,halos[0][ihalo].npart,halos[0][ihalo].ncroco);
        fprintf(stderr,"             [1]:  nHalos=%"PRIu64" -> ",nHalos[1]);
#endif
      halos[1] = (HALOptr) realloc(halos[1], (nHalos[1]+1)*sizeof(HALOS));
      halos[1][nHalos[1]].haloid = halos[0][ihalo].haloid;
      halos[1][nHalos[1]].npart  = halos[0][ihalo].npart;
      halos[1][nHalos[1]].ncroco = halos[0][ihalo].ncroco;
      halos[1][nHalos[1]].mtree  = NULL;  // we do not have any merger tree information for this one (yet)
      halos[1][nHalos[1]].Pid    = (uint64_t *) calloc(halos[0][ihalo].npart, sizeof(uint64_t));
      for(ipart=0; ipart<halos[1][nHalos[1]].npart; ipart++) {
        halos[1][nHalos[1]].Pid[ipart] = halos[0][ihalo].Pid[ipart];
        //if(halos[1][nHalos[1]].Pid[ipart] > PidMax[1]) PidMax[1] = halos[1][nHalos[1]].Pid[ipart]; // not needed as particle_halo_mapping takes MAX(PidMax[0],PidMax[1])
      }
      
      // increment number of halos at [1]
      nHalos[1]++;
#ifdef DEBUG_SNAPSKIPPING
        fprintf(stderr,"%"PRIu64"\n",nHalos[1]);
#endif
    } // if(remi_flag)
  } // for(ihalo)
  
  // have we added additional halos to [1]?
  if(nHalos[1] > nHalos1) {
#ifdef USE_PIDMAP
    // update the PidMap for [1]
    free(PidMap[1]);
    create_PidMap(1);
#endif
    
    // we need to update the particle_halo_maping for [1]
    for(ipart=0; ipart<PidMax_global+1; ipart++) {   // free() old memory first
      if(parts[1][ipart].Hid != NULL){
        free(parts[1][ipart].Hid);
        parts[1][ipart].Hid = NULL;
      }
    }
    free(parts[1]);
    parts[1] = NULL;
    
    fprintf(stderr,"\n");
    particle_halo_mapping(1);                       // now call for re-creation of that mapping
  }
  
  free(remi_flag);
  
  return;
}
