#ifdef AHF2

/***************************************************************************
 *   Includes                                                              *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <inttypes.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

#include "../libutility/utility.h"

#include "../libtree/tree.h"
#include "../libtree/utilities.h"

int iter=0;

//*********************************************************************************************************
// PROTOTYPES
//*********************************************************************************************************
void generate_daughter_halos  (HALO *halos, patch_t *parent_patch, uint64_t ihost, uint64_t *ihalo);
int  qcmpHALOnpart            (const void *halo1, const void *halo2);
void set_gatherRad_like_AHF1  (HALO *halos, uint64_t numHalos);
void read_spatialRef          (patch_t *** patch_tree, uint64_t **n_patches);


//*********************************************************************************************************
// patchtree2halos()
//*********************************************************************************************************
halo_s* patchtree2halos (patch_t **patch_tree, uint64_t *n_patches)
{
  int      initial_level, final_level, ilevel;
  uint64_t ipatch, iparent;
  HALO    *halos;
  uint64_t ihalo, nhalos, ihost;
  halo_s  *return_values;
  patch_t  parent_patch;
  double   x, y, z;
  

#ifdef AHF2_read_spatialRef
  //======================================================================================
  // read spatialRef[][] from file and put it into patch_tree[][]
  //======================================================================================
  read_spatialRef(&patch_tree, &n_patches);
#endif
  
  
  //======================================================================================
  // determine initial and final level
  //======================================================================================
  get_patch_level_range(patch_tree, n_patches, &initial_level, &final_level);

#ifdef DEBUG_AHF2
  fprintf(stderr,"  initial_level=%d final_level=%d\n",initial_level,final_level);
#endif
  
  //======================================================================================
  // first guess for number of haloes = all the end-leaves
  //======================================================================================
  nhalos = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for default(none) private(ilevel,ipatch) shared(patch_tree,n_patches,initial_level,final_level) schedule(dynamic) reduction( + : nhalos )
#endif
  for(ilevel=initial_level; ilevel<final_level; ilevel++) {
    for(ipatch=0; ipatch<n_patches[ilevel]; ipatch++) {
      if(patch_tree[ilevel][ipatch].n_daughters == 0) {
        nhalos++;
      } // if()
    } // ipatch
  } // ilevel

#ifdef VERBOSE
  fprintf(stderr,"  patchtree2halos: number of end-leaves in patch_tree[][] -> nhalos=%"PRIu64"\n",nhalos);
#endif
  fprintf(io.logfile,"  patchtree2halos: number of end-leaves in patch_tree[][] -> nhalos=%"PRIu64"\n",nhalos);
  fflush(io.logfile);
  
  //======================================================================================
  // allocate memory to store all the prospective halos
  //======================================================================================
  halos = (HALO *) calloc(nhalos, sizeof(HALO));
  if(halos == NULL) {
#ifdef VERBOSE
    fprintf(stderr,"  patchtre2halos: cannot allocate memory for halos[] array\n");
#endif
    fprintf(io.logfile,"  patchtre2halos: cannot allocate memory for halos[] array\n");
    fflush(io.logfile);
    exit(0);
  }
  

  //======================================================================================
  // pick each parent patch in initial_level and treat it as a halo
  //======================================================================================
  // use initial_level as starting point
  ilevel      = initial_level;
  ahf.min_ref = initial_level + (int)(log(simu.NGRID_DOM)/log(2));
  
  // accumulate haloes
  ihalo = 0;

  // this loop is intrinsically serial due to the use of 'ihalo' as the counter for the created halos
  // !!!!! NO OpenMP possible here !!!!
  for(iparent=0; iparent<n_patches[ilevel]; iparent++) {
    
    // set parent patch for easier access
    parent_patch = (patch_tree[ilevel][iparent]);
    
    // TODO: deal with possible mergers on the highest level
    if(parent_patch.trunk == NO_SINGLE_TRUNK) {
    }
  
    // this walks down the trunk to get the most refined centre
    get_patch_centre(parent_patch, &x, &y, &z);
    
#ifdef DEBUG_AHF2
    fprintf(stderr,"creating parent halo %"PRIu64" from patch-%"PRIi32"-%"PRIi64" at x=%lf y=%lf z=%lf\n",ihalo,parent_patch.level,parent_patch.id,x,y,z);
#endif
    
    // copy properties over to halos[] array
    halos[ihalo].npart         = parent_patch.Npart;
    halos[ihalo].numNodes      = parent_patch.n_subcubes;
    halos[ihalo].refLev        = parent_patch.level;
    halos[ihalo].spaRes        = 1./(pow(2.0,parent_patch.level));
    halos[ihalo].pos.x         = x;
    halos[ihalo].pos.y         = y;
    halos[ihalo].pos.z         = z;
    halos[ihalo].gatherRad     = MIN(parent_patch.radius, simu.MaxGatherRad/simu.boxsize);
    halos[ihalo].R_vir         = MIN(parent_patch.radius, simu.MaxGatherRad/simu.boxsize);
    
    halos[ihalo].hostHalo      = -1;
    halos[ihalo].hostHaloLevel = parent_patch.level;
    halos[ihalo].numSubStruct  = 0;    // will be calculated as we go along
    halos[ihalo].subStruct     = NULL;
    
    // keep track of the host halo
    ihost = ihalo;
    
    // generate daughter halos[]
    generate_daughter_halos(halos, &parent_patch, ihost, &ihalo);

    // next parent halo
    ihalo++;
  }
  
  // consistency check
  if(ihalo != nhalos) {
    fprintf(stderr,"  patchtree2halos: something went wrong => ihalo=%"PRIu64" nhalos=%"PRIu64"\n",ihalo,nhalos);
    exit(0);
  }
  
#ifdef VERBOSE
  fprintf(stderr,"  patchtree2halos: filled halos[] with %"PRIu64" halos (nhalos=%"PRIu64")\n",ihalo,nhalos);
#endif
  fprintf(io.logfile,"  patchtree2halos: filled halos[] with %"PRIu64" halos (nhalos=%"PRIu64")\n",ihalo,nhalos);
  fflush(io.logfile);
  
#ifdef AHF2_set_gatherRad_like_AHF1
  set_gatherRad_like_AHF1(halos, nhalos);
#endif
  
  //======================================================================================
  // put everything into a halo_s container to be returned
  //======================================================================================
  return_values         = (halo_s *) calloc(1, sizeof(halo_s));
  return_values->halos  = halos;
  return_values->nhalos = nhalos;
  return(return_values);
}


//*********************************************************************************************************
// generate_daughter_halos()
//*********************************************************************************************************
void generate_daughter_halos(HALO *halos, patch_t *parent_patch, uint64_t ihost, uint64_t *ihalo)
{
  uint64_t inewhost;
  int32_t  idaughter;
  double   x, y, z;
  patch_t *daughter_patch;
  
  //===========================================================
  // check for consistency
  //===========================================================
  if(parent_patch->n_daughters > 0 && parent_patch->trunk == END_OF_TRUNK) {
    fprintf(stderr,"generate_daughter_halos(): something wrong here: there are n_daughters=%"PRIu64", but no trunk=%d\n",
            parent_patch->n_daughters,parent_patch->trunk);
    exit(0);
  }
  
  //===========================================================
  // loop over all daughters
  //===========================================================
  for(idaughter=0; idaughter<parent_patch->n_daughters; idaughter++) {

    //===========================================================
    // TODO: deal with mergers
    //===========================================================
    if(parent_patch->trunk == NO_SINGLE_TRUNK) {
#ifdef DEBUG_AHF2
      fprintf(stderr,"generate_daughter_halos(): untreated merger for ihalo=%"PRIu64" ihost=%"PRIu64" trunk=%d\n",*ihalo,ihost,parent_patch->trunk);
#endif
      parent_patch->trunk = 0;
    }
   
    //===========================================================
    // direct access to the daughter patch
    //===========================================================
    daughter_patch = parent_patch->daughters[idaughter];

#ifdef DEBUG_AHF2
    fprintf(stderr,"  *daughter_patch = %ld (idaughter=%"PRIi32", parent_patch->n_daughters=%"PRIu64")\n", *daughter_patch, idaughter, parent_patch->n_daughters);
#endif
    
    //===========================================================
    // do not create the parent halo again
    //===========================================================
    if(idaughter != parent_patch->trunk) {

      // we are inserting halos[] even though they are "just" the daughters
      (*ihalo)++;
      
#ifdef DEBUG_AHF22
      fprintf(stderr,"creating daughter halo %"PRIu64" from patch-%"PRIi32"-%"PRIi64" (whose parent is patch-%"PRIi32"-%"PRIi64")\n",
              *ihalo,daughter_patch->level,daughter_patch->id,parent_patch->level,parent_patch->id);
#endif
      
      // add a pointer to the host
      halos[ihost].numSubStruct++;
      halos[ihost].subStruct = (int *) realloc(halos[ihost].subStruct, halos[ihost].numSubStruct*sizeof(int));
      if(halos[ihost].subStruct == NULL) {
        fprintf(stderr,"generate_daughter_halos(): could not realloc memory for halos[].subStruct\n");
        exit(-1);
      }
      halos[ihost].subStruct[halos[ihost].numSubStruct-1] = (int) (*ihalo);
      
#ifdef DEBUG_AHF22
      fprintf(stderr,"   halo %"PRIu64" is now substructure of ihost=%"PRIu64": halos[%"PRIu64"].numSubStruct=%d \n",
              *ihalo,ihost,ihost,halos[ihost].numSubStruct);
#endif
      
      // this walks down the trunk to get the most refined centre
      get_patch_centre(*daughter_patch, &x, &y, &z);
      
      // copy properties over to halos[] array
      halos[*ihalo].npart          = daughter_patch->Npart;
      halos[*ihalo].numNodes       = daughter_patch->n_subcubes;
      halos[*ihalo].refLev         = daughter_patch->level;
      halos[*ihalo].spaRes         = 1./(pow(2.0,daughter_patch->level));
      halos[*ihalo].pos.x          = x;
      halos[*ihalo].pos.y          = y;
      halos[*ihalo].pos.z          = z;
      halos[*ihalo].R_vir          = MIN(daughter_patch->radius, simu.MaxGatherRad/simu.boxsize);
      halos[*ihalo].gatherRad      = MIN(daughter_patch->radius, simu.MaxGatherRad/simu.boxsize);
      
      halos[*ihalo].hostHalo       = (int)ihost;
      halos[*ihalo].hostHaloLevel  = parent_patch->level;
      halos[*ihalo].numSubStruct   = 0;    // will be calculated as we go along
      halos[*ihalo].subStruct      = NULL;

      // treat daughter as the host to new subhaloes
      inewhost = *ihalo;
    } // if(trunk)
    
    //===========================================================
    // simply follow the trunk refining some of the ihalo values
    //===========================================================
    else {
#ifdef DEBUG_AHF2
      fprintf(stderr,"following the trunk of halo %"PRIu64" from patch-%"PRIi32"-%"PRIi64" (whose parent is patch-%"PRIi32"-%"PRIi64")\n",
              *ihalo,daughter_patch->level,daughter_patch->id,parent_patch->level,parent_patch->id);
#endif
      
#ifdef AHF2_read_spatialRef
      // spatialRef[][] only stored linked particles on each isolated refinement
      halos[ihost].npart     += daughter_patch->Npart;
#endif

      // follow same logic as AHF1
      halos[ihost].numNodes   = daughter_patch->n_subcubes;
      halos[ihost].refLev     = daughter_patch->level;
      halos[ihost].spaRes     = 1./(pow(2.0,daughter_patch->level));
#ifdef DEBUG_AHF2
      fprintf(stderr,"readjusted host: x=%lf y=%lf z=%lf Npart=%"PRIu64" (added %"PRIu64") refLev=%d spaRes=%20.16lf\n",
              halos[ihost].pos.x,halos[ihost].pos.y,halos[ihost].pos.z,halos[ihost].npart,daughter_patch->Npart,halos[ihost].refLev,halos[ihost].spaRes);
#endif

      // treat original halo as the host to new subhaloes
      inewhost = ihost;
    }

    //===========================================================
    // generate daughter halos[] of this daughter
    //===========================================================
    generate_daughter_halos(halos, daughter_patch, inewhost, ihalo);

  } // for(idaughter)
}

/*==============================================================================
 *  compare halo npart (used with qsort)
 *==============================================================================*/
int qcmpHALOnpart(const void *halo1, const void *halo2)
{
	uint64_t n1, n2;
  
	n1 = ((HALO *)halo1)->npart;
	n2 = ((HALO *)halo2)->npart;
  
	return n1 < n2 ? -1 : (n1 > n2 ? 1 : 0);
}

/*==============================================================================
 * the gathering radius is the distance to the closest halo that is more massive than the current halo...
 * we therefore require the halos to be ordered by npart for this part!
 *==============================================================================*/
void set_gatherRad_like_AHF1(HALO *halos, uint64_t numHalos)
{
  double         maxGathRad, dx, dy, dz, tmpRad;
  int64_t        ii, jj, i, j, inumpart, jnumpart;
  long unsigned *idx, *idxtmp;
  double        *fsort;
  
  // sort haloes by npart
  idx    = (long unsigned *)calloc(numHalos, sizeof(long unsigned));
  idxtmp = (long unsigned *)calloc(numHalos, sizeof(long unsigned));
  fsort  = (double *)       calloc(numHalos, sizeof(double));
  for (i = 0; i < numHalos; i++)
    fsort[i] = (double)halos[i].npart;
  indexx(numHalos, fsort-1, idxtmp-1);
  
  // indexx sorts ascending and gives indizes starting at 1
  for (i = 0; i < numHalos; i++)
    idx[numHalos - i - 1] = idxtmp[i] - 1;
  free(idxtmp);
  free(fsort);
  
  // maximal allowed gathering radius (user provided)
	maxGathRad = MIN(simu.MaxGatherRad / simu.boxsize, 1. / 4.);
  
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) private(ii,i,jj,j,inumpart,jnumpart,tmpRad,dx,dy,dz) shared(halos,numHalos,idx,maxGathRad,stderr) default(none)
#endif
	for (ii = numHalos - 1; ii >= 0; ii--) {
    
    // pick halo from ordered list
    i = idx[ii];

    // start with maximum allowed gathering radius
    halos[i].gatherRad = pow2(2.*maxGathRad);
    
    inumpart = halos[i].npart;
    
    // loop over all halos that are more massive...
    for (jj = ii; jj >= 0; jj--) {
      
      // pick halo from ordered list
      j = idx[jj];
      
      jnumpart = halos[j].npart;
      
      if ((i != j) && (jnumpart > inumpart)) {
        // Calculate the distance to this more massive halo
        dx = fabs(halos[i].pos.x - halos[j].pos.x);
        dy = fabs(halos[i].pos.y - halos[j].pos.y);
        dz = fabs(halos[i].pos.z - halos[j].pos.z);
        if (dx > 0.5)		dx = 1.0 - dx;
        if (dy > 0.5) 	dy = 1.0 - dy;
        if (dz > 0.5) 	dz = 1.0 - dz;
        
        tmpRad = pow2(dx) + pow2(dy) + pow2(dz);
        
        // Is it the smallest distance
        if (tmpRad < halos[i].gatherRad) {
          halos[i].gatherRad = tmpRad;
        }
      }
    }
    
    // Finally calculating the gathering radius
    halos[i].gatherRad = (sqrt(halos[i].gatherRad)) * 0.5;
    
		// restrict the gathering radius (Note: R_vir at this stage is closeRefDist or -1)
		halos[i].gatherRad = MAX(halos[i].gatherRad, halos[i].R_vir);
		halos[i].gatherRad = MIN(halos[i].gatherRad, maxGathRad);

#ifdef DEBUG_AHF2
    fprintf(stderr,"set_gatherRad_like_AHF1:  ii=%"PRIi64" halos[%"PRIi64"].gatherRad=%lf\n",ii,i,halos[i].gatherRad);
#endif
  } // for ( ii=numHalos-1; ii>=0; ii-- )
  
  free(idx);
}

//======================================================================================
// read spatialRef[][] from file (as written by AHF1)
//======================================================================================
void read_spatialRef(patch_t *** patch_tree, uint64_t **n_patches)
{
  FILE *fp;
  char     infile[MAXSTRING], fprefix[MAXSTRING], file_no[MAXSTRING];
  int      trunk;
  double   x, y, z;
  int      ilevel;
  uint64_t ipatch;
  
  // AHF language
  int      min_ref, no_grids;
  int      i, j, k, numIsoRef, refLevel, isoRefIndex, numSubStruct;
  double   closeRefDist;
  int      numNodes;
  long     numParts;
  
  fprintf(stderr,"    read_spatialRef():\n");
  fprintf(stderr,"    ==================\n");
  
  // create and fill patch_tree[] with NULL pointers
  *patch_tree = (patch_t **) calloc(NLEVELS, sizeof(patch_t *));
  *n_patches  = (uint64_t *) calloc(NLEVELS, sizeof(uint64_t));
  for(ilevel=0; ilevel<NLEVELS; ilevel++) {
    (*patch_tree)[ilevel] = NULL;
    (*n_patches)[ilevel]  = 0;
  }
  
  
  // open _gridtree file
#ifdef WITH_MPI
	snprintf(fprefix, MAXSTRING, "%s.%04d.", global_io.params->outfile_prefix, global_mpi.rank);
#else
	snprintf(fprefix, MAXSTRING, "%s.", global_io.params->outfile_prefix);
#endif
  
	/* OUTPUT CONVENTION: 3 digits*/
	sprintf(file_no, "z%.3f", fabs(global.z));
	strcat(fprefix, file_no);
  
  sprintf(infile,"%s.AHF_gridtree",fprefix);
  fprintf(stderr,"Reading all patchtree[][] arrays from %s ... ",infile);
  fp = fopen(infile,"r");
  if(fp == NULL) {
    fprintf(stderr,"patchtree2halos: could not open %s\n",infile);
    exit(0);
  }
  
  // read min level and total number of levels (in AHF1 language)
  fscanf(fp,"%d %d", &min_ref, &no_grids);
  
#ifdef DEBUG_AHF2
  fprintf(stderr,"  min_ref=%d no_grids=%d (ahf.min_ref=%d)\n",min_ref,no_grids,ahf.min_ref);
#endif
  
  // ignore all levels before min_ref (patch_tree[] pointers are NULL)
  for(ilevel=min_ref; ilevel<min_ref+no_grids; ilevel++) {
    
    // number of patches on ilevel
    fscanf(fp, "%d %d", &i, &numIsoRef);
    
    // store number of patches
    (*n_patches)[ilevel] = numIsoRef;
    
    if(i != ilevel) {
      fprintf(stderr,"patchtree2halos: level mis-match: i=%d ilevel=%"PRIu64" (%d %d , min_ref=%d)\n",i,ilevel,i,numIsoRef,min_ref);
      exit(0);
    }
    
    // memory for patch_tree[][]
    (*patch_tree)[ilevel] = (patch_t *) calloc(numIsoRef, sizeof(patch_t));
    if((*patch_tree)[ilevel] == NULL) {
      fprintf(stderr,"patchtree2halos: could not allocate memory for patch_tree[%"PRIu64"]: numIsoRef=%d\n",ilevel,numIsoRef);
      exit(0);
    }
    
    // loop over all patches
    for(ipatch=0; ipatch<numIsoRef; ipatch++) {
      
      // read line
      fscanf(fp,"%lf %lf %lf %lf %d %ld %d %d %d",
             &x,
             &y,
             &z,
             &closeRefDist,
             &numNodes,
             &numParts,
             &refLevel,          // of daughter!
             &isoRefIndex,       // of daughter!
             &numSubStruct);
      
      if(isoRefIndex != -1 && refLevel != ilevel+1) {
        fprintf(stderr,"  a) spatialRef[][] not as expected: refLevel=%d ilevel+1=%d (%d %d , min_ref=%d)\n",refLevel,ilevel+1,refLevel,isoRefIndex,min_ref);
        exit(0);
      }
      
      // copy to patch_tree[][]
      (*patch_tree)[ilevel][ipatch].id          = ipatch;
      (*patch_tree)[ilevel][ipatch].pos[0]      = x;
      (*patch_tree)[ilevel][ipatch].pos[1]      = y;
      (*patch_tree)[ilevel][ipatch].pos[2]      = z;
      (*patch_tree)[ilevel][ipatch].Npart       = numParts;
      (*patch_tree)[ilevel][ipatch].n_subcubes  = numNodes;
      (*patch_tree)[ilevel][ipatch].level       = ilevel;
      (*patch_tree)[ilevel][ipatch].radius      = closeRefDist;
      (*patch_tree)[ilevel][ipatch].trunk       = isoRefIndex; // trunk requires more sophistication to be set (see below)
      (*patch_tree)[ilevel][ipatch].n_daughters = numSubStruct;
      
      if(numSubStruct > 0) {
        (*patch_tree)[ilevel][ipatch].daughters   = (patch_t **) calloc(numSubStruct, sizeof(patch_t *));
        (*patch_tree)[ilevel][ipatch].idaughter   = (uint64_t *) calloc(numSubStruct, sizeof(uint64_t));
      }
      else {
        (*patch_tree)[ilevel][ipatch].daughters   = NULL;
        (*patch_tree)[ilevel][ipatch].idaughter   = NULL;
      }
      
      // loop over all substructures
      for(k=0; k<numSubStruct; k++) {
        fscanf(fp,"%d %d", &refLevel, &isoRefIndex);
        if(refLevel != ilevel+1) {
          fprintf(stderr,"  b) spatialRef[][] not as expected: refLevel=%d ilevel+1=%d (%d %d , min_ref=%d)\n",refLevel,ilevel+1,refLevel,isoRefIndex,min_ref);
          exit(0);
        }
        (*patch_tree)[ilevel][ipatch].idaughter[k] = isoRefIndex;
        
        // find the trunk
        if(isoRefIndex == (*patch_tree)[ilevel][ipatch].trunk) {
          trunk = k;
        }
      } //k
      
      // set the trunk properly
      if((*patch_tree)[ilevel][ipatch].n_daughters == 0) {
        (*patch_tree)[ilevel][ipatch].trunk = -1;
      }
      else  {
        (*patch_tree)[ilevel][ipatch].trunk = trunk;
      }
      
#ifdef DEBUG_AHF2
      fprintf(stderr," patch_tree[%ld][%ld].trunk = %d (n_daughters=%d ilevel=%d)\n",(long)ilevel,(long)ipatch,(*patch_tree)[ilevel][ipatch].trunk,(int)(*patch_tree)[ilevel][ipatch].n_daughters,(int)(*patch_tree)[ilevel][ipatch].level);
#endif
      
    } // ipatch
  } // ilevel
  
  fclose(fp);
  
  // we still need to set the correct daughter pointers as patchtree[][] requires 'patch_t*' instead of indizes!
  for(ilevel=min_ref; ilevel<min_ref+no_grids-1; ilevel++) {
    for(ipatch=0; ipatch<(*n_patches)[ilevel]; ipatch++) {
      
      // loop over all daughters
      for(k=0; k<(*patch_tree)[ilevel][ipatch].n_daughters; k++) {
        (*patch_tree)[ilevel][ipatch].daughters[k] = &( (*patch_tree)[ilevel+1][(*patch_tree)[ilevel][ipatch].idaughter[k]] );
      } // k
      
    } // ipatch
  } // ilevel
  
  fprintf(stderr,"done\n");
}

#endif // AHF2

