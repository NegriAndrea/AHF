/*==================================================================================================
 *
 *  ahfHaloHistoryStats:   use AHF_halos and MergerTree's _mtree_idx file to follow individual haloes
 *
 *   - code built upon ahfHaloHistory.c
 *   - add-ons to allow for comparison of certain halo properties through time
 *
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <libgen.h>
#include <ctype.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"
#include "../src/libutility/utility.h"

/*-------------------------------------------------------------------------------------
 *                                      DEFINES
 *-------------------------------------------------------------------------------------*/
#undef VERBOSE

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct AHFhalos {
  uint64_t haloid;
  uint64_t hostHalo;
  uint32_t numSubStruct;
  float    Mvir;
  uint32_t npart;
  float    Xc;
  float    Yc;
  float    Zc;
  float    VXc;
  float    VYc;
  float    VZc;
  float    Rvir;
  float    Rmax;
  float    r2;
  float    mbp_offset;
  float    com_offset;
  float    Vmax;
  float    v_esc;
  float    sigV;
  float    lambda;
  float    lambdaE;
  float    Lx;
  float    Ly;
  float    Lz;
  float    b;
  float    c;
  float    Eax;
  float    Eay;
  float    Eaz;
  float    Ebx;
  float    Eby;
  float    Ebz;
  float    Ecx;
  float    Ecy;
  float    Ecz;
  float    ovdens;
  uint32_t nbins;
  float    fMhires;
  float    Ekin;
  float    Epot;
  float    SurfP;
  float    Phi0;
  float    cNFW;
} halo_t;


/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------
 *                                    PROTOTYPES
 *-------------------------------------------------------------------------------------*/
halo_t  *read_halos           (char *, int *);
int      qcompareHostHaloIDs  (const void *, const void *);
int      qcompareHaloIDs      (const void *, const void *);
int      bcompareHostHaloIDs  (const void *, const void *);
int      bcompareHaloIDs      (const void *, const void *);

/*==================================================================================================
 * main
 *==================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
  char      AHF_halos[MAXSTRING];
  halo_t   *halos, *itmp_halo;
  int       num_halos, num_sub;
  uint64_t  haloid, ihalo, jhalo, ipos;
  uint64_t  nhalos_with_sub, n_mismatch, n_pointer;
  uint64_t  haloid_old;
  
  /*==================================================================*
   *                          USER INTERFACE                          *
   *==================================================================*/
  if(argc<2)
   {
    fprintf(stderr,"usage: %s AHF_halos\n",*argv);
    exit(1);
   }
  strcpy(AHF_halos,  argv[1]);

  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"            check substructure information in AHF_halos file\n");
  fprintf(stderr,"=======================================================================\n");
  
  /*==================================================================*
   *                       DEAL WITH USER DATA                        *
   *==================================================================*/
  halos = read_halos(AHF_halos, &num_halos);
  qsort((void *)halos, num_halos, sizeof(halo_t), qcompareHostHaloIDs);
  fprintf(stderr," o read and sorted %d haloes according to 'hostHalo' from file %s\n",num_halos,AHF_halos);

#ifdef DEBUG
  for(ihalo=0; ihalo<num_halos; ihalo++) {
    if(halos[ihalo].hostHalo != 0)
      fprintf(stderr,"ihalo %d  haloid %"PRIu64" hostHalo %"PRIu64"\n",ihalo,halos[ihalo].haloid,halos[ihalo].hostHalo);
  }
#endif
  
  /*==================================================================*
   *                      CHECK numSubStruct                          *
   *==================================================================*/
  fprintf(stderr,"\n");
  
  // check all haloes
  fprintf(stderr," o checking how many haloes point to a halo with numSubStruct>0\n");
  nhalos_with_sub = 0;
  n_mismatch      = 0;
  for(ihalo=0; ihalo<num_halos; ihalo++) {
    
    // but we are only interested in the ones with substructure
    if(halos[ihalo].numSubStruct > 0) {
      
      nhalos_with_sub++;
      
      haloid = halos[ihalo].haloid;
      
#ifdef VERBOSE
      fprintf(stderr," o haloid %021"PRIu64" should have %08"PRIu32" subhaloes: ", haloid, halos[ihalo].numSubStruct);
#endif
      
      // locate all subhaloes using the hostHalo pointer
      itmp_halo = (halo_t *) bsearch(&haloid, halos, num_halos, sizeof(halo_t), bcompareHostHaloIDs);

      if(itmp_halo == NULL) {
#ifdef VERBOSE
        fprintf(stderr," not a single subhalo points to this host!?\n");
#endif
      }
      else {
        ipos = (uint64_t)(itmp_halo-halos);
        
        // move to beginning of matches (bsearch() is undefined if there is more than one match)
        while(halos[ipos].hostHalo == haloid && ipos>=0) {
          ipos--;
        }
        ipos++;
        
        // move forward counting pointers to hostHalo
        num_sub = 0;
        while(halos[ipos].hostHalo == haloid && ipos<num_halos && ipos>=0) {
          num_sub++;
          ipos++;
        }
#ifdef VERBOSE
        fprintf(stderr," found %08d pointers\n",num_sub);
#endif
        if(num_sub != halos[ihalo].numSubStruct) {
          n_mismatch++;
          fprintf(stderr,"        -> mismatch: haloid = %"PRIu64" numSubStruct = %"PRIu32" num_sub = %d\n",
                  halos[ihalo].haloid,halos[ihalo].numSubStruct,num_sub);
        }
      }
    }
  }
  
  fprintf(stderr,"   -> result: nhalos_with_sub = %"PRIu64" n_mismatch = %"PRIu64"\n",nhalos_with_sub,n_mismatch);
  
  
  /*==================================================================*
   *                     CHECK hostHalo pointer                       *
   *==================================================================*/
  qsort((void *)halos, num_halos, sizeof(halo_t), qcompareHaloIDs);
  fprintf(stderr,"\n o sorted haloes according to 'haloid'\n");
  
  fprintf(stderr,"\n o checking whether all hostHalo pointers point to existing halo\n");
  n_pointer  = 0;
  n_mismatch = 0;
  for(ihalo=0; ihalo<num_halos; ihalo++) {
    if(halos[ihalo].hostHalo > 0) {
      n_pointer++;
      haloid    = halos[ihalo].hostHalo;
      itmp_halo = (halo_t *) bsearch(&haloid, halos, num_halos, sizeof(halo_t), bcompareHaloIDs);
      if(itmp_halo == NULL) {
        fprintf(stderr,"   -> could not find haloid=%"PRIu64" (ihalo=%"PRIu64")\n",haloid,ihalo);
        n_mismatch++;
      }
    }
  }
  fprintf(stderr,"   -> result: n_pointer = %"PRIu64" n_mismatch = %"PRIu64"\n",n_pointer,n_mismatch);
  
  /*==================================================================*
   *                 CHECK FOR DUPLICATE HALOIDS                      *
   *==================================================================*/

  fprintf(stderr,"\n o checking for duplicate haloids\n");
  haloid_old   = halos[0].haloid;
  n_mismatch = 0;
  for(ihalo=1; ihalo<num_halos; ihalo++) {
    if(halos[ihalo].haloid == haloid_old) {
      n_mismatch++;
      fprintf(stderr,"     -> mismatch: ihalo = %"PRIu64" haloid = %"PRIu64", ihalo_old = %"PRIu64" haloid_old = %"PRIu64"\n",
              ihalo,halos[ihalo].haloid,ihalo-1,haloid_old);
    }
    haloid_old = halos[ihalo].haloid;
  }
  fprintf(stderr,"   -> result: n_mismatch = %"PRIu64"\n",n_mismatch);
  
  
  /*==================================================================*
   *                           BYE-BYE                                *
   *==================================================================*/
  // free remaining memory
  if(halos) free(halos);
    
  printf("STOP\n");
  return(1);
}

/*==================================================================================================
 * read_halos
 *==================================================================================================*/
halo_t *read_halos(char *infile, int *num_halos)
{
  halo_t *halos;
  char    line[MAXSTRING];
  
  FILE *fp;
  
  // open file
  fp = fopen(infile,"r");
  assert(fp != NULL);
  
  // ignore header line
  fgets(line,MAXSTRING,fp);
  
  // read until EOF
  fgets(line,MAXSTRING,fp);
  *num_halos = 0;
  halos      = NULL;
  while(!feof(fp)) {
    
    // make room for one more halo
    halos = (halo_t *) realloc(halos, ((*num_halos)+1)*sizeof(halo_t));
    
    // read line
    sscanf(line,
           "%"SCNi64" %"SCNi64" %"SCNi32" %f %"SCNi32" %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %"SCNi32" %f %f %f %f %f %f",
           &(halos[*num_halos].haloid),
           &(halos[*num_halos].hostHalo),
           &(halos[*num_halos].numSubStruct),
           &(halos[*num_halos].Mvir),
           &(halos[*num_halos].npart),
           &(halos[*num_halos].Xc),
           &(halos[*num_halos].Yc),
           &(halos[*num_halos].Zc),
           &(halos[*num_halos].VXc),
           &(halos[*num_halos].VYc),
           &(halos[*num_halos].VZc),
           &(halos[*num_halos].Rvir),
           &(halos[*num_halos].Rmax),
           &(halos[*num_halos].r2),
           &(halos[*num_halos].mbp_offset),
           &(halos[*num_halos].com_offset),
           &(halos[*num_halos].Vmax),
           &(halos[*num_halos].v_esc),
           &(halos[*num_halos].sigV),
           &(halos[*num_halos].lambda),
           &(halos[*num_halos].lambdaE),
           &(halos[*num_halos].Lx),
           &(halos[*num_halos].Ly),
           &(halos[*num_halos].Lz),
           &(halos[*num_halos].b),
           &(halos[*num_halos].c),
           &(halos[*num_halos].Eax),
           &(halos[*num_halos].Eay),
           &(halos[*num_halos].Eaz),
           &(halos[*num_halos].Ebx),
           &(halos[*num_halos].Eby),
           &(halos[*num_halos].Ebz),
           &(halos[*num_halos].Ecx),
           &(halos[*num_halos].Ecy),
           &(halos[*num_halos].Ecz),
           &(halos[*num_halos].ovdens),
           &(halos[*num_halos].nbins),
           &(halos[*num_halos].fMhires),
           &(halos[*num_halos].Ekin),
           &(halos[*num_halos].Epot),
           &(halos[*num_halos].SurfP),
           &(halos[*num_halos].Phi0),
           &(halos[*num_halos].cNFW)
           );
    
    if((int64_t)halos[*num_halos].hostHalo < 0)
      halos[*num_halos].hostHalo = 0;
    
    // increment halo counter
    (*num_halos)++;
    
    // next halo
    fgets(line,MAXSTRING,fp);
  }
  
  // close file
  fclose(fp);
  
  return(halos);
}

/*==============================================================================
 *  compare halo ids (used with qsort)
 *==============================================================================*/
int qcompareHostHaloIDs(const void *halo1, const void *halo2)
{
	uint64_t n1, n2;
  
	n1 = ((halo_t *)halo1)->hostHalo;
	n2 = ((halo_t *)halo2)->hostHalo;
  
	return n1 < n2 ? -1 : (n1 > n2 ? 1 : 0);
}
int qcompareHaloIDs(const void *halo1, const void *halo2)
{
	uint64_t n1, n2;
  
	n1 = ((halo_t *)halo1)->haloid;
	n2 = ((halo_t *)halo2)->haloid;
  
	return n1 < n2 ? -1 : (n1 > n2 ? 1 : 0);
}

/*==============================================================================
 *  compare particle ids (used with bsearch)
 *==============================================================================*/
int bcompareHostHaloIDs(const void *id, const void *halo)
{
	const uint64_t *i       = (const uint64_t *)id;
  const uint64_t hostHalo = ((halo_t *)halo)->hostHalo;
  
	return *i < hostHalo ? -1 : (*i > hostHalo ? 1 : 0);
}

int bcompareHaloIDs(const void *id, const void *halo)
{
	const uint64_t *i       = (const uint64_t *)id;
  const uint64_t haloid   = ((halo_t *)halo)->haloid;
  
	return *i < haloid ? -1 : (*i > haloid ? 1 : 0);
}


