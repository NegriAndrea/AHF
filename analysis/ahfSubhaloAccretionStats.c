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
#define hostHalo_changes(x,y)          (((x) == 0 && (y) != 0) || ((x) != 0 && (y) == 0)) // either becoming a host or becoming a subhalo
#define hostHalo_becomes_hosthalo(x,y) (((x) >  0 && (y) == 0))                           // only becoming a subhalo
#define hostHalo_becomes_subhalo(x,y)  (((x) == 0 && (y)  > 0))                           // only becoming a hosthalo

#define condition(x,y)          hostHalo_becomes_subhalo(x,y)                            // condition used in code below

#define SET_ZRED_AND_ISNAP    {  if(filemap == NULL) { \
z     = get_redshift(AHF_halos[n]); \
isnap = -1; \
} else { \
z     = filemap[n].zred; \
isnap = filemap[n].isnap; }}



//#define FULL_AHF_INFORMATION   // will write the full AHF_halos line into the output file
#define MINIMAL_OUTPUT // only writes the information needed by Noam to re-do the subaccretion analysis excl. the backsplashed haloes

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
#ifdef FULL_AHF_INFORMATION
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
  float    TheRest[28];
#endif
} halo_t;

typedef struct MergerTree {
  uint64_t haloid;
  uint64_t progid;
} mtree_t;

typedef struct {
  double zred;
  int    isnap;
} zred_isnap_map_t;

/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/
char **AHF_halos;
char **mtree_idx;
int  num_files;

uint64_t *haloid;
int  num_haloid;


/*-------------------------------------------------------------------------------------
 *                                    PROTOTYPES
 *-------------------------------------------------------------------------------------*/
void              get_filename_list(char *);
void              get_haloid_list(char *);
zred_isnap_map_t *get_filemap(char *);
void              get_headerline(char *, char*);
void              get_haloline(char *, uint64_t, char *);
int64_t           get_haloidx(char *, uint64_t);
double            get_redshift(char *);
int               cmp_haloid(const void *, const void *);
halo_t  *read_halos                (char *, int *);
mtree_t *read_mtree_idx            (char *, int *);
int      qcompareHaloIDs           (const void *, const void *);
int      qcompareMtreeIDs          (const void *, const void *);
int      bcompareHaloIDs           (const void *, const void *);
int      bcompareMtreeIDs          (const void *, const void *);


/*==================================================================================================
 * main
 *==================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
  FILE *fpout;
  char prefix_list[MAXSTRING], haloid_list[MAXSTRING], outfile[MAXSTRING], haloline[MAXSTRING], zred_list[MAXSTRING];
  zred_isnap_map_t *filemap;
  int i, n;
  int64_t ihalo, ihalo_0, hostHalo_track;
  double z;
  int    isnap;
  
  // the property we aim at tracking
  uint64_t idummy;
  uint64_t hostHalo_now, hostHalo_prev;
  
  // the pointers to the arrays holding the all halos and mtree files
  halo_t  **halos;
  mtree_t **mtree;
  int     *num_halos, *num_mtree;
  halo_t  *itmp_halo;
  uint64_t jhalo;
  mtree_t *itmp_mtree;
  uint64_t jmtree;
  
  /*==================================================================*
   *                          USER INTERFACE                          *
   *==================================================================*/
  if(argc<3)
   {
    fprintf(stderr,"usage: %s haloid_list prefix_list zred_list\n",*argv);
    exit(1);
   }
  strcpy(haloid_list,      argv[1]);
  strcpy(prefix_list,      argv[2]);
  strcpy(zred_list, argv[3]);
  
  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"               follow AHF_halos through _mtree_idx files\n");
  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"will read mtree_idx and AHF_halos files using prefixes found in %s\n",prefix_list);
  fprintf(stderr,"will read haloids to follow from                                %s\n",haloid_list);
  fprintf(stderr,"will read mapping between snapid and redshift from              %s\n\n",zred_list);
  
  /*==================================================================*
   *                          READ USER DATA                          *
   *==================================================================*/
  get_haloid_list(haloid_list);
  get_filename_list(prefix_list);
  filemap = get_filemap(zred_list);  // returns NULL, if zred_list = "-1" (and then zred will be extracted from prefix_list)
  
  // open general outfile
  sprintf(outfile,"ahfSubhaloAccretionStats_%03d.dat",filemap[0].isnap);
  fpout = fopen(outfile,"w");
  assert(fpout != NULL);
#ifndef MINIMAL_OUTPUT
  fprintf(fpout,"#haloid(1) hostHalo(2) numSubStruct(3) Mvir(4) npart(5) Xc(6) Yc(7) Zc(8)\n");
  fprintf(fpout,"#  h/s/b(1) ISNAP(2) z(3) hostHalo(4) ihalo(5)\n");
#endif
  
  /*==================================================================*
   *                            DO THE WORK                           *
   *==================================================================*/
  fprintf(stderr,"\n");
  
  // allocate memory for all data involved
  halos     = (halo_t **)  calloc(num_files,   sizeof(halo_t *));
  num_halos = (int *)      calloc(num_files,   sizeof(int));
  mtree     = (mtree_t **) calloc(num_files,   sizeof(mtree_t *)); // one too many, but who cares
  num_mtree = (int *)      calloc(num_files,   sizeof(int));
  
  // read all files involved
  for(i=0; i<num_files; i++) {
    halos[i] = read_halos(AHF_halos[i], &(num_halos[i]));
    qsort((void *)halos[i], num_halos[i], sizeof(halo_t), qcompareHaloIDs);
    fprintf(stderr,"read and sorted %d haloes from file %s\n",num_halos[i],AHF_halos[i]);
  }
  for(i=0; i<num_files-1; i++) {
    mtree[i]  = read_mtree_idx(mtree_idx[i], &(num_mtree[i]));
    qsort((void *)(mtree[i]), num_mtree[i], sizeof(mtree_t), qcompareMtreeIDs);
    fprintf(stderr,"read and sorted %d entries from mtree file %s\n",num_mtree[i],mtree_idx[i]);
  }
  
  for(i=0; i<num_haloid; i++) {
    ihalo   = haloid[i];
    ihalo_0 = ihalo;    // remember the halo traced (primarily for MINIMAL_OUTPUT)
    fprintf(stderr," + tracing halo %"PRIu64" through %d files: ",ihalo,num_files);
    
    
    //--------------------------------------------
    // initialize halo property
    //get_haloline(AHF_halos[0], ihalo, haloline);
    //sscanf(haloline,"%"SCNi64"%"SCNi64,&idummy, &hostHalo_prev);
    //fprintf(fpout,"%s",haloline);
    itmp_halo = (halo_t *)bsearch(&(ihalo), halos[0], num_halos[0], sizeof(halo_t), bcompareHaloIDs);
    if(itmp_halo == NULL) {
      fprintf(stderr,"could not find %"PRIi64" in %s\n",ihalo,AHF_halos[0]);
      exit(0);
    }
    jhalo         = (uint64_t) (itmp_halo - halos[0]);
    hostHalo_prev = halos[0][jhalo].hostHalo;
#ifndef MINIMAL_OUTPUT
    fprintf(fpout,"%"PRIu64" %"PRIu64" %"PRIu32" %g %"PRIu32" %f %f %f\n",
            halos[0][jhalo].haloid,
            halos[0][jhalo].hostHalo,
            halos[0][jhalo].numSubStruct,
            halos[0][jhalo].Mvir,
            halos[0][jhalo].npart,
            halos[0][jhalo].Xc,
            halos[0][jhalo].Yc,
            halos[0][jhalo].Zc);
#endif
    //--------------------------------------------
    
    
    hostHalo_track = hostHalo_prev;
    
    // loop over all files
    for(n=0; n<num_files-1; n++) {
      
      //--------------------------------------------
      //get_haloline(AHF_halos[n], ihalo, haloline);
      //sscanf(haloline,"%"SCNi64"%"SCNi64,&idummy, &hostHalo_now);
      itmp_halo = (halo_t *)bsearch(&(ihalo), halos[n], num_halos[n], sizeof(halo_t), bcompareHaloIDs);
      if(itmp_halo == NULL) {
        fprintf(stderr,"could not find %"PRIi64" in %s\n",ihalo,AHF_halos[n]);
        exit(0);
      }
      jhalo        = (uint64_t) (itmp_halo - halos[i]);
      hostHalo_now = halos[i][jhalo].hostHalo;
      //--------------------------------------------
      
      // hostHalo ID changed and now points to nowhere
      if(hostHalo_becomes_hosthalo(hostHalo_prev,hostHalo_now)) {
        
        // keep track of the last host it belonged to
        hostHalo_track = hostHalo_prev;
        
#ifndef MINIMAL_OUTPUT
        SET_ZRED_AND_ISNAP;
        //fprintf(fpout," h %03d %8.4lf %021"PRIu64" %021"PRIu64" %021"PRIi64"/",isnap,z,hostHalo_now,hostHalo_prev,hostHalo_track);
        fprintf(fpout,"h %03d %8.4lf %"PRIu64" %"PRIu64"\n",isnap,z,hostHalo_prev,ihalo);
#endif
        
        //--------------------------------------------
        // but update hostHalo_track to the present host halo
        //hostHalo_track = get_haloidx(mtree_idx[n-1], hostHalo_track);
        itmp_mtree = (mtree_t *) bsearch(&(hostHalo_track), (void *)(mtree[n-1]), num_mtree[n-1], sizeof(mtree_t), bcompareMtreeIDs);
        if(itmp_mtree == NULL) {
          fprintf(stderr,"could not find %"PRIi64" in %s",hostHalo_track,mtree_idx[n-1]);
          hostHalo_track = -1;
        }
        else {
          jmtree         = (uint64_t)(itmp_mtree - mtree[n-1]);
          hostHalo_track = mtree[n-1][jmtree].progid;
        }
        //--------------------------------------------
      }
      
      // hostHalo ID changed and now points to a host halo
      if(hostHalo_becomes_subhalo(hostHalo_prev,hostHalo_now)) {
        SET_ZRED_AND_ISNAP;
        
        // does it becomes subhalo of a new host halo
        if(hostHalo_track != hostHalo_now) {
#ifndef MINIMAL_OUTPUT
          //fprintf(fpout," s %03d %8.4lf %021"PRIu64" %021"PRIu64" %021"PRIi64"\n",isnap,z,hostHalo_now,hostHalo_prev,hostHalo_track);
          fprintf(fpout,"s %03d %8.4lf %"PRIu64" %"PRIu64"\n",isnap,z,hostHalo_now,ihalo);
#endif
        } else {
#ifndef MINIMAL_OUTPUT
          //fprintf(fpout," b %03d %8.4lf %021"PRIu64" %021"PRIu64" %021"PRIi64"\n",isnap,z,hostHalo_now,hostHalo_prev,hostHalo_track);
          fprintf(fpout,"b %03d %8.4lf %"PRIu64" %"PRIu64"\n",isnap,z,hostHalo_now,ihalo);
#else
          fprintf(fpout,"%"PRIu64"\n",ihalo_0);
          break;
#endif
        }
      }
      
      // remember hostHalo_now pointer
      hostHalo_prev = hostHalo_now;
      
      // follow last hostHalo backward in time
      if(hostHalo_track > 0) {
        //--------------------------------------------
        //hostHalo_track = get_haloidx(mtree_idx[n], hostHalo_track);
        itmp_mtree = (mtree_t *) bsearch(&(hostHalo_track), (void *)(mtree[n]), num_mtree[n], sizeof(mtree_t), bcompareMtreeIDs);
        if(itmp_mtree == NULL) {
          fprintf(stderr,"could not find %"PRIi64" in %s",hostHalo_track,mtree_idx[n]);
          hostHalo_track = -1;
        }
        else {
          jmtree = (uint64_t)(itmp_mtree - mtree[n]);
          hostHalo_track = mtree[n][jmtree].progid;
        }
        //--------------------------------------------
      }
      
      //--------------------------------------------
      // follow ihalo backwards in time (if possible)
      //ihalo = get_haloidx(mtree_idx[n], ihalo);
      itmp_mtree = (mtree_t *) bsearch(&(ihalo), (void *)(mtree[n]), num_mtree[n], sizeof(mtree_t), bcompareMtreeIDs);
      if(itmp_mtree == NULL) {
        ihalo = -1;
      }
      else {
        jmtree = (uint64_t)(itmp_mtree - mtree[n]);
        ihalo  = mtree[n][jmtree].progid;
      }
      //--------------------------------------------
      
      fprintf(stderr,"%"PRIi64" ",ihalo);
      if(ihalo < 0)
        break;
    }
    fprintf(stderr,"\n");
    
    
    
    
    
    
    // the last file requires special treatment as we do not need to find the progenitor
    if(n == num_files-1 && ihalo > 0) {
      //--------------------------------------------
      //get_haloline(AHF_halos[n], ihalo, haloline);
      //sscanf(haloline,"%"SCNi64"%"SCNi64,&idummy, &hostHalo_now);
      itmp_halo = (halo_t *)bsearch(&(ihalo), halos[n], num_halos[n], sizeof(halo_t), bcompareHaloIDs);
      if(itmp_halo == NULL) {
        fprintf(stderr,"could not find %"PRIi64" in %s\n",ihalo,AHF_halos[n]);
        exit(0);
      }
      jhalo = (uint64_t) (itmp_halo - halos[i]);
      hostHalo_now = halos[i][jhalo].hostHalo;
      //--------------------------------------------
      
      // hostHalo ID changed and now points to nowhere
      if(hostHalo_becomes_hosthalo(hostHalo_prev,hostHalo_now)) {
        
        // keep track of the last host it belonged to
        hostHalo_track = hostHalo_prev;
        
#ifndef MINIMAL_OUTPUT
        SET_ZRED_AND_ISNAP;
        //fprintf(fpout," h %03d %8.4lf %021"PRIu64" %021"PRIu64" %021"PRIi64"/",isnap,z,hostHalo_now,hostHalo_prev,hostHalo_track);
        fprintf(fpout," h %03d %8.4lf END\n",isnap,z);
#endif
        
        //--------------------------------------------
        // but update hostHalo_track to the present host halo
        //hostHalo_track = get_haloidx(mtree_idx[n-1], hostHalo_track);
        itmp_mtree = (mtree_t *) bsearch(&(hostHalo_track), (void *)(mtree[n-1]), num_mtree[n-1], sizeof(mtree_t), bcompareMtreeIDs);
        if(itmp_mtree == NULL) {
          fprintf(stderr,"could not find %"PRIi64" in %s",hostHalo_track,mtree_idx[n-1]);
          hostHalo_track = -1;
        }
        else {
          jmtree = (uint64_t)(itmp_mtree - mtree[n-1]);
          hostHalo_track = mtree[n-1][jmtree].progid;
        }
        //--------------------------------------------
      }
      
      // hostHalo ID changed and now points to a host halo
      if(hostHalo_becomes_subhalo(hostHalo_prev,hostHalo_now)) {
        SET_ZRED_AND_ISNAP;
        
        // does it becomes subhalo of a new host halo
        if(hostHalo_track != hostHalo_now) {
#ifndef MINIMAL_OUTPUT
          //fprintf(fpout," s %03d %8.4lf %021"PRIu64" %021"PRIu64" %021"PRIi64"\n",isnap,z,hostHalo_now,hostHalo_prev,hostHalo_track);
          fprintf(fpout," s %03d %8.4lf END\n",isnap,z);
#endif
        } else {
#ifndef MINIMAL_OUTPUT
          //fprintf(fpout," b %03d %8.4lf %021"PRIu64" %021"PRIu64" %021"PRIi64"\n",isnap,z,hostHalo_now,hostHalo_prev,hostHalo_track);
          fprintf(fpout," b %03d %8.4lf END\n",isnap,z);
#else
          fprintf(fpout,"%"PRIu64"\n",ihalo_0);
#endif
        }
      }
    }
    fflush(fpout);
  }
  
  
  /*==================================================================*
   *                              BYE-BYE                             *
   *==================================================================*/
  // close output file
  fclose(fpout);
  
  free(haloid);
  for(i=0; i<num_files; i++) {
    // free filename strings
    if(AHF_halos[i]) free(AHF_halos[i]);
    if(mtree_idx[i]) free(mtree_idx[i]);
    
    // free data
    if(halos[i]) free(halos[i]);
    if(mtree[i]) free(mtree[i]);
  }
  if(AHF_halos) free(AHF_halos);
    if(mtree_idx) free(mtree_idx);
      if(halos)     free(halos);
        if(mtree)     free(mtree);
          if(filemap)   free(filemap);
            
            printf("STOP\n");
            return(1);
}

/*==================================================================================================
 * get_filenam_list
 *==================================================================================================*/
void get_filename_list(char *prefix_list)
{
  FILE *fp;
  char line[MAXSTRING];
  int  i;
  
  // open file
  fp = fopen(prefix_list,"r");
  assert(fp != NULL);
  
  // figure out number of AHF_halos files
  num_files = 0;
  while(!feof(fp)) {
    fgets(line,MAXSTRING,fp);
    num_files++;
  }
  num_files--; // we counted one too many
  
  // allocate memory
  AHF_halos = (char **) calloc(num_files, sizeof(char *));
  mtree_idx = (char **) calloc(num_files, sizeof(char *)); // one too many, but who cares
  
  // eventually read haloids
  fprintf(stderr," + will use the following files for the tracing:\n");
  rewind(fp);
  for(i=0; i<num_files; i++) {
    fscanf(fp,"%s",line);
    AHF_halos[i] = (char *) calloc(MAXSTRING,sizeof(char));
    mtree_idx[i] = (char *) calloc(MAXSTRING,sizeof(char));
    sprintf(AHF_halos[i],"%s_halos",line);
    sprintf(mtree_idx[i],"%s_mtree_idx",line);
    if(i==num_files-1)
      fprintf(stderr,"     %s\n", AHF_halos[i]);
    else
      fprintf(stderr,"     %s, %s\n", AHF_halos[i], mtree_idx[i]);
  }
  
  // close file
  fclose(fp);
}

/*==================================================================================================
 * get_haloid_list
 *==================================================================================================*/
void get_haloid_list(char *haloid_list)
{
  FILE *fp;
  char line[MAXSTRING];
  int  i;
  
  // open file
  fp = fopen(haloid_list,"r");
  assert(fp != NULL);
  
  // figure out number of haloids
  num_haloid = 0;
  while(!feof(fp)) {
    fgets(line,MAXSTRING,fp);
    num_haloid++;
  }
  num_haloid--; // we counted one too many
  
  // allocate memory
  haloid = (uint64_t *) calloc(num_haloid, sizeof(uint64_t));
  
  // eventually read haloids
  fprintf(stderr," + found %d haloids to be traced:\n",num_haloid);
  rewind(fp);
  for(i=0; i<num_haloid; i++) {
    fgets(line,MAXSTRING,fp);
    sscanf(line,"%"SCNu64,&(haloid[i]));
    if(i == num_haloid-1)
      fprintf(stderr,"     %"PRIu64,haloid[i]);
    else
      fprintf(stderr,"     %"PRIu64", ",haloid[i]);
  }
  fprintf(stderr,"\n");
  
  // close file
  fclose(fp);
  
  // sort haloid[] ascending
  //qsort(haloid, num_haloid, sizeof(uint64_t), &cmp_haloid);
}


/*==================================================================================================
 * get_headerline
 *==================================================================================================*/
void get_headerline(char *halofile, char *line)
{
  FILE *fp;
  int64_t haloid;
  
  fp = fopen(halofile,"r");
  assert(fp != NULL);
  
  strcpy(line,"a");
  while(strncmp(line,"#",1) != 0 && !feof(fp)) {
    fgets(line,MAXSTRING,fp);
  }
  
  if(feof(fp)) {
    strcpy(line,"# could not find any header!?\n");
  }
  
  fclose(fp);
}

/*==================================================================================================
 * get_haloline
 *==================================================================================================*/
void get_haloline(char *halofile, uint64_t id, char *line)
{
  FILE *fp;
  int64_t haloid;
  
  fp = fopen(halofile,"r");
  assert(fp != NULL);
  
  haloid = -1;
  while(haloid != id && !feof(fp)) {
    fgets(line,MAXSTRING,fp);
    if(strncmp(line,"#",1) != 0)
      sscanf(line,"%"SCNi64,&haloid);
    else
      haloid = -1;
  }
  
  if(feof(fp)) {
    strcpy(line,"-1\n");
  }
  
  fclose(fp);
}

/*==================================================================================================
 * get_haloidx
 *==================================================================================================*/
int64_t get_haloidx(char *halofile, uint64_t id)
{
  FILE *fp;
  int64_t haloid, haloidx;
  char line[MAXSTRING];
  
  fp = fopen(halofile,"r");
  assert(fp != NULL);
  
  haloid = -1;
  while(haloid != id && !feof(fp)) {
    fgets(line,MAXSTRING,fp);
    if(strncmp(line,"#",1) != 0)
      sscanf(line,"%"SCNi64,&haloid);
    else
      haloid = -1;
  }
  
  if(feof(fp)) {
    fclose(fp);
    return(-1);
  }
  else {
    sscanf(line,"%"SCNu64"%"SCNu64,&haloid,&haloidx);
    fclose(fp);
    return(haloidx);
  }
  
}

/*==================================================================================================
 * get_redshift
 *==================================================================================================*/
double get_redshift(char *prefix)
{
	double z    = -99.0;
	int    zIdx = strlen(prefix);
  
  // ignore all EOS's
  while(prefix[zIdx] == '\0')
    zIdx--;
  
  // hunt down the 'z'
  while(prefix[zIdx] != 'z')
    zIdx--;
  
	if (prefix[zIdx] == 'z') {
		(void)sscanf(prefix + zIdx, "z%lf", &z);
	}
  
	return z;
}

/*==============================================================================
 *  compare particle ids (used with qsort and bsearch)
 *==============================================================================*/
int cmp_haloid(const void *i1, const void *i2)
{
	uint64_t *haloid1, *haloid2;
  
	haloid1 = (uint64_t *) i1;
	haloid2 = (uint64_t *) i2;
  
	return *haloid1 < *haloid2 ? -1 : (*haloid1 > *haloid2 ? 1 : 0);
}

/*==================================================================================================
 * get_filemap
 *==================================================================================================*/
zred_isnap_map_t *get_filemap(char *zred_list)
{
  FILE *fp;
  char line[MAXSTRING];
  int  num_zred, i_zred;
  zred_isnap_map_t *map;
  double zred;
  int    isnap;
  
  map = NULL;
  
  // we intent to extract the redshift from the filename prefix and not using a mapping file
  if(strncmp(zred_list,"-1",2) == 0) {
    return (NULL);
  }
  
  // open file
  fp = fopen(zred_list,"r");
  
  // count number of mappings
  num_zred = 0;
  while(!feof(fp)) {
    fgets(line,MAXSTRING,fp);
    num_zred++;
  }
  num_zred--;  // we counted one too many
  
  // the number of redshifts should match at least match the number of files
  if(num_zred < num_files) {
    fprintf(stderr,"there are fewer redshifts (%d) in %s than files (%d) to be used during the tracking!?\nABORTING\n",
            num_zred,zred_list,num_files);
    exit(0);
  }
  
  // allocate memory for redshifts
  map = (zred_isnap_map_t *) calloc(num_zred, sizeof(zred_isnap_map_t));
  
  // read redshifts
  rewind(fp);
  for(i_zred=0; i_zred<num_zred; i_zred++) {
    fgets(line,MAXSTRING,fp);
    sscanf(line,"%lf %d",&zred,&isnap);
    map[i_zred].zred  = zred;
    map[i_zred].isnap = isnap;
  }
  
  // close file
  fclose(fp);
  
  // return pointer to zred[] array
  return(map);
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
#ifdef FULL_AHF_INFORMATION
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
#else
    sscanf(line,
           "%"SCNi64" %"SCNi64" %"SCNi32" %f %"SCNi32" %f %f %f",
           &(halos[*num_halos].haloid),
           &(halos[*num_halos].hostHalo),
           &(halos[*num_halos].numSubStruct),
           &(halos[*num_halos].Mvir),
           &(halos[*num_halos].npart),
           &(halos[*num_halos].Xc),
           &(halos[*num_halos].Yc),
           &(halos[*num_halos].Zc)
           );
#endif
    
    // increment halo counter
    (*num_halos)++;
    
    // next halo
    fgets(line,MAXSTRING,fp);
  }
  
  // close file
  fclose(fp);
  
  return(halos);
}

/*==================================================================================================
 * read_mtree_idx
 *==================================================================================================*/
mtree_t *read_mtree_idx(char *infile, int *num_mtree)
{
  mtree_t *mtree;
  char    line[MAXSTRING];
  
  FILE *fp;
  
  // open file
  fp = fopen(infile,"r");
  assert(fp != NULL);
  
  // ignore header line
  fgets(line,MAXSTRING,fp);
  
  // read until EOF
  fgets(line,MAXSTRING,fp);
  *num_mtree = 0;
  mtree      = NULL;
  while(!feof(fp)) {
    
    // make room for one more entry
    mtree = (mtree_t *) realloc(mtree, ((*num_mtree)+1)*sizeof(mtree_t));
    
    // read line
    sscanf(line, "%"SCNi64" %"SCNi64, &(mtree[(*num_mtree)].haloid), &(mtree[(*num_mtree)].progid));
    
    // increment line counter
    (*num_mtree)++;
    
    // next line
    fgets(line,MAXSTRING,fp);
  }
  
  // close file
  fclose(fp);
  
  return(mtree);
}

/*==============================================================================
 *  compare halo ids (used with qsort)
 *==============================================================================*/
int qcompareHaloIDs(const void *halo1, const void *halo2)
{
	uint64_t n1, n2;
  
	n1 = ((halo_t *)halo1)->haloid;
	n2 = ((halo_t *)halo2)->haloid;
  
	return n1 < n2 ? -1 : (n1 > n2 ? 1 : 0);
}

/*==============================================================================
 *  compare halo ids in mtree (used with qsort)
 *==============================================================================*/
int qcompareMtreeIDs(const void *mtree1, const void *mtree2)
{
	uint64_t n1, n2;
  
	n1 = ((mtree_t *)mtree1)->haloid;
	n2 = ((mtree_t *)mtree2)->haloid;
  
	return n1 < n2 ? -1 : (n1 > n2 ? 1 : 0);
}

/*==============================================================================
 *  compare particle ids (used with bsearch)
 *==============================================================================*/
int bcompareMtreeIDs(const void *id, const void *mtree)
{
  const uint64_t haloid = ((mtree_t *)mtree)->haloid;
	const uint64_t *i     = (const uint64_t *)id;
  
	return *i < haloid ? -1 : (*i > haloid ? 1 : 0);
}

/*==============================================================================
 *  compare particle ids (used with bsearch)
 *==============================================================================*/
int bcompareHaloIDs(const void *id, const void *halo)
{
	const uint64_t *i       = (const uint64_t *)id;
  const uint64_t haloid = ((halo_t *)halo)->haloid;
  
	return *i < haloid ? -1 : (*i > haloid ? 1 : 0);
}





