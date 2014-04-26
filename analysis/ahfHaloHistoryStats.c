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
#define condition(x,y) (((x) == 0 && (y) != 0) || ((x) != 0 && (y) == 0))


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
  float    TheRest[28];
} halo_t;

typedef struct MergerTree {
  uint64_t haloid;
  uint64_t progid;
} mtree_t;



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
void     get_filename_list(char *);
void     get_haloid_list(char *);
void     get_headerline(char *, char*);
void     get_haloline(char *, uint64_t, char *);
int64_t  get_haloidx(char *, uint64_t);
double   get_redshift(char *);
double  *get_zred(char *);
int      cmp_haloid(const void *, const void *);


/*==================================================================================================
 * main
 *==================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
  FILE *fpout;
  char prefix_list[MAXSTRING], haloid_list[MAXSTRING], outfile[MAXSTRING], haloline[MAXSTRING], zred_list[MAXSTRING];
  double *zred;
  int i, n;
  int64_t ihalo;
  double z;
  
  // the property we aim at tracking
  uint64_t idummy;
  uint64_t hostHalo_now, hostHalo_prev;
  
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
  zred = get_zred(zred_list); // returns NULL, if zred_list = "-1" (and then zred will be extracted from prefix_list)
  
  // open general outfile
  sprintf(outfile,"ahfHaloHistoryStats.dat");
  fpout = fopen(outfile,"w");
  assert(fpout != NULL);

  /*==================================================================*
   *                            DO THE WORK                           *
   *==================================================================*/
  for(i=0; i<num_haloid; i++) {
    ihalo = haloid[i];
    fprintf(stderr," + tracing halo %"PRIu64" through %d files: ",ihalo,num_files);
    
    // initialize halo property
    get_haloline(AHF_halos[0], ihalo, haloline);
    sscanf(haloline,"%"SCNi64"%"SCNi64,&idummy, &hostHalo_prev);
    fprintf(fpout,"%"PRIu64" %"PRIu64"\n",ihalo,hostHalo_prev);
    
    // loop over all files
    for(n=0; n<num_files-1; n++) {
      get_haloline(AHF_halos[n], ihalo, haloline);
      sscanf(haloline,"%"SCNi64"%"SCNi64,&idummy, &hostHalo_now);
      if(zred == NULL)
        z = get_redshift(AHF_halos[n]);
      else
        z = zred[n];
      
      if(condition(hostHalo_prev,hostHalo_now)) {
        fprintf(fpout," %lf %s",z,haloline);
      }
      hostHalo_prev = hostHalo_now;
      
      
      ihalo = get_haloidx(mtree_idx[n], ihalo);
      fprintf(stderr,"%"PRIi64" ",ihalo);
      if(ihalo < 0)
        break;
    }
    fprintf(stderr,"\n");
    
    // the last file requires special treatment as we do not need to find the progenitor
    if(ihalo >= 0) {
      get_haloline(AHF_halos[n], ihalo, haloline);
      sscanf(haloline,"%"SCNi64"%"SCNi64,&idummy, &hostHalo_now);
      if(zred == NULL)
        z = get_redshift(AHF_halos[n]);
      else
        z = zred[n];

      if(condition(hostHalo_prev,hostHalo_now)) {
        fprintf(fpout," %lf %s",z,haloline);
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
    free(AHF_halos[i]);
    free(mtree_idx[i]);
  }
  free(AHF_halos);
  free(mtree_idx);
  if(zred) free(zred);
  
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
  mtree_idx = (char **) calloc(num_files, sizeof(char *));
  
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
 * get_zred
 *==================================================================================================*/
double *get_zred(char *zred_list)
{
  FILE *fp;
  char line[MAXSTRING];
  int  num_zred, i_zred;
  double *zred;
  
  zred = NULL;
  
  // we intent to extract the redshift from the filename prefix and not using a mapping file
  if(strncmp(zred_list,"-1",2) == 0) {
    return (zred);
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
  zred = (double *) calloc(num_zred, sizeof(double));
  
  // read redshifts
  rewind(fp);
  for(i_zred=0; i_zred<num_zred; i_zred++) {
    fgets(line,MAXSTRING,fp);
    sscanf(line,"%lf",&(zred[i_zred]));
  }
  
  // close file
  fclose(fp);
  
  // return pointer to zred[] array
  return(zred);
}



