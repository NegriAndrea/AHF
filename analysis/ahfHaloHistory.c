/*==================================================================================================
 *
 *  ahfHaloHistory:   use AHF_halos and MergerTree's _mtree_idx file to follow individual haloes
 *
 *  WARNING:
 *       this code is not tuned for performance and works better the fewer haloes are traced!
 *            
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <libgen.h>
#include <ctype.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"
#include "../src/libutility/utility.h"

/*-------------------------------------------------------------------------------------
 *                                     FEATURES
 *-------------------------------------------------------------------------------------*/
#define PROVIDE_LINENUMBER_AS_HALOID

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct my_halofile {
  char *line;
} halofile_t;

typedef struct my_mtreefile {
  uint64_t ida;
  uint64_t idx;
} mtreefile_t;

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
int64_t  get_haloidx(mtreefile_t *, uint64_t, uint64_t);
double   get_redshift(char *);
double  *get_zred(char *);
int      cmp_haloid(const void *, const void *);
uint64_t read_halos(char *halofile, halofile_t **halobuffer);
uint64_t read_mtree_idx(char *, mtreefile_t **);


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
  int i, n, k;
  int64_t ihalo, id;
  double z;
  halofile_t  *halobuffer;
  mtreefile_t **mtreebuffer;
  uint64_t     Nhalos, *Nmtree;
  
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
  
  
  /*==================================================================*
   *                            DO THE WORK                           *
   *==================================================================*/
  // open all halo files putting the header into them
  for(i=0; i<num_haloid; i++) {
    
    // id of halo to be traced
    ihalo = haloid[i];
    
    // (re-)open outfile for this halo
    sprintf(outfile,"halo_%07"PRIu64".dat",ihalo);
    fpout = fopen(outfile,"w");
    if(fpout == NULL) {
      fprintf(stderr,"could not open %s\n",outfile);
      exit(0);
    }
    get_headerline(AHF_halos[0], haloline);
    fprintf(fpout,"#added redshift as first column so please add +1 for all other column numbers\n");
    fprintf(fpout,"%s",haloline);
    fclose(fpout);
  }
  
  // read all mtree_idx files
  mtreebuffer = (mtreefile_t **) calloc(num_files-1, sizeof(mtreefile_t *));
  Nmtree      = (uint64_t *) calloc(num_files-1, sizeof(uint64_t));
  for(n=0; n<num_files-1; n++) {
    Nmtree[n] = read_mtree_idx(mtree_idx[n], &(mtreebuffer[n]));
  }
  
  
  //-----------------------------------------
  // loop over all available AHF_halos files
  //-----------------------------------------
  for(n=0; n<num_files; n++) {
    
    fprintf(stderr,"working with file %s ",AHF_halos[n]);

    // read AHF_halos file into character halobuffer[]
    Nhalos = read_halos(AHF_halos[n], &halobuffer);
    
    // which redshift does this correspond to?
    if(zred == NULL)  z = get_redshift(AHF_halos[n]);
    else              z = zred[n];

    fprintf(stderr,"(Nhalos=%"PRIu64"): ",Nhalos);

    
    //-----------------------------------------
    // loop over all haloes to be traced
    //-----------------------------------------
#ifdef WITH_OPENMP
#pragma omp parallel for private(i,ihalo,outfile,fpout,k,haloline,id) shared(stderr,z,num_haloid,haloid,n,halobuffer,mtreebuffer,Nmtree,Nhalos) schedule(dynamic) default(none)
#endif
    for(i=0; i<num_haloid; i++) {
      
      // id of halo to be traced
      ihalo = haloid[i];
      
      // (re-)open outfile for this halo
      sprintf(outfile,"halo_%07"PRIu64".dat",ihalo);
      fpout = fopen(outfile,"a");
      if(fpout == NULL) {
        fprintf(stderr,"could not open %s\n",outfile);
        exit(0);
      }

      // find haloid corresponding to present AHF_halos file
      for(k=0; k<n; k++) {
        ihalo = get_haloidx(mtreebuffer[k], Nmtree[k], ihalo);
        
        // halo cannot be tracked further
        if(ihalo < 0) {
          break;
        }
      }
      fprintf(stderr,"%016"PRIi64" ",ihalo);

      if(ihalo >= 0) {
        // copy haloline over from halobuffer[]
#ifdef PROVIDE_LINENUMBER_AS_HALOID
        if(ihalo >= Nhalos) {
          fprintf(stderr,"hmmmm: ihalo=%"PRIu64" vs Nhalos=%"PRIu64"\n",ihalo,Nhalos);
        }
        strcpy(haloline,(halobuffer+ihalo)->line);
#else
        for(k=0; k<Nhalos; k++) {
          strcpy(haloline,(halobuffer+k)->line);
          sscanf(haloline,"%"SCNi64,&id);
          if(id == ihalo) {
            break;
          }
        }
#endif
        
        // and write to outfile
        fprintf(fpout,"%lf %s",z,haloline);
        fflush(fpout);
      } // if(ihalo >= 0)
      
      // close outfile for this halo
      fclose(fpout);
    } // for(num_haloid)
    
    // free AHF_halos[n] halobuffer[]
    for(ihalo=0; ihalo<Nhalos; ihalo++)
      free((halobuffer+ihalo)->line);
    free(halobuffer);
    halobuffer = NULL;
    
    fprintf(stderr,"\n");
  } //for(num_files)
  
  /*==================================================================*
   *                              BYE-BYE                             *
   *==================================================================*/
  free(haloid);
  for(i=0; i<num_files; i++) {
    free(AHF_halos[i]);
    free(mtree_idx[i]);
  }
  free(AHF_halos);
  free(mtree_idx);
  if(zred) {
    free(zred);
  }
  for(n=0; n<num_files-1; n++) {
    free(mtreebuffer[n]);
  }
  free(mtreebuffer);
  free(Nmtree);
  
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
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",prefix_list);
    exit(0);
  }

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
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",haloid_list);
    exit(0);
  }
  
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
  qsort(haloid, num_haloid, sizeof(uint64_t), &cmp_haloid);
}


/*==================================================================================================
 * get_headerline
 *==================================================================================================*/
void get_headerline(char *halofile, char *line)
{
  FILE *fp;
  int64_t id;
  
  fp = fopen(halofile,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",halofile);
    exit(0);
  }
  
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
 * get_haloidx_from_file
 *==================================================================================================*/
int64_t get_haloidx_from_file(char *mtreefile, uint64_t id)
{
  FILE *fp;
  int64_t ida, idx;
  char line[MAXSTRING];
  
  fp = fopen(mtreefile,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",mtreefile);
    exit(0);
  }
  
  ida = -1;
  while(ida != id && !feof(fp)) {
    fgets(line,MAXSTRING,fp);
    if(strncmp(line,"#",1) != 0)
      sscanf(line,"%"SCNi64,&ida);
    else
      ida = -1;
  }

  if(feof(fp)) {
    fclose(fp);
    return(-1);
  }
  else {
    sscanf(line,"%"SCNu64"%"SCNu64,&ida,&idx);
    fclose(fp);
    return(idx);
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

/*==================================================================================================
 * read_halos()
 *==================================================================================================*/
uint64_t read_halos(char *halofile, halofile_t **halobuffer)
{
  uint64_t Nhalos, i;
  FILE *fp;
  char  line[MAXSTRING];
  
  // open file
  fp = fopen(halofile,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",halofile);
    exit(0);
  }
  
  // ignore header line
  fgets(line,MAXSTRING,fp);

  // empty halobuffer
  *halobuffer = NULL;
  
  // and read it all in
  Nhalos = 0;
  while(!feof(fp)) {
    *halobuffer = realloc(*halobuffer, (Nhalos+1)*sizeof(halofile_t *));
    fgets(line,MAXSTRING,fp);
    ((*halobuffer)+Nhalos)->line = malloc((strlen(line)+1)*sizeof(char));
    strcpy(((*halobuffer)+Nhalos)->line,line);
    Nhalos++;
  }
  
  // close file
  fclose(fp);
  
#ifdef DEBUG
  for(i=0; i<Nhalos; i++) {
    fprintf(stderr,"%s",((*halobuffer)+i)->line);
  }
#endif
  
  // return number of lines
  return(Nhalos);
}

/*==================================================================================================
 * read_mtree_idx
 *==================================================================================================*/
uint64_t read_mtree_idx(char *mtreefile, mtreefile_t **mtreebuffer)
{
  FILE *fp;
  int64_t Nmtree, ida, idx;
  char line[MAXSTRING];
  
  fp = fopen(mtreefile,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",mtreefile);
    exit(0);
  }

  *mtreebuffer = NULL;
  Nmtree       = 0;
  
  while(!feof(fp)) {
    fgets(line,MAXSTRING,fp);
    if(strncmp(line,"#",1) != 0) {
      *mtreebuffer = (mtreefile_t *) realloc(*mtreebuffer, (Nmtree+1)*sizeof(mtreefile_t));
      sscanf(line,"%"SCNi64"%"SCNi64,&ida,&idx);
      ((*mtreebuffer)+Nmtree)->ida = ida;
      ((*mtreebuffer)+Nmtree)->idx = idx;
      Nmtree++;
    }
  }
  
  fclose(fp);
  
  return(Nmtree);
}

/*==================================================================================================
 * get_haloidx
 *==================================================================================================*/
int64_t get_haloidx(mtreefile_t *mtreebuffer, uint64_t Nmtree, uint64_t ihalo)
{
  uint64_t i, idx;
  
  i = 0;
  while(ihalo != mtreebuffer[i].ida && i < Nmtree) {
    i++;
  }
  
  if(i==Nmtree)
    return(-1);
  else
    return(mtreebuffer[i].idx);
}


