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

typedef struct MergerTree {
  uint64_t haloid;
  uint64_t progid;
} mtree_t;


/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/
char **AHF_halos;
char **mtree_idx;
char **outfile;

/*-------------------------------------------------------------------------------------
 *                                    PROTOTYPES
 *-------------------------------------------------------------------------------------*/
int      get_filename_list         (char *);
double   get_redshift              (char *);
double  *get_zred                  (char *, int);
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
  double *zred;
  int i, n, num_files, num_halos[2], num_mtree;
  uint64_t ihalo, isub, haloid, subhaloid, subhostHalo, hostHalo, subprogid, haloprogid, ipos;
  mtree_t *itmp_mtree;
  halo_t  *itmp_halo;
  halo_t  *host_halo, *sub_halo, *host_prog, *sub_prog;
  double   z, absAngMom;
  
  // user provided input
  char prefix_list[MAXSTRING], zred_list[MAXSTRING];

  // the pointers to the arrays holding the 2 halos and 1 mtree file
  halo_t  *halos[2];
  mtree_t *mtree;
  
  
  /*==================================================================*
   *                          USER INTERFACE                          *
   *==================================================================*/
  if(argc<2)
   {
    fprintf(stderr,"usage: %s prefix_list zred_list\n",*argv);
    exit(1);
   }
  strcpy(prefix_list,      argv[1]);
  strcpy(zred_list, argv[2]);

  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"               follow AHF_halos through _mtree_idx files\n");
  fprintf(stderr,"=======================================================================\n");
  fprintf(stderr,"will read mtree_idx and AHF_halos files using prefixes found in %s\n",prefix_list);
  fprintf(stderr,"will read mapping between snapid and redshift from              %s\n\n",zred_list);
  
  /*==================================================================*
   *                       DEAL WITH USER DATA                        *
   *==================================================================*/
  num_files = get_filename_list(prefix_list);
  zred      = get_zred(zred_list, num_files); // returns NULL, if zred_list = "-1" (and then zred will be extracted from prefix_list)
  
  
  /*==================================================================*
   *                          DO THE WORK                             *
   *==================================================================*/
  fprintf(stderr,"\n");
  
  // start-up
  halos[0] = read_halos(AHF_halos[0], &(num_halos[0]));
  qsort((void *)halos[0], num_halos[0], sizeof(halo_t), qcompareHaloIDs);
  fprintf(stderr,"read and sorted %d haloes from file %s\n",num_halos[0],AHF_halos[i]);

  // loop over all file pairs
  for(i=1; i<num_files-1; i++) {
    
    //------------------
    // READING OF FILES
    //------------------
    
    // read second halos file
    halos[1] = read_halos(AHF_halos[i], &(num_halos[1]));
    qsort((void *)halos[1], num_halos[1], sizeof(halo_t), qcompareHaloIDs);
    fprintf(stderr,"read and sorted %d haloes from file %s\n",num_halos[1],AHF_halos[i]);
    
    // no point continuing when there are no halos
    if(num_halos[1] == 0)
      break;
    
    // read the connecting mtree file
    mtree  = read_mtree_idx(mtree_idx[i-1], &num_mtree);
    qsort((void *)mtree, num_mtree, sizeof(mtree_t), qcompareMtreeIDs);
    fprintf(stderr,"read and sorted %d entries from mtree file %s\n",num_mtree,mtree_idx[i-1]);
    
    //------------------
    // OPEN OUTPUT FILE
    //------------------
    fprintf(stderr," -> writing to %s\n\n",outfile[i]);
    fpout = fopen(outfile[i],"w");
    assert(fpout);
    
    // the header line
//    fprintf(fpout,"# haloid(1) Mvir(2) Rvir(3)"
//                   " haloprogid(4) Mvir(5) Rvir(6) Xc(7) Yc(8) Zc(9)"
//                   " subhaloid(10) Mvir(11) Rvir(12)"
//                   " subhaloprogid(13) Mvir(14) Rvir(15)"
//                   " Xsubhaloprog-Xhaloprog(16) Ysubhaloprog-Yhaloprog(17) Zsubhaloprog-Zhaloprog(18)"
//                   " VXsubhaloprog-VXhaloprog(19) VYsubhaloprog-VYhaloprog(20) VZsubhaloprog-VZhaloprog(21)\n");
    
    fprintf(fpout,
            "# haloid(1) Mvir(2) Rvir(3) Lx(4) Ly(5) Lz(6)"
            " haloprogid(7) Mvir(8) Rvir(9) Xc(10) Yc(11) Zc(12) Lx(13) Ly(14) Lz(15)"
            " subhaloid(16) Mvir(17) Rvir(18) Lx(19) Ly(20) Lz(21)"
            " subhaloprogid(22) Mvir(23) Rvir(24)"
            " Xsubhaloprog-Xhaloprog(25) Ysubhaloprog-Yhaloprog(26) Zsubhaloprog-Zhaloprog(27)"
            " VXsubhaloprog-VXhaloprog(28) VYsubhaloprog-VYhaloprog(29) VZsubhaloprog-VZhaloprog(30) Lx(31) Ly(32) Lz(33)\n");
    
    //------------------
    // WORK-LOOP
    //------------------
    
#ifdef WITH_OPENMP
#pragma omp parallel for private(ihalo, haloid, hostHalo, host_halo, itmp_mtree, itmp_halo, ipos, haloprogid, host_prog, isub, subhaloid, subhostHalo, sub_halo, subprogid, sub_prog) shared(halos, mtree, num_halos, num_mtree, fpout)
#endif
    // loop over all relevant halos in halos[0]
    for(ihalo=0; ihalo<num_halos[0]; ihalo++) {
      
      haloid   = (halos[0]+ihalo)->haloid;
      hostHalo = (halos[0]+ihalo)->hostHalo;
      
      // provide convenient name for haloid
      host_halo = (halos[0]+ihalo);
      
      // we are only interested in host halos
      if(hostHalo == 0) {
        
        // locate haloid in mtree[] array
        itmp_mtree = (mtree_t *)bsearch(&(haloid), mtree, num_mtree, sizeof(mtree_t), bcompareMtreeIDs);
        
        
        // only makes sense to continue if haloid does have a progenitor
        if(itmp_mtree != NULL) {
          ipos       = (uint64_t) (itmp_mtree-mtree);
          haloprogid = mtree[ipos].progid;
          
          // locate haloprogid in halos[1] array
          itmp_halo = (halo_t *)bsearch(&(haloprogid), halos[1], num_halos[1], sizeof(halo_t), bcompareHaloIDs);
          if(itmp_halo == NULL) {
            fprintf(stderr,"could not find haloprogid %"PRIu64" in halos[1] array!?\nABORTING\n",haloprogid);
            exit(0);
          }
          // provide convenient name for haloprogid
          host_prog = itmp_halo;
          
          // hunt down the subhalos of host_halo with haloid
          for(isub=0; isub<num_halos[0]; isub++) {
            subhaloid   = (halos[0]+isub)->haloid;
            subhostHalo = (halos[0]+isub)->hostHalo;
            
            // we found a subhalo
            if(subhostHalo == haloid) {
              
              // locate subhaloid in halos[0] array
              itmp_halo = (halo_t *)bsearch(&(subhaloid), halos[0], num_halos[0], sizeof(halo_t), bcompareHaloIDs);
              if(itmp_halo == NULL) {
                fprintf(stderr,"could not find subhaloid %"PRIu64" in halos[1] array!?\nABORTING\n",subhaloid);
                exit(0);
              }
              // provide convenient name for subhaloid
              sub_halo = itmp_halo;
              
              // locate subhaloid in mtree[]
              itmp_mtree = (mtree_t *)bsearch(&(subhaloid), mtree, num_mtree, sizeof(mtree_t), bcompareMtreeIDs);
              
              // only makes sense to continue if subhaloid does have a progenitor
              if(itmp_mtree != NULL) {
                ipos      = (uint64_t)(itmp_mtree-mtree);
                subprogid = mtree[ipos].progid;
                
                // locate subprogid in halos[1] array
                itmp_halo = (halo_t *)bsearch(&(subprogid), halos[1], num_halos[1], sizeof(halo_t), bcompareHaloIDs);
                if(itmp_halo == NULL) {
                  fprintf(stderr,"could not find subprogid %"PRIu64" in halos[1] array!?\nABORTING\n",subprogid);
                  exit(0);
                }
                ipos = (uint64_t)(itmp_halo-halos[1]);
                
                // provide convenient name for subprogid
                sub_prog  = (halos[1]+ipos);
                
                // we are only interested in it if it is a field halo (and not pointing to the host)
                if( sub_prog->hostHalo == 0 && subprogid != haloprogid ) {
                  
                  //-----------------------------------------------------------------------------------------
                  // WRITE RELEVANT INFORMATION TO FILE
                  //
                  // host_halo->    gives access to host halo from file AHF_halos[0]
                  // host_prog->    gives access to host halo from file AHF_halos[1]
                  //
                  // sub_halo->     gives access to  sub halo from file AHF_halos[0]
                  // sub_prog->     gives access to  sub halo from file AHF_halos[1]
                  //
                  //-----------------------------------------------------------------------------------------
//                  fprintf(fpout,
//                          "%"PRIu64" %g %f %"PRIu64" %g %f %f %f %f %"PRIu64" %g %f %"PRIu64" %g %f %f %f %f %f %f %f\n",
//                          host_halo->haloid,
//                          host_halo->Mvir,
//                          host_halo->Rvir,
//                          //
//                          host_prog->haloid,
//                          host_prog->Mvir,
//                          host_prog->Rvir,
//                          host_prog->Xc,
//                          host_prog->Yc,
//                          host_prog->Zc,
//                          //
//                          sub_halo->haloid,
//                          sub_halo->Mvir,
//                          sub_halo->Rvir,
//                          //
//                          sub_prog->haloid,
//                          sub_prog->Mvir,
//                          sub_prog->Rvir,
//                          sub_prog->Xc  - host_prog->Xc,
//                          sub_prog->Yc  - host_prog->Yc,
//                          sub_prog->Zc  - host_prog->Zc,
//                          sub_prog->VXc - host_prog->VXc,
//                          sub_prog->VYc - host_prog->VYc,
//                          sub_prog->VZc - host_prog->VZc
//                          );
                  fprintf(fpout,
                          "%"PRIu64" %g %f %g %g %g %"PRIu64" %g %f %f %f %f %g %g %g %"PRIu64" %g %f %g %g %g %"PRIu64" %g %f %f %f %f %f %f %f %g %g %g\n",
                          host_halo->haloid,
                          host_halo->Mvir,
                          host_halo->Rvir,
                          host_halo->Lx * sqrt(2.)*host_halo->lambda*host_halo->Mvir*sqrt(Grav*host_halo->Mvir*host_halo->Rvir/1000.),
                          host_halo->Ly * sqrt(2.)*host_halo->lambda*host_halo->Mvir*sqrt(Grav*host_halo->Mvir*host_halo->Rvir/1000.),
                          host_halo->Lz * sqrt(2.)*host_halo->lambda*host_halo->Mvir*sqrt(Grav*host_halo->Mvir*host_halo->Rvir/1000.),
                          //
                          host_prog->haloid,
                          host_prog->Mvir,
                          host_prog->Rvir,
                          host_prog->Xc,
                          host_prog->Yc,
                          host_prog->Zc,
                          host_prog->Lx * sqrt(2.)*host_prog->lambda*host_prog->Mvir*sqrt(Grav*host_prog->Mvir*host_prog->Rvir/1000.),
                          host_prog->Ly * sqrt(2.)*host_prog->lambda*host_prog->Mvir*sqrt(Grav*host_prog->Mvir*host_prog->Rvir/1000.),
                          host_prog->Lz * sqrt(2.)*host_prog->lambda*host_prog->Mvir*sqrt(Grav*host_prog->Mvir*host_prog->Rvir/1000.),
                          //
                          sub_halo->haloid,
                          sub_halo->Mvir,
                          sub_halo->Rvir,
                          sub_halo->Lx * sqrt(2.)*sub_halo->lambda*sub_halo->Mvir*sqrt(Grav*sub_halo->Mvir*sub_halo->Rvir/1000.),
                          sub_halo->Ly * sqrt(2.)*sub_halo->lambda*sub_halo->Mvir*sqrt(Grav*sub_halo->Mvir*sub_halo->Rvir/1000.),
                          sub_halo->Lz * sqrt(2.)*sub_halo->lambda*sub_halo->Mvir*sqrt(Grav*sub_halo->Mvir*sub_halo->Rvir/1000.),
                          //
                          sub_prog->haloid,
                          sub_prog->Mvir,
                          sub_prog->Rvir,
                          sub_prog->Xc  - host_prog->Xc,
                          sub_prog->Yc  - host_prog->Yc,
                          sub_prog->Zc  - host_prog->Zc,
                          sub_prog->VXc - host_prog->VXc,
                          sub_prog->VYc - host_prog->VYc,
                          sub_prog->VZc - host_prog->VZc,
                          sub_prog->Lx * sqrt(2.)*sub_prog->lambda*sub_prog->Mvir*sqrt(Grav*sub_prog->Mvir*sub_prog->Rvir/1000.),
                          sub_prog->Ly * sqrt(2.)*sub_prog->lambda*sub_prog->Mvir*sqrt(Grav*sub_prog->Mvir*sub_prog->Rvir/1000.),
                          sub_prog->Lz * sqrt(2.)*sub_prog->lambda*sub_prog->Mvir*sqrt(Grav*sub_prog->Mvir*sub_prog->Rvir/1000.)
                          );
                  fflush(fpout);
                  
                }
                
              } // if(itmp)
            } // if(hosthaloid==haloid)
          } // for(isub)
        } // if(itmp)
      } // if(hostHalo)
    } // for(ihalo)
    
    //------------------
    // WORK-LOOP
    //------------------
    
    
    
    
    //-------------------
    // CLOSE OUTPUT FILE
    //-------------------
    fclose(fpout);
    
    //-------------
    // FREE MEMORY
    //-------------
    if(mtree) {
      free(mtree);
      mtree = NULL;
    }
    if(halos[0]) {
      free(halos[0]);
      halos[0] = NULL;
    }
    
    //--------------------------------
    // make halos[1] the new halos[0]
    //--------------------------------
    halos[0]     = halos[1];
    num_halos[0] = num_halos[1];
  }
  
  /*==================================================================*
   *                              BYE-BYE                             *
   *==================================================================*/
  // free remaining memory
  if(halos[0]) free(halos[0]);
  if(halos[1]) free(halos[1]);
  if(mtree) free(mtree);
        
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
int get_filename_list(char *prefix_list)
{
  FILE *fp;
  char line[MAXSTRING];
  int  i, num_files;
  
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
  outfile   = (char **) calloc(num_files, sizeof(char *));
  
  // eventually read haloids
  fprintf(stderr," + will use the following files for the tracing:\n");
  rewind(fp);
  for(i=0; i<num_files; i++) {
    fscanf(fp,"%s",line);
    AHF_halos[i] = (char *) calloc(MAXSTRING,sizeof(char));
    mtree_idx[i] = (char *) calloc(MAXSTRING,sizeof(char));
    outfile[i]   = (char *) calloc(MAXSTRING,sizeof(char));
    sprintf(AHF_halos[i],"%s_halos",line);
    sprintf(mtree_idx[i],"%s_mtree_idx",line);
    sprintf(outfile[i],  "%s_subaccretion",basename(line));
    if(i==num_files-1)
      fprintf(stderr,"     %s\n", AHF_halos[i]);
    else
      fprintf(stderr,"     %s, %s\n", AHF_halos[i], mtree_idx[i]);      
  }
  
  // close file
  fclose(fp);
  
  // return number of files
  return(num_files);
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

/*==================================================================================================
 * get_zred
 *==================================================================================================*/
double *get_zred(char *zred_list, int num_files)
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
int bcompareMtreeIDsOLD(const void *mtree, const void *id)
{
  const uint64_t haloid = ((mtree_t *)mtree)->haloid;
	const uint64_t *i     = (const uint64_t *)id;
  
	return *i > haloid ? -1 : (*i < haloid ? 1 : 0);
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
int bcompareHaloIDsOLD(const void *halo, const void *id)
{
  const uint64_t haloid = ((halo_t *)halo)->haloid;
	const uint64_t *i    = (const uint64_t *)id;
  
	return *i > haloid ? -1 : (*i < haloid ? 1 : 0);
}




