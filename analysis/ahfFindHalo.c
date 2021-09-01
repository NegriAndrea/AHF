/*==================================================================================================
 *  FindHalo:   read *.AHF_halos file and finds a certain halo based upon user supplied (X,Y,Z)
 *
 *
 *  input:    - 1x *.AHF_halos file
 *            - (X,Y,Z)
 *
 *  output:   - ID (and position) of halo closest to (X,Y,Z)
 *
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//#define DEBUG

#define MINPART         50000

#define MAXSTRING       2048

#define pow2(x)         ((x)*(x))
#define pow3(x)         ((x)*(x)*(x))

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct HALO *HALOptr;
typedef struct HALO
{   
  long   npart;
  
  /* read_halos() will convert lengths into grid units [0,LGRID-1]! */
  double Xc;
  double Yc;
  double Zc;
  double Vx;
  double Vy;
  double Vz;
  double Rvir;
  
  /* feel free to add whatever halo property you fancy... */
  /* ...but for the time being we are only concerned with position and radius */
  
}HALO;

/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/

HALOptr halo;
long    nhalos;
double  Xhalo, Yhalo, Zhalo, BoxSize;


/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
void read_halos     (char *halofile);

/*==================================================================================================
 * main:
 *
 *
 *==================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
  long   ihalo;
  double dX, dY, dZ;
  double d2, d2min;
  long   id;
  
  if(argc<6)
   {
    fprintf(stderr,"usage: %s .AHF_halos X Y Z BoxSize\n",*argv);
    exit(1);
   }
  
  printf("==================================================\n");
  printf("  Read *.AHF_halos file and identify close pairs\n");
  printf("==================================================\n");
  Xhalo   = (double)atof(argv[2]);
  Yhalo   = (double)atof(argv[3]);
  Zhalo   = (double)atof(argv[4]);
  BoxSize = (double)atof(argv[5]);
  
  /* read *.AHF_halos file (converts everything to internal units, too!) */
  read_halos(argv[1]);
  
  /* simply loop over all haloes picking the one with the minimum distance to (Xhalo,Yhalo,Zhalo) */
  id    = -1;
  d2min = 1E40;
  for(ihalo=0; ihalo<nhalos; ihalo++)
   {
    dX = fabs(halo[ihalo].Xc-Xhalo);
    dY = fabs(halo[ihalo].Yc-Yhalo);
    dZ = fabs(halo[ihalo].Zc-Zhalo);
    if(dX > BoxSize/2.) dX-=BoxSize;
    if(dY > BoxSize/2.) dY-=BoxSize;
    if(dZ > BoxSize/2.) dZ-=BoxSize;
    d2 = (pow2(dX)+pow2(dY)+pow2(dZ));
    if(d2 < d2min)
     {
      id    = ihalo;
      d2min = d2;
     }
   }
  
  /* dump result to screen */
  fprintf(stderr," o closest halo => %ld   %ld   %lf %lf %lf  %lf %lf %lf\n",
          id,halo[id].npart,
          halo[id].Xc,halo[id].Yc,halo[id].Zc,
          halo[id].Vx,halo[id].Vy,halo[id].Vz);
  
  /* be nice */
  free(halo);
  
  printf("STOP\n");
  return(1);
}

/*==============================================================================
 * read in initial halo positions from AHF_halos
 *==============================================================================*/
void read_halos(char *infile)
{
  FILE   *fpin;
  char    dummyline[MAXSTRING];
  long    ihalo, haloID, hostHalo, numSubStruct;
  double  Xc, Yc, Zc, VXc, VYc, VZc, npart, nvpart, Mvir, Rvir;
  
  
  fprintf(stderr," o reading _halos file %s ... ",infile);
  
  /* open AHF_halos file */
  if((fpin=fopen(infile,"r"))==NULL)
   {
    printf("I cannot open %s\n", infile);
    exit(1);
   }
  
  
  /* overread header line */
  fgets(dummyline,MAXSTRING,fpin);
  
  /* how many halos are there? */
  ihalo = 0;
  while(fgets(dummyline,MAXSTRING,fpin) != NULL)
    ihalo++;
  nhalos = ihalo;
  fprintf(stderr,"(nhalos = %ld) ... ",nhalos);
  
  /* rewind file back to start and overread header */
  rewind(fpin);
  fgets(dummyline,MAXSTRING,fpin);
  
  /* allocate memory for halos */
  halo = (HALOptr) calloc(nhalos, sizeof(HALO));
  
  /* eventually read halos file */
  for(ihalo=0; ihalo<nhalos; ihalo++)
   {
    /* read info from file */
    fgets(dummyline,MAXSTRING,fpin);
    
    /* extract information from last read dummyline */
    //sscanf(dummyline,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&npart,&nvpart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Mvir,&Rvir);
    sscanf(dummyline,"%ld %ld %ld %lf %lf %lf %lf %lf %lf %lf",&haloID,&hostHalo,&numSubStruct,&Mvir,&npart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Rvir);
    
    /* transfer info to halo structure */
    halo[ihalo].npart = (long) npart;
    halo[ihalo].Xc    = Xc;
    halo[ihalo].Yc    = Yc;
    halo[ihalo].Zc    = Zc;
    halo[ihalo].Rvir  = Rvir/1000.; /* Rvir in file is given in kpc/h */
    halo[ihalo].Vx    = VXc;
    halo[ihalo].Vy    = VYc;
    halo[ihalo].Vz    = VZc;      
   }
  
  fclose(fpin);
  
  fprintf(stderr,"done\n");
}
