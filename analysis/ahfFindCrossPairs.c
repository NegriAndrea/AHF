/*==================================================================================================
 *  FindCrossPairs:   reads 2 *.AHF_halos files and cross-identifies halo pairs
 *
 *
 *        the code tries to find a counterpart for each halo in .AHF_halos1 in .AHF_halos2
 *
 *
 *
 *  input:    - 2x *.AHF_halos file
 *            - box size of simulation
 *            - fraction of virial radii for "pair definition"
 *            - fraction of masses for "pair definition"
 *
 *  output:   - 1x file containing halo id's of pairs
 *
 *==================================================================================================*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//#define DEBUG

#define Rfrac           1.0
#define Mfrac           2.0
#define MINPART         100

#define MAXSTRING       2048

#define pow2(x)         ((x)*(x))
#define pow3(x)         ((x)*(x)*(x))

/* we use a LGRID^3 grid for a linked list (that grid is a 1d array!) */
#define LGRID           256
#define GridIdx(i,j,k)  ((i)+(j)*LGRID+(k)*LGRID*LGRID)

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
typedef struct PAIR *PAIRptr;
typedef struct PAIR
{
  long id1;
  long id2;
}PAIR;

typedef struct HALO *HALOptr;
typedef struct HALO
{
  /* next-in-chain pointer*/
  HALOptr nic;
  
  long   npart;
  
  /* read_halos() will convert lengths into grid units [0,LGRID-1]! */
  double Xc;
  double Yc;
  double Zc;
  double Vx;
  double Vy;
  double Vz;
  double Rvir;
  double Mvir;
  
  /* feel free to add whatever halo property you fancy... */
  /* ...but for the time being we are only concerned with position and radius */
  
  long    npairs;
  PAIRptr pair;
}HALO;

typedef struct NODE *NODEptr;
typedef struct NODE
{
  /* head-of-chain pointer */
  HALOptr hoc;
}NODE;


/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/

long    nhalos[2];
HALOptr halo[2];
NODEptr node[2];
//PAIRptr pair;
//long    npairs;
double  BoxSize;


/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
void read_halos     (char *halofile, int ifile);
void ll             (int ifile);
void find_pairs     ();
void write_pairs    (char *outprefix);

/*==================================================================================================
 * main:
 *
 *
 *==================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
  
  if(argc<4)
    {
    fprintf(stderr,"usage: %s .AHF_halos1 .AHF_halos2 BoxSize\n",*argv);
    exit(1);
    }
  
  printf("=====================================================\n");
  printf("  Read 2 *.AHF_halos files and identify close pairs\n");
  printf("=====================================================\n");
  BoxSize = (double)atof(argv[3]);
  
  /* read *.AHF_halos file (converts everything to internal units, too!) */
  read_halos(argv[1], 0);
  read_halos(argv[2], 1);
  
  /* organize halos into a linked-list */
  ll(0);
  ll(1);
  
  /* eventually find pairs */
  find_pairs();
  
  /* and write those pairs to file */
  write_pairs(argv[1]);
  
  
  /* be nice and gentle */
  free(halo[0]);
  free(halo[1]);
  free(node[0]);
  free(node[1]);
  
  printf("STOP\n");
  return(1);
}

/*==============================================================================
 * read in initial halo positions from AHF_halos
 *==============================================================================*/
void read_halos(char *infile, int ifile)
{
  FILE   *fpin;
  char    dummyline[MAXSTRING];
  long    ihalo, haloID, hostHalo, numSubStruct;
  double  Xc, Yc, Zc, VXc, VYc, VZc, npart, nvpart, Mvir, Rvir;
  double  x_fac;
  
  
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
  nhalos[ifile] = ihalo;
  fprintf(stderr,"(nhalos = %ld) ... ",nhalos[ifile]);
  
  /* rewind file back to start and overread header */
  rewind(fpin);
  fgets(dummyline,MAXSTRING,fpin);
  
  /* allocate memory for halos */
  halo[ifile] = (HALOptr) calloc(nhalos[ifile], sizeof(HALO));
  
  /* conversion factor to internal grid units [0,L-1] */
  x_fac = (double)(LGRID-1)/BoxSize;
  
  /* eventually read halos file */
  for(ihalo=0; ihalo<nhalos[ifile]; ihalo++)
  {
    /* read info from file */
    fgets(dummyline,MAXSTRING,fpin);
    
    /* extract information from last read dummyline */
     //sscanf(dummyline,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&npart,&nvpart,&Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Mvir,&Rvir);
    sscanf(dummyline,"%ld %ld %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	   &haloID,&hostHalo,&numSubStruct,&Mvir,&npart,
	   &Xc,&Yc,&Zc,&VXc,&VYc,&VZc,&Rvir);

    /* transfer info to halo structure */
    halo[ifile][ihalo].npart = (long) npart;
    halo[ifile][ihalo].Xc    = Xc      * x_fac;
    halo[ifile][ihalo].Yc    = Yc      * x_fac;
    halo[ifile][ihalo].Zc    = Zc      * x_fac;
    halo[ifile][ihalo].Rvir  = Rvir    * x_fac; /*/1000.; /* Rvir in file is given in kpc/h */
    halo[ifile][ihalo].Mvir  = Mvir;
    halo[ifile][ihalo].Vx    = VXc;
    halo[ifile][ihalo].Vy    = VYc;
    halo[ifile][ihalo].Vz    = VZc;
    
    /* initialize linked-list, too */
    halo[ifile][ihalo].nic   = NULL;
    }
  
  fclose(fpin);
  
  fprintf(stderr,"done\n");
}

/*==================================================================================================
 * organize halos into a linked-list on a grid of dimension LGRID^3 
 *==================================================================================================*/
void ll(int ifile)
{
  int     i,j,k;
  long    ihalo;
#ifdef DEBUG
  FILE   *fpout;
  HALOptr cur_halo;
#endif
  
  fprintf(stderr," o generating linked-list for file %5d ... ",ifile);
  
  /* allocate a 1D array consiting of LGRID^3 nodes */
  node[ifile] = (NODEptr) calloc(pow3(LGRID), sizeof(NODE));
  
  /* initialize hoc's */
  for(i=0; i<LGRID; i++)
    for(j=0; j<LGRID-1; j++)
      for(k=0; k<LGRID-1; k++)
        node[ifile][GridIdx(i,j,k)].hoc = NULL;
  
  for(ihalo=0; ihalo<nhalos[ifile]; ihalo++)
    {
    /* which node to attach halo to? */
    i = (int) (halo[ifile][ihalo].Xc+0.5);
    j = (int) (halo[ifile][ihalo].Yc+0.5);
    k = (int) (halo[ifile][ihalo].Zc+0.5);
    
#ifdef DEBUG
    /* double check node id's */
    if(i > LGRID-1 || i < 0) fprintf(stderr,"ll: strange i = %d\n",i);
    if(j > LGRID-1 || j < 0) fprintf(stderr,"ll: strange j = %d\n",j);
    if(k > LGRID-1 || k < 0) fprintf(stderr,"ll: strange k = %d\n",k);
#endif
    
    /* attach to (existing) linked-list */
    if(node[ifile][GridIdx(i,j,k)].hoc == NULL)
      {
      node[ifile][GridIdx(i,j,k)].hoc = (halo[ifile]+ihalo);
      }
    else
      {
      halo[ifile][ihalo].nic          = node[ifile][GridIdx(i,j,k)].hoc;
      node[ifile][GridIdx(i,j,k)].hoc = (halo[ifile]+ihalo);
      }
    }
  
#ifdef DEBUG
  /* check linked-list */
  fpout = fopen("test.geom","w");
  for(i=0; i<LGRID; i++)
    for(j=0; j<LGRID; j++)
      for(k=0; k<LGRID; k++)
        for(cur_halo=node[ifile][GridIdx(i,j,k)].hoc; cur_halo != NULL; cur_halo=cur_halo->nic)
          fprintf(fpout,"s %g %g %g %g 1 0 0\n",
                  cur_halo->Xc,
                  cur_halo->Yc,
                  cur_halo->Zc,
                  cur_halo->Rvir);
  
  fclose(fpout);
#endif
  
  fprintf(stderr,"done\n");
}

/*==================================================================================================
 * find halo pairs 
 *==================================================================================================*/
void find_pairs()
{
  long    ihalo;
  int     i,j,k;
  int     ix, jy, kz;
  int     imin,imax,jmin,jmax,kmin,kmax;
  double  Rmax;
  double  dx, dy, dz;
  HALOptr cur_halo;
  
  fprintf(stderr," o identifying pairs ... ");
  fflush(stderr);
  /* loop over all halos in the first file (ifile=0) */
  for(ihalo=0; ihalo<nhalos[0]; ihalo++)
    {
    /* progress meter */
	if( !(ihalo % (int)(nhalos[0]/10.)) ) {
	    fprintf(stderr,"%4.1f ",(float)ihalo/(float)nhalos[0]);
	    fflush(stderr);
	}

    /* initialize pair counter */
    halo[0][ihalo].npairs = 0;
    halo[0][ihalo].pair   = NULL;
    
    /* locate linked-list nodes to search for pair candidates */
    Rmax = Rfrac *halo[0][ihalo].Rvir;
    imin = (int) (halo[0][ihalo].Xc   - Rmax);
    jmin = (int) (halo[0][ihalo].Yc   - Rmax);
    kmin = (int) (halo[0][ihalo].Zc   - Rmax);
    imax = (int) (halo[0][ihalo].Xc+1 + Rmax);
    jmax = (int) (halo[0][ihalo].Yc+1 + Rmax);
    kmax = (int) (halo[0][ihalo].Zc+1 + Rmax);
    
    for(ix=imin; ix<=imax; ix++)
      for(jy=jmin; jy<=jmax; jy++)
        for(kz=kmin; kz<=kmax; kz++)
          {
          i = (ix+LGRID)%LGRID;
          j = (jy+LGRID)%LGRID;
          k = (kz+LGRID)%LGRID;
          
          /* we can use the same indices for the linked-list in the other file (ifile=1)
           * as both simulations have been covered with the same grid */
          for(cur_halo=node[1][GridIdx(i,j,k)].hoc; cur_halo != NULL; cur_halo=cur_halo->nic)
            {

            /* calculate distance of centres */
            dx = fabs(halo[0][ihalo].Xc-cur_halo->Xc);
            dy = fabs(halo[0][ihalo].Yc-cur_halo->Yc);
            dz = fabs(halo[0][ihalo].Zc-cur_halo->Zc);
            if(dx > BoxSize/2.) dx = dx-BoxSize;
            if(dy > BoxSize/2.) dy = dy-BoxSize;
            if(dz > BoxSize/2.) dz = dz-BoxSize;
            
            /* criterion to be considered as a "pair halo" */
            if(
               (pow2(dx)+pow2(dy)+pow2(dz))/pow2(halo[0][ihalo].Rvir)                    < Rfrac &&  // spatially close
               halo[0][ihalo].Vx*cur_halo->Vx +
               halo[0][ihalo].Vy*cur_halo->Vy +                                                      // flying in the same direction
               halo[0][ihalo].Vz*cur_halo->Vz > 0                                                &&
               fabs(halo[0][ihalo].Rvir -cur_halo->Rvir) /(double)halo[0][ihalo].Rvir    < Rfrac &&  // similar virial radii
               fabs(halo[0][ihalo].Mvir-cur_halo->Mvir)/(double)halo[0][ihalo].Mvir      < Mfrac     // similar masses
               )
              {
              /* we found a pair halo! */
              halo[0][ihalo].npairs++;
              
              /* make room for that new pair */
              halo[0][ihalo].pair = (PAIRptr) realloc(halo[0][ihalo].pair, halo[0][ihalo].npairs*sizeof(PAIR));
              
              /* ...and store the halo id's */
              halo[0][ihalo].pair[halo[0][ihalo].npairs-1].id1 = ihalo;
              halo[0][ihalo].pair[halo[0][ihalo].npairs-1].id2 = cur_halo-halo[1];
              }
            }
          }
    }
  
  fprintf(stderr,"done\n");
}


/*==================================================================================================
 * write halo pairs 
 *==================================================================================================*/
void write_pairs(char *outprefix)
{
  char   outfile[MAXSTRING];
  FILE  *fpout;
  int    ipair;
  long   ihalo, idcmp, id2, npairs;
  double x_fac;
  double  dQ, dQmin;
  
  strcpy(outfile,outprefix);
  strcat(outfile,"-pairs");
  
  fprintf(stderr," o writing pairs to file %s ... ",outfile);
  fpout = fopen(outfile,"w");
  
  /* convert positions back to Mpc/h */
  x_fac = BoxSize/(double)(LGRID-1);
  
  /* count total number of pairs */
  npairs = 0;
  
  for(ihalo=0; ihalo<nhalos[0]; ihalo++)
    {
    if(halo[0][ihalo].npart > MINPART && halo[0][ihalo].npairs > 0)
      {
      npairs++;
      dQmin = 1E40;
      
      for(ipair=0; ipair<halo[0][ihalo].npairs; ipair++)
        {
        /* id of companion */
        id2 = halo[0][ihalo].pair[ipair].id2;
        
        /* whatever difference to define the companion (radius proofed best) */
        dQ = fabs(halo[0][ihalo].Rvir-halo[1][id2].Rvir);
        
        if(dQ < dQmin)
          {
          dQmin = dQ;
          idcmp = id2;
          }
        }
      
      fprintf(fpout,"%12ld %12ld    %12ld %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g  %12ld %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n",
              ihalo,idcmp,
              halo[0][ihalo].npart,halo[0][ihalo].Xc*x_fac,halo[0][ihalo].Yc*x_fac,halo[0][ihalo].Zc*x_fac,
              halo[0][ihalo].Vx,   halo[0][ihalo].Vy,      halo[0][ihalo].Vz,
              halo[0][ihalo].Mvir, halo[0][ihalo].Rvir*x_fac,
              halo[1][idcmp].npart,halo[1][idcmp].Xc*x_fac,halo[1][idcmp].Yc*x_fac,halo[1][idcmp].Zc*x_fac,
              halo[1][idcmp].Vx,   halo[1][idcmp].Vy,      halo[1][idcmp].Vz,
              halo[1][idcmp].Mvir, halo[1][idcmp].Rvir*x_fac);      
      }
    }
  
  fclose(fpout);
  
  fprintf(stderr,"done (npairs=%ld)\n",npairs);
}
