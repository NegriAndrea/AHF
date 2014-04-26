/*==================================================================================================
 *  FindHaloPairs:   read *.AHF_halos file and identify halo pairs
 *
 *
 *  input:    - 1x *.AHF_halos file
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

#define MINPART         50000

#define MAXSTRING       2048

#define pow2(x)         ((x)*(x))
#define pow3(x)         ((x)*(x)*(x))

/* we use a LGRID^3 grid for a linked list (that grid is a 1d array!) */
#define LGRID           256
#define GridIdx(i,j,k)  ((i)+(j)*LGRID+(k)*LGRID*LGRID)

/*-------------------------------------------------------------------------------------
 *                                  THE STRUCTURES
 *-------------------------------------------------------------------------------------*/
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
   
   /* feel free to add whatever halo property you fancy... */
   /* ...but for the time being we are only concerned with position and radius */
   
}HALO;

typedef struct NODE *NODEptr;
typedef struct NODE
{
   /* head-of-chain pointer */
   HALOptr hoc;
}NODE;


typedef struct PAIR *PAIRptr;
typedef struct PAIR
{
   long id1;
   long id2;
}PAIR;


/*-------------------------------------------------------------------------------------
 *                                 GLOBAL VARIABLES
 *-------------------------------------------------------------------------------------*/

HALOptr halo;
NODEptr node;
PAIRptr pair;
long    nhalos;
long    npairs;
double  BoxSize;
double  Rfrac;
double  Mfrac;


/*-------------------------------------------------------------------------------------
 *                                 ALL THOSE FUNCTIONS
 *-------------------------------------------------------------------------------------*/
void read_halos     (char *halofile);
void ll             ();
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
   
   if(argc<5)
     {
      fprintf(stderr,"usage: %s .AHF_halos BoxSize Rfrac Mfrac\n",*argv);
      exit(1);
     }
   
   printf("==================================================\n");
   printf("  Read *.AHF_halos file and identify close pairs\n");
   printf("==================================================\n");
   BoxSize = (double)atof(argv[2]);
   Rfrac   = (double)atof(argv[3]);
   Mfrac   = (double)atof(argv[4]);

   /* read *.AHF_halos file (converts everything to internal units, too!) */
   read_halos(argv[1]);
   
   /* organize halos into a linked-list */
   ll();
   
   /* eventually find pairs */
   find_pairs();
   
   
   /* and write those pairs to file */
   write_pairs(argv[1]);
   
   free(pair);
   free(halo);
   free(node);
   
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
   nhalos = ihalo;
   fprintf(stderr,"(nhalos = %ld) ... ",nhalos);
   
   /* rewind file back to start and overread header */
   rewind(fpin);
   fgets(dummyline,MAXSTRING,fpin);
   
   /* allocate memory for halos */
   halo = (HALOptr) calloc(nhalos, sizeof(HALO));
   
   /* conversion factor to internal grid units [0,L-1] */
   x_fac = (double)(LGRID-1)/BoxSize;
   
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
      halo[ihalo].Xc    = Xc      * x_fac;
      halo[ihalo].Yc    = Yc      * x_fac;
      halo[ihalo].Zc    = Zc      * x_fac;
      halo[ihalo].Rvir  = Rvir    * x_fac/1000.; /* Rvir in file is given in kpc/h */
      halo[ihalo].Vx    = VXc;
      halo[ihalo].Vy    = VYc;
      halo[ihalo].Vz    = VZc;
      
      /* initialize linked-list, too */
      halo[ihalo].nic   = NULL;
     }
   
   fclose(fpin);
   
   fprintf(stderr,"done\n");
}

/*==================================================================================================
 * organize halos into a linked-list on a grid of dimension LGRID^3 
 *==================================================================================================*/
void ll()
{
   int     i,j,k;
   long    ihalo;
#ifdef DEBUG
   FILE   *fpout;
   HALOptr cur_halo;
#endif
   
   fprintf(stderr," o generatig linked-list ... ");
   
   /* allocate a 1D array consiting of LGRID^3 nodes */
   node = (NODEptr) calloc(pow3(LGRID), sizeof(NODE));
   
   /* initialize hoc's */
   for(i=0; i<LGRID; i++)
      for(j=0; j<LGRID-1; j++)
         for(k=0; k<LGRID-1; k++)
            node[GridIdx(i,j,k)].hoc = NULL;
   
   for(ihalo=0; ihalo<nhalos; ihalo++)
     {
      /* which node to attach halo to? */
      i = (int) (halo[ihalo].Xc+0.5);
      j = (int) (halo[ihalo].Yc+0.5);
      k = (int) (halo[ihalo].Zc+0.5);
      
#ifdef DEBUG
      /* double check node id's */
      if(i > LGRID-1 || i < 0) fprintf(stderr,"ll: strange i = %d\n",i);
      if(j > LGRID-1 || j < 0) fprintf(stderr,"ll: strange j = %d\n",j);
      if(k > LGRID-1 || k < 0) fprintf(stderr,"ll: strange k = %d\n",k);
#endif
      
      /* attach to (existing) linked-list */
      if(node[GridIdx(i,j,k)].hoc == NULL)
        {
         node[GridIdx(i,j,k)].hoc = (halo+ihalo);
        }
      else
        {
         halo[ihalo].nic          = node[GridIdx(i,j,k)].hoc;
         node[GridIdx(i,j,k)].hoc = (halo+ihalo);
        }
     }
   
#ifdef DEBUG
   /* check linked-list */
   fpout = fopen("test.geom","w");
   for(i=0; i<LGRID; i++)
      for(j=0; j<LGRID; j++)
         for(k=0; k<LGRID; k++)
            for(cur_halo=node[GridIdx(i,j,k)].hoc; cur_halo != NULL; cur_halo=cur_halo->nic)
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
   HALOptr cur_halo;
   
   fprintf(stderr," o identifying pairs ... ");

   /* initialize pair counter */
   npairs = 0;
   pair   = NULL;
   
   /* loop over all halos */
   for(ihalo=0; ihalo<nhalos; ihalo++)
     {
      /* locate linked-list nodes to search for pair candidates */
      Rmax = Rfrac*halo[ihalo].Rvir;
      imin = (int) (halo[ihalo].Xc   - Rmax);
      jmin = (int) (halo[ihalo].Yc   - Rmax);
      kmin = (int) (halo[ihalo].Zc   - Rmax);
      imax = (int) (halo[ihalo].Xc+1 + Rmax);
      jmax = (int) (halo[ihalo].Yc+1 + Rmax);
      kmax = (int) (halo[ihalo].Zc+1 + Rmax);
      
      for(ix=imin; ix<=imax; ix++)
         for(jy=jmin; jy<=jmax; jy++)
            for(kz=kmin; kz<=kmax; kz++)
              {
               i = (ix+LGRID)%LGRID;
               j = (jy+LGRID)%LGRID;
               k = (kz+LGRID)%LGRID;
               
               
               for(cur_halo=node[GridIdx(i,j,k)].hoc; cur_halo != NULL; cur_halo=cur_halo->nic)
                 {
                  /* crtierion to be considered as a "pair halo" */
                  if(fabs(halo[ihalo].npart-cur_halo->npart)/(double)halo[ihalo].npart < Mfrac &&
                     (ihalo != cur_halo-halo))
                    {
                     /* we found a pair halo! */
                     npairs++;
                     
                     /* make room for that new pair */
                     pair = (PAIRptr) realloc(pair, npairs*sizeof(PAIR));
                     
                     /* ...and store the halo id's */
                     pair[npairs-1].id1 = ihalo;
                     pair[npairs-1].id2 = cur_halo-halo;
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
   long   id1, id2;
   double x_fac;
   
   strcpy(outfile,outprefix);
   strcat(outfile,"-pairs");
   
   fprintf(stderr," o writing %ld pairs to file %s ... ",npairs,outfile);
   fpout = fopen(outfile,"w");
   
   /* convert positions back to Mpc/h */
   x_fac = BoxSize/(double)(LGRID-1);
   
   for(ipair=0; ipair<npairs; ipair++)
     {
      id1 = pair[ipair].id1;
      id2 = pair[ipair].id2;
      if(halo[id1].npart > MINPART)
        {
         fprintf(fpout,"%12ld %12ld    %12ld %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g   %12ld %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n",
                 id1,id2,
                 halo[id1].npart,halo[id1].Xc*x_fac,halo[id1].Yc*x_fac,halo[id1].Zc*x_fac,
                                 halo[id1].Vx,      halo[id1].Vy,      halo[id1].Vz,
                 halo[id2].npart,halo[id2].Xc*x_fac,halo[id2].Yc*x_fac,halo[id2].Zc*x_fac,
                                 halo[id2].Vx,      halo[id2].Vy,      halo[id2].Vz);
        }
     }
   
   fclose(fpout);

   fprintf(stderr,"done\n");
}
