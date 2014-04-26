#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*=============================================================================
 *                                   DEFINES
 *=============================================================================*/
#define MAXSTRING   2048
#define GADGET_SKIP        ReadUInt(icfile,&blklen,SWAPBYTES);
#define SIZEOFGADGETHEADER 256
#define MZERO             (1e-10)
#define X                  0
#define Y                  1
#define Z                  2
#define PI           3.14159265358979323846264338
#define TWOPI       (2*PI)

#define TRUE         1
#define FALSE        0

#define pow2(x)  ((x)*(x))
#define pow3(x)  ((x)*(x)*(x))
#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))

/*=============================================================================
 *                                   COMMONS
 *=============================================================================*/
long  nhaloes;
long *halo;

double       GADGET_LUNIT    = 1.E-3;
double       GADGET_MUNIT    = 1.E+10;
int          SWAPBYTES       = 1;
int          FORMAT          = 2;
int          LGADGET         = 0;

long         IDmin =  4294967296;
long         IDmax = -4294967296; //simply for debugging purposes

unsigned int blklen;
struct info_gadget
{
  int      no_gadget_files;
  int      i_gadget_file;
  long   *(np[6]);
  long     nall;
  
  struct io_gadget_header
 {
  int      np[6];
  double   massarr[6];
  double   expansion;
  double   redshift;
  int      flagsfr;
  int      flagfeedback;
  int      nall[6];
  int      flagcooling;
  int      NumFiles;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     unused[SIZEOFGADGETHEADER- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
 } header; 
  
} gadget;

struct particle_data 
{
  float     Pos[3];       /* particle position   */  
  float     Vel[3];       /* particle velocity   */  
  float     Mass;         /* particle mass       */
  float     u;            /* gas internal energy */
  long      ID;
} *Part;

long *ID;

/*=============================================================================
 *                                   PROTOYPES
 *=============================================================================*/
void   get_haloes(char *);
double get_Rstar(long, long, char *);
void   read_gadget(FILE *);
int ReadFloat          (FILE *fptr,float *n, int swap);
int ReadDouble         (FILE *fptr,double *n,int swap);
int ReadLongLong       (FILE *fptr,long long *n,int swap);
int ReadLong           (FILE *fptr,long *n,int swap);
int ReadUInt           (FILE *fptr,unsigned int *n,int swap);
int ReadInt            (FILE *fptr,int *n,int swap);
int ReadChars          (FILE *fptr,char *s,int n);
int ReadString         (FILE *fptr,char *s,int n);
long get_pid(int i);

/*=============================================================================
 *                                   MAIN
 *=============================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
  long   ihalo, jhalo, fhalo, TOTnhaloes, npart, ipart, pid, pstar;
  int    ptype;
  double pmass;
  double MinStar, MinR, Mstar, Rstar;
  FILE  *fpin, *fpout;
  
  char haloIDfile[MAXSTRING];
  char particlesfile[MAXSTRING];
  char halosfile[MAXSTRING];
  char snapshotfile[MAXSTRING];
  char outfile[MAXSTRING];
  
  if(argc<6)
   {
    fprintf(stderr,"usage: %s particlesfile halosfile snapshotfile haloIDfile outfile\n",*argv);
    exit(1);
   }
  strcpy(particlesfile,argv[1]);
  strcpy(halosfile,argv[2]);
  strcpy(snapshotfile,argv[3]);
  strcpy(haloIDfile,argv[4]);
  strcpy(outfile,argv[5]);
  
  /*------------------
   * read halo list
   *------------------*/
  get_haloes(haloIDfile);
  
  
  /*------------------
   * read GADGET file
   *------------------*/
  fpin = fopen(snapshotfile,"rb");
  gadget.no_gadget_files  = 1;
  gadget.i_gadget_file    = 0;
  
  /* allocate temporary storage for no. of particles arrays */
  gadget.np[0]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
  gadget.np[1]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
  gadget.np[2]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
  gadget.np[3]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
  gadget.np[4]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
  gadget.np[5]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
  
  read_gadget(fpin);
  fclose(fpin);
  
  /* remove temporary storage again */
  free(gadget.np[0]);
  free(gadget.np[1]);
  free(gadget.np[2]);
  free(gadget.np[3]);
  free(gadget.np[4]);
  free(gadget.np[5]);

  
  /*---------------------
   * read particles file
   *---------------------*/
  fpout = fopen(outfile,"w");
  fpin  = fopen(particlesfile,"r");
  fscanf(fpin,"%ld",&TOTnhaloes);
  fhalo = 0;

  for(jhalo=0; jhalo<nhaloes; jhalo++)
   {
    ihalo = halo[jhalo];
    
    fprintf(stderr,"scanning for halo %ld ... ",ihalo);
    
    /* skip until we reach ihalo */
    while(fhalo<ihalo)
     {
      fhalo++;
      fscanf(fpin,"%ld",&npart);
      for(ipart=0; ipart<npart; ipart++)
        fscanf(fpin,"%ld %d %lf",&pid,&ptype,&pmass);
     }
    
    /* this is the halo */
    fhalo++;
    MinStar = 0.;
    MinR    = 0.;
    Mstar   = 0.;
    fscanf(fpin,"%ld",&npart);
    for(ipart=0; ipart<npart; ipart++)
     {
      fscanf(fpin,"%ld %d %lf",&pid,&ptype,&pmass);
      MinR += pmass;
      
      if(ptype==4)
       {
        Mstar  += pmass;
        MinStar = MinR;
        pstar   = pid;
       }
     }
    
    if(Mstar > 0)
      Rstar = get_Rstar(ihalo, pstar, halosfile);
    else
      Rstar = -1.0;
    
    fprintf(fpout, "%ld %g %g %g\n",ihalo,MinStar,Rstar,Mstar);
    fprintf(stderr,"%ld %g %g %g\n",ihalo,MinStar,Rstar,Mstar);
   }
  fclose(fpin);
  fclose(fpout);
}

/*=============================================================================
 *                                   get_haloes
 *=============================================================================*/
void get_haloes(char *haloIDfile)
{
  FILE *fpin;
  long  ihalo, idummy;
  
  fpin = fopen(haloIDfile,"r");
  fscanf(fpin,"%ld",&nhaloes);
  fprintf(stderr,"will look for %ld haloes: ",nhaloes);
  
  halo = (long *) calloc(nhaloes, sizeof(long));
  
  for(ihalo=0; ihalo<nhaloes; ihalo++)
   {
    fscanf(fpin,"%ld %ld",&(halo[ihalo]),&idummy);
    //fprintf(stderr," %ld",halo[ihalo]);
   }
  fprintf(stderr,"\n");
  fclose(fpin);
}

/*=============================================================================
 *                                   get_Rstar
 *=============================================================================*/
double get_Rstar(long ihalo, long pstar, char *halosfile)
{
  FILE *fpin;
  char dummyline[MAXSTRING];
  long fhalo, ipart;
  double dummy, Xc, Yc, Zc, Xp, Yp, Zp, Rvir, Rstar;
  
  /*--------------------------------------
   * find halo centre in *.AHF_halos file
   *--------------------------------------*/
  fpin = fopen(halosfile,"r");
  fgets(dummyline,MAXSTRING,fpin);
  for(fhalo=0; fhalo<=ihalo; fhalo++)
   {
    fgets(dummyline,MAXSTRING,fpin);
    //sscanf(dummyline,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &dummy,&dummy,&Xc,&Yc,&Zc,&dummy,&dummy,&dummy,&dummy,&Rvir);
    sscanf(dummyline,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &dummy,&dummy,&dummy,&dummy,&dummy,&Xc,&Yc,&Zc,&dummy,&dummy,&dummy,&Rvir);
   }
  fclose(fpin);

  // convert to Mpc/h
  Xc   /= 1000.;
  Yc   /= 1000.;
  Zc   /= 1000.;
  Rvir /= 1000.;
  fprintf(stderr,"\nhalo position          = %g %g %g\n",Xc,Yc,Zc);
    
  /*----------------------------------
   * determine star particle position
   *----------------------------------*/
  for(ipart=0; ipart<gadget.nall; ipart++)
   {
    if(Part[ipart].ID == pstar)
      break;
   }
  fprintf(stderr,"farthest star particle = %g %g %g (%ld)\n",Part[ipart].Pos[X],Part[ipart].Pos[Y],Part[ipart].Pos[Z],pstar);
  Rstar = sqrt(pow2(Part[ipart].Pos[X]-Xc)+pow2(Part[ipart].Pos[Y]-Yc)+pow2(Part[ipart].Pos[Z]-Zc));
  fprintf(stderr,"                         %g %g\n",Rvir,Rstar*1000.);
  return Rstar;
}


/*=============================================================================
 *                                READ_GADGET
 *=============================================================================*/
void read_gadget(FILE *icfile)
{
  long unsigned  ipart;
  
  double         tot_mass[6];
  
  int            i,j,k;
  int            no_part;
  int            massflag;
  char           DATA[MAXSTRING];
  float          dummy[3];
  double         x_fac, v_fac, m_fac;
  long           pid, ldummy;
  int            idummy;
  
  /*================= read in GADGET IO header =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading... %s\n",DATA);
   }
  
  GADGET_SKIP;
  
  ReadInt(icfile,&(gadget.header.np[0]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.np[1]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.np[2]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.np[3]),SWAPBYTES);    /* number of particles in current file */
  ReadInt(icfile,&(gadget.header.np[4]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.np[5]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[0]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[1]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[2]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[3]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[4]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[5]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.expansion),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.redshift),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.flagsfr),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.flagfeedback),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[0]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[1]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[2]),SWAPBYTES);  /* total number of particles in simulation */
  ReadInt(icfile,&(gadget.header.nall[3]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[4]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[5]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.flagcooling),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.NumFiles),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.BoxSize),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.Omega0),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.OmegaLambda),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.HubbleParam),SWAPBYTES);
  ReadChars(icfile,&(gadget.header.unused[0]),SIZEOFGADGETHEADER- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8);
  
  GADGET_SKIP;
  /*================= read in GADGET IO header =================*/
  
  
  /* keep track of no. of particles in each GADGET file */
  gadget.np[0][gadget.i_gadget_file] = gadget.header.np[0];
  gadget.np[1][gadget.i_gadget_file] = gadget.header.np[1];
  gadget.np[2][gadget.i_gadget_file] = gadget.header.np[2];
  gadget.np[3][gadget.i_gadget_file] = gadget.header.np[3];
  gadget.np[4][gadget.i_gadget_file] = gadget.header.np[4];
  gadget.np[5][gadget.i_gadget_file] = gadget.header.np[5];
  
  /* conversion factors to Mpc/h, km/sec, Msun/h */
  x_fac  = GADGET_LUNIT;
  v_fac  = sqrt(gadget.header.expansion);
  m_fac  = GADGET_MUNIT;
  
  /* count total no. of particles in current file (and set massflag) */
  massflag    = 0;
  no_part     = 0;
  gadget.nall = 0;
  for(i=0;i<6;i++) 
   {
    no_part     += gadget.header.np[i];
    gadget.nall += gadget.header.nall[i];
    if(gadget.header.massarr[i] < MZERO && gadget.header.np[i] > 0)
      massflag=1;  
   }  
  
  /* be verbose */
  fprintf(stderr,"expansion factor: %lf\n",             gadget.header.expansion);
  fprintf(stderr,"redshift:         %lf\n",             gadget.header.redshift);
  fprintf(stderr,"boxsize:          %lf (%lf Mpc/h)\n", gadget.header.BoxSize,gadget.header.BoxSize*GADGET_LUNIT);
  fprintf(stderr,"omega0:           %lf\n",             gadget.header.Omega0);
  fprintf(stderr,"lambda0:          %lf\n",             gadget.header.OmegaLambda);
  fprintf(stderr,"HubbleParam:      %lf\n\n",           gadget.header.HubbleParam);
  
  fprintf(stderr,"gas:    np[0]=%9d\t nall[0]=%9d\t massarr[0]=%g\n",gadget.header.np[0],gadget.header.nall[0],gadget.header.massarr[0]); 
  fprintf(stderr,"halo:   np[1]=%9d\t nall[1]=%9d\t massarr[1]=%g\n",gadget.header.np[1],gadget.header.nall[1],gadget.header.massarr[1]); 
  fprintf(stderr,"disk:   np[2]=%9d\t nall[2]=%9d\t massarr[2]=%g\n",gadget.header.np[2],gadget.header.nall[2],gadget.header.massarr[2]); 
  fprintf(stderr,"bulge:  np[3]=%9d\t nall[3]=%9d\t massarr[3]=%g\n",gadget.header.np[3],gadget.header.nall[3],gadget.header.massarr[3]); 
  fprintf(stderr,"stars:  np[4]=%9d\t nall[4]=%9d\t massarr[4]=%g\n",gadget.header.np[4],gadget.header.nall[4],gadget.header.massarr[4]); 
  fprintf(stderr,"bndry:  np[5]=%9d\t nall[5]=%9d\t massarr[5]=%g\n",gadget.header.np[5],gadget.header.nall[5],gadget.header.massarr[5]); 
  
  fprintf(stderr,"\n-> reading %d particles from  GADGET file #%d/%d...\n\n", no_part, gadget.i_gadget_file+1, gadget.no_gadget_files);
  
  /* allocate particle array (only once when reading the first file, of course!) */
  if(gadget.i_gadget_file == 0)
   {
    fprintf(stderr,"-> allocating %f GB of RAM for particles\n\n",(float)(gadget.nall*sizeof(struct particle_data))/1024./1024./1024.);
    if(!(Part=(struct particle_data *) calloc(gadget.nall, sizeof(struct particle_data))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
   }
  
  /*================= read in GADGET particles =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
   }
  else
   {
    fprintf(stderr,"reading ");
   }
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
  
  for(i=0;i<no_part;i++)
   {    
     /* read */
     ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
     ReadFloat(icfile,&(dummy[1]),SWAPBYTES);
     ReadFloat(icfile,&(dummy[2]),SWAPBYTES);
     
     /* get proper position in Part[] array */
     pid = get_pid(i);
     
     /* storage and conversion to comoving physical units */
     Part[pid].Pos[0] = dummy[0] * x_fac;
     Part[pid].Pos[1] = dummy[1] * x_fac;
     Part[pid].Pos[2] = dummy[2] * x_fac;      
   }
  fprintf(stderr,"Pos[X]=%12.6g Pos[Y]=%12.6g Pos[Z]=%12.6g ... ",Part[no_part-1].Pos[X],Part[no_part-1].Pos[Y],Part[no_part-1].Pos[Z]);
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  /*================= read in GADGET particles =================*/
  
  
  
  /*================= read in GADGET velocities =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;  
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
   }
  else
   {
    fprintf(stderr,"reading ");
   }
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
  
  for(i=0;i<no_part;i++)
   {
    /* read */
    ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
    ReadFloat(icfile,&(dummy[1]),SWAPBYTES);
    ReadFloat(icfile,&(dummy[2]),SWAPBYTES);
    
    /* get proper position in Part[] array */
    pid = get_pid(i);
    
    /* storage and conversion to comoving physical units */
    Part[pid].Vel[0] = dummy[0] * v_fac;
    Part[pid].Vel[1] = dummy[1] * v_fac;
    Part[pid].Vel[2] = dummy[2] * v_fac; 
   }
  fprintf(stderr,"Vel[X]=%12.6g Vel[Y]=%12.6g Vel[Z]=%12.6g ... ",Part[no_part-1].Vel[X],Part[no_part-1].Vel[Y],Part[no_part-1].Vel[Z]);
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  /*================= read in GADGET velocities =================*/
  
  
  /*================= read in GADGET id's =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
   }
  else
   {
    fprintf(stderr,"reading ");
   }
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
  
  for(i=0;i<no_part;i++)
   {
    /* get proper position in Part[] array */
    pid = get_pid(i);
    
    if(LGADGET)
     {
      ReadLong(icfile,&ldummy,SWAPBYTES);
      Part[pid].ID = ldummy;
     }
    else
     {
      ReadInt(icfile,&idummy,SWAPBYTES);
      Part[pid].ID = (long) idummy;
     }
    
    /* check the ID range of the "halo" particles */
    if(gadget.header.np[0] <= i && i < gadget.header.np[0]+gadget.header.np[1])
     {
      if(Part[pid].ID > IDmax) IDmax = Part[pid].ID; 
      if(Part[pid].ID < IDmin) IDmin = Part[pid].ID; 
     }
   }
  
  fprintf(stderr,"ID=%12ld ...  ",Part[no_part-1].ID);
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  /*================= read in GADGET id's =================*/
  
  
  k = 0;
  /* massflag == 1 indicates that massarr[i] = 0 and hence need to read in particle masses */
  if(massflag==1) 
   {
    /*================= read in GADGET individual particle masses =================*/
    if(FORMAT == 2)
     {
      GADGET_SKIP; 
      fread(DATA,sizeof(char),4,icfile);
      DATA[4] = '\0';
      GADGET_SKIP;
      
      GADGET_SKIP;
      fprintf(stderr,"reading %s",DATA);
     }
    else
     {
      fprintf(stderr,"reading ");
     }
    
    GADGET_SKIP;
    fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
    
    for(i=0;i<6;i++)
     {
      tot_mass[i] = 0.;
      if (gadget.header.np[i] > 0 && gadget.header.massarr[i] < MZERO  ) 
       {
        
        fprintf(stderr,"  %d    ",i);
        
        for(j=0; j<gadget.header.np[i]; j++)
         {
          /* read */
          ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
          
          /* get proper position in Part[] array */
          pid = get_pid(k);
          
          /* store */
          Part[pid].Mass  = dummy[0];
          tot_mass[i]    += dummy[0];
          k++;
         }
       }
      else
       {
        /* simply copy appropriate massarr[i] to particles */
        for(j=0; j<gadget.header.np[i]; j++) 
         {
          /* get proper position in Part[] array */
          pid = get_pid(k);
          
          /* store */
          Part[pid].Mass = gadget.header.massarr[i];
          k++;
         }
        tot_mass[i] = gadget.header.np[i]*gadget.header.massarr[i];
       }
     }
    
    GADGET_SKIP;
    fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
    /*================= read in GADGET individual particle masses =================*/
   } 
  
  /* simply copy appropriate massarr[i] to particles */
  else 
   {
    k=0;
    for(i=0;i<6;i++)
     {
      for(j=0;j<gadget.header.np[i];j++) 
       {
        /* get proper position in Part[] array */
        pid = get_pid(k);
        
        /* store */
        Part[pid].Mass = gadget.header.massarr[i];
        k++;
       }
      tot_mass[i] = gadget.header.np[i]*gadget.header.massarr[i];
     }
   }
  
  /* convert masses to Msun/h */
  k=0;
  for(i=0;i<6;i++)
   {
    for(j=0;j<gadget.header.np[i];j++) 
     {
      /* get proper position in Part[] array */
      pid = get_pid(k);
      
      Part[pid].Mass *= GADGET_MUNIT;
      k++;
     }
   }
  
  
  /*================= read in GADGET gas particle energies =================*/
  if(gadget.header.np[0] > 0) 
   {      
     if(FORMAT == 2)
      {
       GADGET_SKIP;
       
       fread(DATA,sizeof(char),4,icfile);
       DATA[4] = '\0';
       GADGET_SKIP;
       
       GADGET_SKIP;
       fprintf(stderr,"reading %s",DATA);
      }
     else
      {
       fprintf(stderr,"reading ");
      }
     
     GADGET_SKIP; 
     fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
     
     for(i=0; i<gadget.header.np[0]; i++)
      {
       /* store */
       ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
       
       /* get proper position in Part[] array */
       pid = get_pid(i);
       
       /* store additional gas particle property */
       Part[pid].u = dummy[0];         
      }
     
     GADGET_SKIP;
     fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
   } 
  /*================= read in GADGET gas particle energies =================*/
  
  
  /* be verbose */
  fprintf(stderr,"\n");
  if(gadget.header.np[0] > 0) fprintf(stderr,"    gas:    tot_mass[0]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[0]*GADGET_MUNIT,tot_mass[0]/(double)gadget.header.np[0]*GADGET_MUNIT);
  if(gadget.header.np[1] > 0) fprintf(stderr,"    halo:   tot_mass[1]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[1]*GADGET_MUNIT,tot_mass[1]/(double)gadget.header.np[1]*GADGET_MUNIT);
  if(gadget.header.np[2] > 0) fprintf(stderr,"    disk:   tot_mass[2]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[2]*GADGET_MUNIT,tot_mass[2]/(double)gadget.header.np[2]*GADGET_MUNIT);
  if(gadget.header.np[3] > 0) fprintf(stderr,"    bulge:  tot_mass[3]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[3]*GADGET_MUNIT,tot_mass[3]/(double)gadget.header.np[3]*GADGET_MUNIT);
  if(gadget.header.np[4] > 0) fprintf(stderr,"    stars:  tot_mass[4]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[4]*GADGET_MUNIT,tot_mass[4]/(double)gadget.header.np[4]*GADGET_MUNIT);
  if(gadget.header.np[5] > 0) fprintf(stderr,"    bndry:  tot_mass[5]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[5]*GADGET_MUNIT,tot_mass[5]/(double)gadget.header.np[5]*GADGET_MUNIT);
              
  fprintf(stderr,"===================================================================\n");
}

/*
 Read a string of n characters
 */
int ReadString(FILE *fptr,char *s,int n)
{
  int i,c;
  
  if(sizeof(char) != 1)
   {
    fprintf(stderr,"ReadString: sizeof(char)=%ld and not 1\n",sizeof(char));
    exit(0);
   }
  
  s[0] = '\0';
  for (i=0;i<n;i++) {
    c = fgetc(fptr);
    if (c == EOF)
      return(FALSE);
    s[i] = c;
    s[i+1] = '\0';
  }
  return(TRUE);
}

/*
 Read an array of n characters
 NOTE: the difference to ReadString() is that we do not '\0'-terminate the array
 */
int ReadChars(FILE *fptr,char *s,int n)
{
  int i,c;
  
  if(sizeof(char) != 1)
   {
    fprintf(stderr,"ReadChars: sizeof(char)=%ld and not 1\n",sizeof(char));
    exit(0);
   }
  
  s[0] = '\0';
  for (i=0;i<n;i++) {
    c = fgetc(fptr);
    if (c == EOF)
      return(FALSE);
    s[i] = c;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped integer
 */
int ReadInt(FILE *fptr,int *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
   {
    fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
   }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped unsigned integer
 */
int ReadUInt(FILE *fptr,unsigned int *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
   {
    fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
   }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped long integer
 */
int ReadLong(FILE *fptr,long *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(long) == 4)
   {
    if (fread(n,4,1,fptr) != 1)
      return(FALSE);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
    }
   }
  else if(sizeof(long) == 8)
   {
    if (fread(n,8,1,fptr) != 1)
      return(FALSE);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[7];
      cptr[7] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[6];
      cptr[6] = tmp;
      tmp     = cptr[2];
      cptr[2] = cptr[5];
      cptr[5] = tmp;
      tmp     = cptr[3];
      cptr[3] = cptr[4];
      cptr[4] = tmp;
    }
   }
  else
   {
    fprintf(stderr,"ReadLong: something wrong...cannot read long\n");
    exit(0);
   }
  
  
  
  return(TRUE);
}

/*
 Read a possibly byte swapped long long integer
 */
int ReadLongLong(FILE *fptr,long long *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if (fread(n,8,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[7];
    cptr[7] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[6];
    cptr[6] = tmp;
    tmp     = cptr[2];
    cptr[2] = cptr[5];
    cptr[5] = tmp;
    tmp     = cptr[3];
    cptr[3] = cptr[4];
    cptr[4] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped double precision number
 Assume IEEE
 */
int ReadDouble(FILE *fptr,double *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(double) != 8)
   {
    fprintf(stderr,"ReadDouble: sizeof(double)=%ld and not 8\n",sizeof(double));
    exit(0);
   }
  
  if (fread(n,8,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[7];
    cptr[7] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[6];
    cptr[6] = tmp;
    tmp     = cptr[2];
    cptr[2] = cptr[5];
    cptr[5] = tmp;
    tmp     = cptr[3];
    cptr[3] = cptr[4];
    cptr[4] = tmp;
  }
  
  return(TRUE);
}

/*
 Read a possibly byte swapped floating point number
 Assume IEEE format
 */
int ReadFloat(FILE *fptr,float *n, int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(float) != 4)
   {
    fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float));
    exit(0);
   }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) 
   {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
   }
  return(TRUE);
}

/*=============================================================================
 *                        get proper position in Part[] array
 *=============================================================================*/
long get_pid(int i)
{
  long pid;
  long itype, ifile;
  
  pid = 0;
  for(ifile=0; ifile<gadget.i_gadget_file; ifile++)
    for(itype=0; itype<6; itype++)
      pid += gadget.np[itype][ifile];
  pid += i;
  
  return(pid);
}


