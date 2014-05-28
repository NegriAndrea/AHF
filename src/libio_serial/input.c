#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "io_serial.h"
#include "../libutility/utility.h"

#define Xc   11.4264
#define Yc    8.5032
#define Zc    8.3975
#define Rvir  0.01

#ifdef MLAPM
/*=========================================================================
 * definitions for read_mlapm()
 *=========================================================================*/
#define MLAPM_HEADERSTRING  100      /* this value must agree with the value in param.h */
//#define MLAPM_FILLHEADER (508-MLAPM_HEADERSTRING-4-5*8-4-6*8-20*8)
#define MLAPM_FILLHEADER (508 - MLAPM_HEADERSTRING*sizeof(char) - 1*sizeof(long) - 5*sizeof(double) - 1*sizeof(int) - 6*sizeof(double) - 10*sizeof(float))

struct mlapm_header
{
   char          header[MLAPM_HEADERSTRING];
   long unsigned no_part;
   double        boxsize;
   double        omega0;
   double        lambda0;
   double        a_initial;
   double        a_current;
   int           no_timestep;
   double        K_initial;
   double        U_initial;
   double        K_current;
   double        U_current;
   double        Eintegral;
   double        Econst;
   float         numerics[10];
   char          dummy[MLAPM_FILLHEADER];
} mlapm_header;
#endif /* MLAPM */

struct ascii_header
{
   char          header[MAXSTRING];
   long unsigned no_part;
   double        omega0;
   double        lambda0;
   double        omegab;
   double        gamma;
   double        H_frac;
   double        T_init;
   double        B_init;
   double        boxsize;
   double        a_initial;
   double        a_current;
   double        z_initial;
   double        z_current;
   int           no_timestep;
   int           multi_mass;
   int           double_precision;
} ascii_header;

#ifdef MARE_NOSTRUM
#define MN_SKIP  ReadFloat(icfile,&fdummy,1);
struct MN_header
{
   float  boxsize;
   float  aexpn;
   float  Om0;
   float  Oml0;
   float  xm_DM;
   float  xm_gas;
   float  boxsize_cl;
   float  xc1;
   float  yc1;
   float  zc1;
   int    lkl1;
   float  xm1;
   float  vdisp1;
   float  xc;
   float  yc;
   float  zc;
   int    i_sub;
   float  xcnew;
   float  ycnew;
   float  zcnew;
   float  r_vir;
   double xm_vir;
   int    i_DM_vir;
   int    i_gas_vir;
   float  bfrac_vir;
   float  over_dens_vir;
   
   int    n_cl_DM;
   int    n_cl_gas;
   
} MN_header;
#endif /* MARE_NOSTRUM */

#ifdef GADGET
/*=========================================================================
* definitions for read_gadget()
*=========================================================================*/
#ifdef LGADGET
long long    blklen;
#define GADGET_SKIP  ReadLongLong (icfile,&blklen,SWAPBYTES);
//#define GADGET_SKIP  fread(&blklen,sizeof(long long),1,icfile)
#else
unsigned int blklen;
#define GADGET_SKIP  ReadUInt     (icfile,&blklen,SWAPBYTES);
//#define GADGET_SKIP  fread(&blklen,sizeof(int),1,icfile)
#endif

struct particle_data 
{
   float     Pos[3];       /* particle position   */  
   float     Vel[3];       /* particle velocity   */  
   float     Mass;         /* particle mass       */
   float     u;            /* gas internal energy */
#ifdef LGADGET
   long long ID;           /* unique particle identifier */  
#else
   int       ID;
#endif
} *P_gadget;
#endif

#ifdef ART
/*=========================================================================*/
/* definitions for read_ART()*/
/* for binary fortran files one needs to skip one float before actually reading */
int dustbin;
#define ART_SKIP(fp) ReadInt(fp,&dustbin,1);
//long dustbin;
//#define ART_SKIP(fp) ReadLong(fp,&dustbin,1);
struct CONTROL
{
   char    header_string[45];
   float   aexpn;
   float   aexp0;
   float   amplt;
   float   astep;
   int     istep;
   float   partw;
   float   tintg;
   float   ekin;
   float   ekin1;
   float   ekin2;
   float   au0;
   float   aeu0;
   int     nrowc;
   int     ngridc;
   int     nspecies;
   int     nseed;
   float   Om0;
   float   Oml0;
   float   hubble;
   float   wp5;
   float   Ocurv;
   float   extras[100];
} ART_header;
#endif

#ifdef DEVA
#define itstar   -1
#define itdark    0
#define itgas     1
#define itgas2    2
#define nmetals   16

#ifdef DEVA_DOUBLE
#define DEVAFLOUBLE 8
typedef double devaflouble;
void ReadFlouble(FILE *icfile, double *ddummy, int SWAPBYTES)
{
  ReadDouble(icfile, ddummy, SWAPBYTES);
}
#else /* DEVA_DOUBLE */
#define DEVAFLOUBLE 4
typedef float  devaflouble;
void ReadFlouble(FILE *icfile, float *fdummy, int SWAPBYTES)
{
  ReadFloat(icfile, fdummy, SWAPBYTES);
}
#endif /* DEVA_DOUBLE */

#ifdef DEVA_DIRECT_ACCESS
#define DEVA_SKIP {;} 
#else /* DEVA_DIRECT_ACCESS */
unsigned int deva_dummy;
#define DEVA_SKIP ReadUInt(icfile,&deva_dummy,SWAPBYTES); 
#endif /* DEVA_DIRECT_ACCESS */

#ifdef BYTESWAP
int SWAPBYTES = TRUE;
#else
int SWAPBYTES = FALSE;
#endif

#ifdef DEVA2
struct ibuf_1
{
  int   itime;
  int   itstop;
  int   itdump;
  int   iout;
  int   nsformed;
  int   nsdead;
  int   irun;
  int   nobj;
  int   ngas;
  int   ndark;
  int   L;
  float CHEMEVOL;
  float ANOTHER;
  float COOL;
  float REFINEMENT;
  float HYDRO;
  float GRAVITY;
  int   ISOLATED;
  float EXPAND;
  float COMOVING;
  float STARFORM;
  float GRADH;
  int   INITIALCOND;
  int   nstar;
  int   iseed1;
  int   ispec;
  int   indxsp;
  int   n_neigh;
  int   lastbar;
  
  int   fill[100-29];
} ibuf1;

struct ibuf_2
{
  devaflouble time;
  devaflouble atime;
  devaflouble htime;
  devaflouble dtime;
  devaflouble E_init;
  devaflouble E_kin;
  devaflouble E_ther;
  devaflouble E_pot;
  devaflouble Radiation;
  devaflouble Esum;
  devaflouble Rsum;
  devaflouble cpu;
  devaflouble time_end;
  devaflouble tout;
  devaflouble padding;
  devaflouble Tlost;
  devaflouble Qlost;
  devaflouble Ulost;
  devaflouble delta_min;
  devaflouble delta_max;
  devaflouble T_min;
  devaflouble avisc;
  devaflouble bvisc;
  devaflouble eta2;
  devaflouble rho_star;
  devaflouble c_star;
  devaflouble rmtot;
  devaflouble rmsep;
  devaflouble dnthres;
  devaflouble sft0;
  devaflouble sftmin;
  devaflouble sftmax;
  devaflouble h100;
  devaflouble box100;
  devaflouble rmgas;
  devaflouble rmdark;
  devaflouble omega0;
  devaflouble xlambda0;
  devaflouble h0t0;
  devaflouble omegab0;
  devaflouble sigma80;
  devaflouble ztime0;
  devaflouble e0;
  
  devaflouble fill[100-43];
} ibuf2;

#else /* DEVA2 */
struct ibuf_1
{
  int   itime;
  int   itstop;
  int   itdump;
  int   iout;
  float time;
  float atime;
  float htime;
  float dtime;
  float E_init;
  float E_kin;
  float E_ther;
  float E_pot;
  float Radiation;
  float Esum;
  float Rsum;
  float cpu;
  float time_end;
  float tout;
  int   icdump;     //ifree1
  float padding;
  float Tlost;
  float Qlost;
  float Ulost;
  float delta_min;
  float delta_max;
  float T_min;
  float avisc;
  float bvisc;
  float eta2;
  float rho_star;
  float c_star;
  float pid;        //xnul
  float rmtot;
  int   nsformed;
  int   nsdead;
  float rmsep;

  int   fill[100-36];
} ibuf1;

struct ibuf_2
{
  int   irun;
  int   nobj;
  int   ngas;
  int   ndark;
  int   L;
  int   intl;    //CHEMEVOL
  int   nlmx;    //ANOTHER
  float perr;    //xfree1
  float dtnorm;  //xfree2
  float sft0;
  float sftmin;
  float sftmax;
  float h100;
  float box100;
  float zmet0;   //xnul2
  float spc0;    //rmgas
  float rmbary;  //rmdark
  float rmgas;   //xfree3
  float rmdark;  //xfree4
  float rmnorm;  //xfree5
  float tstart;  //xfree6
  float omega0;
  float xlambda0;
  float h0t0;
  int   mcool;   //COOL
  int   mref;    //REFINEMENT
  int   mhyd;    //HYDRO
  int   mgrav;   //GRAVITY
  int   miso;    //ISOLATED
  int   mexp;    //EXPAND
  int   mcom;    //COMOVING
  int   mstf;    //STARFORM
  int   mgrh;    //GRADH
  int   nstar;
  int   iseed1;
  float omegab0;
  int   ispec;
  int   indxsp;
  float sigma80;
  float ztime0;
  float e0;
  int   ifree2;
  int   n_neigh;
  int   lastbar;
  int   nidsv;

  int   fill[100-45];
} ibuf2;

struct ibuf_0
{
  float time_out[100];
} ibuf;
#endif /* DEVA2 */

#endif /* DEVA */

//========================================================================
//
//                              AMIGA
//
//========================================================================
void read_amiga(FILE *icfile)
{
  partptr        cur_part;
  long unsigned  ipart;
  double         dpos, dmom, dweight;
  float          fpos, fmom, fweight;
  int            i, j, k;
  int            SWAPBYTES;
  
  fprintf(stderr,"\n===================================================================\n");
#ifdef BYTESWAP
  fprintf(stderr,"           start reading BYTESWAPed AMIGA particles\n\n");
#else
  fprintf(stderr,"                 start reading AMIGA particles\n\n");
#endif

  /* read AMIGA header */
  read_amiga_header(icfile, &io, &SWAPBYTES);
  
  /* allocate memory for particles */
  io.no_part  = io.header.no_part;
  io.fst_part = c_part(io.no_part);
  io.no_gas   = 0;
  io.fst_gas  = NULL;
  io.no_stars = 0;
  io.fst_star = NULL;
  
  /* loop over all particles */
  for(cur_part = io.fst_part; cur_part < io.fst_part+io.no_part; cur_part++)
    {
      /* read positions and velocities */
      for(j = X; j <= Z; j++)
        {
#ifdef DOUBLE
          if(io.header.double_precision == 0)
            {
              ReadFloat(icfile,&fpos,SWAPBYTES);
              ReadFloat(icfile,&fmom,SWAPBYTES);
              dpos = (double)fpos;
              dmom = (double)fmom;
            }
          else
            {
              ReadDouble(icfile,&dpos,SWAPBYTES);
              ReadDouble(icfile,&dmom,SWAPBYTES);
            }
#else
          if(io.header.double_precision == 1)
            {
              ReadDouble(icfile,&dpos,SWAPBYTES);
              ReadDouble(icfile,&dmom,SWAPBYTES);
            }
          else
            {
              ReadFloat(icfile,&fpos,SWAPBYTES);
              ReadFloat(icfile,&fmom,SWAPBYTES);
              dpos = (double)fpos;
              dmom = (double)fmom;
            }
#endif
          /* brute force cast to whatever pos[] and mom[] are... */
          cur_part->pos[j] = f1mod(dpos + 1.0, 1.0);
          cur_part->mom[j] = dmom;
        }
      
      /* read masses (if required) */
#ifdef MULTIMASS
#ifdef DOUBLE
      if(io.header.double_precision == 0)
        {
          ReadFloat(icfile,&fweight,SWAPBYTES);
          dweight = (double)fweight;
        }
      else
        {
          ReadDouble(icfile,&dweight,SWAPBYTES);
        }
#else // DOUBLE
      if(io.header.double_precision == 1)
        {
          ReadDouble(icfile,&dweight,SWAPBYTES);
        }
      else
        {
          ReadFloat(icfile,&fweight,SWAPBYTES);
          dweight = (double)fweight;
        }
#endif // DOUBLE
      
      cur_part->weight  = dweight; 
#endif  // MULTIMASS
      
            
#ifdef SWAP_XY
      {
        double tmp_x, tmp_vx;
        
        tmp_x  = cur_part->pos[X];
        tmp_vx = cur_part->mom[X];
        
        cur_part->pos[X] = cur_part->pos[Y];
        cur_part->mom[X] = cur_part->mom[Y];
        
        cur_part->pos[Y] = tmp_x;
        cur_part->mom[Y] = tmp_vx;
      }
#endif
      
#ifdef SWAP_XZ
      {
        double tmp_x, tmp_vx;
        
        tmp_x  = cur_part->pos[X];
        tmp_vx = cur_part->mom[X];
        
        cur_part->pos[X] = cur_part->pos[Z];
        cur_part->mom[X] = cur_part->mom[Z];
        
        cur_part->pos[Z] = tmp_x;
        cur_part->mom[Z] = tmp_vx;
      }
#endif
    }
  

  
  
  
  /*=====================================================================================
   * perform some simple tests whether or not the data is compliant with the DEFINEFLAGS
   *             (this check has to be done here and >>not<< earlier!!!)
   *=====================================================================================*/
  sanity_check();
  
  
  
  
  
  
  /*--------------------------------------------------------------------
   *                           misc stuff
   *--------------------------------------------------------------------*/
  
#ifdef GRIDPARTICLES
  {
    int    n1dim;
    int    int_dummy;
    double shift;
    
    n1dim     = (long unsigned) pow((double)io.no_part,0.333333333333333333);
    int_dummy = -123456;
    ipart     = 0;
    
    /* do not place the particles exactly onto the cells or cell boundaries */
    shift = ran3(&int_dummy)/(double)n1dim;
    
    for(i=0; i<n1dim; i++)
      for(j=0; j<n1dim; j++)
        for(k=0; k<n1dim; k++)
          {
            cur_part = io.fst_part + ipart;
            
            cur_part->pos[X] = f1mod((double)i/(double)n1dim + shift + 1.0, 1.0);
            cur_part->pos[Y] = f1mod((double)j/(double)n1dim + shift + 1.0, 1.0);
            cur_part->pos[Z] = f1mod((double)k/(double)n1dim + shift + 1.0, 1.0);
            
            ipart++;
          }
  }
#endif
  
#ifdef SHIFTPARTICLES
  {
    double shift;
    
    shift = 0.5*1/64.; /* one-cell */
    
    for(ipart=0; ipart<io.no_part; ipart++)
      {
        cur_part = io.fst_part + ipart;
        
        cur_part->pos[X] = f1mod(cur_part->pos[X] + shift + 1.0, 1.0);
        cur_part->pos[Y] = f1mod(cur_part->pos[Y] + shift + 1.0, 1.0);
        cur_part->pos[Z] = f1mod(cur_part->pos[Z] + shift + 1.0, 1.0);      
      }
  }
#endif
  
  fprintf(stderr,"\n                  finished reading AMIGA particles\n");
  fprintf(stderr,"====================================================================\n\n");
}

#ifdef MLAPM
//========================================================================
//
//                              MLAPM
//
//========================================================================
/*========================================================================
 * read data from MLAPM initial conditions file
 *
 * NOTE: this routine is tailored to work with the STH2 data only!
 *       the STH2 simulations were run with mlapm-v4.6 and for
 *       later versions the header format changed!
 *========================================================================*/
void read_mlapm(FILE *icfile)
{
   partptr        cur_part;
   int            j, ino_part;
   long           lno_part;
   long unsigned  ipart, ndummy;
   flouble        pos, mom;
   int            i, idummy;
   double         dummy;
   int            SWAPBYTES;
   double         x_fac, v_fac, m_fac;
   flouble        weight;
   
#ifdef BYTESWAP
   SWAPBYTES = TRUE;
#else
   SWAPBYTES = FALSE;
#endif
   
   fprintf(stderr,"\n===================================================================\n");
#ifdef BYTESWAP
   fprintf(stderr,"            start reading BYTESWAPed MLAPM particles\n");
#else
   fprintf(stderr,"                 start reading MLAPM particles\n\n");
#endif
      
   fprintf(stderr,"sizeof(mlapm_header): %d = 508?\n",sizeof(mlapm_header));
   
   /* read in IO header */
   ReadChars(icfile,mlapm_header.header,MLAPM_HEADERSTRING);
#ifdef LONG4
   ReadInt(icfile,&ino_part,SWAPBYTES);
   mlapm_header.no_part = (long unsigned) ino_part;
#else
   ReadLong(icfile,&lno_part,SWAPBYTES);
   mlapm_header.no_part = (long unsigned) lno_part;
#endif
   ReadDouble(icfile,&(mlapm_header.boxsize),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.omega0),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.lambda0),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.a_initial),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.a_current),SWAPBYTES);
   ReadInt(icfile,&(mlapm_header.no_timestep),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.K_initial),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.U_initial),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.K_current),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.U_current),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.Eintegral),SWAPBYTES);
   ReadDouble(icfile,&(mlapm_header.Econst),SWAPBYTES);
   for (i=0;i<10;i++)
      ReadFloat(icfile,&(mlapm_header.numerics[i]),SWAPBYTES);
#ifdef LONG4
   ReadChars(icfile,mlapm_header.dummy,MLAPM_FILLHEADER+4);
#else
   ReadChars(icfile,mlapm_header.dummy,MLAPM_FILLHEADER);
#endif
   
   
#ifdef VERBOSE
   /* dump MLAPM header (*before* reading particles) */
   fprintf(stderr,"%s\n",mlapm_header.header);
   fprintf(stderr,"mlapm_header.no_part               = %ld\n",mlapm_header.no_part);
   fprintf(stderr,"mlapm_header.no_timestep           = %d\n",mlapm_header.no_timestep);
   fprintf(stderr,"mlapm_header.boxsize               = %g\n",mlapm_header.boxsize);
   fprintf(stderr,"mlapm_header.omega0                = %g\n",mlapm_header.omega0);
   fprintf(stderr,"mlapm_header.lambda0               = %g\n",mlapm_header.lambda0);   
   fprintf(stderr,"mlapm_header.pmass                 = %g\n",mlapm_header.numerics[9]);
   fprintf(stderr,"mlapm_header.cur_reflevel          = %g\n",mlapm_header.numerics[7]);
   fprintf(stderr,"mlapm_header.cur_frcres            = %g\n",mlapm_header.numerics[8]);
   fprintf(stderr,"mlapm_header.a_initial             = %g\n",mlapm_header.a_initial);
   fprintf(stderr,"mlapm_header.a_current             = %g\n",mlapm_header.a_current);
   fprintf(stderr,"mlapm_header.K_initial             = %g\n",mlapm_header.K_initial);
   fprintf(stderr,"mlapm_header.K_current             = %g\n",mlapm_header.K_current);
   fprintf(stderr,"mlapm_header.U_initial             = %g\n",mlapm_header.U_initial);
   fprintf(stderr,"mlapm_header.U_current             = %g\n",mlapm_header.U_current);
   fprintf(stderr,"mlapm_header.Eintegral             = %g\n",mlapm_header.Eintegral);
   fprintf(stderr,"mlapm_header.Econst                = %g\n",mlapm_header.Econst);
#endif /* VERBOSE */
   
   /*transfer everything from MLAPM header to AMIGA header */
   strcpy(io.header.header,mlapm_header.header);
   io.header.no_part     = mlapm_header.no_part;
   io.header.boxsize     = mlapm_header.boxsize;
   io.header.omega0      = mlapm_header.omega0;
   io.header.lambda0     = mlapm_header.lambda0;
   io.header.omegab      = 0.0;
   io.header.gamma       = 0.0;
   io.header.H_frac      = 0.0;
   io.header.T_init      = 0.0;
   io.header.B_init      = 0.0;
   io.header.pmass       = mlapm_header.numerics[9];
   io.header.t_unit      = 1./H0;
   io.header.a_initial   = mlapm_header.a_initial;
   io.header.a_current   = mlapm_header.a_current;
   io.header.timestep    = 0.0;
   io.header.no_timestep = mlapm_header.no_timestep;
   io.header.cur_reflevel= mlapm_header.numerics[7];
   io.header.cur_frcres  = mlapm_header.numerics[8];
   io.header.K_initial   = mlapm_header.K_initial;
   io.header.K_current   = mlapm_header.K_current;
   io.header.U_initial   = mlapm_header.U_initial;
   io.header.U_current   = mlapm_header.U_current;
   io.header.Eintegral   = mlapm_header.Eintegral;
   io.header.Econst      = mlapm_header.Econst;
   
   /* allocate memory */
   io.no_part  = io.header.no_part;
   io.fst_part = c_part(io.no_part);
   io.no_gas   = 0;
   io.fst_gas  = NULL;
   io.no_stars = 0;
   io.fst_star = NULL;

   /* we also convert MLAPM data temporarily to proper physical units...
    * ...as ic_unit_conversion() will determine no_vpart, etc. */
   x_fac = mlapm_header.boxsize;
   v_fac = H0*mlapm_header.boxsize/mlapm_header.a_current;
   m_fac = mlapm_header.numerics[9];
   
   /* loop over all particles */
   for(cur_part = io.fst_part; cur_part < io.fst_part+io.no_part; cur_part++)
     {
      for(j = X; j <= Z; j++)
        {
#ifdef DOUBLE
         ReadDouble(icfile,&pos,SWAPBYTES);
         ReadDouble(icfile,&mom,SWAPBYTES);
#else
         ReadFloat(icfile,&pos,SWAPBYTES);
         ReadFloat(icfile,&mom,SWAPBYTES);
#endif /* DOUBLE */
         cur_part->pos[j] = pos * x_fac;
         cur_part->mom[j] = mom * v_fac;
         
        }
#ifdef MULTIMASS
#ifdef DOUBLE
      ReadDouble(icfile,&weight,SWAPBYTES);
#else
      ReadFloat(icfile,&weight,SWAPBYTES);
#endif /* DOUBLE */  
      cur_part->weight = weight * m_fac; 
#endif /* MULTIMASS */      
     }   
   
   fprintf(stderr,"\n                    finished reading MLAPM\n");
   fprintf(stderr,"===================================================================\n");
}
#endif /* MLAPM */

//========================================================================
//
//                              ASCII
//
//========================================================================
void read_ascii(FILE *icfile)
{
   partptr        cur_part;
   int            j, ino_part;
   long unsigned  ipart, ndummy;
   long           no_species;
   int            i;
   double         dummy;
   char           dummyline[MAXSTRING], cdummy;
   double         r[3], v[3], mass;
      
   fprintf(stderr,"\n==========================================================\n");
   fprintf(stderr,"             start reading ASCII particles\n\n ");
   
   /*============================================================*
    * ordering of ASCII header information
    * ------------------------------------
    * #header_string
    * #no_part
    * #boxsize
    * #omega0
    * #lambda0
    * #omegab
    * #gamma
    * #H_frac
    * #T_init
    * #B_init
    * #z_initial
    * #z_currrent
    * #no_timestep
    *===========================================================*/
   
   fgets(dummyline,MAXSTRING,icfile);
   strcpy(ascii_header.header,dummyline);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %ld",&cdummy, &ascii_header.no_part);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.boxsize);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.omega0);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.lambda0);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.omegab);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.gamma);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.H_frac);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.T_init);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.B_init);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.z_initial);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %lf",&cdummy, &ascii_header.z_current);
   fgets(dummyline,MAXSTRING,icfile);
   sscanf(dummyline,"%c %d",&cdummy,  &ascii_header.no_timestep);
  
   ascii_header.a_initial = 1./(1.+ascii_header.z_initial);
   ascii_header.a_current = 1./(1.+ascii_header.z_current);
   
   
   /*transfer everything from ASCII header to AMIGA header */
   strcpy(io.header.header,ascii_header.header);
   io.header.no_part     = ascii_header.no_part;
   io.header.boxsize     = ascii_header.boxsize;
   io.header.omega0      = ascii_header.omega0;
   io.header.lambda0     = ascii_header.lambda0;
   io.header.omegab      = ascii_header.omegab;
   io.header.gamma       = ascii_header.gamma;
   io.header.H_frac      = ascii_header.H_frac;
   io.header.T_init      = ascii_header.T_init;
   io.header.B_init      = ascii_header.B_init;
   io.header.a_initial   = ascii_header.a_initial;
   io.header.a_current   = ascii_header.a_current;
   io.header.timestep    = 0.0;
   io.header.no_timestep = ascii_header.no_timestep;
   io.header.cur_reflevel= 0.0;
   io.header.cur_frcres  = 0.0;
   io.header.K_initial   = 0.0;
   io.header.K_current   = 0.0;
   io.header.U_initial   = 0.0;
   io.header.U_current   = 0.0;
   io.header.Eintegral   = 0.0;
   io.header.Econst      = 0.0;
   io.header.pmass       = 0.0;    /* will be determined in ic_unit_conversion() */
   io.header.t_unit      = 0.0;    /* will be determined in ic_unit_conversion() */
   
   /* allocate memory for particles */
   io.no_part  = io.header.no_part;
   io.fst_part = c_part(io.no_part);
   io.no_gas   = 0;
   io.fst_gas  = NULL;
   io.no_stars = 0;
   io.fst_star = NULL;
   
  /* loop over all particle to be read in... */
   for(cur_part = io.fst_part; cur_part < io.fst_part+io.no_part; cur_part++)
     {
	   /* read line */
      fgets(dummyline,MAXSTRING,icfile);      
#ifdef ASCII_ONLYPOSITIONS
      sscanf(dummyline,"%lf %lf %lf",r,r+1,r+2);
      v[X] = 0.0;
      v[Y] = 0.0;
      v[Z] = 0.0;
#else
#ifdef MULTIMASS
      sscanf(dummyline,"%lf %lf %lf %lf %lf %lf %lf",r,r+1,r+2,v,v+1,v+2,&mass);
#else
      sscanf(dummyline,"%lf %lf %lf %lf %lf %lf",r,r+1,r+2,v,v+1,v+2);
#endif
#endif // ASCII_ONLYPOSITIONS
      
//      if(cur_part == io.fst_part)
//         fprintf(stderr,"%g %g %g    %g %g %g\n",r[0],r[1],r[2],v[0],v[1],v[2]);
      

      /* transfer to AMIGA particles */
      cur_part->pos[X] = r[X];
      cur_part->pos[Y] = r[Y];
      cur_part->pos[Z] = r[Z];
      
      cur_part->mom[X] = v[X];
      cur_part->mom[Y] = v[Y];
      cur_part->mom[Z] = v[Z];
      
      
#ifdef MULTIMASS
      cur_part->weight = mass; 
#endif
     }
   
   fprintf(stderr,"\n             finished reading ASCII particles\n");
   fprintf(stderr,"==========================================================\n");
}   

#ifdef MARE_NOSTRUM
//========================================================================
//
//                           MareNostrum subboxes
//
//========================================================================
/*==============================================================================
 * read data from Stefan's MareNostrum sub-boxes (requires MULTIMASS)
 *==============================================================================*/
void read_mare_nostrum(FILE *icfile)
{
  partptr        cur_part;
  gasptr         cur_gas;
  long unsigned  npart, ipart, ndummy;
  long           n_cl_DM;
  int            i, k, i1, i6, i7, i14, SWAPBYTES, nspecies;
  float         *x_DM, *y_DM, *z_DM, *vx_DM, *vy_DM, *vz_DM;
  float         *x_gas, *y_gas, *z_gas, *vx_gas, *vy_gas, *vz_gas, *e_gas, *rho_gas;
  float          fdummy;
  FILE          *ftmp;
  
#ifdef MARENOSTRUM_SKIP_DM
  int temp_n_cl_DM;
  float dummy;
#endif
  
#ifdef DEBUG_MARE_NOSTRUM
  ftmp = fopen("test.ascii","w");
#endif
  
  
  fprintf(stderr,"\n==================================================\n");
  fprintf(stderr,"        start reading MareNostrum sub-boxes\n\n");
  
  
#ifdef BYTESWAP
  SWAPBYTES = TRUE;
#else
  SWAPBYTES = FALSE;
#endif
  
  /*-------------------
   *   read header
   *-------------------*/
  fprintf(stderr,"\nheader information:\n");
  fprintf(stderr,"-------------------\n");
  
  MN_SKIP;
  
  ReadFloat  (icfile,&(MN_header.boxsize),SWAPBYTES);
  fprintf(stderr,"MN_header.boxsize    = %g\n",MN_header.boxsize);
  ReadFloat  (icfile,&(MN_header.aexpn),SWAPBYTES);
  fprintf(stderr,"MN_header.aexp       = %g (z=%g)\n",MN_header.aexpn,1./MN_header.aexpn-1.);
  ReadFloat  (icfile,&(MN_header.Om0),SWAPBYTES);
  fprintf(stderr,"MN_header.Om0        = %g\n",MN_header.Om0);
  ReadFloat  (icfile,&(MN_header.Oml0),SWAPBYTES);
  fprintf(stderr,"MN_header.Oml0       = %g\n",MN_header.Oml0);
  ReadFloat  (icfile,&(MN_header.xm_DM),SWAPBYTES);
  fprintf(stderr,"MN_header.xm_DM      = %g\n",MN_header.xm_DM);
  ReadFloat  (icfile,&(MN_header.xm_gas),SWAPBYTES);
  fprintf(stderr,"MN_header.xm_gas     = %g\n",MN_header.xm_gas);
  
  MN_SKIP;
  MN_SKIP;
  
  ReadFloat  (icfile,&(MN_header.boxsize_cl),SWAPBYTES);
  fprintf(stderr,"MN_header.boxsize_cl = %g\n",MN_header.boxsize_cl);
  
  MN_SKIP;
  MN_SKIP;
  
  ReadFloat  (icfile,&(MN_header.xc1),SWAPBYTES);
  fprintf(stderr,"MN_header.xc1        = %g\n",MN_header.xc1);
  ReadFloat  (icfile,&(MN_header.yc1),SWAPBYTES);
  fprintf(stderr,"MN_header.yc1        = %g\n",MN_header.yc1);
  ReadFloat  (icfile,&(MN_header.zc1),SWAPBYTES);
  fprintf(stderr,"MN_header.zc1        = %g\n",MN_header.zc1);
  ReadInt    (icfile,&(MN_header.lkl1),SWAPBYTES);
  fprintf(stderr,"MN_header.lkl1       = %d\n",MN_header.lkl1);
  ReadFloat  (icfile,&(MN_header.xm1),SWAPBYTES);
  fprintf(stderr,"MN_header.xm1        = %g\n",MN_header.xm1);
  ReadFloat  (icfile,&(MN_header.vdisp1),SWAPBYTES);
  fprintf(stderr,"MN_header.vdips1     = %g\n",MN_header.vdisp1);
  ReadFloat  (icfile,&(MN_header.xc),SWAPBYTES);
  fprintf(stderr,"MN_header.xc         = %g\n",MN_header.xc);
  ReadFloat  (icfile,&(MN_header.yc),SWAPBYTES);
  fprintf(stderr,"MN_header.yc         = %g\n",MN_header.yc);
  ReadFloat  (icfile,&(MN_header.zc),SWAPBYTES);
  fprintf(stderr,"MN_header.zc         = %g\n",MN_header.zc);
  ReadInt    (icfile,&(MN_header.i_sub),SWAPBYTES);
  fprintf(stderr,"MN_header.i_sub      = %d\n",MN_header.i_sub);
  
  MN_SKIP;
  MN_SKIP;
  
  ReadFloat  (icfile,&(MN_header.xcnew),SWAPBYTES);
  fprintf(stderr,"MN_header.xcnew      = %g\n",MN_header.xcnew);
  ReadFloat  (icfile,&(MN_header.ycnew),SWAPBYTES);
  fprintf(stderr,"MN_header.ycnew      = %g\n",MN_header.ycnew);
  ReadFloat  (icfile,&(MN_header.zcnew),SWAPBYTES);
  fprintf(stderr,"MN_header.zcnew      = %g\n",MN_header.zcnew);
  ReadFloat  (icfile,&(MN_header.r_vir),SWAPBYTES);
  fprintf(stderr,"MN_header.r_vir      = %g\n",MN_header.r_vir);
  ReadDouble (icfile,&(MN_header.xm_vir),SWAPBYTES);
  fprintf(stderr,"MN_header.xm_vir     = %g\n",MN_header.xm_vir);
  ReadInt    (icfile,&(MN_header.i_DM_vir),SWAPBYTES);
  fprintf(stderr,"MN_header.i_DM_vir   = %d\n",MN_header.i_DM_vir);
  ReadInt    (icfile,&(MN_header.i_gas_vir),SWAPBYTES);
  fprintf(stderr,"MN_header.i_gas_vir  = %d\n",MN_header.i_gas_vir);
  ReadFloat  (icfile,&(MN_header.bfrac_vir),SWAPBYTES);
  fprintf(stderr,"MN_header.bfrac_vir  = %g\n",MN_header.bfrac_vir);
  ReadFloat  (icfile,&(MN_header.over_dens_vir),SWAPBYTES);
  fprintf(stderr,"MN_header.ovdens_vir = %g\n",MN_header.over_dens_vir);
  
  MN_SKIP;
  MN_SKIP;
  
  
  /*-------------------
   * read DM particles
   *-------------------*/
  ReadInt    (icfile,&k,SWAPBYTES);
  ReadInt    (icfile,&(MN_header.n_cl_DM),SWAPBYTES);
  ReadInt    (icfile,&i1,SWAPBYTES);
  ReadInt    (icfile,&i6,SWAPBYTES);
  MN_SKIP;
  MN_SKIP;
  
#ifdef MARENOSTRUM_SKIP_DM
  temp_n_cl_DM = MN_header.n_cl_DM;
  MN_header.n_cl_DM = 1;
#endif
  
  fprintf(stderr,"\no reading %d DM particles ... (%d %d %d)",MN_header.n_cl_DM,k,i1,i6);
  
  /* allocate temporary storage to hold DM particles */
  x_DM  = (float *) calloc(MN_header.n_cl_DM, sizeof(float));
  y_DM  = (float *) calloc(MN_header.n_cl_DM, sizeof(float));
  z_DM  = (float *) calloc(MN_header.n_cl_DM, sizeof(float));
  vx_DM = (float *) calloc(MN_header.n_cl_DM, sizeof(float));
  vy_DM = (float *) calloc(MN_header.n_cl_DM, sizeof(float));
  vz_DM = (float *) calloc(MN_header.n_cl_DM, sizeof(float));
  
#ifdef MARENOSTRUM_SKIP_DM
  fprintf(stderr, "\n*********************************************\n");
  fprintf(stderr, "*********************************************\n");
  fprintf(stderr, "** read_mare_nostrum: SKIPING DM PARTICLES **\n");
  fprintf(stderr, "*********************************************\n");
  fprintf(stderr, "*********************************************\n");
  for(i=0; i<temp_n_cl_DM; i++)
    {
      ReadFloat (icfile,&dummy,SWAPBYTES);
      ReadFloat (icfile,&dummy,SWAPBYTES);
      ReadFloat (icfile,&dummy,SWAPBYTES);
      ReadFloat (icfile,&dummy,SWAPBYTES);
      ReadFloat (icfile,&dummy,SWAPBYTES);
      ReadFloat (icfile,&dummy,SWAPBYTES);
    }
  x_DM[0] = MN_header.xc1;
  y_DM[0] = MN_header.yc1;
  z_DM[0] = MN_header.zc1;
  vx_DM[0] = 0.0;
  vy_DM[0] = 0.0;
  vz_DM[0] = 0.0;
#else
  for(i=0; i<MN_header.n_cl_DM; i++)
    {
      ReadFloat (icfile,&x_DM[i],SWAPBYTES);
      ReadFloat (icfile,&y_DM[i],SWAPBYTES);
      ReadFloat (icfile,&z_DM[i],SWAPBYTES);
      ReadFloat (icfile,&vx_DM[i],SWAPBYTES);
      ReadFloat (icfile,&vy_DM[i],SWAPBYTES);
      ReadFloat (icfile,&vz_DM[i],SWAPBYTES);
    }
#endif
  fprintf(stderr,"done\n");
  MN_SKIP;
  MN_SKIP;
  
  
  /*--------------------
   * read gas particles
   *--------------------*/
  ReadInt    (icfile,&k,SWAPBYTES);
  ReadInt    (icfile,&(MN_header.n_cl_gas),SWAPBYTES);
  ReadInt    (icfile,&i7,SWAPBYTES);
  ReadInt    (icfile,&i14,SWAPBYTES);
  MN_SKIP;
  MN_SKIP;
  
  fprintf(stderr,"o reading %d gas particles ... (%d %d %d)",MN_header.n_cl_gas,k,i7,i14);
  
  /* allocate temporary storage to hold gas particles */
  x_gas   = (float *) calloc(MN_header.n_cl_gas, sizeof(float));
  y_gas   = (float *) calloc(MN_header.n_cl_gas, sizeof(float));
  z_gas   = (float *) calloc(MN_header.n_cl_gas, sizeof(float));
  vx_gas  = (float *) calloc(MN_header.n_cl_gas, sizeof(float));
  vy_gas  = (float *) calloc(MN_header.n_cl_gas, sizeof(float));
  vz_gas  = (float *) calloc(MN_header.n_cl_gas, sizeof(float));
  e_gas   = (float *) calloc(MN_header.n_cl_gas, sizeof(float));
  rho_gas = (float *) calloc(MN_header.n_cl_gas, sizeof(float));
  
  for(i=0; i<MN_header.n_cl_gas; i++)
    {
      ReadFloat (icfile,&x_gas[i],SWAPBYTES);
      ReadFloat (icfile,&y_gas[i],SWAPBYTES);
      ReadFloat (icfile,&z_gas[i],SWAPBYTES);
      ReadFloat (icfile,&vx_gas[i],SWAPBYTES);
      ReadFloat (icfile,&vy_gas[i],SWAPBYTES);
      ReadFloat (icfile,&vz_gas[i],SWAPBYTES);
      ReadFloat (icfile,&e_gas[i],SWAPBYTES);
      ReadFloat (icfile,&rho_gas[i],SWAPBYTES);
    }
  fprintf(stderr,"done\n\n");
  
  
  /*------------------------------
   * transfer to AMIGA structures
   *------------------------------*/
#ifdef MARE_NOSTRUM_DM_ONLY
  nspecies           = 1;
  npart              = (long) (MN_header.n_cl_DM);
  MN_header.n_cl_gas = 0;
#else /* MARE_NOTRUM_DM_ONLY */
#ifdef MARE_NOSTRUM_GAS_ONLY
  nspecies           = 1;
  npart              = (long) (MN_header.n_cl_gas);
  MN_header.n_cl_DM  = 0;
#else /* MARE_NOSTRUM_GAS_ONLY */
  /* STANDARD call that uses both DM and GAS */
  nspecies           = 2;
  npart              = (long) (MN_header.n_cl_DM + MN_header.n_cl_gas);
#endif /* MARE_NOSTRUM_GAS_ONLY */
#endif /* MARE_NOTRUM_DM_ONLY */
  
  /* allocate memory */
  io.fst_part       = c_part(npart);
  io.no_part        = npart;
  
  io.fst_gas    = (gasptr) calloc(MN_header.n_cl_gas, sizeof(gasptr));
  io.no_gas     = MN_header.n_cl_gas;
  io.offset_gas = MN_header.n_cl_DM;
  
  io.fst_star   = NULL;
  io.no_stars   = 0;

  /* DM particles */
  fprintf(stderr,"o transfering DM particles             0 - %12d to AMIGA structure... ",MN_header.n_cl_DM-1);
  
  for(ipart=0; ipart<(long)MN_header.n_cl_DM; ipart++)
    {
      cur_part =  io.fst_part + ipart;
      
      cur_part->pos[X] = x_DM[ipart]+MN_header.xc1;
      cur_part->pos[Y] = y_DM[ipart]+MN_header.yc1;
      cur_part->pos[Z] = z_DM[ipart]+MN_header.zc1;
      
      cur_part->mom[X] = vx_DM[ipart];
      cur_part->mom[Y] = vy_DM[ipart];
      cur_part->mom[Z] = vz_DM[ipart];
      
      cur_part->weight = MN_header.xm_DM;
      
#ifdef DEBUG_MARE_NOSTRUM
      fprintf(ftmp,"%g %g %g   %g %g %g\n",
              x_DM[ipart],y_DM[ipart],z_DM[ipart],
              vx_DM[ipart],vy_DM[ipart],vz_DM[ipart]);
#endif
    }
  fprintf(stderr,"done\n");
  
  /* gas particles */
  fprintf(stderr,"o transfering gas particles %12d - %12d to AMIGA structure... ",MN_header.n_cl_DM,npart-1);
  for(ipart=MN_header.n_cl_DM; ipart<npart; ipart++)
    {
      cur_part =  io.fst_part + ipart;
      cur_gas  =  io.fst_gas  + (ipart-MN_header.n_cl_DM);
      
      cur_part->pos[X] = x_gas[ipart-MN_header.n_cl_DM]+MN_header.xc1;
      cur_part->pos[Y] = y_gas[ipart-MN_header.n_cl_DM]+MN_header.yc1;
      cur_part->pos[Z] = z_gas[ipart-MN_header.n_cl_DM]+MN_header.zc1;
      
      cur_part->mom[X] = vx_gas[ipart-MN_header.n_cl_DM];
      cur_part->mom[Y] = vy_gas[ipart-MN_header.n_cl_DM];
      cur_part->mom[Z] = vz_gas[ipart-MN_header.n_cl_DM];
      
      cur_part->weight = MN_header.xm_gas;
      
      cur_gas->u       = e_gas[ipart-MN_header.n_cl_DM];
      //cur_gas->u       = 0.0;
      
#ifdef DEBUG_MARE_NOSTRUM
      fprintf(ftmp,"%g %g %g   %g %g %g\n",
              x_gas[ipart-MN_header.n_cl_DM],y_gas[ipart-MN_header.n_cl_DM],z_gas[ipart-MN_header.n_cl_DM],
              vx_gas[ipart-MN_header.n_cl_DM],vy_gas[ipart-MN_header.n_cl_DM],vz_gas[ipart-MN_header.n_cl_DM]);
#endif
    }
  fprintf(stderr,"done\n\n");
  
  /* transfer everything from MN_header header to AMIGA header */
  strcpy(io.header.header,"MareNostrum-SubBoxAnalysis");
  io.header.no_species  = nspecies;
  io.header.no_part     = io.no_part;
  io.header.boxsize     = MN_header.boxsize;
  io.header.omega0      = MN_header.Om0;
  io.header.lambda0     = MN_header.Oml0;
  io.header.omegab      = 0.0;
  io.header.gamma       = 0.0;
  io.header.H_frac      = 0.0;
  io.header.T_init      = 0.0;
  io.header.B_init      = 0.0;
  io.header.pmass       = 0.0;    /* will be determined in ic_unit_conversion() */
  io.header.t_unit      = 0.0;    /* will be determined in ic_unit_conversion() */
  io.header.a_initial   = 0.01;
  io.header.a_current   = MN_header.aexpn;
  io.header.timestep    = 0.0;
  io.header.no_timestep = 0;
  io.header.cur_reflevel= 0.0;
  io.header.cur_frcres  = 0.0;
  io.header.K_initial   = 0.0;
  io.header.K_current   = 0.0;
  io.header.U_initial   = 0.0;
  io.header.U_current   = 0.0;
  io.header.Eintegral   = 0.0;
  io.header.Econst      = 0.0;
  
  fprintf(stderr,"\n      finished reading MareNostrum sub-boxes\n");
  fprintf(stderr,"==================================================\n");
  
#ifdef DEBUG_MARE_NOSTRUM
  fclose(ftmp);
#endif
  
  /* free temporary arrays */
  free(x_DM);
  free(y_DM);
  free(z_DM);
  free(vx_DM);
  free(vy_DM);
  free(vz_DM);
  
  free(x_gas);
  free(y_gas);
  free(z_gas);
  free(vx_gas);
  free(vy_gas);
  free(vz_gas);
  free(e_gas);
  free(rho_gas);
}   
#endif /* MARE_NOSTRUM */

#ifdef GADGET
//========================================================================
//
//                              GADGET
//
//========================================================================
/*==============================================================================
 * read data from initial conditions file provided in GADGET format
 * kindly provided by Chris Power, August 2004
 * adapted for multiple GADGET files and stars/gas by A. Knebe, August 2005
 *==============================================================================*/
 void skim_gadget(FILE *icfile)
{
   int    i,j,k;            /* Dummy variables */
   char   DATA[MAXSTRING];
   
   int    no_part;    /* Number of particles in simulation */
   float  dummy[3];   /* Dummy float variable, for reading (x,y,z,vx,vy,vz,m) */
   
   double x_fac, v_fac, m_fac;
   
   int    massflag;           /* Are the masses contained in the header? */
   int    SWAPBYTES, version;
   
  /* determine file version */
  version = check_gadgetversion(icfile);
  
#ifdef BYTESWAP
   SWAPBYTES = TRUE;
#else
   SWAPBYTES = FALSE;
#endif
   
  /*================= read in GADGET IO header =================*/
  if(version==2){
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
  }
  
   GADGET_SKIP;
   
#ifdef BYTESWAP
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
#else
   fread(&gadget.header,sizeof(gadget.header),1,icfile);
#endif
   fprintf(stderr,"expansion factor: %lf\n",gadget.header.expansion);
   fprintf(stderr,"redshift:         %lf\n",gadget.header.redshift);
   fprintf(stderr,"boxsize:          %lf (%lf Mpc/h)\n",gadget.header.BoxSize,gadget.header.BoxSize*GADGET_LUNIT);
   fprintf(stderr,"omega0:           %lf\n",gadget.header.Omega0);
   fprintf(stderr,"lambda0:          %lf\n",gadget.header.OmegaLambda);
   fprintf(stderr,"HubbleParam:      %lf\n",gadget.header.HubbleParam);
   fprintf(stderr,"gas:    np[0]=%9d\t nall[0]=%9d\t massarr[0]=%g\n",gadget.header.np[0],gadget.header.nall[0],gadget.header.massarr[0]); 
   fprintf(stderr,"halo:   np[1]=%9d\t nall[1]=%9d\t massarr[1]=%g\n",gadget.header.np[1],gadget.header.nall[1],gadget.header.massarr[1]); 
   fprintf(stderr,"disk:   np[2]=%9d\t nall[2]=%9d\t massarr[2]=%g\n",gadget.header.np[2],gadget.header.nall[2],gadget.header.massarr[2]); 
   fprintf(stderr,"bulge:  np[3]=%9d\t nall[3]=%9d\t massarr[3]=%g\n",gadget.header.np[3],gadget.header.nall[3],gadget.header.massarr[3]); 
   fprintf(stderr,"stars:  np[4]=%9d\t nall[4]=%9d\t massarr[4]=%g\n",gadget.header.np[4],gadget.header.nall[4],gadget.header.massarr[4]); 
   fprintf(stderr,"bndry:  np[5]=%9d\t nall[5]=%9d\t massarr[5]=%g\n",gadget.header.np[5],gadget.header.nall[5],gadget.header.massarr[5]); 
   
   GADGET_SKIP;
   
   /*================= read in GADGET IO header =================*/
   
   /* do we need to read a mass array? (also count no_part in current file...) */
   massflag   = 0;
   no_part    = 0;
   for(i=0;i<6;i++) 
     {
      no_part += gadget.header.np[i];
      if(gadget.header.massarr[i] < MZERO && gadget.header.np[i] > 0)
         massflag=1;  
     }
   
   /* allocate particle array */
   if(!(P_gadget=malloc(no_part*sizeof(struct particle_data))))
     {
      fprintf(io.logfile,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
   
   gadget.np[0][gadget.i_gadget_file] = gadget.header.np[0];
   gadget.np[1][gadget.i_gadget_file] = gadget.header.np[1];
   gadget.np[2][gadget.i_gadget_file] = gadget.header.np[2];
   gadget.np[3][gadget.i_gadget_file] = gadget.header.np[3];
   gadget.np[4][gadget.i_gadget_file] = gadget.header.np[4];
   gadget.np[5][gadget.i_gadget_file] = gadget.header.np[5];
   
   /*================= read in GADGET particles =================*/
  if(version==2){
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
  }
   
   GADGET_SKIP;
   
   for(i=0;i<no_part;i++)
     {
#ifdef BYTESWAP
      ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
      ReadFloat(icfile,&(dummy[1]),SWAPBYTES);
      ReadFloat(icfile,&(dummy[2]),SWAPBYTES);
#else
      fread(&dummy[0],sizeof(float),3,icfile);
#endif      
     }
   
   GADGET_SKIP;
   
   /*================= read in GADGET particles =================*/
   
   
   
   /*================= read in GADGET velocities =================*/
  if(version==2){
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
  }
   
   GADGET_SKIP;   
   
   for(i=0;i<no_part;i++)
     {
#ifdef BYTESWAP
      ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
      ReadFloat(icfile,&(dummy[1]),SWAPBYTES);
      ReadFloat(icfile,&(dummy[2]),SWAPBYTES);
#else
      fread(&dummy[0],sizeof(float),3,icfile);
#endif 
     }
   
   GADGET_SKIP;
   
   /*================= read in GADGET velocities =================*/
   
   
   /*================= read in GADGET id's =================*/   
  if(version==2){
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
  }
   
   GADGET_SKIP;
   
   /* we need to read in *all* ID's as the masses are behind the ID's in the file... */
#ifdef LGADGET      
   for(i=0;i<no_part;i++)
     {
#ifdef BYTESWAP
      ReadLongLong(icfile,&(P_gadget[i].ID),SWAPBYTES);
#else
      fread(&P_gadget[i].ID,sizeof(long long),1,icfile);
#endif
     }
#else /* LGADGET */
   for(i=0;i<no_part;i++)
     {
#ifdef BYTESWAP
      ReadInt(icfile,&(P_gadget[i].ID),SWAPBYTES);
#else
      fread(&P_gadget[i].ID,sizeof(int),1,icfile);
#endif
     }
#endif
   
   GADGET_SKIP;
   
   /*================= read in GADGET id's =================*/
   
   
   k = 0;
   
   /* massflag == 1 indicates that massarr[i] = 0 and hence need to read in particle masses */
   if(massflag==1) 
     {
      /*================= read in GADGET individual particle masses =================*/
      if(version==2){
        GADGET_SKIP;
        
        fread(DATA,sizeof(char),4,icfile);
        DATA[4] = '\0';
        GADGET_SKIP;
        
        GADGET_SKIP;
      }
      
      GADGET_SKIP; 
      
      /* we are only interested in the masses of the halo species [1] */
      for(i=0;i<2;i++)
        {
         if (gadget.header.np[i] > 0 && gadget.header.massarr[i] < MZERO  ) 
           {
            
            fprintf(stderr,"\n reading masses from file for species [%d] (massarr=%g), ",i,gadget.header.massarr[i]);
            
            for(j=0; j<gadget.header.np[i]; j++)
              {
#ifdef BYTESWAP
               ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
#else
               fread(&dummy[0],sizeof(float),1,icfile);
#endif
               P_gadget[k].Mass  = dummy[0];
               
               k++;
              }
           }
         else
           {
            /* simply copy appropriate massarr[i] to particles */
            for(j=0; j<gadget.header.np[i]; j++) 
              {
               P_gadget[k].Mass = gadget.header.massarr[i];
               k++;
              }
           }
        }
      
      GADGET_SKIP;
      
      /*================= read in GADGET individual particle masses =================*/
     } 
   
   /* simply copy appropriate massarr[i] to particles */
   else 
     {
      k=0;
      
      /* we are only interested in the masses of the halo species [1] */
      for(i=0;i<2;i++)
        {
         for(j=0;j<gadget.header.np[i];j++) 
           {
            P_gadget[k].Mass = gadget.header.massarr[i];
            k++;
           }
        }
     }
   
   /* adjust IDmin, IDmax, and DMmmin using current halo paricles */
   for(k=gadget.header.np[0]; k<gadget.header.np[0]+gadget.header.np[1]; k++)
     {
      if(P_gadget[k].Mass < gadget.mmin)
         gadget.mmin  = P_gadget[k].Mass;
      
      if(P_gadget[k].ID   < gadget.IDmin)
         gadget.IDmin = P_gadget[k].ID;
      
      if(P_gadget[k].ID   > gadget.IDmax)
         gadget.IDmax = P_gadget[k].ID;
     }
   
   gadget.nall1 = gadget.header.nall[1];
   
   
   /* now it's time to free GADGET's P array */
   free(P_gadget);
}


/*==========================================
* read_gadget():   actual reading routine
*==========================================*/
void read_gadget(FILE *icfile)
{
   partptr        cur_part;
   gasptr         cur_gas;
   
   long unsigned  ipart;
#ifdef LGADGET
   long long      ID, IDmin, IDmax;
#else
   long           ID, IDmin, IDmax;
#endif
   
   double         tot_mass[6];
   double         DMmmin;
   
   int            i,j,k;
   int            no_part, ppart;
   int            massflag;
   int            SWAPBYTES;
   char           DATA[MAXSTRING];
   float          dummy[3];
   double         x_fac, v_fac, m_fac, u_fac;
  int version;

  /* determine file version */
  version = check_gadgetversion(icfile);

  
   fprintf(stderr,"\n===================================================================\n");
#ifdef BYTESWAP
   SWAPBYTES = TRUE;
   fprintf(stderr,"           start reading BYTESWAPed GADGET particles\n");
#else
   SWAPBYTES = FALSE;
   fprintf(stderr,"                 start reading GADGET particles\n");
#endif
   
   /*================= read in GADGET IO header =================*/
  if(version==2){
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
  }
   
   GADGET_SKIP;
   
#ifdef BYTESWAP
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
#else
   fread(&gadget.header,sizeof(gadget.header),1,icfile);
#endif
   
   GADGET_SKIP;
   
   fprintf(stderr,"\nFinished reading GADGET header\n\n");
   /*================= read in GADGET IO header =================*/
   
   
   /* keep track of no. of particles in each GADGET file */
   gadget.np[0][gadget.i_gadget_file] = gadget.header.np[0];
   gadget.np[1][gadget.i_gadget_file] = gadget.header.np[1];
   gadget.np[2][gadget.i_gadget_file] = gadget.header.np[2];
   gadget.np[3][gadget.i_gadget_file] = gadget.header.np[3];
   gadget.np[4][gadget.i_gadget_file] = gadget.header.np[4];
   gadget.np[5][gadget.i_gadget_file] = gadget.header.np[5];

   
   /* count total no. of particles in current file (and set massflag) */
   massflag   = 0;
   no_part    = 0;
   for(i=0;i<6;i++) 
     {
      no_part += gadget.header.np[i];
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
   
   /* allocate particle array */
   if(!(P_gadget=(struct particle_data *) calloc(no_part, sizeof(struct particle_data))))
     {
      fprintf(io.logfile,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
   
   /*================= read in GADGET particles =================*/
  if(version==2){
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
  }
  else{
    fprintf(stderr,"reading ");
  }
   
   GADGET_SKIP;
   fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
      
   for(i=0;i<no_part;i++)
     {
#ifdef BYTESWAP
      ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
      ReadFloat(icfile,&(dummy[1]),SWAPBYTES);
      ReadFloat(icfile,&(dummy[2]),SWAPBYTES);
#else
      fread(&dummy[0],sizeof(float),3,icfile);
#endif      
      P_gadget[i].Pos[0]=dummy[0];
      P_gadget[i].Pos[1]=dummy[1];
      P_gadget[i].Pos[2]=dummy[2];      
     }
   fprintf(stderr,"Pos[X]=%12.6g Pos[Y]=%12.6g Pos[Z]=%12.6g ... ",P_gadget[no_part-1].Pos[X],P_gadget[no_part-1].Pos[Y],P_gadget[no_part-1].Pos[Z]);
   
   GADGET_SKIP;
   fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
   
   /*================= read in GADGET particles =================*/
   
   
   
   /*================= read in GADGET velocities =================*/
  if(version==2){
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
  }
  else {
    fprintf(stderr,"reading ");
  }
   
   GADGET_SKIP;
   fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
   
   for(i=0;i<no_part;i++)
     {
#ifdef BYTESWAP
      ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
      ReadFloat(icfile,&(dummy[1]),SWAPBYTES);
      ReadFloat(icfile,&(dummy[2]),SWAPBYTES);
#else
      fread(&dummy[0],sizeof(float),3,icfile);
#endif 
      P_gadget[i].Vel[0]=dummy[0];
      P_gadget[i].Vel[1]=dummy[1];
      P_gadget[i].Vel[2]=dummy[2]; 
     }
   fprintf(stderr,"Vel[X]=%12.6g Vel[Y]=%12.6g Vel[Z]=%12.6g ... ",P_gadget[no_part-1].Vel[X],P_gadget[no_part-1].Vel[Y],P_gadget[no_part-1].Vel[Z]);
   
   GADGET_SKIP;
   fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
   
   /*================= read in GADGET velocities =================*/
   
   
   /*================= read in GADGET id's =================*/
  if(version==2){
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
  }
  else {
    fprintf(stderr,"reading ");
  }
   
   GADGET_SKIP;
   fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
   
   for(i=0;i<no_part;i++)
     {
#ifdef LGADGET
#ifdef BYTESWAP
      ReadLongLong(icfile,&(P_gadget[i].ID),SWAPBYTES);
#else
      fread(&P_gadget[i].ID,sizeof(long long),1,icfile);
#endif
#else /* LGADGET */
#ifdef BYTESWAP
      ReadInt(icfile,&(P_gadget[i].ID),SWAPBYTES);
#else
      fread(&P_gadget[i].ID,sizeof(int),1,icfile);
#endif
#endif /* LGADGET */
     }
   
   fprintf(stderr,"ID=%12ld ...  ",P_gadget[no_part-1].ID);

   GADGET_SKIP;
   fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
   /*================= read in GADGET id's =================*/
      
   
   k = 0;
   /* massflag == 1 indicates that massarr[i] = 0 and hence need to read in particle masses */
   if(massflag==1) 
     {
      /*================= read in GADGET individual particle masses =================*/
      if(version==2){
        GADGET_SKIP;
        
        fread(DATA,sizeof(char),4,icfile);
        DATA[4] = '\0';
        GADGET_SKIP;
        
        GADGET_SKIP;
        fprintf(stderr,"reading %s",DATA);
      }
      else {
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
#ifdef BYTESWAP
               ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
#else
               fread(&dummy[0],sizeof(float),1,icfile);
#endif
               P_gadget[k].Mass  = dummy[0];
               tot_mass[i]      += dummy[0];

               k++;
              }
           }
         else
           {
            /* simply copy appropriate massarr[i] to particles */
            for(j=0; j<gadget.header.np[i]; j++) 
              {
               P_gadget[k].Mass = gadget.header.massarr[i];
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
            P_gadget[k].Mass = gadget.header.massarr[i];
            k++;
           }
         tot_mass[i] = gadget.header.np[i]*gadget.header.massarr[i];
        }
     }
   
   /*================= read in GADGET gas particle energies =================*/
   if(gadget.header.np[0] > 0) 
     {      
       if(version==2){
         GADGET_SKIP;
         
         fread(DATA,sizeof(char),4,icfile);
         DATA[4] = '\0';
         GADGET_SKIP;
         
         GADGET_SKIP;
         fprintf(stderr,"reading %s",DATA);
       }
       else {
         fprintf(stderr,"reading ");
       }
      
      GADGET_SKIP; 
      fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
      
      for(i=0; i<gadget.header.np[0]; i++)
        {
#ifdef BYTESWAP
         ReadFloat(icfile,&(dummy[0]),SWAPBYTES);
#else
         fread(&dummy[0],sizeof(float),1,icfile);
#endif
         /* store additional gas particle property */
         P_gadget[i].u = dummy[0];         
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
   
   fprintf(stderr,"\n<- finished reading GADGET particles\n");

   
   /* single GADGET file => determine DMmmin, IDmin, IDmax */
   if(gadget.no_gadget_files == 1)
     {
      if(gadget.header.massarr[1] > MZERO)
         DMmmin = gadget.header.massarr[1];
      else
        {
         DMmmin = 1e30;
         for(k=gadget.header.np[0]; k<gadget.header.np[0]+gadget.header.np[1]; k++)
           {
            if(P_gadget[k].Mass < DMmmin)
               DMmmin = P_gadget[k].Mass;
           }
        }
      gadget.mmin = DMmmin*GADGET_MUNIT;
      
#ifdef LGADGET
      IDmin = (long long) pow(2,63);
#else
      IDmin = 2147483647;
#endif
      IDmax = 0;
      for(i=gadget.header.np[0]; i<gadget.header.np[0]+gadget.header.np[1]; i++)
        {
         if(P_gadget[i].ID < IDmin)
            IDmin = P_gadget[i].ID;
         
         if(P_gadget[i].ID > IDmax)
            IDmax = P_gadget[i].ID;
        }         
      
      if(IDmax-IDmin+1 != gadget.header.nall[1])
        {
         fprintf(stderr,"\nPROBLEM:  halo partice ID's are not consecutive from IDmin (%ld) to IDmax (%ld) [%d]\n\n",IDmin,IDmax,gadget.header.nall[1]);
         exit(0);
        }
      
      gadget.IDmin = IDmin;
      gadget.IDmax = IDmax;

       fprintf(stderr,"\nGADGET file: mmin=%g, IDmin=%ld, IDmax=%ld\n", gadget.mmin, gadget.IDmin, gadget.IDmax);
     }
   else
     {
      /* for multiple GADGET files skim_gadget() already determined mmin, IDmin, IDmax! */
       fprintf(stderr,"\nGADGET file: mmin=%g, IDmin=%ld, IDmax=%ld (as already determined by skim_gadget())\n", gadget.mmin, gadget.IDmin, gadget.IDmax);
     }
   
      
   
   /*====================================== AMIGA header ======================================*/
#ifdef GADGET_GAS_ONLY
   io.no_part      = gadget.header.nall[0];
   
   io.no_gas       = gadget.header.nall[0];
   io.offset_gas   = 0;

   io.no_stars     = 0;
#else
#ifdef GADGET_STARS_ONLY
   io.no_part      = gadget.header.nall[4];
   
   io.no_gas       = 0;
   io.offset_gas   = 0;
   
   io.no_stars     = gadget.header.nall[4];
   io.offset_stars = 0;
#else
   /* STANDARD:  we use all particles found in the input file... */
   io.no_part = 0;
   for(i=0;i<6;i++)
      io.no_part  += gadget.header.nall[i];

   io.no_gas       = gadget.header.nall[0];
   io.offset_gas   = gadget.IDmax-gadget.IDmin+1;

   io.no_stars     = gadget.header.nall[4];
   io.offset_stars = gadget.IDmax-gadget.IDmin+1 + io.no_gas + gadget.header.nall[2] + gadget.header.nall[3];
   
#endif /* STARS */
#endif /* GAS */
   
   /* either allocate memory for AMIGA particles... */
   if(gadget.i_gadget_file == 0)
     {
      io.fst_part = c_part(io.no_part);
      io.fst_gas  = c_gas(io.no_gas);
      io.fst_star = NULL;              /* we only store the star-mass and hence do not need a fst_star[] array!*/
      
      /* cur_part will be calculated based upon P_gadget[].ID */
      cur_gas     = io.fst_gas;
      /* cur_star is not in use! */

     }
   /* ...or attach to existing memory */
   else
     {
      // cur_part will be calculated based upon P_gadget[].ID !!!
      
      cur_gas     = gadget.lst_gas + 1;
      /* cur_star is not in use! */
     }
   
   
   strcpy(io.header.header,"GADGET simulation");
   io.header.no_part     = io.no_part;
   io.header.boxsize     = gadget.header.BoxSize * GADGET_LUNIT;
   io.header.omega0      = gadget.header.Omega0;
   io.header.lambda0     = gadget.header.OmegaLambda;
   io.header.omegab      = 0.0;                               /* omegab could in principle be calculated! */
   io.header.gamma       = 0.0;
   io.header.H_frac      = 0.0;
   io.header.T_init      = 0.0;
   io.header.B_init      = 0.0;
   io.header.pmass       = 0.0;                               /* will be determined in ic_unit_conversion() */
   io.header.t_unit      = 0.0;                               /* will be determined in ic_unit_conversion() */
   io.header.a_current   = 1.0/(1.0+gadget.header.redshift);  
   io.header.a_initial   = 1.0/(1.0+gadget.header.redshift);  /* only correct for IC data */
   io.header.timestep    = 0.0;
   io.header.no_timestep = 0;
   io.header.cur_reflevel= 0.0;
   io.header.cur_frcres  = 0.0;
   io.header.K_initial   = 0.0;
   io.header.K_current   = 0.0;
   io.header.U_initial   = 0.0;
   io.header.U_current   = 0.0;
   io.header.Eintegral   = 0.0;
   io.header.Econst      = 0.0;
   
   
   /*====================================== AMIGA particles ======================================*/
   /* convert to Mpc/h, km/sec, Msun/h */
   x_fac  = GADGET_LUNIT;
   v_fac  = sqrt(io.header.a_current);
   m_fac  = GADGET_MUNIT;
   u_fac  = 1.0;

   
   /*
    we are now doing some nasty things to the particles ID's ;-)
    
    actually not really nasty, but we are shifting the halo particles to the front, i.e.
    the new block of particles stored as io.fst_part -> io.fst_part+io.no_part
    will look as follows:
    
    [halo][gas][disk][bulge][stars][bndry]  (instead of [gas][halo][disk][bulge][stars][bndry])
    
    this allows us to readily use the value of cur_part in a loop over all particles, i.e.
    for(ipart=0; ipart<io.no_part; ipart++)
    cur_part = io.fst_part + ipart;
    to uniquely identify halo particles across files!
    
    or in other words: we do not need to store halo ID's as we know that 
    io.fst_part -> io.fst_part + gadget.header.np[1]-1
    will always be arranged in the same manner...
    
    */
   
   /*--------------------------
    * loop over all particles
    *--------------------------*/
#ifdef GADGET_DM_ONLY
  ppart = gadget.header.np[0];
  i     = 1;
#else
#ifdef GADGET_GAS_ONLY
   ppart = 0;
   i     = 0;
#else
#ifdef GADGET_STARS_ONLY
   ppart = gadget.header.np[0]+gadget.header.np[1]+gadget.header.np[2]+gadget.header.np[3];
   i     = 4;
#else
   ppart = 0;
   for(i=0; i<6; i++)
#endif /* GADGET_STARS_ONLY */
#endif /* GADGET_GAS_ONLY */
#endif /* GADGET_DM_ONLY */
     {
      for(j=0; j<gadget.header.np[i]; j++)
        {
         /* gas properties */
         if(i == 0)
           {
             // we are not doing any unit conversion later on 
             // and hence cur_gas->u should already be here
             // set to internal AMIGA units!
             // NOTE: this is the case for GADGET files!!!!!
            cur_gas->u = (flouble) P_gadget[ppart].u * u_fac;
            cur_gas++;
           }
         
         /* always use same halo particle ID's to allow for identification across multiple output files */
         if(i == 1)
           {
            /* halo particle ID's will now be in the range [0, nall[1]-1] */
            ID       = (long) P_gadget[ppart].ID - gadget.IDmin;
            cur_part = io.fst_part + ID;
           }
         
         /* manually assign ID's for gas, stars, etc. (not consistent across different snapshots!) */
         else
           {
            /* jump past all halo particle ID's */
            ipart = gadget.IDmax-gadget.IDmin+1;
            
            /* jump past all other species ID's */
            for(k=0; k<i; k++)
               /* do not count nall[1] twice! */
               if(k != 1)
                  ipart += gadget.header.nall[k];
            
            /* take care of particles across P-GADGET output files */
            for(k=0; k<gadget.i_gadget_file; k++)
               ipart += gadget.np[i][k];
            
            /* ...and finally assign an ID to the current particle */
            ipart += j;
            
            cur_part = io.fst_part + ipart;
           }

         /* position and velocity */
         for(k = X; k <= Z; k++)
           {
            cur_part->pos[k] = (flouble) P_gadget[ppart].Pos[k] * x_fac;
            cur_part->mom[k] = (flouble) P_gadget[ppart].Vel[k] * v_fac;
           }
         
         /* mass */
         cur_part->weight    = (flouble) P_gadget[ppart].Mass   * m_fac;
          
         /* GADGET ID */
         cur_part->id        =           P_gadget[ppart].ID;
         
#ifdef GAS_PARTICLES
         if(i == 0) {
           cur_part->u       = P_gadget[ppart].u * u_fac;
         }
         else if (i==4) {
           cur_part->u       = PSTAR;
         }
         else {
           cur_part->u       = PDM;
         }
#endif
         
         /* move to next particle in P_gadget[] array */
         ppart++;         
        } /* j<np[i] */

     } /* i=0,5 */
   

   /* remember last particle */
   if(gadget.i_gadget_file == 0)
     {
      /* cur_part will be calculated based upon P_gadget[].ID */
      gadget.lst_gas = cur_gas;
      /* cur_star is not in use! */
      
     }

   fprintf(stderr,"\n              finished reading GADGET particles\n");
   fprintf(stderr,"===================================================================\n");

   /* now it's time to free GADGET's P array */
   free(P_gadget);
}
#endif /* GADGET */

#ifdef ART
//========================================================================
//
//                              ART
//
//========================================================================
void read_ARTbinaries(FILE *fPMcrd, FILE *fPMcrs)
{
  int   i;
  char  c;
  
  float *recdat;
  float xpar, ypar, zpar;
  float vxpar,vypar,vzpar;
  float *wspecies;
  int   *lspecies;
  int   nspecies;
  long   npart;
  double nvpart, pmass, partw, boxsize;
  
  int   NROW, NGRID, NPAGE, NRECL;
  
  double x_fac, v_fac, m_fac;
  
  int   Npages, N_in_last, N_particles, IROW, In_page, IN, Icount, Ipart, iL;
  float w;
  
  partptr cur_part;
  
  /* only PMcrd.DAT requires a SKIP */
  ART_SKIP(fPMcrd);
  
  /*=============================
   * read ART header information
   *=============================*/
  for(i=0; i<45; i++)
    {
      fread(&c, sizeof(char), 1, fPMcrd);
      strcpy(&(ART_header.header_string[i]), &c);
    }
  fread(&ART_header.aexpn,    sizeof(float), 1, fPMcrd);
  fread(&ART_header.aexp0,    sizeof(float), 1, fPMcrd);
  fread(&ART_header.amplt,    sizeof(float), 1, fPMcrd);
  fread(&ART_header.astep,    sizeof(float), 1, fPMcrd);
  fread(&ART_header.istep,    sizeof(int),   1, fPMcrd);
  fread(&ART_header.partw,    sizeof(float), 1, fPMcrd);
  fread(&ART_header.tintg,    sizeof(float), 1, fPMcrd);
  fread(&ART_header.ekin,     sizeof(float), 1, fPMcrd);
  fread(&ART_header.ekin1,    sizeof(float), 1, fPMcrd);
  fread(&ART_header.ekin2,    sizeof(float), 1, fPMcrd);
  fread(&ART_header.au0,      sizeof(float), 1, fPMcrd);
  fread(&ART_header.aeu0,     sizeof(float), 1, fPMcrd);
  fread(&ART_header.nrowc,    sizeof(int),   1, fPMcrd);
  fread(&ART_header.ngridc,   sizeof(int),   1, fPMcrd);
  fread(&ART_header.nspecies, sizeof(int),   1, fPMcrd);
  fread(&ART_header.nseed,    sizeof(int),   1, fPMcrd);
  fread(&ART_header.Om0,      sizeof(float), 1, fPMcrd);
  fread(&ART_header.Oml0,     sizeof(float), 1, fPMcrd);
  fread(&ART_header.hubble,   sizeof(float), 1, fPMcrd);
  fread(&ART_header.wp5,      sizeof(float), 1, fPMcrd);
  fread(&ART_header.Ocurv,    sizeof(float), 1, fPMcrd);
  for(i=0; i<100; i++)
    fread(&(ART_header.extras[i]), sizeof(float), 1, fPMcrd);
  
#ifdef BYTESWAP
  sexchange(&ART_header.aexpn,    sizeof(float));
  sexchange(&ART_header.aexp0,    sizeof(float));
  sexchange(&ART_header.amplt,    sizeof(float));
  sexchange(&ART_header.astep,    sizeof(float));
  sexchange(&ART_header.istep,    sizeof(int));
  sexchange(&ART_header.partw,    sizeof(float));
  sexchange(&ART_header.tintg,    sizeof(float));
  sexchange(&ART_header.ekin,     sizeof(float));
  sexchange(&ART_header.ekin1,    sizeof(float));
  sexchange(&ART_header.ekin2,    sizeof(float));
  sexchange(&ART_header.au0,      sizeof(float));
  sexchange(&ART_header.aeu0,     sizeof(float));
  sexchange(&ART_header.nrowc,    sizeof(int));  
  sexchange(&ART_header.ngridc,   sizeof(int));  
  sexchange(&ART_header.nspecies, sizeof(int));  
  sexchange(&ART_header.nseed,    sizeof(int));  
  sexchange(&ART_header.Om0,      sizeof(float));
  sexchange(&ART_header.Oml0,     sizeof(float));
  sexchange(&ART_header.hubble,   sizeof(float));
  sexchange(&ART_header.wp5,      sizeof(float));
  sexchange(&ART_header.Ocurv,    sizeof(float));
  for(i=0; i<100; i++)
    sexchange(&ART_header.extras[i], sizeof(float));
#endif
  
  /*==========================
   * dump some info to screen
   *==========================*/
  fprintf(stderr,"%s\n",ART_header.header_string);
  fprintf(stderr,"aexpn    =   %g\n",ART_header.aexpn);
  fprintf(stderr,"aexp0    =   %g\n",ART_header.aexp0);
  fprintf(stderr,"astep    =   %g\n",ART_header.astep);
  fprintf(stderr,"istep    =   %d\n",ART_header.istep);
  fprintf(stderr,"Om0      =   %g\n",ART_header.Om0);
  fprintf(stderr,"Oml0     =   %g\n",ART_header.Oml0);
  fprintf(stderr,"hubble   =   %g\n",ART_header.hubble);
  fprintf(stderr,"nspecies =   %d\n",ART_header.nspecies);
  //in older versions of the PM code the box size is not included in the file
  //ART_header.extras[99]=100;
  fprintf(stderr,"boxsize  =   %g\n",ART_header.extras[99]);
  fprintf(stderr,"pmass    =   %g\n",ART_header.extras[59-1]);
  fprintf(stderr,"NROW     =   %d\n",ART_header.nrowc);
  fprintf(stderr,"NGRID    =   %d\n",ART_header.ngridc);
  fprintf(stderr,"nspecies =   %d\n",ART_header.nspecies);
    
#ifndef MULTIMASS
  if(ART_header.nspecies > 1)
    {
      fprintf(stderr,"Please compile read_art with -DMULTIMASS\n");
      exit(0);
    }
#endif
  
  /*=====================
   * declare some things
   *=====================*/
  fprintf(stderr,"\nderived values:\n");
  NGRID    = ART_header.ngridc;
  NROW     = ART_header.nrowc;
  NPAGE    = pow2(NROW);
  NRECL    = 6*NPAGE;
  nspecies = ART_header.nspecies;
  partw    = ART_header.partw;
  boxsize  = ART_header.extras[100-1];
  if(fabs(ART_header.extras[59-1]) > MACHINE_ZERO)
    {
     pmass = ART_header.extras[59-1];
    }
  else
    {
     pmass = ART_header.Om0*rhoc0*pow3(ART_header.extras[99]/(double)ART_header.ngridc);
     fprintf(stderr,"ART file does not contain mass unit!\n");
     fprintf(stderr,"-> will use pmass = %g for initial conversion of ART masses to Msun/h \n",pmass);
    }
  
  /* conversion to Mpc/h, km/sec, Msun/h */
  x_fac    = 1.0/(double)NGRID * boxsize;
  v_fac    = 1.0/(double)NGRID * H0*boxsize / ART_header.aexpn;
  m_fac    = pmass;
  
  
  wspecies =       &(ART_header.extras[0]);
  lspecies = (int*)&(ART_header.extras[10]);
  
  if(nspecies == 0)
    {
      N_particles = pow3(NROW);
      Npages      = NROW;
      N_in_last   = NPAGE;
    }
  else
    {
      N_particles =  lspecies[nspecies-1];
      Npages      = (N_particles -1)/NPAGE +1;
      N_in_last   =  N_particles -NPAGE*(Npages-1);
    }
  
  fprintf(stderr,"no_part  =   %d\n\n", N_particles);
  
  fprintf(stderr,"x_fac    =   %g\n",x_fac);
  fprintf(stderr,"v_fac    =   %g\n",v_fac);
  fprintf(stderr,"m_fac    =   %g\n\n",m_fac);
  
  /*==================
   * read actual data
   *==================*/
  io.fst_part  = c_part(N_particles);
  io.no_part   = N_particles;
  io.fst_gas   = NULL;
  io.no_gas    = 0;
  io.fst_star  = NULL;
  io.no_stars  = 0;
  
  recdat       = (float *) malloc(NRECL*sizeof(float));
  
  cur_part     = io.fst_part;
  for(IROW=0; IROW<Npages; IROW++)
    {
      In_page = NPAGE;
      
      if(IROW == Npages-1)
        In_page = N_in_last;
      
      fprintf(stderr," o reading page %10d:  %10d   %10d\n",IROW,In_page,Npages);
      
      iL = NPAGE*IROW;
      fread(&recdat[0], sizeof(float), NRECL, fPMcrs);
      
#ifdef BYTESWAP
      for(i=0; i<NRECL; i++)
        sexchange(&recdat[i], sizeof(float));
#endif
      
      for(IN=0; IN<In_page; IN++)
        {
          xpar  = recdat[        IN];
          ypar  = recdat[  NPAGE+IN];
          zpar  = recdat[2*NPAGE+IN];
          vxpar = recdat[3*NPAGE+IN];
          vypar = recdat[4*NPAGE+IN];
          vzpar = recdat[5*NPAGE+IN];
          
#ifdef ART_DEBUG
          fprintf(stderr,"%d  %12.8g %12.8g %12.8g  %12.8g %12.8g %12.8g (%12.8g)\n",
                  IN,xpar,ypar,zpar,vxpar,vypar,vzpar,recdat[IN]);
#endif
          /* conversion to Mpc/h, km/sec, and Msun/h */
          cur_part->pos[X] = x_fac * (xpar-1.);
          cur_part->pos[Y] = x_fac * (ypar-1.);
          cur_part->pos[Z] = x_fac * (zpar-1.);
          cur_part->mom[X] = v_fac * (vxpar);
          cur_part->mom[Y] = v_fac * (vypar);
          cur_part->mom[Z] = v_fac * (vzpar);
          
          if(nspecies == 0)
            w = partw;
          else
            {
              Ipart = IN+iL;
              for(i=0; i<nspecies; i++)
                {
                  if(Ipart < lspecies[i])
                    {
                      w = wspecies[i];
                      break;
                    }
                }
            }
          
#ifdef MULTIMASS
          cur_part->weight = w * m_fac;
#endif
          
          Icount++;
          cur_part++;
        }
    }
  
  fclose(fPMcrs);
  
  free(recdat);
  
  
  /*===================
   * fill io.header
   *===================*/
  strcpy(io.header.header,ART_header.header_string);
  io.header.no_part     = N_particles;
  io.header.boxsize     = boxsize;
  io.header.omega0      = ART_header.Om0;
  io.header.lambda0     = ART_header.Oml0;
//  io.header.omegab      = 0.0;
//  io.header.gamma       = 0.0;
//  io.header.H_frac      = 0.0;  // set if needed at all
//  io.header.T_init      = 0.0;
//  io.header.B_init      = 0.0;
  io.header.pmass       = 0.0;    /* will be determined in ic_unit_conversion() */
  io.header.t_unit      = 0.0;    /* will be determined in ic_unit_conversion() */
  io.header.a_initial   = ART_header.aexp0;
  io.header.a_current   = ART_header.aexpn; 
  io.header.timestep    = 0.0;
  io.header.no_timestep = ART_header.istep;
  io.header.K_initial   = 0.0;
  io.header.U_initial   = 0.0;
  io.header.K_current   = 0.0;
  io.header.U_current   = 0.0;
  io.header.Eintegral   = 0.0;
  io.header.Econst      = 0.0;  

  /* will be re-calculated anyways... */
  io.header.no_species  = nspecies;
}

void read_art(FILE *fPMcrd)
{
 char   PMcrs[MAXSTRING];
 FILE  *fPMcrs;
 int    i, iM, slen;
  
 /* fPMcrd points to the PMcrd.DAT file */
 
 /* opening the PMcrs.DAT file */
 strcpy(PMcrs, io.icfile_name);
 slen = strlen(PMcrs);
 //fprintf(stderr,"%d\n",slen);
 for(i=slen-1; i>0; i--)
    if(PMcrs[i] == 'M')
      {
       iM = i;
       break;
      }
 if(iM==0)
   {
    fprintf(stderr,"could not locate the 'M' in %s\nEXIT\n",PMcrs);
    exit(0);
   }
       
 for(i=slen; i>=iM+4; i--)
    PMcrs[i] = PMcrs[i-1];
 
 PMcrs[iM+3]   = 's';
 PMcrs[iM+4]   = '0';
 PMcrs[slen+1] = '\0';
 
 
 if((fPMcrs = fopen(PMcrs,"rb")) == NULL)
   {
   fprintf(stderr,"read_ART: could not open %s\n",PMcrs);
   exit(0);
   }
 
  
 fprintf(stderr,"=============================================================================\n");
#ifdef BYTESWAP
 fprintf(stderr,"                      reading BYTESWAPed ART particles \n");
#else
 fprintf(stderr,"                           reading ART particles \n\n");
#endif

 read_ARTbinaries(fPMcrd, fPMcrs);
 
 fprintf(stderr,"\n                        finished reading ART particles\n");
 fprintf(stderr,"=============================================================================\n");
 
}


#endif /* ART */

#ifdef TIPSY
//========================================================================
//
//                              TIPSY
//
//========================================================================
void read_tipsy(FILE *icfile)
{
  partptr        cur_part;
  long unsigned  ipart;
  
  struct gas_particle  *gas_particles,  *gp;
  struct dark_particle *dark_particles, *dp;
  struct star_particle *star_particles, *sp;
  struct tipsy_dump     header_tipsy;
  
  long unsigned int  j, i, nx;
  int                no_timestep;
  double             box, omega0, lambda0, DMminmass;
  double             z_initial, z_current, a_initial, a_current;
  
  FILE *fout;
  
  int   SWAPBYTES;
#ifdef BYTESWAP
  SWAPBYTES = TRUE;
#else
  SWAPBYTES = FALSE;
#endif
  
  /**********************************************************************************
   * Reading in the Tipsy data  */
  
  fprintf(stderr,"\n=============================================================================\n");
#ifdef BYTESWAP
  fprintf(stderr,"               start reading BYTESWAPed TIPSY particles\n");
#else
  fprintf(stderr,"                       start reading TIPSY particles\n");
#endif
  
  
  if(fread(&header_tipsy,sizeof(header_tipsy),1,icfile) == 0)
    {
      fprintf(io.logfile,"\ncould not read TIPSY header\n");
      exit(0);
    }
  
  /* fix byte order */
#ifdef BYTESWAP
  sexchange(&header_tipsy.time,     sizeof(double));
  sexchange(&header_tipsy.nbodies,  sizeof(int));
  sexchange(&header_tipsy.ndim,     sizeof(int));
  sexchange(&header_tipsy.nsph,     sizeof(int));
  sexchange(&header_tipsy.ndark,    sizeof(int));
  sexchange(&header_tipsy.nstar,    sizeof(int));
#endif
  fprintf(stderr,"Read the TIPSY header...\n");
  fprintf(stderr,"time    = %lf\n",header_tipsy.time);
  fprintf(stderr,"nbodies = %ld\n",header_tipsy.nbodies);
  fprintf(stderr,"ndim    = %ld\n",header_tipsy.ndim);
  fprintf(stderr,"nsph    = %ld\n",header_tipsy.nsph);
  fprintf(stderr,"ndark   = %ld\n",header_tipsy.ndark);
  fprintf(stderr,"nstar   = %ld\n",header_tipsy.nstar);
  
  
  /********************************************************************************
   * Reading in the Tipsy data */
  if(header_tipsy.nsph != 0) 
    {
      gas_particles = (struct gas_particle *) malloc(header_tipsy.nsph*sizeof(*gas_particles));
      if(gas_particles == NULL) 
        {
          fprintf(io.logfile,"\n<sorry, no memory for gas particles, master>\n") ;
          exit(0);
        }
    }
  if(header_tipsy.ndark != 0) 
    {
      dark_particles = (struct dark_particle *) malloc(header_tipsy.ndark*sizeof(*dark_particles));
      if(dark_particles == NULL) 
        {
          fprintf(io.logfile,"\n<sorry, no memory for dark particles, master>\n") ;
          exit(0);
        }
    }
  if(header_tipsy.nstar != 0) 
    {
      star_particles = (struct star_particle *) malloc(header_tipsy.nstar*sizeof(*star_particles));
      if(star_particles == NULL) 
        {
          fprintf(io.logfile,"\n<sorry, no memory for star particles, master>\n") ;
          exit(0);
        }
    }
  
  
  /*******************************************************************************/
  fread((char *)gas_particles,sizeof(struct gas_particle),header_tipsy.nsph,icfile);
  /*******************************************************************************/
  
  /*********************************************************************************
   * Reading in the DM particles and swapping their bytes */
  DMminmass = 1E50;
  for(i=0; i<header_tipsy.ndark; i++) 
    {						
      fread(&dark_particles[i], sizeof(struct dark_particle), 1, icfile);
      
      /* Changing the sex of the particles */
#ifdef BYTESWAP
      sexchange(&dark_particles[i].mass, sizeof(float));
      for(j=0 ; j<3 ; j++) 
        {
          sexchange(&dark_particles[i].pos[j], sizeof(float));
          sexchange(&dark_particles[i].vel[j], sizeof(float));
        }
#endif
            
      if (dark_particles[i].mass < DMminmass) DMminmass = dark_particles[i].mass;
      
      if ( i == 0 )
        fprintf(stderr,"[1] %g %g %g   %g %g %g   %g\n",
                dark_particles[i].pos[0],dark_particles[i].pos[1],dark_particles[i].pos[2],
                dark_particles[i].vel[0],dark_particles[i].vel[1],dark_particles[i].vel[2],
                dark_particles[i].mass);
    }
  
  /*******************************************************************************/
  fread((char *)star_particles,sizeof(struct star_particle),header_tipsy.nstar,icfile);
  /*******************************************************************************/
  
  
  
  /************************************************************
   * filling in standard variables 
   */
  
  a_current = header_tipsy.time;
  z_current = fabs(((double)1.0/(a_current)) -(double)1.0);
  fprintf(stderr,"z_current = %f\n",z_current);
  
  /**********************************************************
   * Fill in io.header */
  strcpy(io.header.header,"TIPSY simulation");
  io.header.no_part        = header_tipsy.ndark;
  io.header.boxsize        = tipsy_boxsize;
  io.header.omega0         = tipsy_omega0;
  io.header.lambda0        = tipsy_lambda0;
  io.header.omegab         = 0.0;
  io.header.gamma          = 0.0;
  io.header.H_frac         = 0.0;
  io.header.T_init         = 0.0;
  io.header.B_init         = 0.0;
  io.header.pmass          = 0.0;       /* will be determined in ic_unit_conversion() */
  io.header.t_unit         = 0.0;       /* will be determined in ic_unit_conversion() */
  io.header.a_initial      = (double)1.0/(tipsy_initalz + (double)1.0);
  io.header.a_current	   = a_current;
  io.header.no_timestep    = tipsy_currentimeno;
  io.header.cur_reflevel   = 0;
  io.header.cur_frcres     = 0.0;
  io.header.K_initial      = 0.0;
  io.header.K_current      = 0.0;
  io.header.U_initial      = 0.0;
  io.header.U_current      = 0.0;
  io.header.Eintegral      = 0.0;
  io.header.Econst         = 0.0;
  
  /* This is unit conversion */
  /* 
   The positions should be in Mpc/h co-moving
   The velocities should be in km/s co-moving  (and hence the factor a_current!)
   The masses should be in Msun/h [already done above]
   */
   
  
  /******************************************************************************/
  /******************************************************************************/
  
  /* FIX units here [Justin] */
  io.no_part  = io.header.no_part;
  io.fst_part = c_part(io.no_part);
  io.no_stars = 0;
  io.fst_star = NULL;
  io.no_gas   = 0;
  io.fst_gas  = NULL;
  
  cur_part    = io.fst_part;
  for(ipart = 0; ipart < io.no_part; ipart++) {
    
    for(j = X; j <= Z; j++) 
      {	  
        /* TIPSY_LUNIT :: gets it to Mpc/h co-moving */
        cur_part->pos[j] = (flouble)(dark_particles[ipart].pos[j] * tipsy_boxsize);
        
        /* TIPSY_VUNIT :: gets it to km/s co-moving  */
        cur_part->mom[j] = (flouble)(dark_particles[ipart].vel[j] * TIPSY_VUNIT * a_current);
      }
    
#ifdef MULTIMASS
     /* TIPSY_MUNIT :: gets it to Msun/h */
     cur_part->weight    = (flouble)(dark_particles[ipart].mass   * TIPSY_MUNIT);     
#endif
    
    if ( ipart == 0 )
      fprintf(stderr,"[2] %g %g %g %g %g %g\n",cur_part->pos[0],cur_part->pos[1],cur_part->pos[2],cur_part->mom[0],cur_part->mom[1],cur_part->mom[2]);
    
    cur_part++;
  }
  
  /**********************************************************
   * Free memory */
  
  if(header_tipsy.nsph != 0) {
    free(gas_particles);
  }
  if(header_tipsy.ndark != 0) {
    free(dark_particles);
  }
  if(header_tipsy.nstar != 0) {
    free(star_particles);
  }
  
  
  fprintf(stderr,"\n                  finished reading TIPSY particles\n");
  fprintf(stderr,"=============================================================================\n");
}
#endif /* TIPSY */

#ifdef DEVA
//========================================================================
//
//                              DEVA
//
//========================================================================

//============================
// the actual reading routine
//============================
#ifdef DEVA2
void read_deva(FILE *icfile)
{
  // particles (nobj)
  devaflouble *rm, *xpos, *ypos, *zpos, *xvel, *yvel, *zvel, zmet; //metalicity ignored for the time being
  int         *itype;
  
  // gas+star properties (nbar=nstar+ngas)
  devaflouble *h, *e, *dn, *t_star, *t_chem;
  
  // variables used for reading
  int            idummy, *iibuf;
  devaflouble    ddummy, *dibuf;
  
  long      iloop, ipart, idark, igas, istar, *idx_hull;
  double    nvpart;
  partptr   cur_part;
  gasptr    cur_gas;
  
  // conversion factors
  double x_fac, v_fac, m_fac, e_fac;
  
  // the qhull parameter
  FILE  *fqhull;
  double center[3], Ab[6][4];
  int    nobj_hull, ndark_hull, ngas_hull, nstar_hull, nbar_hull;
  double dr[3];
  int    iuse, i, j;
  double rm_max;
  
#ifdef DEBUG_DEVA
  int   iseed = 123456;
  FILE *fpout;
#endif
    
#ifdef DEVA2_QHULL_FILE
  fprintf(stderr,"=> start reading DEVA file (using qhull: %s)\n", DEVA2_QHULL_FILE);
  
  //==============================
  // read qhull information
  //==============================
  fqhull = fopen(DEVA2_QHULL_FILE,"r");
  fscanf(fqhull,"%lf %lf %lf", &(center[0]), &(center[1]), &(center[2]));
  for(iloop=0; iloop<6; iloop++)
    fscanf(fqhull, "%lf %lf %lf %lf ", &(Ab[iloop][0]), &(Ab[iloop][1]), &(Ab[iloop][2]), &(Ab[iloop][3]));
  fclose(fqhull);
#endif /* DEVA2_QHULL_FILE */
  
  //==============================
  // read DEVA header information
  //==============================
  //first header block
  iibuf = (int *) &(ibuf1.itime);
  for(iloop=0; iloop<100; iloop++)
    {
      DEVA_SKIP;
      ReadInt(icfile, &idummy, SWAPBYTES);
      DEVA_SKIP;
      memcpy((iibuf+iloop), &idummy, 4);
    }
  //second header block
  dibuf = (devaflouble *) &(ibuf2.time);
  for(iloop=0; iloop<100; iloop++)
    {
      DEVA_SKIP;
      ReadFlouble(icfile, &ddummy, SWAPBYTES);
      DEVA_SKIP;
      memcpy((dibuf+iloop), &ddummy, DEVAFLOUBLE);
    }
  
  fprintf(stderr,"\n");
  fprintf(stderr,"DEVA header information:\n");
  fprintf(stderr,"------------------------\n");
  fprintf(stderr,"z        = %16.8f\n",1/ibuf2.atime-1.);
  fprintf(stderr,"atime    = %16.8f\n",ibuf2.atime);
  fprintf(stderr,"omega0   = %16.8f\n",ibuf2.omega0);
  fprintf(stderr,"xlambda0 = %16.8f\n",ibuf2.xlambda0);
  fprintf(stderr,"omegab0  = %16.8f\n",ibuf2.omegab0);
  fprintf(stderr,"sigma80  = %16.8f\n",ibuf2.sigma80);
  fprintf(stderr,"box100   = %16.8f\n",ibuf2.box100);
  fprintf(stderr,"h100     = %16.8f\n",ibuf2.h100);
  fprintf(stderr,"nobj     = %12d\n",  ibuf1.nobj);
  fprintf(stderr,"ndark    = %12d         rmdark   = %16.8f\n",ibuf1.ndark,ibuf2.rmdark);
  fprintf(stderr,"ngas     = %12d         rmgas    = %16.8f\n",ibuf1.ngas,ibuf2.rmgas);
  fprintf(stderr,"nstar    = %12d         rmtot    = %16.8f\n",ibuf1.nstar,ibuf2.rmtot);
  
  // conversion factors for positions and velocities
  x_fac = ibuf2.box100;
  v_fac = ibuf2.atime*H0*ibuf2.box100/(ibuf2.h0t0);
  e_fac = pow2(H0*ibuf2.box100/(ibuf2.h0t0));
  m_fac = 1.0; // will be determined below!!!
  
  //================
  // allocate memory
  //=================
  // ... for particles
  rm     = (devaflouble *) calloc(ibuf1.nobj, sizeof(devaflouble));
  xpos   = (devaflouble *) calloc(ibuf1.nobj, sizeof(devaflouble));
  ypos   = (devaflouble *) calloc(ibuf1.nobj, sizeof(devaflouble));
  zpos   = (devaflouble *) calloc(ibuf1.nobj, sizeof(devaflouble));
  xvel   = (devaflouble *) calloc(ibuf1.nobj, sizeof(devaflouble));
  yvel   = (devaflouble *) calloc(ibuf1.nobj, sizeof(devaflouble));
  zvel   = (devaflouble *) calloc(ibuf1.nobj, sizeof(devaflouble));
  itype  = (int *)    calloc(ibuf1.nobj, sizeof(int));
  
  // ... for gas + star properties
  h      = (devaflouble *) calloc(ibuf1.ngas+ibuf1.nstar, sizeof(devaflouble));
  e      = (devaflouble *) calloc(ibuf1.ngas+ibuf1.nstar, sizeof(devaflouble));
  dn     = (devaflouble *) calloc(ibuf1.ngas+ibuf1.nstar, sizeof(devaflouble));
  t_star = (devaflouble *) calloc(ibuf1.ngas+ibuf1.nstar, sizeof(devaflouble));
  t_chem = (devaflouble *) calloc(ibuf1.ngas+ibuf1.nstar, sizeof(devaflouble));
  
  //=======================================
  // read particles (and determine m_fac!)
  //=======================================
  //...masses
  nvpart = 0.0;
  rm_max = -1E10;
  DEVA_SKIP;
  for(iloop=0; iloop<ibuf1.nobj; iloop++)
    {
      ReadFlouble(icfile, &(rm[iloop]), SWAPBYTES);
      
      // total mass in internal units (needed to calculate m_fac conversion factor)
      nvpart += rm[iloop];
      
      // record maximum particle mass
      if(rm[iloop] > rm_max) rm_max = rm[iloop];
    }
  DEVA_SKIP;
  
  // conversion factor for masses
  m_fac = ibuf2.omega0*rhoc0*ibuf2.box100*ibuf2.box100*ibuf2.box100/nvpart;
  
  fprintf(stderr,"pmass    = %16.8g Msun/h\n",m_fac);
  fprintf(stderr,"rm_max   = %16.8g Msun/h\n",rm_max*m_fac);
  
#ifdef DEVA_L80
  rm_max *= 10;
#endif
  
  //...positions
  DEVA_SKIP;
  for(iloop=0; iloop<ibuf1.nobj; iloop++)
    {
      ReadFlouble(icfile, &(xpos[iloop]), SWAPBYTES);
      ReadFlouble(icfile, &(ypos[iloop]), SWAPBYTES);
      ReadFlouble(icfile, &(zpos[iloop]), SWAPBYTES);
    }
  DEVA_SKIP;
  
  //...velocities
  DEVA_SKIP;
  for(iloop=0; iloop<ibuf1.nobj; iloop++)
    {
      ReadFlouble(icfile, &(xvel[iloop]), SWAPBYTES);
      ReadFlouble(icfile, &(yvel[iloop]), SWAPBYTES);
      ReadFlouble(icfile, &(zvel[iloop]), SWAPBYTES);
    }
  DEVA_SKIP;
  
  //...itypes
  nobj_hull   = 0;
  ndark_hull  = 0;
  ngas_hull   = 0;
  nstar_hull  = 0;
  idx_hull    = (long unsigned *) calloc(1, sizeof(long unsigned));
  DEVA_SKIP;
  for(iloop=0; iloop<ibuf1.nobj; iloop++)
    {
      ReadInt(icfile, &(itype[iloop]), SWAPBYTES);
      
#ifdef DEVA2_QHULL_FILE
      // check whether particle lies inside hull based upon information in .qhull file
      dr[X] = fabs(xpos[iloop]-center[X]);
      dr[Y] = fabs(ypos[iloop]-center[Y]);
      dr[Z] = fabs(zpos[iloop]-center[Z]);
      if(dr[X] > 0.5) dr[X] = 1.0-dr[X];
      if(dr[Y] > 0.5) dr[Y] = 1.0-dr[Y];
      if(dr[Z] > 0.5) dr[Z] = 1.0-dr[Z];
      
      iuse = TRUE;
      for(i=0; i<5; i++)
        {
          for(j=0; j<3; j++)
            {
              if(Ab[i][j] * dr[j] > -Ab[i][3])
                iuse = FALSE;
            }
        }
      
      if(iuse)
        {
          // store particile's iloop-ID
          idx_hull            = (long unsigned *) realloc(idx_hull, (nobj_hull+1)*sizeof(long unsigned));
          idx_hull[nobj_hull] = iloop;
          
          // keep track of how many (and what type of) particles in the hull
          nobj_hull++;
          if(itype[iloop] == itdark)
            ndark_hull++;
          else if(itype[iloop] == itgas)
            ngas_hull++;
          else if(itype[iloop] == itstar)
            nstar_hull++;
        }
#else /* DEVA2_QHULL_FILE */
      // throw away all "tidal" dark matter particles
      iuse = TRUE;
      
      // 1. check: DM particle?
      if(itype[iloop] == itdark)
        {
          // 2. check: tidal particle?
          if(fabs(rm[iloop]-rm_max) < ZERO)
            {
              iuse = FALSE;
            }
        }
      
      if(iuse)
        {
          // store particile's iloop-ID
          idx_hull            = (long unsigned *) realloc(idx_hull, (nobj_hull+1)*sizeof(long unsigned));
          idx_hull[nobj_hull] = iloop;
          
          // keep track of how many (and what type of) particles in the hull
          nobj_hull++;
          if(itype[iloop] == itdark)
            ndark_hull++;
          else if(itype[iloop] == itgas)
            ngas_hull++;
          else if(itype[iloop] == itstar)
            nstar_hull++;
        }
#endif /* DEVA2_QHULL_FILE */
    }
  DEVA_SKIP;
  
  
#ifdef DEBUG_DEVA
  fpout = fopen("test.ascii","w");
  for(iloop=0; iloop<ibuf1.nobj; iloop++)
    if(ran3(&iseed) < 0.01)
      fprintf(fpout,"%16.8f %16.8f %16.8f    %16.8f %16.8f %16.8f    %16.8f    %12d\n",
              xpos[iloop], ypos[iloop], zpos[iloop],
              xvel[iloop], yvel[iloop], zvel[iloop],
              rm[iloop],itype[iloop]);
  fclose(fpout);
#endif
  
  // do not read gas properties for initial conditions
  if(!ibuf1.INITIALCOND)
    {
      //...softening lengths
      DEVA_SKIP;
      for(iloop=0; iloop<ibuf1.ngas+ibuf1.nstar; iloop++)
        ReadFlouble(icfile, &(h[iloop]), SWAPBYTES);
      DEVA_SKIP;
      
      //...thermal energies
      DEVA_SKIP;
      for(iloop=0; iloop<ibuf1.ngas+ibuf1.nstar; iloop++)
        ReadFlouble(icfile, &(e[iloop]), SWAPBYTES);
      DEVA_SKIP;
      
#ifdef DEVA_DETAILS
      //...densities
      DEVA_SKIP;
      for(iloop=0; iloop<ibuf1.ngas+ibuf1.nstar; iloop++)
        ReadFlouble(icfile, &(dn[iloop]), SWAPBYTES);
      DEVA_SKIP;
      
      //...star formation time
      DEVA_SKIP;
      for(iloop=0; iloop<ibuf1.ngas+ibuf1.nstar; iloop++)
        ReadFlouble(icfile, &(t_star[iloop]), SWAPBYTES);
      DEVA_SKIP;
      
      //...whatever
      DEVA_SKIP;
      for(iloop=0; iloop<ibuf1.ngas+ibuf1.nstar; iloop++)
        ReadFlouble(icfile, &(t_chem[iloop]), SWAPBYTES);
      DEVA_SKIP;
      
      //...metalicities (CURRENTLY IGNORED)
      DEVA_SKIP;
      for(iloop=0; iloop<nbar; iloop++)
        {
        }
      DEVA_SKIP;
#endif
    }
  
  fprintf(stderr,"\n");
  fprintf(stderr,"qhull information:\n");
  fprintf(stderr,"------------------\n");
  fprintf(stderr,"nobj_hull  = %12ld\n",nobj_hull);
  fprintf(stderr,"ndark_hull = %12ld\n",ndark_hull);
  fprintf(stderr,"ngas_hull  = %12ld\n",ngas_hull);
  fprintf(stderr,"nstar_hull = %12ld\n",nstar_hull);
  
  
  //=======================================
  // initialize AMIGA header
  //=======================================
  /*transfer everything from ASCII header to AMIGA header */
  strcpy(io.header.header,"DEVAsimulation");
  io.header.no_part     = nobj_hull;
  io.header.no_vpart    = nvpart;
  io.header.boxsize     = ibuf2.box100;
  io.header.omega0      = ibuf2.omega0;
  io.header.lambda0     = ibuf2.xlambda0;
  io.header.omegab      = ibuf2.omegab0;
  io.header.gamma       = 0.0;
  io.header.H_frac      = 0.0;
  io.header.T_init      = 0.0;
  io.header.B_init      = 0.0;
  io.header.a_initial   = ibuf2.atime;
  io.header.a_current   = ibuf2.atime;
  io.header.timestep    = 0.0;
  io.header.no_timestep = ibuf1.itime;
  io.header.cur_reflevel= 0.0;
  io.header.cur_frcres  = 0.0;
  io.header.K_initial   = 0.0;
  io.header.K_current   = 0.0;
  io.header.U_initial   = 0.0;
  io.header.U_current   = 0.0;
  io.header.Eintegral   = 0.0;
  io.header.Econst      = 0.0;
  io.header.pmass       = 0.0;    /* will be determined in ic_unit_conversion() */
  io.header.t_unit      = 0.0;    /* will be determined in ic_unit_conversion() */
  
  /* allocate memory for particles */
  io.no_part      = nobj_hull;
  io.fst_part     = c_part(io.no_part);
  
  io.no_gas       = ngas_hull;
  io.offset_gas   = ndark_hull;
  io.fst_gas      = c_gas(io.no_gas);
  
  io.no_stars     = nstar_hull;
  io.offset_stars = ndark_hull+ngas_hull;
  io.fst_star     = NULL; // stars are not treated separately in AHF
  
  
  /* eventually transfer and convert to physical units */
  idark = 0;
  igas  = 0;
  istar = 0;
  for(iloop=0; iloop<io.header.no_part; iloop++)
    {
      // select only the hull particles
      ipart = idx_hull[iloop];
      
      // new arrangement in RAM:  | DM | gas | stars |
      if(itype[ipart] == itdark)
        {
          cur_part = io.fst_part                   + idark;
          idark++;
        }
      else if(itype[ipart] == itgas)
        {
          cur_gas  = io.fst_gas  + igas;
          cur_part = io.fst_part + io.offset_gas   + igas;
          igas++;
          
          /* store energy in AMIGA units!!! (there is no conversion later on!) */
          cur_gas->u = e_fac * e[ipart];
        }
      else
        {
          cur_part = io.fst_part + io.offset_stars + istar;
          istar++;
        }
      
      /* this switches off the re-ordering */
      //cur_part = io.fst_part + iloop;
      
      cur_part->id     = iloop;
      cur_part->itype  = (long) itype[ipart];
      cur_part->pos[X] = x_fac * (flouble) xpos[ipart];
      cur_part->pos[Y] = x_fac * (flouble) ypos[ipart];
      cur_part->pos[Z] = x_fac * (flouble) zpos[ipart];
      cur_part->mom[X] = v_fac * (flouble) xvel[ipart];
      cur_part->mom[Y] = v_fac * (flouble) yvel[ipart];
      cur_part->mom[Z] = v_fac * (flouble) zvel[ipart];
      cur_part->weight = m_fac * (flouble) rm[ipart];      
    }
  
  // free all temporary arrays
  free(xpos);
  free(ypos);
  free(zpos);
  free(xvel);
  free(yvel);
  free(zvel);
  free(rm);
  free(itype);
  free(h);
  free(e);
  free(dn);
  free(t_star);
  free(t_chem);
  free(idx_hull);
  
  fprintf(stderr,"<= finished reading DEVA file\n");
  
}
#else /* DEVA2 */
void read_deva(FILE *icfile)
{
  // particles
  float *rm, *xpos, *ypos, *zpos, *xvel, *yvel, *zvel;
  int    *itype;
  
  // gas+star properties
  float *h, *e, *dn, *t_star, *f_g;
  
  // misc variables
  int       idummy, idummy1, idummy2;
  int      *iibuf, *iibuf1, *iibuf2;
  
  int       iseed;
  long      iloop, idark, igas, istar;
  double    nvpart;
  partptr   cur_part;
  gasptr    cur_gas;
  
  // conversion factors
  double x_fac, v_fac, m_fac, e_fac;
  
  // testing
  FILE *fpout;
  
  iseed = 123456;
  
  //==============================
  // read DEVA header information
  //==============================
  iibuf  = &(ibuf.time_out[0]);
  iibuf1 = &(ibuf1.itime);
  iibuf2 = &(ibuf2.irun);
  for(iloop=0; iloop<100; iloop++)
    {
      // just read 3x4bytes into some temporary memory
      DEVA_SKIP;
      ReadInt(icfile, &idummy,  SWAPBYTES);
      ReadInt(icfile, &idummy1, SWAPBYTES);
      ReadInt(icfile, &idummy2, SWAPBYTES);
      DEVA_SKIP;
      
      // copy those 3x4bytes over to where they belong
      memcpy((iibuf +iloop), &idummy,  4);
      memcpy((iibuf1+iloop), &idummy1, 4);
      memcpy((iibuf2+iloop), &idummy2, 4);
    }

  fprintf(stderr,"\nDEVA header information:\n");
  fprintf(stderr,"------------------------\n");
  fprintf(stderr,"z        = %16.8f\n",1/ibuf1.atime-1.);
  fprintf(stderr,"atime    = %16.8f\n",ibuf1.atime);
  fprintf(stderr,"omega0   = %16.8f\n",ibuf2.omega0);
  fprintf(stderr,"xlambda0 = %16.8f\n",ibuf2.xlambda0);
  fprintf(stderr,"box100   = %16.8f\n",ibuf2.box100);
  fprintf(stderr,"h100     = %16.8f\n",ibuf2.h100);
  fprintf(stderr,"nobj     = %12d\n",  ibuf2.nobj);
  fprintf(stderr,"ndark    = %12d         rmdark   = %16.8f\n",ibuf2.ndark,ibuf2.rmdark);
  fprintf(stderr,"ngas     = %12d         rmgas    = %16.8f\n",ibuf2.ngas,ibuf2.rmgas);
  fprintf(stderr,"nstar    = %12d         rmbary   = %16.8f\n",ibuf2.nstar,ibuf2.rmbary);

  // conversion factors for positions and velocities
  x_fac = ibuf2.box100;
  v_fac = ibuf1.atime*H0*ibuf2.box100/(ibuf2.h0t0);
  e_fac = pow2(H0*ibuf2.box100/(ibuf2.h0t0));
  m_fac = 1.0; // will be determined below!!!
  
  //================
  // allocate memory
  //=================
  // ... for particles
  rm    = (float *) calloc(ibuf2.nobj, sizeof(float));
  xpos  = (float *) calloc(ibuf2.nobj, sizeof(float));
  ypos  = (float *) calloc(ibuf2.nobj, sizeof(float));
  zpos  = (float *) calloc(ibuf2.nobj, sizeof(float));
  xvel  = (float *) calloc(ibuf2.nobj, sizeof(float));
  yvel  = (float *) calloc(ibuf2.nobj, sizeof(float));
  zvel  = (float *) calloc(ibuf2.nobj, sizeof(float));
  itype = (int   *) calloc(ibuf2.nobj, sizeof(int));
  
  // ... for gas + star properties
  h     = (float *) calloc(ibuf2.nobj-ibuf2.ndark, sizeof(float));
  e     = (float *) calloc(ibuf2.nobj-ibuf2.ndark, sizeof(float));
  dn    = (float *) calloc(ibuf2.nobj-ibuf2.ndark, sizeof(float));
  t_star= (float *) calloc(ibuf2.nobj-ibuf2.ndark, sizeof(float));
  f_g   = (float *) calloc(ibuf2.nobj-ibuf2.ndark, sizeof(float));
  
  //=======================================
  // read particles (and determine m_fac!)
  //=======================================
  nvpart = 0.0;
  for(iloop=0; iloop<ibuf2.nobj; iloop++)
    {
      DEVA_SKIP;
      ReadFloat(icfile, &(rm[iloop]),    SWAPBYTES);
      ReadFloat(icfile, &(xpos[iloop]),  SWAPBYTES);
      ReadFloat(icfile, &(ypos[iloop]),  SWAPBYTES);
      ReadFloat(icfile, &(zpos[iloop]),  SWAPBYTES);
      ReadFloat(icfile, &(xvel[iloop]),  SWAPBYTES);
      ReadFloat(icfile, &(yvel[iloop]),  SWAPBYTES);
      ReadFloat(icfile, &(zvel[iloop]),  SWAPBYTES);
      ReadInt  (icfile, &(itype[iloop]), SWAPBYTES);
      DEVA_SKIP;
      
      // total mass in internal units (needed to calculate m_fac conversion factor)
      nvpart += rm[iloop];
      
//      fprintf(stderr,"%g %g %g %g %g %g %g %d\n",rm[iloop],xpos[iloop],ypos[iloop],zpos[iloop],xvel[iloop],yvel[iloop],zvel[iloop],itype[iloop]);
//      exit(0);
    }
  
//  fprintf(stderr,"%g %g %g %g\n",ibuf2.omega0,rhoc0,ibuf2.box100,nvpart);
  
  // conversion factor for masses
  m_fac = ibuf2.omega0*rhoc0*ibuf2.box100*ibuf2.box100*ibuf2.box100/nvpart;
  
  fprintf(stderr,"pmass    = %16.8g Msun/h\n",m_fac);
  
  // read gas+star properties
  for(iloop=0; iloop<ibuf2.nobj-ibuf2.ndark; iloop++)
    {
      DEVA_SKIP;
      ReadFloat(icfile, &(h[iloop]),      SWAPBYTES);
      ReadFloat(icfile, &(e[iloop]),      SWAPBYTES);
      ReadFloat(icfile, &(dn[iloop]),     SWAPBYTES);
      ReadFloat(icfile, &(t_star[iloop]), SWAPBYTES);
      ReadFloat(icfile, &(f_g[iloop]),    SWAPBYTES);
      DEVA_SKIP;
    }
  
#ifdef DEBUG_DEVA
  fpout = fopen("test.ascii","w");
  for(iloop=0; iloop<ibuf2.nobj; iloop++)
    if(ran3(&iseed) < 0.1)
      fprintf(fpout,"%16.8f %16.8f %16.8f    %16.8f %16.8f %16.8f    %16.8f    %12d\n",
              xpos[iloop], ypos[iloop], zpos[iloop],
              xvel[iloop], yvel[iloop], zvel[iloop],
              rm[iloop],itype[iloop]);
  fclose(fpout);
#endif
  
  //==========================
  // transfer to AMIGA arrays
  //==========================
  /*transfer everything from ASCII header to AMIGA header */
  strcpy(io.header.header,"DEVAsimulation");
  io.header.no_part     = ibuf2.nobj;
  io.header.no_vpart    = nvpart;
  io.header.boxsize     = ibuf2.box100;
  io.header.omega0      = ibuf2.omega0;
  io.header.lambda0     = ibuf2.xlambda0;
  io.header.omegab      = ibuf2.rmbary*ibuf2.omega0;
  io.header.gamma       = 0.0;
  io.header.H_frac      = 0.0;
  io.header.T_init      = 0.0;
  io.header.B_init      = 0.0;
  io.header.a_initial   = ibuf1.atime;
  io.header.a_current   = ibuf1.atime;
  io.header.timestep    = 0.0;
  io.header.no_timestep = ibuf1.itime;
  io.header.cur_reflevel= 0.0;
  io.header.cur_frcres  = 0.0;
  io.header.K_initial   = 0.0;
  io.header.K_current   = 0.0;
  io.header.U_initial   = 0.0;
  io.header.U_current   = 0.0;
  io.header.Eintegral   = 0.0;
  io.header.Econst      = 0.0;
  io.header.pmass       = 0.0;    /* will be determined in ic_unit_conversion() */
  io.header.t_unit      = 0.0;    /* will be determined in ic_unit_conversion() */
  
  /* allocate memory for particles */
  io.no_part      = io.header.no_part;
  io.fst_part     = c_part(io.no_part);
  
  io.no_gas       = ibuf2.ngas;
  io.offset_gas   = ibuf2.ndark;
  io.fst_gas      = c_gas(io.no_gas);
  
  io.no_stars     = ibuf2.nstar;
  io.offset_stars = ibuf2.ndark+ibuf2.ngas;
  io.fst_star     = NULL; // stars are not treated separately in AHF

  /* eventually transfer and convert to physical units */
  idark = 0;
  igas  = 0;
  istar = 0;
  for(iloop=0; iloop<io.header.no_part; iloop++)
    {
      // new arrangement in RAM:  | DM | gas | stars |
      if(itype[iloop] == itdark)
        {
          cur_part = io.fst_part                   + idark;
          idark++;
        }
      else if(itype[iloop] == itgas || itype[iloop] == itgas2)
        {
          cur_gas  = io.fst_gas  + igas;
          cur_part = io.fst_part + io.offset_gas   + igas;
          igas++;
          
          /* store energy */
          cur_gas->u = e_fac * e[iloop];
        }
      else
        {
          cur_part = io.fst_part + io.offset_stars + istar;
          istar++;
        }
      cur_part->id     = iloop;
      cur_part->itype  = (long) itype[iloop];
      cur_part->pos[X] = x_fac * (flouble) xpos[iloop];
      cur_part->pos[Y] = x_fac * (flouble) ypos[iloop];
      cur_part->pos[Z] = x_fac * (flouble) zpos[iloop];
      cur_part->mom[X] = v_fac * (flouble) xvel[iloop];
      cur_part->mom[Y] = v_fac * (flouble) yvel[iloop];
      cur_part->mom[Z] = v_fac * (flouble) zvel[iloop];
      cur_part->weight = m_fac * (flouble) rm[iloop];
    }
  
  // free all temporary arrays
  free(xpos);
  free(ypos);
  free(zpos);
  free(xvel);
  free(yvel);
  free(zvel);
  free(rm);
  free(itype);
  free(h);
  free(e);
  free(dn);
  free(t_star);
  free(f_g);
}
#endif /* DEVA2 */
#endif /* DEVA */


/*==============================================================================================
* input data from IC/restart file, i.e. simply makes the appropriate call to read_data()
* -> this is a wrapper to allow for multiple GADGET files...
*=============================================================================================*/
void input(char *infile_name)
{
   FILE *icfile;
   fprintf(stderr,"\n*******************************************************************************\n");
   fprintf(stderr,"                               starting input()\n");
   fprintf(stderr,"*******************************************************************************\n\n");
   
#ifdef GADGET
   char    gadget_file[MAXSTRING];
   int     no_gadget_files, i_gadget_file;

   /*=================================
    * there are multiple GADGET files
    *=================================*/
   if((icfile = fopen(infile_name,"rb")) == NULL)
     {
      /* maybe there are multiple GADGET files ... count them! */
      no_gadget_files = 0;
      i_gadget_file   = 0;
      sprintf(gadget_file,"%s.%d",infile_name,i_gadget_file);
      while((icfile = fopen(gadget_file,"rb")) != NULL)
        {
         no_gadget_files++;
         i_gadget_file++;
         sprintf(gadget_file,"%s.%d",infile_name,i_gadget_file);
        }
      
      if(no_gadget_files > 1)
        {
          fprintf(stderr,"\nreading GADGET data from %d files\n",no_gadget_files);
          gadget.no_gadget_files  = no_gadget_files;

          /* allocate temporary storage for no. of particles arrays */
          gadget.np[0]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
          gadget.np[1]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
          gadget.np[2]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
          gadget.np[3]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
          gadget.np[4]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
          gadget.np[5]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
          
          /* open muli-GADGET files one after the other to determine mmin, IDmin, IDmax */
         fprintf(stderr,"  -> determining mmin and IDmin:\n ");
         gadget.IDmin = 1000000000;
         gadget.IDmax = 0;
         gadget.mmin  = 1e30;
         for(i_gadget_file=0; i_gadget_file<no_gadget_files; i_gadget_file++)
           {
            sprintf(gadget_file,"%s.%d",infile_name,i_gadget_file);
            icfile               = fopen(gadget_file,"rb");
            gadget.i_gadget_file = i_gadget_file;
            
            fprintf(stderr,"\n reading file #%d:\n",i_gadget_file+1);
            skim_gadget(icfile);
            
            fclose(icfile);
           }
#ifdef LGADGET
         fprintf(stderr,"\n     mmin = %g    IDmin = %ld   IDmax = %ld\n",gadget.mmin,gadget.IDmin, gadget.IDmax);
#else
         fprintf(stderr,"\n     mmin = %g    IDmin = %d   IDmax = %d\n",gadget.mmin,gadget.IDmin, gadget.IDmax);
#endif
         if(gadget.IDmax-gadget.IDmin+1 != gadget.nall1)
           {
            fprintf(stderr,"\nPROBLEM:  halo partice ID's are not consecutive from IDmin to IDmax (%d)\n",gadget.nall1);
            fprintf(stderr,"            -> you appear not be analysing a decent cosmological simulation?!\n\n");
            fprintf(stderr,"            -> the authors of AHF will be happy to resolve this issue with you, please get in touch...\n\n");
            exit(0);
           }
         
         /* read multi-GADGET files one by one */
         for(i_gadget_file=0; i_gadget_file<no_gadget_files; i_gadget_file++)
           {
            sprintf(gadget_file,"%s.%d",infile_name,i_gadget_file);
            fprintf(stderr,"\n => reading %s\n",gadget_file);
            icfile = fopen(gadget_file,"rb");
            
            /* tell read_gadget() which file we are using at the moment */
            gadget.i_gadget_file = i_gadget_file;
            
            /* read files... */
            if(i_gadget_file < no_gadget_files-1)
              /* simply call read_gadget() for the first (no_gadget_files-1) files */
               read_gadget(icfile);
            else
              /* call read_data() for the last file toensure ic_unit_conversion() etc. */
               read_data(icfile);
            fclose(icfile);
           } 
          
          /* free temporary storage again */
          free(gadget.np[0]);
          free(gadget.np[1]);
          free(gadget.np[2]);
          free(gadget.np[3]);
          free(gadget.np[4]);
          free(gadget.np[5]);
        }
      else
        {
         /* there are no multi-GADGET files */
         fprintf(stderr,"\n\ninput: could not open file with IC's  %s\n",infile_name);
         exit(1);
        }
     }
   else
     {
      /*===============================
      * there is only one GADGET file
      *===============================*/
      gadget.no_gadget_files  = 1;
      gadget.i_gadget_file    = 0;
      
       /* allocate temporary storage for no. of particles arrays */
      gadget.np[0]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
      gadget.np[1]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
      gadget.np[2]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
      gadget.np[3]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
      gadget.np[4]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
      gadget.np[5]     = (int *) calloc(gadget.no_gadget_files, sizeof(int *));
      
      read_data(icfile);
      fclose(icfile);

       /* remove temporary storage again */
       free(gadget.np[0]);
       free(gadget.np[1]);
       free(gadget.np[2]);
       free(gadget.np[3]);
       free(gadget.np[4]);
       free(gadget.np[5]);
     }
      
#else /* GADGET */
   /*===================================================================
      * THIS IS THE STANDARD CALL TO READ A SINGLE FILE VIA read_data()
      *       (what format will be decided within read_data()!)
      *=================================================================*/
   if((icfile = fopen(infile_name,"rb")) == NULL)
     {
      fprintf(stderr,"\n\ninput: could not open file with IC's  %s\n",infile_name);
      exit(1);
     }  

   read_data(icfile);
   
   fclose(icfile);
   
#endif /* GADGET */
   
   fprintf(stderr,"\n*******************************************************************************\n");
   fprintf(stderr,"                  finished input() -> no_timestep = %10d\n",io.header.no_timestep);
   fprintf(stderr,"*******************************************************************************\n\n");
}

/*==============================================================================================
 * read data from IC/restart file, i.e. make appropriate call to read_XYZ()
 *
 *
 * the only global structure that will be initialized by read_data() is
 * --------------------------------------------------------------------
 * io.
 *
 * io.fst_part   -> will give access to all particles
 * io.header.    -> structure that contains all relevant information about simultion
 *
 * => this is sufficient information to initialize all other global structures 
 *
 *
 *
 * read_data() itself further calculates 
 *    no_species
 *    min_weight
 *    max_weight
 *    no_vpart 
 *    -> these values will be checked against the values found in io.header!
 *
 *
 *=============================================================================================*/
void read_data(FILE *icfile)
{  
#ifdef DEVA
  read_deva(icfile);
#else
#ifdef MLAPM
  read_mlapm(icfile);
#else /* MLAPM not defined ... but maybe GADGET ?! */
#ifdef GADGET
  read_gadget(icfile);
#else /* GADGET not defined ... but maybe TIPSY ?! */
#ifdef TIPSY
  read_tipsy(icfile);
#else /* TIPSY not defined ... but maybe ART ?! */
#ifdef ART
  read_art(icfile);
#else /* ART not defined ... but maybe MARE_NOSTRUM ?! */
#ifdef MARE_NOSTRUM
  read_mare_nostrum(icfile);
#else /* MARE_NOSTRUM not defined ... but maybe ASCII ?! */
#ifdef ASCII
  read_ascii(icfile);
#else /* ASCII not defined ... hence read AMIGA*/
  read_amiga(icfile);
#endif /* ASCII */
#endif /* MARE_NOSTRUM */
#endif /* ART */
#endif /* TIPSY */
#endif /* GADGET */
#endif /* MLAPM */
#endif /* DEVA */
  
  
  /*======================================================================
   * ic_unit_conversion:
   * -------------------
   * 
   * 1. determine...
   *
   *    - length unit      -> boxsize
   *
   *    - mass unit        -> most common particle mass
   *
   *    - time unit        -> NO_EXPANSION: t_unit given by 4 PI G = 1
   *                             EXPANSION: t_unit = 1/H0
   *
   * 2. convert input data to internal units
   *
   *  
   *
   * => AMIGA expects...
   *
   *    pos[]  = Mpc/h 
   *    mom[]  = km/sec
   *    weight = Msun/h
   *
   *    ...prior to unit conversion
   *
   *
   *
   *
   * => this routine initializes:
   *
   *    io.header.pmass
   *    io.header.t_unit
   *    io.header.no_vpart
   *    io.header.min_weight
   *    io.header.max_weight
   *    io.header.no_species
   *
   *=======================================================================*/
#if (defined ART || defined GADGET || defined TIPSY || defined ASCII || defined MARE_NOSTRUM || defined MLAPM || defined DEVA)
  ic_unit_conversion();
#endif

 
  /*===================================================================
   * fill remaining io.header values (other than parameters!)
   *===================================================================*/  
  io.header.version           = VERSION;
  io.header.build             = BUILD;
  
#ifdef MULTIMASS
  io.header.multi_mass        = 1;
#else
  io.header.multi_mass        = 0;
#endif
#ifdef DOUBLE
  io.header.double_precision  = 1;
#else
  io.header.double_precision  = 0;
#endif
#ifdef HYDRO
  io.header.hydro             = 1;
#else
  io.header.hydro             = 0;
#endif
#ifdef MHD
  io.header.magneto           = 1;
#else
  io.header.magneto           = 0;
#endif
  
    
  
  fprintf(stderr,"\n\nsummary of io.header.\n");
  fprintf(stderr,"=====================\n");
  fprintf(stderr,"%s\n",                io.header.header);
  fprintf(stderr,"  no_part     = %ld\n", io.no_part);
  fprintf(stderr,"  boxsize     = %f\n",  io.header.boxsize);
  fprintf(stderr,"  omega0      = %f\n",  io.header.omega0);
  fprintf(stderr,"  omegab      = %f\n",  io.header.omegab);
  fprintf(stderr,"  lambda0     = %f\n",  io.header.lambda0);
  fprintf(stderr,"  gamma       = %f\n",  io.header.gamma);
  fprintf(stderr,"  H_frac      = %f\n",  io.header.H_frac);
  fprintf(stderr,"  T_init      = %f\n",  io.header.T_init);
  fprintf(stderr,"  z_initial   = %f\n",  1.0/io.header.a_initial - 1.0);
  fprintf(stderr,"  z_current   = %f\n",  1.0/io.header.a_current - 1.0);
  fprintf(stderr,"  timestep    = %f\n",  io.header.timestep);
  fprintf(stderr,"  no_timestep = %d\n",  io.header.no_timestep);
  fprintf(stderr,"  pmass       = %f\n",  io.header.pmass);
  fprintf(stderr,"  t_unit      = %f\n",  io.header.t_unit);
  fprintf(stderr,"  no_vpart    = %g\n",  io.header.no_vpart);
  fprintf(stderr,"  no_gas      = %ld\n", io.no_gas);
  fprintf(stderr,"  no_stars    = %ld\n", io.no_stars);
  fprintf(stderr,"  no_species  = %ld\n", io.header.no_species);
  fprintf(stderr,"  min_weight  = %g\n",  io.header.min_weight);
  fprintf(stderr,"  max_weight  = %g\n",  io.header.max_weight);
  fprintf(stderr,"  med_weight  = %g\n",  io.header.med_weight);
  fprintf(stderr,"  multi_mass  = %d\n",  io.header.multi_mass);
  fprintf(stderr,"  double_prcsn= %d\n",  io.header.double_precision);
  fprintf(stderr,"  hydro       = %d\n",  io.header.hydro);
  fprintf(stderr,"  magneto     = %d\n",  io.header.magneto);
  fprintf(stderr,"  version     = %3.1f\n",io.header.version);
  fprintf(stderr,"  build       = %d\n",  io.header.build);


  return;
}

