#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "utility.h"
#include "../libio_serial/io_serial.h"

/* relevant parameters for c2f_slope() */
#define C2F_MIDCENTRE_NODELIMITER
//#define C2F_MIDCENTRE
//#define C2F_REDUCED

/* L1DIM_LENGTH = no. of char's for l1dim in filename (cf. write_filename) */
#define L1DIM_LENGTH 6


/*======================================================================
 *  structure holdig information about timing on each grid level
 * (will be calculated from information stored in each grid structure!)
 *======================================================================*/
struct  total_timing
{
  time_t potential;             /* deriving the potenital by GS sweeps      */
  time_t density;               /* deriving the proper densities            */
  time_t DK;                    /* drifting and kicking particles           */
  time_t grid;                  /* everything related to grid hierarchy     */
  time_t hydro;                 /* time spent for hydro-solver              */
  time_t ahf;
}grid_timing;


char *terminate_amiga = TERMINATE_AMIGA;

/* prepare filename for debugging purposes */
void write_filename(char *f_name, char *prefix, unsigned l1dim)
{
   int  slen, i;                  /* string length                  */
   char file_no[10];              /* file number                    */
   
   /* prepare filename */
   f_name = strcpy(f_name, prefix);
   
   itoa_(l1dim, file_no);
   slen = strlen(file_no);
   if (slen < L1DIM_LENGTH)
     {
      file_no[L1DIM_LENGTH] = '\0';
      
      for(i=L1DIM_LENGTH-1; i >= (L1DIM_LENGTH-slen); i--)
        {
         file_no[i] = file_no[i-(L1DIM_LENGTH-slen)];
        }
      for(i=0; i < (L1DIM_LENGTH-slen); i++)
        {
         file_no[i] = '0';
        }
     }
   f_name = strcat(f_name, file_no);
   f_name = strcat(f_name, ".DAT");
   
}

/*===========================================================================
 * calculate Laplace operator acting on potential at a given node 
 *===========================================================================*/
double Laplace_pot(nptr tsc_nodes[3][3][3], double spacing2)
{
#ifndef AHFlean
  double Lpot;
  
  Lpot = (  (double)tsc_nodes[1][1][2]->pot + (double)tsc_nodes[1][1][0]->pot
          + (double)tsc_nodes[1][2][1]->pot + (double)tsc_nodes[1][0][1]->pot
          + (double)tsc_nodes[2][1][1]->pot + (double)tsc_nodes[0][1][1]->pot
          - 6. * (double)tsc_nodes[1][1][1]->pot) / spacing2;
  
  return Lpot;
#else
   fprintf(stderr, "%s should not be called when AHFlean is defined.\n",
           __func__);
   exit(1);
   return 0.0;
#endif /* AHFlean*/
}

/*===========================================================================
 * calculate Laplace operator acting on temp[1] at a given node 
 *===========================================================================*/
double Laplace_temp1(nptr tsc_nodes[3][3][3], double spacing2)
{
#ifndef AHFlean
  double Ltemp1;
  
  Ltemp1 = (  (double)tsc_nodes[1][1][2]->force.temp[1] + (double)tsc_nodes[1][1][0]->force.temp[1]
            + (double)tsc_nodes[1][2][1]->force.temp[1] + (double)tsc_nodes[1][0][1]->force.temp[1]
            + (double)tsc_nodes[2][1][1]->force.temp[1] + (double)tsc_nodes[0][1][1]->force.temp[1]
            - 6. * (double)tsc_nodes[1][1][1]->force.temp[1]) / spacing2;
  
  return Ltemp1;
#else
   fprintf(stderr, "%s should not be called when AHFlean is defined.\n",
           __func__);
   exit(1);
   return 0.0;
#endif /* AHFlean*/
}


/* f1mod:  x modulo y for double numbers */
double f1mod(double x, double y)
{
  /* the following part is fully tuned for use with AMIGA !!!!! */
   if(x >= 2.0)
      return(x-2.0);
   else if(x >= 1.0)
      return(x-1.0);
   else
      return(x);   
}


/*==============================================================================
*  get eigenvalues of inertia tensor
*==============================================================================*/
void get_axes(double itensor[3][3], double *axis1, double *axis2, double *axis3)
{
   int           n, i, j, nrot;
   unsigned long idx[4];
   double        a[4][4], d[4], v[4][4], tmp[4];
   
   n = NDIM;
   
   for(i=0; i<4; i++)
      for(j=0; j<4; j++)
         a[i][j] = 0.0;
   
   a[1][1] = itensor[0][0];
   a[2][1] = itensor[1][0];
   a[3][1] = itensor[2][0];
   a[1][2] = itensor[0][1];
   a[2][2] = itensor[1][1];
   a[3][2] = itensor[2][1];
   a[1][3] = itensor[0][2];
   a[2][3] = itensor[1][2];
   a[3][3] = itensor[2][2];
   
   jacobi(a, n, d, v, &nrot);
   
   for(i=1; i<=n; i++)
      tmp[i] = d[i];
   indexx((unsigned long)n, tmp, idx);
   
   *axis1 = d[idx[3]];
   *axis2 = d[idx[2]];
   *axis3 = d[idx[1]];

   itensor[0][0] = v[1][idx[3]];
   itensor[1][0] = v[2][idx[3]];
   itensor[2][0] = v[3][idx[3]];
   
   itensor[0][1] = v[1][idx[2]];
   itensor[1][1] = v[2][idx[2]];
   itensor[2][1] = v[3][idx[2]];
   
   itensor[0][2] = v[1][idx[1]];
   itensor[1][2] = v[2][idx[1]];
   itensor[2][2] = v[3][idx[1]];
}

/*==============================================================================
*  inverts the usage of idx[] returned by indexx sorting routine
*==============================================================================*/
int idx_inv(unsigned long *idx, int numHalos, int i)
{
  int j;
  
  if(i<0)
    return(-1);
  
  for(j=0; j<numHalos; j++)
   {
    if(idx[j] == i)
      break;
   }
  
  return(j);
}

/*
 ************************************************************
 ************************************************************
 *  Calculates the minima and the maxima
 */
MINMAX MinMax(double x,double xmin,double xmax) {
   
   MINMAX tmpMinMax;		
   
   if (x < xmin)
      tmpMinMax.min = x;
   else
      tmpMinMax.min = xmin;
   
   if (x > xmax)
      tmpMinMax.max = x;
   else
      tmpMinMax.max = xmax;
   
   
   return(tmpMinMax);
   
}
/*
 ************************************************************
 ************************************************************
 *  Calculates the minima and the maxima for the periodic refinements
 */
MINMAX MinMaxBound(double div, double x,double xmin,double xmax) {
   
   MINMAX tmpMinMax;		
   
   
   if ( x < div ) {
      
      if (x > xmax)
         tmpMinMax.max = x;
      else
         tmpMinMax.max = xmax;
      
      tmpMinMax.min = xmin;
      
   } else {
      
      if (x < xmin)
         tmpMinMax.min = x;
      else
         tmpMinMax.min = xmin;
      
      tmpMinMax.max = xmax;
						
   }
   
   return(tmpMinMax);
   
}

/*==============================================================================
 * initialize binning for density profiles
 *==============================================================================*/
void binning_parameter(HALO halo, int *nbins, double *dist_min, double *dist_max)
{
  partptr cur_part;
  double  dX, dY, dZ, Xc, Yc, Zc;
  long    npart, binpart;
  
  /* how many bins should be used for cur_halo */
  *nbins   = (int) (6.2*(log10((double)halo.npart))-3.5);
  *nbins  *= AHF_NBIN_MULTIPLIER;
  *nbins   = MAX(2,*nbins);
  binpart  = (double)(halo.npart-1)/(double)*nbins;
  
  /* halo position in AMIGA units */
  Xc = halo.pos.x;
  Yc = halo.pos.y;
  Zc = halo.pos.z;
  
  /* minimum distance: we are not using the innermost particle but move out to 1/10 of AHF_MINPART */
//  npart = 0;
//  while(npart < halo.npart && npart < (int)((double)simu.AHF_MINPART/10.))
//   {
//    cur_part = global.fst_part + halo.ipart[npart];
//    npart++;
//   }
//  dX = fabs(cur_part->pos[X] - Xc);
//  dY = fabs(cur_part->pos[Y] - Yc);
//  dZ = fabs(cur_part->pos[Z] - Zc);
//  if(dX > 0.5) dX = 1.0-dX;
//  if(dY > 0.5) dY = 1.0-dY;
//  if(dZ > 0.5) dZ = 1.0-dZ;
//  *dist_min  = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
  
  
  /* minimum distance: search for first particle with dist != 0 */
  npart     = (long)floor(((double)simu.AHF_MINPART/10.)+0.5);
  *dist_min = -1.0;
  while(npart < halo.npart-1 && *dist_min < MACHINE_ZERO)
   {
    cur_part = global.fst_part + halo.ipart[npart];
    dX = fabs(cur_part->pos[X] - Xc);
    dY = fabs(cur_part->pos[Y] - Yc);
    dZ = fabs(cur_part->pos[Z] - Zc);
    if(dX > 0.5) dX = 1.0-dX;
    if(dY > 0.5) dY = 1.0-dY;
    if(dZ > 0.5) dZ = 1.0-dZ;
    *dist_min  = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
    npart++;
   }
  
  /* maximum distance */
  cur_part = global.fst_part + halo.ipart[halo.npart-1];
  dX = fabs(cur_part->pos[X] - Xc);
  dY = fabs(cur_part->pos[Y] - Yc);
  dZ = fabs(cur_part->pos[Z] - Zc);
  if(dX > 0.5) dX = 1.0-dX;
  if(dY > 0.5) dY = 1.0-dY;
  if(dZ > 0.5) dZ = 1.0-dZ;
  *dist_max  = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
  
  /* try to capture the strange instance that half of the particles are right in the centre */
  if(*dist_min < MACHINE_ZERO) {
    *dist_min = *dist_max/2.;
  }
}


/*==============================================================================
 * read AMIGA header 
 *==============================================================================*/
void read_amiga_header(FILE *infile, info_io *io, int *SWAPBYTES)
{
   int       file_sizeof_long, machine_sizeof_long, one;
   int       idummy;
   long      ldummy;
   long long lldummy;
      
#ifdef AMIGA_ONE_FORMAT
   fread(&one, sizeof(int), 1, infile);
   if(one != 1)
     {
      *SWAPBYTES = TRUE;
      fprintf(stderr," => start reading BYTESWAPed AMIGA header ('one' format assuming sizeof(long)=4!) ... ");
     }
   else
     {
      *SWAPBYTES = FALSE;
      fprintf(stderr," => start reading AMIGA header ('one' format assuming sizeof(long)=4!) ... ");
     }

   /* don't bother with this check... */
   machine_sizeof_long = file_sizeof_long = 4;
#else
   /* do we need to do BYTESWAP (simultaneously determine the size of a long...) */
   fread(&(file_sizeof_long), sizeof(int), 1, infile);
   if(file_sizeof_long > 8)
     {
      *SWAPBYTES = TRUE;
      sexchange(&file_sizeof_long, sizeof(int));
      fprintf(stderr," => start reading BYTESWAPed AMIGA header ... ");
     }
   else
     {
      *SWAPBYTES = FALSE;
      fprintf(stderr," => start reading AMIGA header ... ");
     }
   machine_sizeof_long = sizeof(long);
   
   if(machine_sizeof_long != file_sizeof_long)
      fprintf(stderr,"(sizeof(long) mismatch: %d vs. %d!) ",machine_sizeof_long,file_sizeof_long);
#endif
   
   if(*SWAPBYTES == TRUE || machine_sizeof_long != file_sizeof_long)
     {
      /* read in IO header */
      ReadChars(infile,io->header.header,HEADERSTRING);
      ReadInt(infile,&(io->header.multi_mass),*SWAPBYTES);
      ReadInt(infile,&(io->header.double_precision),*SWAPBYTES);
      
      /* read 2x long unsigned */
      if(machine_sizeof_long == 8 && file_sizeof_long == 4)
        {
         ReadInt(infile,&idummy,*SWAPBYTES);
         io->header.no_part    = (long unsigned) idummy;
         ReadInt(infile,&idummy,*SWAPBYTES);
         io->header.no_species = (long unsigned) idummy;
        }
      else if(machine_sizeof_long == 4 && file_sizeof_long == 8)
        {
         ReadLongLong(infile,&lldummy,*SWAPBYTES);
         io->header.no_part    = (long unsigned) lldummy;
         ReadLongLong(infile,&lldummy,*SWAPBYTES);
         io->header.no_species = (long unsigned) lldummy;
        }
      else if((machine_sizeof_long == 8 && file_sizeof_long == 8) ||
              (machine_sizeof_long == 4 && file_sizeof_long == 4))
        {
         ReadLong(infile,&ldummy,*SWAPBYTES);
         io->header.no_part    = (long unsigned) ldummy;
         ReadLong(infile,&ldummy,*SWAPBYTES);
         io->header.no_species = (long unsigned) ldummy;
        }
      else
        {
         fprintf(stderr,"read_amiga_header: machine_sizeof_long=%d vs. file_sizeof_long=%d\n",
                 machine_sizeof_long, file_sizeof_long);
         exit(0);
        }
      
      ReadDouble(infile,&(io->header.no_vpart),*SWAPBYTES);
      ReadDouble(infile,&(io->header.timestep),*SWAPBYTES);
      ReadInt(infile,&(io->header.no_timestep),*SWAPBYTES);
      if(machine_sizeof_long == 4 && file_sizeof_long == 8)
         ReadInt(infile,&idummy,*SWAPBYTES);
      ReadDouble(infile,&(io->header.boxsize),*SWAPBYTES);
      ReadDouble(infile,&(io->header.omega0),*SWAPBYTES);
      ReadDouble(infile,&(io->header.lambda0),*SWAPBYTES);
      ReadDouble(infile,&(io->header.pmass),*SWAPBYTES);
      ReadDouble(infile,&(io->header.cur_reflevel),*SWAPBYTES);
      ReadDouble(infile,&(io->header.cur_frcres),*SWAPBYTES);
      ReadDouble(infile,&(io->header.a_initial),*SWAPBYTES);
      ReadDouble(infile,&(io->header.a_current),*SWAPBYTES);
      ReadDouble(infile,&(io->header.K_initial),*SWAPBYTES);
      ReadDouble(infile,&(io->header.K_current),*SWAPBYTES);
      ReadDouble(infile,&(io->header.U_initial),*SWAPBYTES);
      ReadDouble(infile,&(io->header.U_current),*SWAPBYTES);
      ReadDouble(infile,&(io->header.Eintegral),*SWAPBYTES);
      ReadDouble(infile,&(io->header.Econst),*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.paramNSTEPS),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramNGRID_DOM),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramNth_dom),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramNth_ref),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramE_UPDATE),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramCELLFRAC_MAX),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramCELLFRAC_MIN),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramCA_CRIT),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramMAX_L1DIM),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramDOMSWEEPS),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramREFSWEEPS),*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.paramAHF_MINPART),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramAHF_VTUNE),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramAHF_RISE),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramAHF_SLOPE),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramAHF_MAXNRISE),*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.min_weight),*SWAPBYTES);
      ReadDouble(infile,&(io->header.max_weight),*SWAPBYTES);
      ReadDouble(infile,&(io->header.t_unit),*SWAPBYTES);
      ReadDouble(infile,&(io->header.B_init),*SWAPBYTES);
      ReadDouble(infile,&(io->header.param_dummy5),*SWAPBYTES);
      ReadDouble(infile,&(io->header.param_dummy6),*SWAPBYTES);
      ReadDouble(infile,&(io->header.param_dummy7),*SWAPBYTES);
      ReadDouble(infile,&(io->header.param_dummy8),*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.version),*SWAPBYTES);
      ReadInt(infile,&(io->header.build),*SWAPBYTES);
      if(machine_sizeof_long == 4 && file_sizeof_long == 8)
         ReadInt(infile,&idummy,*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.omegab),*SWAPBYTES);
      ReadDouble(infile,&(io->header.gamma), *SWAPBYTES);
      ReadDouble(infile,&(io->header.H_frac),*SWAPBYTES);
      ReadDouble(infile,&(io->header.T_init),*SWAPBYTES);

      ReadInt(infile,&(io->header.hydro),*SWAPBYTES);
      if(machine_sizeof_long == 4 && file_sizeof_long == 8)
         ReadInt(infile,&idummy,*SWAPBYTES);
      ReadInt(infile,&(io->header.magneto),*SWAPBYTES);
      if(machine_sizeof_long == 4 && file_sizeof_long == 8)
         ReadInt(infile,&idummy,*SWAPBYTES);

      /* take care of this "long issue" when reading the FILLHEADER stuff */
      if(machine_sizeof_long == 8 && file_sizeof_long == 4)
        {
         /* we read 2x4 bytes fewer than FILLHEADER... */
         ReadChars(infile,io->header.dummy,FILLHEADER +4+4);
        }
      else if(machine_sizeof_long == 4 && file_sizeof_long == 8)
        {
         /* we read 2x4 bytes more than FILLHEADER... */
         ReadChars(infile,io->header.dummy,FILLHEADER -4-4);
        }
      else if((machine_sizeof_long == 8 && file_sizeof_long == 8) ||
              (machine_sizeof_long == 4 && file_sizeof_long == 4))
        {
         ReadChars(infile,io->header.dummy,FILLHEADER);
        }
     } 
   else /* SWAPBYTES */
     {
      /* read header as complete structure */
      if(fread(&(io->header), sizeof(io->header), 1, infile) != 1)
        {
         fprintf(stderr,"\n\ninput: could not read AMIGA IO header\n");
         fflush(stderr);
         fclose(stderr);
         exit(1);
        }
     } /* SWAPBYTES */
  fprintf(stderr,"done <=\n");
}


/*===================================================================================
* get_c2fslope:   calculate the slope f'(x_i) needed with c2f_pot
*
*             f(x)    = f(x_i) + f'(x_i) * (x-x_i)
*
*===================================================================================*/
void get_c2fslope(double func[3][3][3], double slope[NDIM])
{
   double slope_left, slope_right;
   
   /* not that the spacing of the sampling of func is 1 */
   
#ifdef C2F_MIDCENTRE_NODELIMITER
   
   /* mid-centred slope (no delimiter!) */
   slope[X] = (func[1][1][2] - func[1][1][0]) / 2.;
   slope[Y] = (func[1][2][1] - func[1][0][1]) / 2.;
   slope[Z] = (func[2][1][1] - func[0][1][1]) / 2.;
   
#endif /* C2F_MIDCENTRE_NODELIMITER */
   
#ifdef C2F_MIDCENTRE
   
   /* mide-centred slope */
   slope[X] = (func[1][1][2] - func[1][1][0]) / 2.;
   slope[Y] = (func[1][2][1] - func[1][0][1]) / 2.;
   slope[Z] = (func[2][1][1] - func[0][1][1]) / 2.;
   
   /* delimit the slope */
   
   slope_left  = func[1][1][1] - func[1][1][0];
   slope_right = func[1][1][2] - func[1][1][1];
   
   if(slope_left*slope_right < 0.)
      slope[X] = 0.0;
   
   slope_left  = func[1][1][1] - func[1][0][1];
   slope_right = func[1][2][1] - func[1][1][1];
   
   if(slope_left*slope_right < 0.)
      slope[Y] = 0.0;
   
   slope_left  = func[1][1][1] - func[0][1][1];
   slope_right = func[2][1][1] - func[1][1][1];
   
   if(slope_left*slope_right < 0.)
      slope[Z] = 0.0;
   
#endif /* C2F_MIDCENTRE */
   
#ifdef C2F_REDUCED
   
   /* delimit the slope */
   
   slope[X] = 0.0;
   slope[Y] = 0.0;
   slope[Z] = 0.0;
   
   slope_left  = func[1][1][1] - func[1][1][0];
   slope_right = func[1][1][2] - func[1][1][1];
   
   if(slope_left*slope_right > MACHINE_ZERO)
      slope[X] = 2*slope_left*slope_right / (slope_left+slope_right);
   
   slope_left  = func[1][1][1] - func[1][0][1];
   slope_right = func[1][2][1] - func[1][1][1];
   
   if(slope_left*slope_right > MACHINE_ZERO)
      slope[Y] = 2*slope_left*slope_right / (slope_left+slope_right);
   
   slope_left  = func[1][1][1] - func[0][1][1];
   slope_right = func[2][1][1] - func[1][1][1];
   
   if(slope_left*slope_right > MACHINE_ZERO)
      slope[Z] = 2*slope_left*slope_right / (slope_left+slope_right);
   
#endif /* C2F_REDUCED */
   
}
  

/*========================================================================
 * determine...
 *  no_species              no. of different particle species
 *  no_vpart                total mass in box
 *  min_weight              min. particle mass
 *  max_weight              max. particle mass
 *  most_common_weight      most common particle mass
 *
 * -> initialize...
 *  io.header.no_species
 *  io.header.no_vpart
 *  io.header.min_weight
 *  io.header.max_weight
 *  io.header.med_weight
 *  ...and return the most common particle weight
 *========================================================================*/
double init_header_masses()
{
   double  *wspecies, cur_weight;
   long    *nspecies, nmax, imax, ipart;
   int      no_species;
   int      ispecies;
   boolean  old_species;
   double   no_vpart;
   double   min_weight;
   double   max_weight;
   double   K, T, C; // for Kahan summation
   partptr  cur_part;
   double   most_common_weight, med_weight, med_norm;
   double   npart, inv_mean_dens, pmass;
   
   fprintf(stderr,"     init_header_masses():\n");

#ifdef MULTIMASS
   no_species             = 1;
   wspecies               = (double*)calloc(no_species, sizeof(double));
   nspecies               = (long*)  calloc(no_species, sizeof(long));
   
   /* we determine particle weights in units of first particle mass */
   pmass                  = io.fst_part->weight;
   min_weight             = 1.0;
   max_weight             = 1.0;
   no_vpart               = 1.0;
   wspecies[no_species-1] = 1.0;
   nspecies[no_species-1] = 1;
   
   /* for Kahan summation of no_vpart */
   C = 0.0;
   
   for(cur_part=io.fst_part+1; cur_part<io.fst_part+io.no_part; cur_part++)
     {
      /* compare particle weight in units of first particle mass */
      cur_weight  = cur_part->weight/pmass;
             
      /* no_vpart += cur_weight ala Kahan summation */
      K = cur_weight-C;
      T = no_vpart + K;
      C = (T-no_vpart) - K;
      no_vpart = T;

      old_species = FALSE;
      for(ispecies=0; ispecies<no_species; ispecies++)
        {
         if(fabs(cur_weight - wspecies[ispecies])/wspecies[ispecies] < 0.1)
           {
            old_species         = TRUE;
            nspecies[ispecies] += 1;
           }
        }
      
      if(old_species == FALSE)
        {
         no_species++;
          
         if(cur_weight > max_weight) max_weight = cur_weight;
         if(cur_weight < min_weight) min_weight = cur_weight;
         
         wspecies               = (double*) realloc(wspecies, no_species*sizeof(double));
         nspecies               = (long*)   realloc(nspecies, no_species*sizeof(long));

         wspecies[no_species-1] = cur_weight;
         nspecies[no_species-1] = 1;
        }
     }
   
   fprintf(stderr,"     -> the following species have been found...\n");
   med_weight = 0.0;
   med_norm = 0.0;
   for(ispecies=0; ispecies<no_species; ispecies++)
     {
       fprintf(stderr,"       %12d      %16.8g [Msun/h]     %16ld\n",ispecies,wspecies[ispecies]*pmass,nspecies[ispecies]);
       med_weight += nspecies[ispecies]*wspecies[ispecies];
       med_norm   += (double)nspecies[ispecies];
     }
   med_weight /= med_norm;
   fprintf(stderr,"\n");
   
   
   /* determine most common particle weight */
   imax = 0;
   nmax = nspecies[0];
   for(ispecies=1; ispecies<no_species; ispecies++)
     {
      if(nspecies[ispecies] > nmax)
        {
         nmax = nspecies[ispecies];
         imax = ispecies;
        }
     }
   most_common_weight = wspecies[imax];

   free(wspecies); 
   free(nspecies);    
   
#else /* MULTIMASS */
   
#ifdef NO_EXPANSION
   fprintf(stderr,"     YOU ARE USING A NON-COSMOLOGICAL SETTING WITHOUT PROVIDING PARTICLE MASSES\n");
   fprintf(stderr,"      -> not implemented yet...exiting\n");
   exit(0);
#else /* NO_EXPANSION */
   /* use cosmological background density to determine pmass */
   npart                  = (double)io.no_part;
   inv_mean_dens          = pow3(io.header.boxsize)/npart;
   pmass                  = io.header.omega0*rhoc0*inv_mean_dens;
   
   /* these values are in internal units */
   no_species             = 1;
   no_vpart               = npart;
   most_common_weight     = 1.0;
   min_weight             = 1.0;
   max_weight             = 1.0;
   med_weight             = 1.0;
   most_common_weight     = 1.0;

   fprintf(stderr,"     you are using a cosmological setting without providing particle masses\n");
   fprintf(stderr,"      -> will use pmass = %16.8g as particle mass\n\n",pmass);
#endif /* NO_EXPANSION */
   
#endif /* MULTIMASS */
   
   
   
   
   
   /* dump no_vpart, no_species, min_weight, and max_weight */
   fprintf(stderr,"     no_species          = %16d\n",             no_species);
   fprintf(stderr,"     pmass               = %16.8g  Msun/h (%16.8g)\n",               pmass, io.header.pmass);
   fprintf(stderr,"     no_vpart            = %16.8g  Msun/h\n",   no_vpart           * pmass);
   fprintf(stderr,"     most_common_weight  = %16.8g  Msun/h\n\n", most_common_weight * pmass);
   fprintf(stderr,"     min_weight          = %16.8g  Msun/h\n",   min_weight         * pmass);
   fprintf(stderr,"     max_weight          = %16.8g  Msun/h\n",   max_weight         * pmass);
   fprintf(stderr,"     med_weight          = %16.8g  Msun/h\n\n", med_weight         * pmass);
   
   /* store newly calculated values */
   io.header.no_species = no_species;
   io.header.no_vpart   = no_vpart   * pmass;
   io.header.min_weight = min_weight * pmass;    // these values are in Msun/h
   io.header.max_weight = max_weight * pmass;
  
   // TODO: what is the right choice for med_weight anyways?!
   io.header.med_weight = med_weight * pmass; //most_common_weight * pmass;

   most_common_weight   = most_common_weight * pmass;    // this will be used as the mass unit in AMIGA

  return (most_common_weight);
}

/*========================================================================
 * convert pos[], mom[], and weight to internal units
 *========================================================================*/
void ic_unit_conversion()
{
   partptr       cur_part;
   long          ipart;
   double        x_fac, v_fac, m_fac;
   double        xmin, xmax, ymin, ymax, zmin, zmax;
   double        rho_mean;
   
   fprintf(stderr,"\n=================================================================\n");
   fprintf(stderr,"                    ic_unit_conversion()\n");
   fprintf(stderr,"=================================================================\n");

   /*=====================================================================
    * 1. CHECK POSITIONS RANGE
    *=====================================================================*/
   fprintf(stderr," 1. position range consistency check:\n");
   fprintf(stderr," ------------------------------------\n");
   xmax = -1E40; ymax = -1E40; zmax = -1E40;
   xmin = +1E40; ymin = +1E40; zmin = +1E40;
   
   /* there is no MIN,MAX reduction for OpenMP in C :-( */
   for(cur_part=io.fst_part; cur_part<io.fst_part+io.no_part; cur_part++)
     {
      if(cur_part->pos[X] > xmax) xmax = cur_part->pos[X];
      if(cur_part->pos[Y] > ymax) ymax = cur_part->pos[Y];
      if(cur_part->pos[Z] > zmax) zmax = cur_part->pos[Z];
      
      if(cur_part->pos[X] < xmin) xmin = cur_part->pos[X];
      if(cur_part->pos[Y] < ymin) ymin = cur_part->pos[Y];
      if(cur_part->pos[Z] < zmin) zmin = cur_part->pos[Z];
     }
   
   if(xmin >= 0.0               && ymin >= 0.0               && zmin >= 0.0 && 
      xmax <= io.header.boxsize && ymax <= io.header.boxsize && zmax <= io.header.boxsize)
     {
      fprintf(stderr,"   -> everything's fine!\n");
     }
   
   /* negative coordinates? */
   if(xmin < 0.0 || ymin < 0.0 || zmin < 0.0)
     {
      fprintf(stderr,"   -> negative coordinates:\n");
      fprintf(stderr,"      boxsize = %16.8g\n", io.header.boxsize);
      fprintf(stderr,"      min     = %16.8g    %16.8g    %16.8g\n",xmin,ymin,zmin);
      fprintf(stderr,"      max     = %16.8g    %16.8g    %16.8g\n",xmax,ymax,zmax);

      /* is it possible to simply shift the positions? */
      if(xmax-xmin <= io.header.boxsize || ymax-ymin <= io.header.boxsize || zmax-zmin <= io.header.boxsize)
        {
         fprintf(stderr,"     -> shifting coordinates to [0,%g]\n",io.header.boxsize);
#ifdef WITH_OPENMP
#pragma omp parallel private(cur_part) shared(xmin, xmax, ymin, ymax, zmin, zmax, io)
#pragma omp for schedule(static)
#endif
         for(ipart=0; ipart<io.no_part; ipart++)
           {
            cur_part         = io.fst_part + ipart;
            cur_part->pos[X] = fmod(cur_part->pos[X]+fabs(xmin), io.header.boxsize);
            cur_part->pos[Y] = fmod(cur_part->pos[Y]+fabs(ymin), io.header.boxsize);
            cur_part->pos[Z] = fmod(cur_part->pos[Z]+fabs(zmin), io.header.boxsize);
           }
        }
      else
        {
         fprintf(stderr,"     -> and coordinates exceeding allowed range:\n");
         fprintf(stderr,"        boxsize = %16.8g\n", io.header.boxsize);
         fprintf(stderr,"        min     = %16.8g    %16.8g    %16.8g\n",xmin,ymin,zmin);
         fprintf(stderr,"        max     = %16.8g    %16.8g    %16.8g\n",xmax,ymax,zmax);

         fprintf(stderr,"        not fixed yet...exiting\n");
         exit(0);
        }
     }
   
   /* coordinates exceeding boxsize limit? */
   if(xmax > io.header.boxsize || ymax > io.header.boxsize || zmax > io.header.boxsize)
     {
      fprintf(stderr,"   -> coordinates exceeding allowed range:\n");
      fprintf(stderr,"      boxsize = %16.8g\n", io.header.boxsize);
      fprintf(stderr,"      min     = %16.8g    %16.8g    %16.8g\n",xmin,ymin,zmin);
      fprintf(stderr,"      max     = %16.8g    %16.8g    %16.8g\n",xmax,ymax,zmax);
      
      fprintf(stderr,"      not fixed yet...exiting\n");
      exit(0);
     }
   
   
   
   
   
   /*==================================================================
    * 2. DETERMINE MOST APPROPRIATE MASS UNIT
    *==================================================================*/
   fprintf(stderr,"\n 2. determining most appropriate mass unit:\n");
   fprintf(stderr," ------------------------------------------\n");
   io.header.pmass = init_header_masses();
#ifdef GADGET_PMASS
   io.header.pmass = GADGET_MUNIT;
   fprintf(stderr,"   -> will use GADGET mass unit %g Msun/h as mass unit for AMIGA\n",io.header.pmass);
#endif
#ifdef MANUAL_PMASS
   io.header.pmass = 42570.70209705;   // set the mass unit to whatever you fancy...
#endif
   
   /* AMIGA expects these values to be in internal units! */
   io.header.no_vpart   /= io.header.pmass;
   io.header.min_weight /= io.header.pmass;
   io.header.max_weight /= io.header.pmass;
   io.header.med_weight /= io.header.pmass;
   
   fprintf(stderr,"   -> will use %g Msun/h as mass unit for AMIGA\n",io.header.pmass);

   
   
   /*==================================================================
    * 3. CHOOSE TIME UNIT
    *==================================================================*/
   fprintf(stderr,"\n 3. setting the time unit:\n");
   fprintf(stderr," -------------------------\n");
#ifdef NO_EXPANSION
   rho_mean         = io.header.no_vpart*io.header.pmass/pow3(io.header.boxsize);
   io.header.t_unit = 1./(4.*PI*Grav*rho_mean);
   io.header.t_unit = sqrt(io.header.t_unit);
   
   fprintf(stderr,"   -> non-cosmological setup, will use %16.8g as time unit for AMIGA\n",io.header.t_unit);
#else
   io.header.t_unit = 1./H0;

   fprintf(stderr,"   -> cosmological setup, will use 1/H0 = %g as time unit for AMIGA\n",io.header.t_unit);
#endif
   
   
   
   
   
   /*==================================================================
    *                        UNIT CONVERSION
    *==================================================================*/
   x_fac = 1./io.header.boxsize;
   v_fac = io.header.a_current*io.header.t_unit/io.header.boxsize;
   m_fac = 1./io.header.pmass;
   
#ifdef WITH_OPENMP
#pragma omp parallel private(cur_part) shared(x_fac, v_fac, m_fac, io)
#pragma omp for schedule(static)
#endif
   for(ipart=0; ipart<io.no_part; ipart++)
     {
      cur_part          = io.fst_part + ipart;
      cur_part->pos[X]  = f1mod(cur_part->pos[X]*x_fac + 1.0, 1.0);
      cur_part->pos[Y]  = f1mod(cur_part->pos[Y]*x_fac + 1.0, 1.0);
      cur_part->pos[Z]  = f1mod(cur_part->pos[Z]*x_fac + 1.0, 1.0);
      cur_part->mom[X] *= v_fac;
      cur_part->mom[Y] *= v_fac;
      cur_part->mom[Z] *= v_fac;
#ifdef MULTIMASS
      cur_part->weight *= m_fac;
#endif
     }
   fprintf(stderr,"=================================================================\n");
   fprintf(stderr,"                finished ic_unit_conversion()\n");
   fprintf(stderr,"=================================================================\n");
   
}

/*========================================================================
 * here we simply check whether the DEFINFLAGS in combination with the
 * io.header parameters are actually meaningful...
 *========================================================================*/
void sanity_check()
{
   double pmass, rho_mean;
   
   fprintf(stderr," => performing sanity check of header information:\n");
   
   /* dump AMIGA header (*before* reading particles) */
   fprintf(stderr,"  io.header (as read from input file):\n");
   fprintf(stderr,"  ------------------------------------\n");
   fprintf(stderr,"   %s\n",io.header.header);
   fprintf(stderr,"   header.multi_mass            = %d\n",io.header.multi_mass);
   fprintf(stderr,"   header.double_precision      = %d\n",io.header.double_precision);
   fprintf(stderr,"   header.hydro                 = %d\n",io.header.hydro);
   fprintf(stderr,"   header.magneto               = %d\n",io.header.magneto);
   fprintf(stderr,"   header.no_part               = %ld\n",io.header.no_part);
   fprintf(stderr,"   header.no_species            = %ld\n",io.header.no_species);
   fprintf(stderr,"   header.min_weight            = %g\n",io.header.min_weight);
   fprintf(stderr,"   header.max_weight            = %g\n",io.header.max_weight);
   fprintf(stderr,"   header.med_weight            = %g\n",io.header.med_weight);
   fprintf(stderr,"   header.no_vpart              = %g\n",io.header.no_vpart);
   fprintf(stderr,"   header.timestep              = %g\n",io.header.timestep);
   fprintf(stderr,"   header.no_timestep           = %d\n",io.header.no_timestep);
   fprintf(stderr,"   header.t_unit                = %g\n",io.header.t_unit);
   fprintf(stderr,"   header.pmass                 = %g\n",io.header.pmass);
   fprintf(stderr,"   header.boxsize               = %g\n",io.header.boxsize);
   fprintf(stderr,"   header.omega0                = %g\n",io.header.omega0);
   fprintf(stderr,"   header.omegab                = %g\n",io.header.omegab);
   fprintf(stderr,"   header.lambda0               = %g\n",io.header.lambda0);   
   fprintf(stderr,"   header.gamma                 = %g\n",io.header.gamma);
   fprintf(stderr,"   header.H_frac                = %g\n",io.header.H_frac);
   fprintf(stderr,"   header.T_init                = %g\n",io.header.T_init);
   fprintf(stderr,"   header.cur_reflevel          = %g\n",io.header.cur_reflevel);
   fprintf(stderr,"   header.cur_frcres            = %g\n",io.header.cur_frcres);
   fprintf(stderr,"   header.a_initial             = %g\n",io.header.a_initial);
   fprintf(stderr,"   header.a_current             = %g\n",io.header.a_current);
   fprintf(stderr,"   header.K_initial             = %g\n",io.header.K_initial);
   fprintf(stderr,"   header.K_current             = %g\n",io.header.K_current);
   fprintf(stderr,"   header.U_initial             = %g\n",io.header.U_initial);
   fprintf(stderr,"   header.U_current             = %g\n",io.header.U_current);
   fprintf(stderr,"   header.Eintegral             = %g\n",io.header.Eintegral);
   fprintf(stderr,"   header.Econst                = %g\n",io.header.Econst);
   fprintf(stderr,"   header.paramNSTEPS           = %g\n",io.header.paramNSTEPS);
   fprintf(stderr,"   header.paramNGRID_DOM        = %g\n",io.header.paramNGRID_DOM);   
   fprintf(stderr,"   header.paramNth_dom          = %g\n",io.header.paramNth_dom);
   fprintf(stderr,"   header.paramNth_ref          = %g\n",io.header.paramNth_ref);
   fprintf(stderr,"   header.paramE_UPDATE         = %g\n",io.header.paramE_UPDATE);
   fprintf(stderr,"   header.paramCELLFRAC_MAX     = %g\n",io.header.paramCELLFRAC_MAX);
   fprintf(stderr,"   header.paramCELLFRAC_MIN     = %g\n",io.header.paramCELLFRAC_MIN);
   fprintf(stderr,"   header.paramCA_CRIT          = %g\n",io.header.paramCA_CRIT);
   fprintf(stderr,"   header.paramMAX_L1DIM        = %g\n",io.header.paramMAX_L1DIM);
   fprintf(stderr,"   header.paramDOMSWEEPS        = %g\n",io.header.paramDOMSWEEPS);
   fprintf(stderr,"   header.paramREFSWEEPS        = %g\n",io.header.paramREFSWEEPS);
   fprintf(stderr,"   header.paramAHF_MINPART      = %g\n",io.header.paramAHF_MINPART);
   fprintf(stderr,"   header.paramAHF_VTUNE        = %g\n",io.header.paramAHF_VTUNE);
   fprintf(stderr,"   header.paramAHF_RISE         = %g\n",io.header.paramAHF_RISE);
   fprintf(stderr,"   header.paramAHF_SLOPE        = %g\n",io.header.paramAHF_SLOPE);
   fprintf(stderr,"   header.paramAHF_MAXNRISE     = %g\n",io.header.paramAHF_MAXNRISE);
   
   /*=========
    *  UNITS
    *=========*/
   /* check mass unit */
   if(fabs(io.header.pmass) < ZERO)
     {
      fprintf(stderr," o mass unit not set  ");
      
#ifdef NO_EXPANSION
      fprintf(stderr,"-> non-cosmological setup ... exiting!\n");
      exit(0);
#else
      /* we brute force use the cosmological setting! */
      io.header.pmass = io.header.omega0*rhoc0*pow3(io.header.boxsize)/io.header.no_vpart;
      fprintf(stderr,"-> using cosmological value pmass = %g\n",io.header.pmass);
#endif
     }
   
   /* check time unit */
   if(fabs(io.header.t_unit) < ZERO)
     {
      fprintf(stderr," o time unit not set  ");
      
#ifdef NO_EXPANSION
      rho_mean         = io.header.no_vpart*io.header.pmass/pow3(io.header.boxsize);
      io.header.t_unit = 1./(4.*PI*Grav*rho_mean);
      io.header.t_unit = sqrt(io.header.t_unit);
      fprintf(stderr,"-> using non-cosmological value t_unit = %g\n",io.header.t_unit);
#else
      /* we brute force use the cosmological setting! */
      io.header.t_unit = 1./H0;   
      fprintf(stderr,"-> using cosmological value 1/H0  = %g\n",io.header.t_unit);
#endif
     }
   

   /*====================
    *  io.header masses
    *====================*/
   if(fabs(io.header.no_vpart)   < ZERO || 
      fabs(io.header.no_species) < ZERO || 
      fabs(io.header.min_weight) < ZERO || 
      fabs(io.header.max_weight) < ZERO ||
      fabs(io.header.med_weight) < ZERO)
     {
      fprintf(stderr," o io.header masses not initialized ...\n");
      
      init_header_masses();   
      
       /* init_header_masses() returns the values in physical units Msun/h */
      io.header.no_vpart   /= io.header.pmass;
      io.header.min_weight /= io.header.pmass;
      io.header.max_weight /= io.header.pmass;
      io.header.med_weight /= io.header.pmass;
     }
   

   
   /*=============
    *  MULTIMASS
    *=============*/
#ifndef MULTIMASS
   if(io.header.multi_mass == 1)
     {
      fprintf(stderr,"\n=================================================================\n");      
      fprintf(stderr,"      your are trying to run a multi-mass simulation:\n");      
      fprintf(stderr,"          please recompile AMIGA with -DMULTIMASS\n");
      fprintf(stderr,"=================================================================\n");  
      exit(0);
     }
#endif
#ifdef MULTIMASS
   if(io.header.multi_mass == 0)
     {
      fprintf(stderr,"\n=================================================================\n");      
      fprintf(stderr,"     your are trying to run a sinlge-mass simulation:\n");      
      fprintf(stderr,"        please recompile AMIGA without -DMULTIMASS\n");
      fprintf(stderr,"=================================================================\n");      
      exit(0);
     }
#endif
   
   /*============
    *   DOUBLE
    *============*/
#ifdef DOUBLE
   if(io.header.double_precision != 1)
     {
      fprintf(stderr,"\n input file is single precision but simulation will be run in double precision\n");
      fprintf(stderr,"   => will upcast\n");      
     }
   
   io.header.double_precision = 1;
#else
   if(io.header.double_precision == 1)
     {
      fprintf(stderr,"\n input file is double precision but simulation will be run in single precision\n");
      fprintf(stderr,"   => will downcast\n");      
     }
   
   io.header.double_precision = 0;
#endif
   
   fprintf(stderr," <= finished sanity_check()\n");
}


/**
 * cmp_sfckey_part: compares the sfc keys of two particles, used
 * for qsort 
 */
extern int
cmp_sfckey_part(const void *p1, const void *p2) 
{
	if (((partptr)p1)->sfckey < ((partptr)p2)->sfckey)
		return -1;

	if (((partptr)p1)->sfckey > ((partptr)p2)->sfckey)
		return 1;

	return 0;
}


/*==============================================================================
 * here we decide what to write into the logfile to keep the user up-to-date...
 *==============================================================================*/
void write_logfile(double timecounter, double timestep, int no_timestep)
{
  
  gridls *for_grid;          /* for looping over all grids          */
  double  RAM;               /* how much RAM has been used                    */
  double  PDgrowth;
  long    dom_nopart;        /* number of particles linked to domain grid     */
  double  da_a;
  double  total_time_step;
  
  /* write current timestep into logfile */
  fprintf(io.logfile, "\n");
  fprintf(io.logfile, "supercomoving T = %f\n", timecounter);
  fprintf(io.logfile, "scale factor  a = %f\n", calc_super_a(timecounter));
  fprintf(io.logfile, "redshift      z = %f\n", 1.0/calc_super_a(timecounter)-1.0);
  fflush(io.logfile);
  
  /* reset time counter */
  grid_timing.potential = 0;
  grid_timing.density   = 0;
  grid_timing.DK        = 0;
  grid_timing.grid      = 0;
  grid_timing.hydro     = 0;
  
  
  /* get number of particles attached to domain grid */
  if(global.fin_l1dim > global.dom_grid->l1dim)
   {
    dom_nopart = global.no_part;
    for(for_grid=(global.dom_grid+1); for_grid->l1dim<global.fin_l1dim; for_grid++)
      dom_nopart -= for_grid->size.no_part;
    dom_nopart -= for_grid->size.no_part;
    
    global.dom_grid->size.no_part = dom_nopart;
   }
  
  /* start accumulating RAM used during this step */
  RAM = simu.no_halos*sizeof(HALO) * bytes2GB;
  
  
  fprintf(io.logfile,"\n");
  fprintf(io.logfile,"grid information\n");
  fprintf(io.logfile,"----------------\n");
  /* write grid information to logfile */
  for(for_grid = global.dom_grid; for_grid->l1dim < global.fin_l1dim; for_grid++)
   {
    fprintf(io.logfile,"GRID %12ld: nodes=%12lu (%8.3g GB) npart=%12lu - TIME: pot=%6ld dens=%6ld DK=%6ld grid=%6ld hydro=%6ld - SWEEPS: %ld %g\n",
            for_grid->l1dim,
            for_grid->size.no_nodes, for_grid->size.no_nodes * global.bytes_node * bytes2GB,
            for_grid->size.no_part,
            for_grid->time.potential,
            for_grid->time.density,
            for_grid->time.DK,
            for_grid->time.grid,
            for_grid->time.hydro,
            for_grid->no_sweeps,
            for_grid->cur_resid);
    
    grid_timing.potential += for_grid->time.potential;
    grid_timing.density   += for_grid->time.density;
    grid_timing.DK        += for_grid->time.DK;
    grid_timing.grid      += for_grid->time.grid;
    grid_timing.hydro     += for_grid->time.hydro;
    RAM                   += for_grid->size.no_nodes * global.bytes_node * bytes2GB;
    
    for_grid->time.potential = 0;
    for_grid->time.density   = 0;
    for_grid->time.DK        = 0;
    for_grid->time.grid      = 0;
    for_grid->time.hydro     = 0;
   }
  
  /* treat grid with for_grid->l1dim==global.fin_l1dim separately
   * in order to avoid incrementing for_grid++ too far!!! */
  
  /* write grid information to logfile */
  /*-----------------------------------*/
  fprintf(io.logfile,"GRID %12ld: nodes=%12lu (%8.3g GB) npart=%12lu - TIME: pot=%6ld dens=%6ld grid=%6ld\n",
          for_grid->l1dim,
          for_grid->size.no_nodes, for_grid->size.no_nodes * global.bytes_node * bytes2GB,
          for_grid->size.no_part,
          for_grid->time.potential,
          for_grid->time.density,
          for_grid->time.grid);
  
  grid_timing.potential += for_grid->time.potential;
  grid_timing.density   += for_grid->time.density;
  grid_timing.DK        += for_grid->time.DK;
  grid_timing.grid      += for_grid->time.grid;
  grid_timing.hydro     += for_grid->time.hydro;
  RAM                   += for_grid->size.no_nodes * global.bytes_node * bytes2GB;
  
  for_grid->time.potential = 0;
  for_grid->time.density   = 0;
  for_grid->time.DK        = 0;
  for_grid->time.grid      = 0;
  for_grid->time.hydro     = 0;  
  fprintf(io.logfile, "                                                                                   %6ld      %6ld    %6ld\n", grid_timing.potential, grid_timing.density, grid_timing.grid);
  
  
  /* write detailed breakdown of timing */
  /*------------------------------------*/
  fprintf(io.logfile,"detailed timing information (in sec.)\n");
  fprintf(io.logfile,"-------------------------------------\n");
  fprintf(io.logfile,"io           = %ld\n",timing.io);
  fprintf(io.logfile,"      - startrun     = %ld\n",timing.startrun);
#ifdef AHFptfocus
  fprintf(io.logfile,"      - ptfocus      = %ld\n",timing.ptfocus);
#endif
#ifdef AHFrfocus
  fprintf(io.logfile,"      - rfocus       = %ld\n",timing.rfocus);
#endif
#ifdef WITH_MPI
  fprintf(io.logfile,"      - loadbalance  = %ld\n",timing.loadbalance);
  fprintf(io.logfile,"      - distribution = %ld\n",timing.distribution);
#else
  fprintf(io.logfile,"      - sfckey       = %ld\n",timing.sfckey);
#endif
  fprintf(io.logfile,"gendomgrids  = %ld\n",timing.gendomgrids);
  fprintf(io.logfile,"ll           = %ld\n",timing.ll);
  fprintf(io.logfile,"genrefgrids  = %ld\n",timing.genrefgrids);
  fprintf(io.logfile,"densrecovery = %ld\n",timing.densrecovery);
#ifdef AHFpotcentre
  fprintf(io.logfile,"potcentre    = %ld\n",timing.potcentre);
#endif
  fprintf(io.logfile,"ahf_gridinfo = %ld\n",timing.ahf_gridinfo);
  fprintf(io.logfile,"ahf_halos    = %ld\n",timing.ahf_halos);
  fprintf(io.logfile,"      - RefCentre                   = %ld\n",timing.RefCentre);
  fprintf(io.logfile,"      - analyseRef                  = %ld\n",timing.analyseRef);
  fprintf(io.logfile,"      - spatialRef2halos            = %ld\n",timing.spatialRef2halos);
  fprintf(io.logfile,"      - ahf_halos_sfc_constructHalo = %ld\n",timing.ahf_halos_sfc_constructHalo);
  fprintf(io.logfile,"      - I/O                         = %ld\n",timing.ahf_io);
  
  
  /* write summary information */
  /*---------------------------*/
  /* finish logfile entries for this step... */
  total_time_step = (grid_timing.potential +
                     grid_timing.density   +
                     grid_timing.DK        +
                     grid_timing.grid      +
                     grid_timing.hydro     +
                     ahf.time              +
                     energy.time           +
                     PkSpectrum.time        )     /3600.;
  total_time_step   += timing.io/3600.;
  global.total_time += total_time_step;
  
  /* add particles to RAM */
#ifndef WITH_MPI
  RAM += global.no_part * global.bytes_part * bytes2GB;
#else
  RAM += global_info.no_part * global.bytes_part * bytes2GB;
#endif /* WITH_MPI */
  
  fprintf(io.logfile,"\n");
  fprintf(io.logfile,"summary information\n");
  fprintf(io.logfile,"-------------------\n");
  fprintf(io.logfile, "force resolution    ~ %8.2g kpc/h\n", 3000.*simu.boxsize/(double)for_grid->l1dim);
  if(ahf.time > 0)
    fprintf(io.logfile, "time for AHF        = %8ld seconds (%8.3g hours)\n",ahf.time,ahf.time/3600.);
  fprintf(io.logfile, "total time          = %8.0f seconds (%8.3g hours)\n", total_time_step * 3600., total_time_step);
  fprintf(io.logfile, "cumulative time     = %8.3g hours (%g days)\n", global.total_time, global.total_time/24.);
  fprintf(io.logfile, "memory during step  = %8.3g GB\n",RAM);
#	ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(0),
	               "total time (MPI measure) =  %f seconds (%f hours)",
	               (global_mpi.stop-global_mpi.start),
	               (global_mpi.stop-global_mpi.start)/3600.);
#	endif
  
  fflush(io.logfile);
  
  ahf.time        = 0;
  energy.time     = 0;
  PkSpectrum.time = 0;
}


/*==============================================================================
 * just dump all relevant parameter into a file >prefix<_parameter
 *==============================================================================*/
void write_parameterfile()
{
  FILE  *fpparam;
  FILE  *codeinfo;
  double dummy;
  char   param_file[MAXSTRING], dummyline[MAXSTRING];
  double a, a3, omega, lambda, ovlim, rho_crit, rho_b, rho_vir, Hubble;
  
#ifdef WITH_MPI
  /* make sure to only write one parameter file when running on multile CPUs */
  if(global_mpi.rank == 0)
#endif
    if(global_io.params->outfile_prefix != NULL)
     {
      strcpy(param_file, global_io.params->outfile_prefix);
      strcat(param_file,".parameter");
      
#ifdef RESTART_PARAMETER_FILE
      /* check if there is already a paramater file */
      if( (fpparam = fopen(param_file,"r")) != NULL)
       {
        fclose(fpparam);
        fprintf(stderr," NOTE: parameter file %s already exists\n",param_file);
        strcat(param_file,"-restart");
        fprintf(stderr,"       writing parameter to %s instead!\n",param_file);
       }
#endif
      
      if( (fpparam = fopen(param_file,"w")) == NULL)
       {
        fprintf(io.logfile," NOTE: could not open %s to dump parameter\n",param_file);
       }
      else
       {
        /* write >>all<< relevant paramater to file */
        WRITEAHFLOGO(fpparam);
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"AHF related parameter:\n");
        fprintf(fpparam,"----------------------\n");
        fprintf(fpparam,"AHF_MINPART                \t\t%d\n",          simu.AHF_MINPART);
        fprintf(fpparam,"AHF_MINPART_GAS            \t\t%d\n",           AHF_MINPART_GAS);
        fprintf(fpparam,"AHF_MINPART_STARS          \t\t%d\n",         AHF_MINPART_STARS);
        fprintf(fpparam,"AHF_VTUNE                  \t\t%g\n",            simu.AHF_VTUNE);
        fprintf(fpparam,"AHF_RISE                   \t\t%g\n",                  AHF_RISE);
        fprintf(fpparam,"AHF_SLOPE                  \t\t%g\n",                 AHF_SLOPE);
        fprintf(fpparam,"AHF_MAXNRISE               \t\t%d\n",              AHF_MAXNRISE);
        fprintf(fpparam,"AHF_MAX_GATHER_RAD         \t\t%g\n",         simu.MaxGatherRad);
        fprintf(fpparam,"AHF_MIN_REF_OFFSET         \t\t%d\n",        AHF_MIN_REF_OFFSET);
        fprintf(fpparam,"AHF_NBIN_MULTIPLIER        \t\t%d\n",       AHF_NBIN_MULTIPLIER);
        fprintf(fpparam,"AHF_HIRES_DM_WEIGHT        \t\t%g\n",       AHF_HIRES_DM_WEIGHT);
        fprintf(fpparam,"AHF_HOSTHALOLEVEL          \t\t%d\n",         AHF_HOSTHALOLEVEL);
        fprintf(fpparam,"PGAS                       \t\t%g\n",                      PGAS);
        fprintf(fpparam,"PDM                        \t\t%g\n",                       PDM);
        fprintf(fpparam,"PSTAR                      \t\t%g\n",                     PSTAR);
        fprintf(fpparam,"PDMbndry                   \t\t%g\n",                  PDMbndry);
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"simulation related values:\n");
        fprintf(fpparam,"--------------------------\n");
        a        = global.a;
        a3       = pow3(a);
        omega    = calc_omega(a);
        lambda   = calc_lambda(a);
        ovlim    = calc_virial(a);
        Hubble   = calc_Hubble(a);          /* in km/sec/Mpc */
        rho_crit = a3 * calc_rho_crit(a);   /* comoving(!) critical density   */
        rho_b    = omega * rho_crit;        /* comoving(!) background density */
        rho_vir  = a3*calc_rho_vir(a);      /* comoving(!) normalisation density */
        fprintf(fpparam,"a                          \t\t%g\n", a);
        fprintf(fpparam,"z                          \t\t%g\n", 1./a-1.);
        fprintf(fpparam,"Omega(z)                   \t\t%g\n", omega);
        fprintf(fpparam,"OmegaL(z)                  \t\t%g\n", lambda);
        fprintf(fpparam,"rho_crit(z)                \t\t%g\n", rho_crit);
        fprintf(fpparam,"rho_back(z)                \t\t%g\n", rho_b);
        fprintf(fpparam,"rho_vir(z)                 \t\t%g\n", rho_vir);
        fprintf(fpparam,"Delta_vir(z)               \t\t%g\n", ovlim);
        fprintf(fpparam,"Hubble(z)                  \t\t%g\n", Hubble);
        
        if(fabs(simu.GADGET_m2Msunh) > ZERO)
         {
          fprintf(fpparam,"\n");
          fprintf(fpparam,"GADGET related parameter:\n");
          fprintf(fpparam,"-------------------------\n");
          fprintf(fpparam,"GADGET_MUNIT               \t\t%g\n",       simu.GADGET_m2Msunh);
          fprintf(fpparam,"GADGET_LUNIT               \t\t%g\n",        simu.GADGET_l2Mpch);
         }
        if((codeinfo = fopen("tipsy.info","r")))
         {
          fprintf(fpparam,"\n");
          fprintf(fpparam,"TIPSY related parameter:\n");
          fprintf(fpparam,"------------------------\n");
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"TIPSY_OMEGA0                \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"TIPSY_LAMBDA0               \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"TIPSY_BOXSIZE               \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"TIPSY_VUNIT                 \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"TIPSY_MUNIT                 \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"TIPSY_EUNIT                 \t\t%g\n",       dummy);
          fclose(codeinfo);
         }
        if((codeinfo = fopen("art.info","r")))
         {
          fprintf(fpparam,"\n");
          fprintf(fpparam,"ART related parameter:\n");
          fprintf(fpparam,"----------------------\n");
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"ART_BOXSIZE                \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"ART_MUNIT                  \t\t%g\n",       dummy);
          fclose(codeinfo);
         }
        if((codeinfo = fopen("cubep3m.info","r")))
         {
          fprintf(fpparam,"\n");
          fprintf(fpparam,"CUBEP3M related parameter:\n");
          fprintf(fpparam,"--------------------------\n");
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"CUBEP3M_OMEGA0                \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"CUBEP3M_LAMBDA0               \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"CUBEP3M_BOXSIZE               \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"CUBEP3M_NGRID                 \t\t%g\n",       dummy);
          fgets(dummyline,MAXSTRING,codeinfo);
          sscanf(dummyline,"%lf", &dummy);
          fprintf(fpparam,"CUBEP3M_NODES_DIM             \t\t%g\n",       dummy);
          fclose(codeinfo);
         }
       
        fprintf(fpparam,"\n");
        fprintf(fpparam,"general parameter:\n");
        fprintf(fpparam,"------------------\n");
        fprintf(fpparam,"NGRID_DOM                  \t\t%d\n",            simu.NGRID_DOM);
        fprintf(fpparam,"NGRID_MAX                  \t\t%d\n",            simu.NGRID_MAX);
        fprintf(fpparam,"MIN_NNODES                 \t\t%d\n",                MIN_NNODES);
        fprintf(fpparam,"MAXTIME                    \t\t%d\n",                   MAXTIME);
        fprintf(fpparam,"MAXSTRING                  \t\t%d\n",                 MAXSTRING);
        fprintf(fpparam,"AMIGAHEADER                \t\t%d\n",               AMIGAHEADER);
        fprintf(fpparam,"HEADERSTRING               \t\t%d\n",              HEADERSTRING);
        fprintf(fpparam,"HEADERSIZE                 \t\t%d\n",           (int)HEADERSIZE);
        fprintf(fpparam,"FILLHEADER                 \t\t%d\n",           (int)FILLHEADER);
        fprintf(fpparam,"CRITMULTI                  \t\t%g\n",                 CRITMULTI);
        fprintf(fpparam,"NP_RATIO                   \t\t%g\n",                  NP_RATIO);
        fprintf(fpparam,"ZERO                       \t\t%g\n",                      ZERO);
        fprintf(fpparam,"MACHINE_ZERO               \t\t%g\n",              MACHINE_ZERO);
        fprintf(fpparam,"MZERO                      \t\t%g\n",                     MZERO);
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"gravity solver related parameter:\n");
        fprintf(fpparam,"---------------------------------\n");
        fprintf(fpparam,"DOMSWEEPS                  \t\t%d\n",                 DOMSWEEPS);
        fprintf(fpparam,"REFSWEEPS                  \t\t%d\n",                 REFSWEEPS);
        fprintf(fpparam,"W_SOR                      \t\t%g\n",                     W_SOR);
        fprintf(fpparam,"ETA                        \t\t%g\n",                       ETA);
        fprintf(fpparam,"CONVCRIT                   \t\t%g\n",                  CONVCRIT);
        fprintf(fpparam,"DOMCORRECT                 \t\t%g\n",                DOMCORRECT);
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"MPI related parameter:\n");
        fprintf(fpparam,"----------------------\n");
        fprintf(fpparam,"BITS_PER_DIMENSION         \t\t%d\n",        BITS_PER_DIMENSION);
        fprintf(fpparam,"LOADBALANCE_DOMAIN_LEVEL   \t\t%d\n",             simu.lb_level);
        fprintf(fpparam,"MAX_SEND_PARTICLES         \t\t%d\n",        MAX_SEND_PARTICLES);
        fprintf(fpparam,"VERBOSITY                  \t\t%d\n",                 VERBOSITY);
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"physical constants:\n");
        fprintf(fpparam,"-------------------\n");
        fprintf(fpparam,"Gyr                        \t\t%g\n",                       Gyr);
        fprintf(fpparam,"Mpc                        \t\t%g\n",                       Mpc);
        fprintf(fpparam,"H0                         \t\t%g\n",                        H0);
        fprintf(fpparam,"rhoc0                      \t\t%g\n",                     rhoc0);
        fprintf(fpparam,"Grav                       \t\t%g\n",                      Grav);
        fprintf(fpparam,"cH0                        \t\t%g\n",                       cH0);
        fprintf(fpparam,"kB_per_mp                  \t\t%g\n",                 kB_per_mp);
        fprintf(fpparam,"kBoltzman                  \t\t%g\n",                 kBoltzman);
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"*=====================*\n");
        fprintf(fpparam,"     DEFINEFLAGS:\n");
        fprintf(fpparam,"*=====================*\n");
        fprintf(fpparam,"\n");
        fprintf(fpparam,"AHF related flags:\n");
        fprintf(fpparam,"------------------\n");
#ifdef AHF
        fprintf(fpparam,"AHF                        \t\t1\n");
#else
        fprintf(fpparam,"AHF                        \t\t0\n");
#endif
#ifdef AHFlean
        fprintf(fpparam,"AHFlean                    \t\t1\n");
#else
        fprintf(fpparam,"AHFlean                    \t\t0\n");
#endif
#ifdef AHFfast
        fprintf(fpparam,"AHFfast                    \t\t1\n");
#else
        fprintf(fpparam,"AHFfast                    \t\t0\n");
#endif
#ifdef AHFnoHubbleDrag
        fprintf(fpparam,"AHFnoHubbleDrag            \t\t1\n");
#else
        fprintf(fpparam,"AHFnoHubbleDrag            \t\t0\n");
#endif
#ifdef DARK_ENERGY
        fprintf(fpparam,"DARK_ENERGY                \t\t1\n");
#else
        fprintf(fpparam,"DARK_ENERGY                \t\t0\n");
#endif
#ifdef AHFdisks
        fprintf(fpparam,"AHFdisks                   \t\t1\n");
#else
        fprintf(fpparam,"AHFdisks                   \t\t0\n");
#endif
#ifdef AHFdensrecovery
        fprintf(fpparam,"AHFdensrecovery            \t\t1\n");
#else
        fprintf(fpparam,"AHFdensrecovery            \t\t0\n");
#endif
#ifdef AHFmaxdenscentre
        fprintf(fpparam,"AHFmaxdenscentre           \t\t1\n");
#else
        fprintf(fpparam,"AHFmaxdenscentre           \t\t0\n");
#endif
#ifdef AHFptfocus
        fprintf(fpparam,"AHFptfocus=                 \t\t%d\n",AHFptfocus);
#else
        fprintf(fpparam,"AHFptfocus=                 \t\t-\n");
#endif
#ifdef AHFrfocus
        fprintf(fpparam,"AHFrfocus=                  \t\tX=%g Y=%g Z=%g R=%g [Mpc/h]\n",AHFrfocusX,AHFrfocusY,AHFrfocusZ,AHFrfocusR);
#else
        fprintf(fpparam,"AHFrfocus=                  \t\t-\n");
#endif
#ifdef AHFnoremunbound
        fprintf(fpparam,"AHFnoremunbound             \t\t1\n");
#else
        fprintf(fpparam,"AHFnoremunbound             \t\t0\n");
#endif
#ifdef AHFundoPositionShiftAndScale
        fprintf(fpparam,"AHFundoPositionShiftAndScale\t\t1\n");
#else
        fprintf(fpparam,"AHFundoPositionShiftAndScale\t\t0\n");
#endif
#ifdef AHFvmbp
        fprintf(fpparam,"AHFvmbp                     \t\t1\n");
#else
        fprintf(fpparam,"AHFvmbp                     \t\t0\n");
#endif
#ifdef AHFaddDMonlyproperties
        fprintf(fpparam,"AHFaddDMonlyproperties      \t\t1\n");
#else
        fprintf(fpparam,"AHFaddDMonlyproperties      \t\t0\n");
#endif
#ifdef AHFignore_ugas
        fprintf(fpparam,"AHFignore_ugas              \t\t1\n");
#else
        fprintf(fpparam,"AHFignore_ugas              \t\t0\n");
#endif
#ifdef AHFgeom
        fprintf(fpparam,"AHFgeom                    \t\t1\n");
#else
        fprintf(fpparam,"AHFgeom                    \t\t0\n");
#endif
#ifdef AHFmaxdenscentre
        fprintf(fpparam,"AHFmaxdenscentre           \t\t1\n");
#else
        fprintf(fpparam,"AHFmaxdenscentre           \t\t0\n");
#endif
#ifdef AHFgeomcentre
        fprintf(fpparam,"AHFgeomcentre              \t\t1\n");
#else
        fprintf(fpparam,"AHFgeomcentre              \t\t0\n");
#endif
#ifdef AHFcomcentre
        fprintf(fpparam,"AHFcomcentre               \t\t1\n");
#else
        fprintf(fpparam,"AHFcomcentre               \t\t0\n");
#endif
#ifdef AHFpotcentre
        fprintf(fpparam,"AHFpotcentre               \t\t1\n");
#else
        fprintf(fpparam,"AHFpotcentre               \t\t0\n");
#endif
#ifdef AHFphi_infty
        fprintf(fpparam,"AHFphi_infty               \t\t1\n");
#else
        fprintf(fpparam,"AHFphi_infty               \t\t0\n");
#endif
#ifdef AHFsplinefit
        fprintf(fpparam,"AHFsplinefit               \t\t1\n");
#else
        fprintf(fpparam,"AHFsplinefit               \t\t0\n");
#endif
#ifdef AHFprofilerise
        fprintf(fpparam,"AHFprofilerise             \t\t1\n");
#else
        fprintf(fpparam,"AHFprofilerise             \t\t0\n");
#endif
#ifdef AHFreducedinertiatensor
        fprintf(fpparam,"AHFreducedinertiatensor    \t\t1\n");
#else
        fprintf(fpparam,"AHFreducedinertiatensor    \t\t0\n");
#endif
#ifdef AHFcentrefile
        fprintf(fpparam,"AHFcentrefile              \t\t1\n");
#else
        fprintf(fpparam,"AHFcentrefile              \t\t0\n");
#endif
#ifdef AHFcentrefileBASIC
        fprintf(fpparam,"AHFcentrefileBASIC         \t\t1\n");
#else
        fprintf(fpparam,"AHFcentrefileBASIC         \t\t0\n");
#endif
#ifdef DPhalos
        fprintf(fpparam,"DPhalos                    \t\t1\n");
#else
        fprintf(fpparam,"DPhalos                    \t\t0\n");
#endif
#ifdef AHFshellshape
        fprintf(fpparam,"AHFshellshape              \t\t1\n");
#else
        fprintf(fpparam,"AHFshellshape              \t\t0\n");
#endif
#ifdef AHFsplit_only
        fprintf(fpparam,"AHFsplit_only              \t\t1\n");
#else
        fprintf(fpparam,"AHFsplit_only              \t\t0\n");
#endif
#ifdef AHFexciseSubhaloStars
        fprintf(fpparam,"AHFexciseSubhaloStars      \t\t1\n");
#else
        fprintf(fpparam,"AHFexciseSubhaloStars      \t\t0\n");
#endif
#ifdef AHFnewHaloIDs
        fprintf(fpparam,"AHFnewHaloIDs              \t\t1\n");
#else
        fprintf(fpparam,"AHFnewHaloIDs              \t\t0\n");
#endif
#ifdef AHFrestart
        fprintf(fpparam,"AHFrestart                 \t\t1\n");
#else
        fprintf(fpparam,"AHFrestart                 \t\t0\n");
#endif
#ifdef AHFbinary
        fprintf(fpparam,"AHFbinary                  \t\t1\n");
#else
        fprintf(fpparam,"AHFbinary                  \t\t0\n");
#endif
#ifdef AHFsubstructure
        fprintf(fpparam,"AHFsubstructure            \t\t1\n");
#else
        fprintf(fpparam,"AHFsubstructure            \t\t0\n");
#endif
#ifdef AHFdmonly_Rmax_r2
        fprintf(fpparam,"AHFdmonly_Rmax_r2          \t\t1\n");
#else
        fprintf(fpparam,"AHFdmonly_Rmax_r2          \t\t0\n");
#endif
#ifdef AHFparticle_Rmax_r2
        fprintf(fpparam,"AHFparticle_Rmax_r2        \t\t1\n");
#else
        fprintf(fpparam,"AHFparticle_Rmax_r2        \t\t0\n");
#endif
#ifdef AHFgridinfofile
        fprintf(fpparam,"AHFgridinfofile            \t\t1\n");
#else
        fprintf(fpparam,"AHFgridinfofile            \t\t0\n");
#endif
#ifdef AHFspatialReffile
        fprintf(fpparam,"AHFspatialReffile          \t\t1\n");
#else
        fprintf(fpparam,"AHFspatialReffile          \t\t0\n");
#endif
#ifdef AHFnewCloseRefDist
        fprintf(fpparam,"AHFnewCloseRefDist         \t\t1\n");
#else
        fprintf(fpparam,"AHFnewCloseRefDist         \t\t0\n");
#endif
#ifdef AHFmaxGatherRadTest
        fprintf(fpparam,"AHFmaxGatherRadTest        \t\t1\n");
#else
        fprintf(fpparam,"AHFmaxGatherRadTest        \t\t0\n");
#endif
#ifdef AHF_LRSI
        fprintf(fpparam,"AHF_LRSI                   \t\t1\n");
#else
        fprintf(fpparam,"AHF_LRSI                   \t\t0\n");
#endif
#ifdef PARDAU_DISTANCE
        fprintf(fpparam,"PARDAU_DISTANCE            \t\t1\n");
#else
        fprintf(fpparam,"PARDAU_DISTANCE            \t\t0\n");
#endif
#ifdef PARDAU_NODES
        fprintf(fpparam,"PARDAU_NODES               \t\t1\n");
#else
        fprintf(fpparam,"PARDAU_NODES               \t\t0\n");
#endif
#ifdef PARDAU_PARTS
        fprintf(fpparam,"PARDAU_PARTS               \t\t1\n");
#else
        fprintf(fpparam,"PARDAU_PARTS               \t\t0\n");
#endif
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"various flags:\n");
        fprintf(fpparam,"--------------\n");
#ifdef WITH_MPI
        fprintf(fpparam,"WITH_MPI                   \t\t1\n");
#else
        fprintf(fpparam,"WITH_MPI                   \t\t0\n");
#endif
#ifdef BYTESWAP
        fprintf(fpparam,"BYTESWAP                   \t\t1\n");
#else
        fprintf(fpparam,"BYTESWAP                   \t\t0\n");
#endif
#ifdef MULTIMASS
        fprintf(fpparam,"MULTIMASS                  \t\t1\n");
#else
        fprintf(fpparam,"MULTIMASS                  \t\t0\n");
#endif
#ifdef GAS_PARTICLES
        fprintf(fpparam,"GAS_PARTICLES              \t\t1\n");
#else
        fprintf(fpparam,"GAS_PARTICLES              \t\t0\n");
#endif
#ifdef VERBOSE
        fprintf(fpparam,"VERBOSE                    \t\t1\n");
#else
        fprintf(fpparam,"VERBOSE                    \t\t0\n");
#endif
#ifdef VERBOSE2
        fprintf(fpparam,"VERBOSE2                   \t\t1\n");
#else
        fprintf(fpparam,"VERBOSE2                   \t\t0\n");
#endif
#ifdef DOUBLE
        fprintf(fpparam,"DOUBLE                     \t\t1\n");
#else
        fprintf(fpparam,"DOUBLE                     \t\t0\n");
#endif
        fprintf(fpparam,"\n");
        fprintf(fpparam,"grid related flags:\n");
        fprintf(fpparam,"-------------------\n");
#ifdef REFINE_BARYONIC_MASS
        fprintf(fpparam,"REFINE_BARYONIC_MASS       \t\t1\n");
#else
        fprintf(fpparam,"REFINE_BARYONIC_MASS       \t\t0\n");
#endif
#ifdef TSC
        fprintf(fpparam,"TSC                        \t\t1\n");
#else
        fprintf(fpparam,"TSC                        \t\t0\n");
#endif
#ifdef CIC
        fprintf(fpparam,"CIC                        \t\t1\n");
#else
        fprintf(fpparam,"CIC                        \t\t0\n");
#endif
#ifdef NGP
        fprintf(fpparam,"NGP                        \t\t1\n");
#else
        fprintf(fpparam,"NGP                        \t\t0\n");
#endif
        
        
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"gravity solver related flags:\n");
        fprintf(fpparam,"-----------------------------\n");
#ifdef FIXED_SWEEPS
        fprintf(fpparam,"FIXED_SWEEPS               \t\t1\n");
#else
        fprintf(fpparam,"FIXED_SWEEPS               \t\t0\n");
#endif
#ifdef NO_SOR
        fprintf(fpparam,"NO_SOR                     \t\t1\n");
#else
        fprintf(fpparam,"NO_SOR                     \t\t0\n");
#endif
#ifdef SWAP_XY
        fprintf(fpparam,"SWAP_XY                    \t\t1\n");
#else
        fprintf(fpparam,"SWAP_XY                    \t\t0\n");
#endif
#ifdef SWAP_XZ
        fprintf(fpparam,"SWAP_XZ                    \t\t1\n");
#else
        fprintf(fpparam,"SWAP_XZ                    \t\t0\n");
#endif
        
        
        
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"GADGET related flags:\n");
        fprintf(fpparam,"---------------------\n");
#ifdef GADGET
        fprintf(fpparam,"GADGET                     \t\t1\n");
#else
        fprintf(fpparam,"GADGET                     \t\t0\n");
#endif
#ifdef GADGET_GAS_ONLY
        fprintf(fpparam,"GADGET_GAS_ONLY            \t\t1\n");
#else
        fprintf(fpparam,"GADGET_GAS_ONLY            \t\t0\n");
#endif
#ifdef GADGET_STARS_ONLY
        fprintf(fpparam,"GADGET_STARS_ONLY          \t\t1\n");
#else
        fprintf(fpparam,"GADGET_STARS_ONLY          \t\t0\n");
#endif
#ifdef GADGET_LUNIT_KPC
        fprintf(fpparam,"GADGET_LUNIT_KPC           \t\t1\n");
#else
        fprintf(fpparam,"GADGET_LUNIT_KPC           \t\t0\n");
#endif
#ifdef METALHACK
        fprintf(fpparam,"METALHACK                  \t\t1\n");
#else
        fprintf(fpparam,"METALHACK                  \t\t0\n");
#endif
        
        
        fprintf(fpparam,"\n");
        fprintf(fpparam,"misc. flags:\n");
        fprintf(fpparam,"------------\n");
#ifdef NCPUREADING_EQ_NFILES
        fprintf(fpparam,"NCPUREADING_EQ_NFILES      \t\t1\n");
#else
        fprintf(fpparam,"NCPUREADING_EQ_NFILES      \t\t0\n");
#endif
#ifdef BCASTHEADER
        fprintf(fpparam,"BCASTHEADER                \t\t1\n");
#else
        fprintf(fpparam,"BCASTHEADER                \t\t0\n");
#endif
#ifdef FOPENCLOSE
        fprintf(fpparam,"FOPENCLOSE                 \t\t1\n");
#else
        fprintf(fpparam,"FOPENCLOSE                 \t\t0\n");
#endif
#ifdef CHECK_RLIMIT_NOFILE
        fprintf(fpparam,"CHECK_RLIMIT_NOFILE        \t\t1\n");
#else
        fprintf(fpparam,"CHECK_RLIMIT_NOFILE        \t\t0\n");
#endif
#ifdef REF_TEST
        fprintf(fpparam,"REF_TEST                   \t\t1\n");
#else
        fprintf(fpparam,"REF_TEST                   \t\t0\n");
#endif
#ifdef AHFDEBUG
        fprintf(fpparam,"AHFDEBUG                   \t\t1\n");
#else
        fprintf(fpparam,"AHFDEBUG                   \t\t0\n");
#endif
        
#ifdef ART
        fprintf(fpparam,"ART                        \t\t1\n");
#else
        fprintf(fpparam,"ART                        \t\t0\n");
#endif
#ifdef TIPSY
        fprintf(fpparam,"TIPSY                      \t\t1\n");
#else
        fprintf(fpparam,"TIPSY                      \t\t0\n");
#endif
#ifdef TIPSY_ZOOMDATA
        fprintf(fpparam,"TIPSY_ZOOMDATA             \t\t1\n");
#else
        fprintf(fpparam,"TIPSY_ZOOMDATA             \t\t0\n");
#endif
#ifdef ASCII
        fprintf(fpparam,"ASCII                      \t\t1\n");
#else
        fprintf(fpparam,"ASCII                      \t\t0\n");
#endif
#ifdef MLAPM
        fprintf(fpparam,"MLAPM                      \t\t1\n");
#else
        fprintf(fpparam,"MLAPM                      \t\t0\n");
#endif
#ifdef AMIGA_ONE_FORMAT
        fprintf(fpparam,"AMIGA_ONE_FORMAT           \t\t1\n");
#else
        fprintf(fpparam,"AMIGA_ONE_FORMAT           \t\t0\n");
#endif
#ifdef ASCII_ONLYPOSITIONS
        fprintf(fpparam,"ASCII_ONLYPOSITIONS        \t\t1\n");
#else
        fprintf(fpparam,"ASCII_ONLYPOSITIONS        \t\t0\n");
#endif
        
        fclose(fpparam);
       }
      
      
     } 
  
}   

/*==============================================================================
 * small routine calculating cNFW according to Eq.(9) in Prada et al. (2012)
 * Note: we are only searching for concentrations in the unique range
 *       c = [2.3,100]   corresponding to approx.   V2_ratio = [1,5.9]
 *==============================================================================*/
double cNFWroot(double c, double V2_ratio)
{
  return(0.216*c/(log(1+c)-c/(1+c)) - V2_ratio);
}

double calc_cNFW(double V2_max, double V2_vir)
{
  double V2_ratio, a, b, c, r1;
  
  V2_ratio = V2_max/V2_vir;
  
  /* we cannot find a root in these cases */
  if (V2_ratio <= 1 || V2_ratio > 5.9)
    return(-1);
  
  /* the most simple bi-section root-finding method on intervall [2.2,100] */
  a=2.2;
  b=100;
  /* we only care about the concentration up to the 3rd decimal */
  while (b-a > 1e-3)
   {
    c = (a+b)/2;
    
    if (cNFWroot(a, V2_ratio)*cNFWroot(c, V2_ratio) > 0)
      a = c;
    else
      b = c;  
   }
  
#ifdef NEW_cNFW
  a=2.0;
  b=100.0;
  c = 5.0;  // starting at 5 is a better choice than 50 (thanks, Julian Onions, for that tip)
  while (b-a > 1e-5) {
    r1 = nfwroot(c) - ratio;
    if (r1 < 0)
      a = c;
    else
      b = c;
    c = (a+b)/2.0;
  }
#endif

  c = (a+b)/2.0;
  
  return(c);

}

/*==============================================================================
 * construct a unique halo ID out of some halo properties
 *==============================================================================*/
uint64_t getHaloID(HALO *halos, int i)
{
  uint64_t ID;
  
#ifdef AHFmixHaloIDandSnapID
  uint64_t snapid, haloid;
  
  haloid = (uint64_t)i;  
  snapid = (uint64_t)simu.isnap;
  
  /* combine */
  ID     = (snapid << 54) | haloid;
#else
  HALO    *halo;
  uint64_t Nid, Xid, Yid, Zid;
  double   Xc,Yc,Zc;
  
  /* access halo from halos[] array */
  halo = halos+i;
  
  /* Nid shall be encoded in the first 4 bits  = left-shift by 60 */
  Nid = (uint64_t) log10(halo->npart);
  Nid = Nid << 60;
  
  /* encode position using 20 bits each */
  Xc  = halo->pos.x;
  Xid = (uint64_t) (Xc*(double)(1LU<<20));
  Xid = Xid << 40;
  
  Yc  = halo->pos.y;
  Yid = (uint64_t) (Yc*(double)(1LU<<20));
  Yid = Yid << 20;
  
  Zc  = halo->pos.z;
  Zid = (uint64_t) (Zc*(double)(1LU<<20));
  Zid = Zid << 0;
  
  /* combine all IDs into one single number */
  ID = Nid;
  ID = ID | Xid;
  ID = ID | Yid;
  ID = ID | Zid;
#endif
  
  return(ID);
}

/*==============================================================================
 * construct a unique halo ID out of some halo properties
 *==============================================================================*/
uint64_t getSussing2013ID(int isnap, int ihalo)
{
  uint64_t ID;
  
  ID = (uint64_t) ( ((uint64_t)1e12)*isnap + (ihalo+1) );
  
  return(ID);
}



/*==============================================================================
 * return the version of the GADGET file format
 *==============================================================================*/
#ifdef LGADGET
long long    blklen;
#define GADGET_SKIP  ReadLongLong (fp,&blklen,FALSE);
#else
unsigned int blklen;
#define GADGET_SKIP  ReadUInt     (fp,&blklen,FALSE);
#endif
int check_gadgetversion(FILE *fp)
{
  char string[MAXSTRING];
  int version;
#ifdef LGADGET
#else
#endif
  
  rewind(fp);
  
  GADGET_SKIP;
  fread(string,sizeof(char),4,fp);
  if(strncmp("HEAD",string,4) == 0)
    version = 2;
  else
    version = 1;

  rewind(fp);
  
  return version;
}

/*==============================================================================
 * calculate the mininum level for the generation of the patch-tree
 *==============================================================================*/
int ptree_min_level(double Nthreshold)
{
  double a, a3, omega, ovlim, rho_vir, refine_mass, refine_ovdens, refine_vol, refine_len, fl1dim;
  int    min_level, max_level, i;

#ifdef VERBOSE
  fprintf(stderr,"\nptree_min_level: determining minimum level for patch-tree:\n");
#endif
  fprintf(io.logfile,"\nptree_min_level: determining minimum level for patch-tree:\n");
  
  // cosmology related stuff
	a        = global.a;
	a3       = pow3(a);
	omega    = calc_omega(a);
	ovlim    = calc_virial(a);
  rho_vir  = a3 * calc_rho_vir(a); // comoving(!) density used to normalize densities
  
  // get the number of the grid satisfying virial overdensity criterion
  refine_mass  = simu.pmass*simu.med_weight;
  
  
  // loop over all possible levels in the particle tree
  min_level = 0;
  max_level = (int)(log(global_io.params->NGRID_MAX)/log(2.0));
  
  
  for ( i=0; i <= max_level ; i++ )
   {
    fl1dim        = pow(2.0, (double)i);
    refine_len    = (double)(simu.boxsize/((double)(fl1dim)));
    refine_vol    = pow3(refine_len);
    refine_ovdens = (Nthreshold*refine_mass/refine_vol) / rho_vir;
#ifdef VERBOSE
    fprintf(stderr,"    level = %16d l1dim = %16.0f refine_ovdens = %16.4g ovlim = %16.4g\n", i,fl1dim,refine_ovdens,ovlim);
#endif
    fprintf(io.logfile,"    level = %16d l1dim = %16.0f refine_ovdens = %16.4g ovlim = %16.4g\n", i,fl1dim,refine_ovdens,ovlim);
    if ( refine_ovdens <  ovlim )
      min_level = i;
   }
 
  // store maximum overdensity possible with current refinement hierarchy
  global.max_ovdens = refine_ovdens;
  
  /* allow some leverage via arc/param.h */
  min_level += AHF_MIN_REF_OFFSET;
  
  fprintf(io.logfile,"\n   Cosmology:\n",a);
  fprintf(io.logfile,"   ==========\n",a);
  fprintf(io.logfile,"    a                = %lf\n",a);
  fprintf(io.logfile,"    omega            = %lf\n",omega);
  fprintf(io.logfile,"    ovlim            = %lf\n",ovlim);
  fprintf(io.logfile,"    rho_vir          = %g\n",rho_vir);
  fprintf(io.logfile,"    simu.pmass       = %g\n",simu.pmass);
  fprintf(io.logfile,"    simu.med_weight  = %lf\n",simu.med_weight);
  fprintf(io.logfile,"    refine_mass      = %g\n",refine_mass);
  fprintf(io.logfile,"    min_level        = %d\n",min_level);
  fprintf(io.logfile,"    max_level        = %d\n\n",max_level);
  fflush(io.logfile);
  
  return(min_level);
}





