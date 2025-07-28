#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif


#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libio_serial/io_serial.h"
#include "../src/libutility/utility.h"

/*==============================
 * switch on/off various things
 *==============================*/
#define TSC            // use TSC or NGP mass assignment? (no CIC implemented!)
//#define GAS_ONLY   // use only the gas particles

/*===================
 * COMMON structures
 *===================*/
info_io     io;
info_global global;
gridls     *cur_grid;

uint64_t NGRID;
uint64_t NGRID3;
FILE  *fpout;
double RhoMean, R_sphereReal;

/*===================
 * functions
 *===================*/
void  statistic (flouble *delta, double R_sphere, int N_sphere);
void  assign    (flouble *delta);
void  add_gas   (flouble *delta);

#define INDX(i,j,k)  ( (uint64_t)k + (uint64_t)j*(uint64_t)NGRID + (uint64_t)i*pow2((uint64_t)NGRID) )

/*===================
 * main()
 *===================*/
int main(argc,argv)
int argc;
char **argv;
{
  flouble      *delta;
  long unsigned i,j,k,np;
  int           no_timestep, ilen;
  double        vmean, vmin, vmax, vtemp;
  double        x_fac, v_fac;
  char          indata[MAXSTRING], header[HEADERSTRING];
  FILE         *fpin;
  double        R_sphere, R_min, R_max;
  int           N_sphere, N_rad;
  int           SWAPBYTES;
  
  printf("============================================================================\n");
  printf(" Read simulation data and calculate mass variance in sphere of given radius\n");
  printf("============================================================================\n");
  printf("Please give name of input  file:              ");
  scanf("%s", indata);
  printf("%s\n",indata);
  printf("please give radius of sphere:    ");
  scanf("%lf", &R_sphereReal);
  printf("%g\n",R_sphereReal);
  if(R_sphereReal < 0)
  {
    printf("please give number of radii:     ");
    scanf("%d", &N_rad);
    printf("%d\n",N_rad);
    printf("please give min. radius of sphere:    ");
    scanf("%lf", &R_min);
    printf("%g\n",R_min);
    printf("please give max. radius of sphere:    ");
    scanf("%lf", &R_max);
    printf("%g\n",R_max);
  }
  printf("please give number of spheres:   ");
  scanf("%d", &N_sphere);
  printf("%d\n",N_sphere);
  printf("please give NGRID:               ");
  scanf("%d", &NGRID);
  printf("%d\n",NGRID);
  NGRID3 = (uint64_t)NGRID*(uint64_t)NGRID*(uint64_t)NGRID;
  
  /*========================
   * read DM particles file
   *========================*/
  input(indata);
  
  /*==========================
   * assign particles to grid
   *==========================*/
  fprintf(stderr,"o allocating %g GB of memory for delta[] (NGRID^3=%"PRIu64")... ",((double)NGRID3)*sizeof(flouble)/1024./1024./1024.,NGRID3);
  delta = (flouble *) calloc(NGRID3, sizeof(flouble));
  fprintf(stderr,"done\n");
  if(delta == NULL) {
    fprintf(stderr,"not enough memory for dens[] (needed=%g GB)\n",((double)NGRID3)*sizeof(flouble)/1024./1024./1024.);
    exit(0);
  }
  fprintf(stderr,"o assigning particles to grid ... ");
  assign(delta);
  RhoMean = io.header.no_vpart / (double)NGRID3;
  fprintf(stderr,"done\n");
  
  /*==========================
   * calculate sigma(r)
   *==========================*/
  fpout = fopen("SigmaR.DAT","w");
  fprintf(fpout,"#N_spheres=%d, NGRID=%d\n",N_sphere,NGRID);
  fprintf(fpout,"#R(1) sigmaR(2) <M>(3,code units)\n");
  
  fprintf(stderr,"o calculating sigma(r) ... \n");
  
  if(R_sphereReal < 0)
  {
    for(i=0; i<N_rad; i++)
    {
      /* size of "i-sphere" */
      R_sphereReal = (double)i * (R_max-R_min)/((double)N_rad-1.0) + R_min;   //0.5 * (double)i/(double)N_rad * io.header.boxsize;
      
      /* scale Sphere to internal units */
      R_sphere = R_sphereReal * (double)NGRID / io.header.boxsize;
      
      statistic(delta, R_sphere, N_sphere);
    }
  }
  else
  {
    /* scale Sphere to internal units */
    R_sphere = R_sphereReal * (double)NGRID / io.header.boxsize;
    
    statistic(delta, R_sphere, N_sphere);
  }
  fclose(fpout);
}


/*===========================================================================
 * do you own stuff with simulation data given now in proper units...
 *===========================================================================*/
void statistic(flouble *delta, double R_sphere, int N_sphere)
{
  int    n,i,j,k,ii,jj,kk,imin,jmin,kmin,imax,jmax,kmax,icount;
  double D2, D_n, Rx, Ry, Rz, dr2, sigma_R, R_sphere2;
  double Nv_mean, Nv_in_sphere, sigN2;
  int    ISEED=-123456;
  uint64_t idx;
  
  
  /* how many particles do we expect in current sphere */
  Nv_mean   = RhoMean * 4.*PI/3.*pow3(R_sphere);
  
  R_sphere2 = pow2(R_sphere);
  
  /* loop over N_sphere Spheres chosen at random positions */
  sigN2 = 0.0;
#ifdef WITH_OPENMP
#pragma omp default(none) parallel reduction(+:sigN2) firstprivate(NGRID, NGRID3, ISEED, R_sphere, R_sphere2, Nv_mean) private(n, dr2, Rx, Ry, Rz, icount, Nv_in_sphere, idx, ii, jj, kk, i, j, k, imin, jmin, kmin, imax, jmax, kmax) shared(delta, N_sphere)
#pragma omp for schedule(static)
#endif
  for(n=0; n<N_sphere; n++)
  {
    /* random centre position [0,NGRID] */
    Rx = (double) (NGRID) * ran3(&ISEED);
    Ry = (double) (NGRID) * ran3(&ISEED);
    Rz = (double) (NGRID) * ran3(&ISEED);
    
    icount       = 0;
    Nv_in_sphere = 0.0;
    
    /* possible range of grid cells */
    imin = (int) (Rx-R_sphere);
    jmin = (int) (Ry-R_sphere);
    kmin = (int) (Rz-R_sphere);
    imax = (int) (Rx+R_sphere)+1;
    jmax = (int) (Ry+R_sphere)+1;
    kmax = (int) (Rz+R_sphere)+1;
    
    for(i=imin; i<=imax; i++)
      for(j=jmin; j<=jmax; j++)
        for(k=kmin; k<=kmax; k++)
        {
          dr2 = (pow2(Rx-i)+pow2(Ry-j)+pow2(Rz-k));
          if(dr2 <= R_sphere2)
          {
            icount++;
            
            /* take care of periodic boundary conditions */
            ii = (int) fmod((double)(i+NGRID),(double)NGRID);
            jj = (int) fmod((double)(j+NGRID),(double)NGRID);
            kk = (int) fmod((double)(k+NGRID),(double)NGRID);
            
            idx = INDX(ii,jj,kk);
            Nv_in_sphere += delta[idx];
          }
        }
    
    sigN2 += pow2(Nv_in_sphere - Nv_mean)/pow2(Nv_mean);
  }
  
  sigma_R = sqrt(sigN2 / (double)N_sphere);
  
  fprintf(stderr,"     sigma_%g = %g (<M>=%g)\n",R_sphereReal,sigma_R,Nv_mean);
  fprintf(fpout,"%g %g %g\n",R_sphereReal,sigma_R,Nv_mean);
  fflush(fpout);
}


/*====================================================================
 * assign:   assign particles to the grid points using a TSC scheme
 *====================================================================*/
void assign(flouble *dens)
{
  partptr       cur_part;
  uint64_t      ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1;
  long unsigned ipart;
  double        rrx, rry, rrz;
  double        hx, hy, hz;
  double        hx0, hy0, hz0, hxp1, hyp1, hzp1, hxm1, hym1, hzm1;
  double        weight, f_b;
  
  /* baryon fraction */
  f_b = 0.0;
  
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(NGRID) private(iz, iy, ix) shared(dens)
#pragma omp for schedule(static)
#endif
  for(ix=0; ix<NGRID; ix++)
    for(iy=0; iy<NGRID; iy++)
      for(iz=0; iz<NGRID; iz++) {
        dens[INDX(ix,iy,iz)]=0.0;
      }
  
#ifndef GAS_ONLY
  for(ipart=0; ipart<io.header.no_part; ipart++)
  {
    cur_part   = io.fst_part + ipart;
    
    /* particle mass */
#ifdef MULTIMASS
    weight     = (1.-f_b) * cur_part->weight;
#else
    weight     = (1.-f_b) * 1.0;
#endif
    /* coordinates in grid units [0,NGRID] */
    rrx = cur_part->pos[X] * (float)(NGRID);
    rry = cur_part->pos[Y] * (float)(NGRID);
    rrz = cur_part->pos[Z] * (float)(NGRID);
    
    /* index of nearest grid point [0,NGRID] */
    ix  = (int)(rrx+0.5);
    iy  = (int)(rry+0.5);
    iz  = (int)(rrz+0.5);
    
    /* distance to nearest grid point */
    hx  = rrx - (double)ix;
    hy  = rry - (double)iy;
    hz  = rrz - (double)iz;
    
    /* keep track of peridoc boundaries -> [0,NGRID-1] ; NGRID=0 */
    ix=(uint64_t)fmod(ix,NGRID);
    iy=(uint64_t)fmod(iy,NGRID);
    iz=(uint64_t)fmod(iz,NGRID);
    
#ifdef TSC
    /* calculate TSC weights */
    hx0=0.75 - hx*hx;
    hxp1=0.5* pow2(0.5 + hx);
    hxm1=0.5* pow2(0.5 - hx);
    hy0=0.75 - hy*hy;
    hyp1=0.5* pow2(0.5 + hy);
    hym1=0.5* pow2(0.5 - hy);
    hz0= 0.75 - hz*hz;
    hzp1=0.5* pow2(0.5 + hz);
    hzm1=0.5* pow2(0.5 - hz);
    
    /* keep track of peridoc boundaries */
    ixp1=(uint64_t)fmod(ix+1,NGRID);
    iyp1=(uint64_t)fmod(iy+1,NGRID);
    izp1=(uint64_t)fmod(iz+1,NGRID);
    ixm1=(uint64_t)fmod(ix-1+NGRID,NGRID);
    iym1=(uint64_t)fmod(iy-1+NGRID,NGRID);
    izm1=(uint64_t)fmod(iz-1+NGRID,NGRID);
    
    //     fprintf(stderr,"%ld %ld %ld   %ld %ld %ld   %ld %ld %ld\n",ix,iy,iz,ixp1,iyp1,izp1,ixm1,iym1,izm1);
    
    /* assign particle according to weights to 27 neighboring nodes */
    dens[INDX(ixm1,iym1,izm1)] = dens[INDX(ixm1,iym1,izm1)]+ hxm1*hym1 *hzm1 * weight;
    dens[INDX(ix,  iym1,izm1)] = dens[INDX(ix,  iym1,izm1)]+ hx0 *hym1 *hzm1 * weight;
    dens[INDX(ixp1,iym1,izm1)] = dens[INDX(ixp1,iym1,izm1)]+ hxp1*hym1 *hzm1 * weight;
    dens[INDX(ixm1,iy,  izm1)] = dens[INDX(ixm1,iy,  izm1)]+ hxm1*hy0  *hzm1 * weight;
    dens[INDX(ix,  iy,  izm1)] = dens[INDX(ix,  iy,  izm1)]+ hx0 *hy0  *hzm1 * weight;
    dens[INDX(ixp1,iy,  izm1)] = dens[INDX(ixp1,iy,  izm1)]+ hxp1*hy0  *hzm1 * weight;
    dens[INDX(ixm1,iyp1,izm1)] = dens[INDX(ixm1,iyp1,izm1)]+ hxm1*hyp1 *hzm1 * weight;
    dens[INDX(ix,  iyp1,izm1)] = dens[INDX(ix,  iyp1,izm1)]+ hx0 *hyp1 *hzm1 * weight;
    dens[INDX(ixp1,iyp1,izm1)] = dens[INDX(ixp1,iyp1,izm1)]+ hxp1*hyp1 *hzm1 * weight;
    dens[INDX(ixm1,iym1,iz)]   = dens[INDX(ixm1,iym1,iz)]+ hxm1*hym1 *hz0 * weight;
    dens[INDX(ix,  iym1,iz)]   = dens[INDX(ix,  iym1,iz)]+ hx0 *hym1 *hz0 * weight;
    dens[INDX(ixp1,iym1,iz)]   = dens[INDX(ixp1,iym1,iz)]+ hxp1*hym1 *hz0 * weight;
    dens[INDX(ixm1,iy,  iz)]   = dens[INDX(ixm1,iy,  iz)]+ hxm1*hy0  *hz0 * weight;
    dens[INDX(ix,  iy,  iz)]   = dens[INDX(ix,  iy,  iz)]+ hx0 *hy0  *hz0 * weight;
    dens[INDX(ixp1,iy,  iz)]   = dens[INDX(ixp1,iy,  iz)]+ hxp1*hy0  *hz0 * weight;
    dens[INDX(ixm1,iyp1,iz)]   = dens[INDX(ixm1,iyp1,iz)]+ hxm1*hyp1 *hz0 * weight;
    dens[INDX(ix,  iyp1,iz)]   = dens[INDX(ix,  iyp1,iz)]+ hx0 *hyp1 *hz0 * weight;
    dens[INDX(ixp1,iyp1,iz)]   = dens[INDX(ixp1,iyp1,iz)]+ hxp1*hyp1 *hz0 * weight;
    dens[INDX(ixm1,iym1,izp1)] = dens[INDX(ixm1,iym1,izp1)]+ hxm1*hym1 *hzp1 * weight;
    dens[INDX(ix,  iym1,izp1)] = dens[INDX(ix,  iym1,izp1)]+ hx0 *hym1 *hzp1 * weight;
    dens[INDX(ixp1,iym1,izp1)] = dens[INDX(ixp1,iym1,izp1)]+ hxp1*hym1 *hzp1 * weight;
    dens[INDX(ixm1,iy,  izp1)] = dens[INDX(ixm1,iy,  izp1)]+ hxm1*hy0  *hzp1 * weight;
    dens[INDX(ix,  iy,  izp1)] = dens[INDX(ix,  iy,  izp1)]+ hx0 *hy0  *hzp1 * weight;
    dens[INDX(ixp1,iy,  izp1)] = dens[INDX(ixp1,iy,  izp1)]+ hxp1*hy0  *hzp1 * weight;
    dens[INDX(ixm1,iyp1,izp1)] = dens[INDX(ixm1,iyp1,izp1)]+ hxm1*hyp1 *hzp1 * weight;
    dens[INDX(ix,  iyp1,izp1)] = dens[INDX(ix,  iyp1,izp1)]+ hx0 *hyp1 *hzp1 * weight;
    dens[INDX(ixp1,iyp1,izp1)] = dens[INDX(ixp1,iyp1,izp1)]+ hxp1*hyp1 *hzp1 * weight;
#else /* TSC */
    /* simple NGP mass assignment */
    dens[INDX(ix,iy,iz)] += weight;
#endif /* TSC */
  }
  
#endif /* GAS_ONLY */
}
