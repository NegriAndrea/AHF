#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/***********************************************************************
 *                          STRUCTURES  ETC.
 ***********************************************************************/
#define NDIM      3
#define X         0
#define Y         1
#define Z         2
#define ZERO      1e-6
#define MACHINE_ZERO 5e-16

#define AHF_MINPART 20
#define AHF_VTUNE   1.5

#define MULTIMASS
#define VERBOSE
#define DEBUG

#define pow2(x)         ((x)*(x))
#define SWAP(a,b,temp)  temp=(a);(a)=(b);(b)=temp;

//simulation particulars
#define atime  0.998100
#define H0     100.069

//halo particulars
#define Nhalo  203
#define Xhalo  23.275803
#define Yhalo  33.153991
#define Zhalo  27.488322


typedef struct particle *partptr;
typedef struct particle
{
   float   pos[NDIM];
   float   mom[NDIM];
   float   weight;
} part;  

typedef struct {
 double x,y,z;
} XYZ;

typedef struct {   
 long unsigned  npart;
 long unsigned *ipart;
 
 XYZ     pos;
 XYZ     vel;
 double  M_vir;
 double  R_vir;
 double  Phi0;
} HALO;


typedef struct {
   double boxsize;
   double pmass;
   double t_unit;
   double a;
   double Hubble;
   double Grav;
} info_simu;

/***********************************************************************
 *                          COMMON VARIABLES
 ***********************************************************************/
HALO      halo;
info_simu simu;
double    r_fac, x_fac, v_fac, m_fac, phi_fac;
partptr   fst_part;


/***********************************************************************
 *                              FUNCTIONS
 ***********************************************************************/
void read_halo     (char*);
void write_halo    (char*);
void rem_unbound   ();
void indexx        (unsigned long n, double arr[], unsigned long indx[]);


/***********************************************************************
 *                                MAIN
 ***********************************************************************/
int main(argc,argv)
int argc;
char **argv;
{  
   char infile[2048];
   
   if(argc<2)
     {
      fprintf(stderr,"usage: %s halofile\n",*argv);
      exit(1);
     }
   strcpy(infile,argv[1]);

   /* initialize simulation */
   simu.boxsize = 64.;
   simu.pmass   = 1678784.125000;
   simu.t_unit  = 0.01;
   simu.Grav    = 4.3006485e-9;
 
   simu.a       = atime;
   simu.Hubble  = H0;
   
   halo.npart = Nhalo;
   halo.pos.x = Xhalo / simu.boxsize;
   halo.pos.y = Yhalo / simu.boxsize;
   halo.pos.z = Zhalo / simu.boxsize;
      

   /* calculate conversion factors */
   x_fac   = simu.boxsize;
   r_fac   = simu.boxsize*simu.a;
   v_fac   = simu.boxsize/simu.t_unit/simu.a;
   m_fac   = simu.pmass;
   phi_fac = simu.Grav*simu.pmass/(simu.boxsize*simu.a);

   /* read particles from file */
   read_halo(infile);
   
   rem_unbound();
   
   write_halo(infile);
}


/***********************************************************************
 *                            write_halo()
 ***********************************************************************/
void write_halo(char *infile)
{
   char          outfile[2048];
   FILE         *fpout;
   long unsigned jpart;
   partptr       cur_part;
   
   strcpy(outfile,infile);
   strcat(outfile,"-remunbound");
   fpout = fopen(outfile,"w");
   
   for(jpart=0; jpart<halo.npart; jpart++)
     {
      cur_part = fst_part + halo.ipart[jpart];
      
      fprintf(fpout,"%g %g %g    %g %g %g    %g\n",
              (cur_part->pos[X]), (cur_part->pos[Y]), (cur_part->pos[Z]),
              (cur_part->mom[X]), (cur_part->mom[Y]), (cur_part->mom[Z]),
              (cur_part->weight));
     }
   fclose(fpout);
}

/***********************************************************************
 *                            read_halo()
 ***********************************************************************/
void read_halo(char *infile)
{
   FILE         *fpin;
   long unsigned jpart;
   partptr       cur_part;
   double        Xp, Yp, Zp, Vxp, Vyp, Vzp, Mp;
   
   fst_part   = (partptr)         calloc(halo.npart, sizeof(part));
   halo.ipart = (long unsigned *) calloc(halo.npart, sizeof(long unsigned));
   
   fpin = fopen(infile,"r");
   for(jpart=0; jpart<halo.npart; jpart++)
     {
      cur_part          = fst_part + jpart;
      halo.ipart[jpart] = cur_part - fst_part;
      
      fscanf(fpin,"%lf %lf %lf %lf %lf %lf %lf",&Xp,&Yp,&Zp,&Vxp,&Vyp,&Vzp,&Mp);
      
      cur_part->pos[X] = Xp;
      cur_part->pos[Y] = Yp;
      cur_part->pos[Z] = Zp;
      cur_part->mom[X] = Vxp;
      cur_part->mom[Y] = Vyp;
      cur_part->mom[Z] = Vzp;
      cur_part->weight = Mp;      
     }
   fclose(fpin);
}

/***********************************************************************
 *                            rem_unbound()
 ***********************************************************************/
void rem_unbound()
{
   
   long unsigned  nremove, npart_old;
   double         Xc, Yc, Zc, VXc, VYc, VZc;
   double         Xp, Yp, Zp, VXp, VYp, VZp;
   double         dX, dY, dZ, dVX, dVY, dVZ;
   double         v2_tune, weight;
   double         Phi, Phi0, Phi_infty, M_r, M_vir, M_vel, vel2, v_esc2, d_prev, R_vir;
   double         I_now, I_mid, I_prev, dr;
   partptr        cur_part, pre_part, host_part, tmp_part;
   double         x,y,z,dx,dy,dz,distA,distB,dist;
   double         xx,yy,zz;
   int            niter;
   int 		      hostHalo;
   double        *mom2;
   long           jpart;
   long unsigned *bound_ipart, bound_npart, *idx, no_vbulk;
   
   
#ifdef VERBOSE
   fprintf(stderr,"    rem_unbound:      npart=%12ld -> ",halo.npart);
   fflush(stderr);
#endif
   
   /* velocity tune parameter */
   v2_tune = pow2(AHF_VTUNE);
   
  /* how many central particles to use for the initial bulk velocity guess */
  no_vbulk = (long unsigned) (AHF_MINPART/2); 

  /* remember initial number of particles */
   /*  npart_old = halos[i].npart; */
   
   /* halo centre in AMIGA units */
   Xc = halo.pos.x;
   Yc = halo.pos.y;
   Zc = halo.pos.z;
   
   
   nremove = 4;
   niter   = 0;
   /*----------------------------------------------------------------------------
      * remove particles until all particles are bound or halo mass gets too small
      *----------------------------------------------------------------------------*/
   while( (nremove > 3))
     {      
      /* iteration counter */
      niter++;
      
      /*---------------------------------------------------------------
      * determine Phi0  (the zero point of the potential is infinity)
      *---------------------------------------------------------------*/
      
      /* reset values */
      I_prev  = 0.0;
      d_prev  = 0.0;
      M_r     = 0.0;
      Phi0    = 0.0;
      
      /************************************************************/
      /* loop over all sorted particles from inside out */
      /* calculate Phi0 */
      for(jpart=0; jpart<halo.npart; jpart++)
        {
         /* access particle */
         cur_part = fst_part + halo.ipart[jpart];
         
#ifdef MULTIMASS
         weight = (double)cur_part->weight;
#else
         weight = (double)1.0;
#endif
         
         /* cumulative mass */
         M_r += weight;
         
         /* particle position */
         Xp  = (double)cur_part->pos[X];
         Yp  = (double)cur_part->pos[Y];
         Zp  = (double)cur_part->pos[Z];
         
         /* put particle into halo rest frame */
         dX  = fabs(Xp - Xc);
         dY  = fabs(Yp - Yc);
         dZ  = fabs(Zp - Zc);
         
         /* take care of periodic boundary conditions */
         if(dX >  0.5) dX -= 1.0;
         if(dY >  0.5) dY -= 1.0;
         if(dZ >  0.5) dZ -= 1.0;
         
         /* finally calculate distance (conversion to dist_phys=a*dist via phi_fac!)  */
         dist = sqrt( pow2(dX) + pow2(dY) + pow2(dZ) );
      
         /* accumulate Phi0 */
         if( dist > MACHINE_ZERO)
           {
            
            /* mid-point integration */
            I_now = M_r/pow2(dist);
            I_mid = (I_now + I_prev)/2.;
            dr    =  dist - d_prev;
            
            /* accumulate Phi0 */
            Phi0 += I_mid * dr;
           }
         
         d_prev   = dist;
         I_prev   = I_now;
         
        } /* particle loop for Phi0 determination */
      
      /* finally calculate Phi0 */
      Phi0 += M_r/dist;
      
      /* remember Phi0 as it will be used when calculating halos[].prof.Epot */
      halo.Phi0 = Phi0;

#ifdef DEBUG
      fprintf(stderr,"\nPhi0 = %g ", Phi0*m_fac*phi_fac);
#endif

      /*--------------------------------------------------------------------
       * determine Phi             ( v_esc^2 = 2 * |Phi| )
       *--------------------------------------------------------------------*/
      /* reset values */
      nremove = 0;
      I_prev  = 0.0;
      d_prev  = 0.0;
      M_r     = 0.0;
      M_vir   = 0.0;
      Phi     = 0.0;

      /* initial bulk velocity of halo */
       if(niter == 1)
         {
           mom2 = (double *)        calloc(no_vbulk + 1, sizeof(double));
           idx  = (long unsigned *) calloc(no_vbulk + 1, sizeof(long unsigned));
           
           for(jpart=0; jpart<no_vbulk; jpart++)
             {
               cur_part       = fst_part + halo.ipart[jpart];
               mom2[jpart+1]  = pow2(cur_part->mom[X])+pow2(cur_part->mom[Y])+pow2(cur_part->mom[Z]);
             }
           indexx((long)(no_vbulk), mom2, idx);
           
           /* use the median of the innermost particle's velocity */
           jpart = idx[(int)(no_vbulk/2)] - 1;
           
           free(mom2);
           free(idx);
         }
       else
         {
           /* use the most bound particle's velocity */
           jpart = 0;
         }
       
       fprintf(stderr,"niter = %d  jpart = %ld\n", niter, jpart);       
       
      cur_part = fst_part + halo.ipart[jpart];
#ifdef MULTIMASS
      weight   = (double)cur_part->weight;
#else
      weight   = (double)1.0;
#endif
      M_vel    = weight;
      VXc      = weight*cur_part->mom[X];
      VYc      = weight*cur_part->mom[Y];
      VZc      = weight*cur_part->mom[Z];
       
      /************************************************************/
      /* loop over all sorted particles from inside out */
      bound_npart = 0;
      bound_ipart = (long unsigned *) calloc(1, sizeof(long unsigned));    // some realloc()'s do not like NULL pointers...
      for(jpart=0; jpart<halo.npart; jpart++)
        {
         /* access particle */
         cur_part = fst_part + halo.ipart[jpart];
         
#ifdef MULTIMASS
         weight = (double)cur_part->weight;
#else
         weight = (double)1.0;
#endif
         /* cumulative mass */
         M_r += weight;
         
         /* particle position */
         Xp  = (double)cur_part->pos[X];
         Yp  = (double)cur_part->pos[Y];
         Zp  = (double)cur_part->pos[Z];
         
         /* put particle into halo rest frame :: *no* fabs() this time ! */
         dX  = (Xp - Xc);
         dY  = (Yp - Yc);
         dZ  = (Zp - Zc);
         
         /* take care of periodic boundary conditions */
         if(dX >  0.5) dX -= 1.0;
         if(dY >  0.5) dY -= 1.0;
         if(dZ >  0.5) dZ -= 1.0;
         if(dX < -0.5) dX += 1.0;
         if(dY < -0.5) dY += 1.0;
         if(dZ < -0.5) dZ += 1.0;
         
         /* finally calculate distance (conversion to dist_phys=a*dist via phi_fac!) */
         dist = sqrt( pow2(dX) + pow2(dY) + pow2(dZ) );
         
         /* get potential escape velocity */
         if( dist > MACHINE_ZERO )
           {  	  
            /* mid-point integration */
            I_now = M_r/pow2(dist);
            I_mid = (I_now + I_prev)/2.;
            dr    =  dist - d_prev;
            
            /* accumulate potential */
            Phi += I_mid * dr;
            
            /* get escape velocity */
            v_esc2 = (2*fabs(Phi-Phi0)*phi_fac);            
           }
         
         /* potential="inf" for dist=0 -> set v_esc manually */
         else
           {
            v_esc2 = 1e30;
           }
         
         /* get particle velocity in halo rest frame velocity plus Hubble flow */
         VXp = (double)cur_part->mom[X];
         VYp = (double)cur_part->mom[Y];
         VZp = (double)cur_part->mom[Z];
         
         dVX = (VXp - VXc/M_vel) * v_fac + simu.Hubble * dX * r_fac;
         dVY = (VYp - VYc/M_vel) * v_fac + simu.Hubble * dY * r_fac;
         dVZ = (VZp - VZc/M_vel) * v_fac + simu.Hubble * dZ * r_fac;
         
         /* absolute velocity of current particle */
         vel2  = pow2(dVX) + pow2(dVY) + pow2(dVZ);
         
#ifdef DEBUG
         fprintf(stderr,"\ndist=%16.8g   vel=%16.8g   v_esc=%16.8g    v_bulk=%16.8g    v_raw=%16.8g",
                 dist*x_fac*1000.,sqrt(vel2),sqrt(v_esc2),
                 sqrt(pow2(VXc/M_r)+pow2(VYc/M_r)+pow2(VZc/M_r))* v_fac,
                 sqrt(pow2(VXp)+pow2(VYp)+pow2(VZp))* v_fac);
#endif
         
         /* unbound particle? */
         if( vel2 > v2_tune*v_esc2)
           {
            /* count number of particles to be removed */
            nremove++;
           }         
         else 
           {
            /* let particle contribute to (interior) bulk velocity */
            M_vel += weight;
            VXc   += weight*cur_part->mom[X];
            VYc   += weight*cur_part->mom[Y];
            VZc   += weight*cur_part->mom[Z];
               
            /* store bound particles temporarily in bound_ipart[] */
            bound_npart++;
            bound_ipart = (long unsigned *) realloc(bound_ipart, bound_npart*sizeof(long unsigned));
            bound_ipart[bound_npart-1] = cur_part - fst_part;
            
            /* accumulate virial values (NOTE: for v_esc2 we are accumulating M_r and not M_vir!) */
            M_vir += weight;
            R_vir  = dist;
           }
         
         I_prev   = I_now;
         d_prev   = dist;
         
        } /* particle loop */
    
      /* double-check new number of bound particles */
      if(bound_npart != (halo.npart - nremove))
         fprintf(stderr,"rem_unbound: better check the unbinding procedure! bound_part=%ld vs. halo.npart-nremove=%ld\n",
                 bound_npart, halo.npart-nremove);
         
         
      /* update number of particles in halo */
      free(halo.ipart);
      halo.ipart = bound_ipart;
      halo.npart = bound_npart;
      halo.M_vir = M_vir;
      halo.R_vir = R_vir;

       nremove = 0;
       
       
     } /* while( (nremove > n)  ) */
   
#ifdef VERBOSE
         fprintf(stderr,"%12ld (Rvir=%g kpc/h)\n",halo.npart, halo.R_vir*x_fac*1000.);
         fflush(stderr);
#endif   
}





#define NR_END 1
#define M 7
#define NSTACK 50
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
double *vector(long nl, long nh)
/* allocate a vector with subscript range v[nl..nh] */
{
  double *v;
  
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
}
void free_vector(double *v, long nl, long nh)
/* free a vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
}
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void indexx(unsigned long n, double arr[], unsigned long indx[])
{
  unsigned long i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  double a;
  
  istack=ivector(1,NSTACK);
  for (j=1;j<=n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        indxt=indx[j];
        a=arr[indxt];
        for (i=j-1;i>=1;i--) {
          if (arr[indx[i]] <= a) break;
          indx[i+1]=indx[i];
        }
        indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1],itemp);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP(indx[l+1],indx[ir],itemp)
      }
      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l],indx[ir],itemp)
      }
      if (arr[indx[l+1]] > arr[indx[l]]) {
        SWAP(indx[l+1],indx[l],itemp)
      }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
        do i++; while (arr[indx[i]] < a);
        do j--; while (arr[indx[j]] > a);
        if (j < i) break;
        SWAP(indx[i],indx[j],itemp)
      }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
  free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */
