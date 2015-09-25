#include <math.h>
#include <stdio.h>
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
#ifdef DARK_ENERGY
#include "darkenergy.h"
#endif

/*=============================================================================
*
* This file conaints everything related to the underlying cosmological model
*
*=============================================================================*/


/*==========================================================================
 * calc_rho_crit: calculate rho_crit(a)
 *==========================================================================*/
double calc_rho_crit(double a)
{
  double rho_crit;
  
//  double a2, a3, omega0, lambda0;
//  a2      = pow2(a);
//	a3      = pow3(a);
//  omega0  = simu.omega0;
//  lambda0 = simu.lambda0;
  
  rho_crit = 3.*pow2(calc_Hubble(a))/8./PI/Grav;
  
  return (rho_crit);
}

/*==========================================================================
 * calc_rho_vir: density used for halo edge definition
 *==========================================================================*/
double calc_rho_vir(double a)
{
  if(simu.UseRhoBack == TRUE)
    return (calc_omega(a)*calc_rho_crit(a));
  else
    return (calc_rho_crit(a));
}

/*==========================================================================
* calc_omega: calculate omega(a)
*==========================================================================*/
double calc_omega(double a)
{
   double omega, omega0, lambda0;
   
#ifdef DARK_ENERGY
	omega = OmegaM_DE(a);
#else
  
   omega0  = simu.omega0;
   lambda0 = simu.lambda0;
   
   omega   = omega0 / (a + omega0*(1.-a) + lambda0*(a*a*a-a));
#endif
  
   return (omega);
}
/*==========================================================================
* calc_lambda: calculate lambda(a)
*==========================================================================*/
double calc_lambda(double a)
{
   double omega, omega0, lambda0, lambda;
   
#ifdef DARK_ENERGY
 /* Assuming a flat universe here. Otherwise we need to load in another DE user provided table */ 
   lambda = 1 - OmegaM_DE(a);
#else
   omega0  = simu.omega0;
   lambda0 = simu.lambda0;
   
   lambda  = a*a*a*lambda0 / (a + omega0*(1.-a) + lambda0*(a*a*a-a));
#endif
  
   return (lambda);
}

/*==========================================================================
* calc_Hubble: calculate H(a)
*==========================================================================*/
double calc_Hubble(double a)
{
   double Hubble;
   
#ifdef NO_EXPANSION
   Hubble = 0.0;
#else /* NO_EXPANSION */
#ifdef DARK_ENERGY
   Hubble = H0 * Hubble_DE(a); 
#else
   Hubble = H0 * sqrt(simu.lambda0 * (1.-1./pow2(a)) + simu.omega0 * (1./pow3(a)-1./pow2(a)) + 1./pow2(a));
#endif /* DARK_ENERGY */
#endif /* NO_EXPANSION */
   
   return (Hubble);
}

/*==========================================================================
* calc_virial: calculates the virial overdensity
*==========================================================================*/
double calc_virial(double a)
{
  double virial, age, omega_ta, eta, reduce, t1, t2;
  
  /* Check for a manual overdensity and return it if specified. Nothing else to do here then. */
  if(simu.UserDvir > 0)
    virial = (double)(simu.UserDvir);
  else
   {
#ifdef DARK_ENERGY
	/* Calculation of STH collapse is NOT implemented for Dark Energy models. Print a warning message, here or somewhere else */
	fprintf(stderr, "Warning: the calculation of the virial overdensity in dynamical dark energy cosmologies is NOT implemented in AHF.\n");
#else
    /* Figure out how old the universe is */
    age      = calc_t(a);
    
    /*  Figure out how much overdensity we need to reach maximum expansion by
     half the age of the universe.  The overdense expansion factor is defined
     to be 1 at maximum. */
    omega_ta = collapse(age/(double)2., (double)1.e-8);
    
    
    /* Figure out how far an object collapses to Virial equilibrium. */
    eta      = (double)2.*simu.lambda0/omega_ta; ///pow3(a); this obscure a^3 factor will prevent virial to approach the SCDM value of 178 at high redshift; not sure why it was there in the first place!
    
    if (eta < ZERO)
     {
      reduce = 0.5;
     }
    else
     {
      t1      = 2.*eta;
      t2      = -(2. + eta);
      reduce  = cubic(t1, (double)0., t2, (double)1.);
     }
    
    /* as found in M.A.K. Gross' original version */
    virial = pow3(a) * omega_ta/simu.omega0/pow3(reduce);

    /* check whether we want virial in terms of RhoCrit */
    if(simu.UseRhoBack == FALSE)
      virial *= calc_omega(a);
    
#endif
   }

  return (virial);
}
/*=========================================================================
* collapse:
c
c     Subroutine COLLAPSE performs a brute-force search using bracketing
c     to find the value of the cosmic curvature that gives turnaround at the
c     specified expansion parameter.  The expansion parameter is defined to
c     be 1 at turnaround.
c
c     Note that for a constant cosmological constant, the age at a given
c     expansion factor decreases monotonically with omegam (before turnaround;
c     after, there is more than one age at a given expansion factor).
c
c     The second argument is the desired fractional accuracy.
c
*=========================================================================*/
#define omax0 ((double)1.e10)
#define omin0 ((double)1.)
double collapse(double age0, double acc)
{
   double age, omega0_safe, omega_ta, omax, omin;
   double dtda_ta(double);
   
   /* this routine (ab-)uses COMMON variable simu.omega0 */
   omega0_safe = simu.omega0;
   
   age  = (double)-1.;
   omax = omax0;
   omin = omin0;
   
   while( fabs(age-age0) > acc*age0 && omax-omin > acc*simu.omega0)
     {
      simu.omega0 = (double)0.5*(omax+omin);
      
      age = INTEGRATE(dtda_ta, (double)0., (double)1., (double)0.1, acc);
      
      if (age > age0)
         omin = simu.omega0;
      else
         omax = simu.omega0;
      
      
     }
   
#ifdef VERBOSE
   if (omax == omax0 || omin == omin0)
     {
      fprintf(stderr,"WARNING: presumed bounds for omega are inadequate in COLLAPSE.\n");
      fprintf(stderr,"WARNING: omax=%g,omax0=%g,omin=%g,omin0=%g\n",(omax),(omax0),(omin),(omin0));
     }
#endif
  
//   fprintf(stderr, "age0 = %f, age = %f\n", age0, age);
//   fprintf(stderr, "omega0 = %f\n", simu.omega0);
//   
//   FILE *fp;
//   fp = fopen ("por.dat", "w");
//   for(simu.omega0 = 0; simu.omega0 <= 20; simu.omega0+=.1)
//   {
//     age = INTEGRATE(dtda_ta, (double)0., (double)1., (double)0.1, acc);
//     fprintf(fp, "%f  %f\n", simu.omega0, age);
//   }
  
   /* get back original value */
   omega_ta    = simu.omega0;
   simu.omega0 = omega0_safe;
   
   return (omega_ta);
}

/*==========================================================================
* dtda_ta:
c
c     Find dt/da given a, from the Friedmann Equation.  This is exact for
c     any isotropic-metric cosmology consistent with General Relativity.
c     Here, "t" is understood to be in units of the inverse Hubble constant
c     (i.e. "t" = H0*t).
c
c     For this version, define the curvature so that turnaround occurs for
c     a = 1.  Note that omegam+omegal+omegak sum to ZERO, not one, in that
c     case.
c
c     Definitions for parameters are as in Peebles 1993, eqn (5.53).
c
*==========================================================================*/
double dtda_ta(double a)
{
   if(a > ZERO)
      /* we drop the factor (H0*simu.t_unit) as for cosmological settings we use t_unit=1/H0 */
      return(1.0/sqrt(simu.lambda0*(a*a-1.) + simu.omega0/a - simu.omega0));
   else
      return(0.0);
}

/*==========================================================================
* calc_q: calculate the deceleration parameter q
*==========================================================================*/
double calc_q(double a)
{
   double omega, lambda;
   
   omega  = calc_omega(a);
   lambda = calc_lambda(a);
   
   return (omega/2. - lambda);
}

/*==========================================================================
* calc_a_dot: calculate the time derivative of a
*==========================================================================*/
double calc_a_dot(double a)
{
#ifdef NO_EXPANSION
   return (0.0);
#else
   /* we drop the factor (H0*simu.t_unit) as for cosmological settings we use t_unit=1/H0 */
   return sqrt(1.-simu.omega0 + simu.omega0/a + simu.lambda0*(a*a-1));
#endif
}

/*==========================================================================
* calc_a_dotdot: calculate the second time derivative of a
*==========================================================================*/
double calc_a_ddot(double a)
{
   double q, omega, lambda, a_ddot, a_dot;
   
#ifdef NO_EXPANSION
   return (0.0);
#else
   q      = calc_q(a);
   omega  = calc_omega(a);
   lambda = calc_lambda(a);
   a_dot  = calc_a_dot(a);
   
   /* we drop the factor (H0*simu.t_unit)^2 as for cosmological settings we use t_unit=1/H0 */
   a_ddot = pow2(a_dot)/a * (lambda-omega/2.);
   
   return (a_ddot);
#endif
}

/*==========================================================================
 * calc_t: calculate t(a) via integration
 *==========================================================================*/
double calc_t(double a)
{
   double  t;
   double  dtda(double);
   
   t = (double) INTEGRATE(dtda, (double)0.0, a, (double)0.1*a, (double)1.e-5);
   
   return(t);
}

/*==========================================================================
 * calc_super_t: calculate super_t(a) via integration
 *==========================================================================*/
double calc_super_t(double a)
{
   double  super_t;
   double  dsuper_tda(double);
   
   super_t = (double) INTEGRATE(dsuper_tda, a, (double)1.0, (double)0.1*a, (double)1.e-5);
   
   return(-super_t);
}

/*==========================================================================
 * calc_a: calculate a(t) via interpolation
 *==========================================================================*/
double calc_a(double t)
{
   int    IORD=2;
   int    j,jl,ju,jm;
   double a,da;
   double *t_array, *a_array;
   
   t_array = &(simu.timeline.t[0]);
   a_array = &(simu.timeline.a[0]);
   
   /* find correct index within timeline */
   jl = 0;
   ju = MAXTIME+1;
   while(ju-jl > 1)
     {
      jm = (ju+jl)/2;
      if(t > t_array[jm]) 
         jl=jm;
      else
         ju=jm;
     }
   j=jl;
   
   if(j >= MAXTIME-IORD)
      j = MAXTIME-IORD-1;
   
   /* interpolate using IORD points around t_array[j] and a_array[j] */
   polint(&t_array[j],&a_array[j],IORD,t,&a,&da);
   
   return(a);
}


/*==========================================================================
* calc_super_a: calculate a(super_t) via interpolation
*==========================================================================*/
double calc_super_a(double super_t)
{
   int    IORD=2;
   int    j,jl,ju,jm;
   double a,da;
   double *t_array, *a_array;
   
#ifdef NO_EXPANSION
   
   return (1.0);
   
#else /* NO_EXPANSION */
   
   t_array = &(simu.timeline.super_t[0]);
   a_array = &(simu.timeline.a[0]);
   
   /* find correct index within timeline */
   jl = 0;
   ju = MAXTIME+1;
   
   while(ju-jl > 1)
     {
      jm = (ju+jl)/2;
      if(super_t > t_array[jm]) 
         jl=jm;
      else
         ju=jm;
     }
   j=jl;
   
   if(j >= MAXTIME-IORD)
      j = MAXTIME-IORD-1;
   
   /* interpolate using IORD points around t_array[j] and a_array[j] */
   polint(&t_array[j],&a_array[j],IORD,super_t,&a,&da);
   
   return(a);
#endif
}


/*===========================================================================
 * da/dt  as given by Friedman equation
 c
 c     Find dt/da given a, from the Friedmann Equation.  This is exact for
 c     any isotropic-metric cosmology consistent with General Relativity.
 c     Here, "t" is understood to be in units of the inverse Hubble constant
 c     (i.e. "t" = H0*t).
 c
 c     Definitions for parameters are as in Peebles 1993, eqn (5.53).
 c
 *===========================================================================*/
double dtda(double a)
{
   if(a > ZERO)
      /* we drop the factor (H0*simu.t_unit) as for cosmological settings we use t_unit=1/H0 */
      return(1.0/sqrt(simu.lambda0*(pow2(a)-1.) + simu.omega0/a - simu.omega0 + 1.));
   else
      return(0.0);
}


/*===========================================================================
 * da/dsuper_t  as given by supercomoving Friedman equation
 *
 *===========================================================================*/
double dsuper_tda(double a)
{
   if(a > ZERO)
      return( dtda(a)/pow2(a) );
   else
      return(0.0);
}


/*==========================================================================
 * calc_growth: calculates the linear growth factor (without Omega_r)
 *==========================================================================*/
double dtda3(double a)
{
   double t1;
   
   t1 = dtda(a);
   
   return(pow3(t1));
}
double calc_growth(double a)
{
   double t1, t2;
   
   t1 = INTEGRATE(dtda3, 0.0, a, a/10., 1e-8);
   t2 = calc_a_dot(a)/a;
   
   return(2.5 * simu.omega0 * t1 * t2);
}

/*==========================================================================
 * calc_growthr: calculates the linear growth factor
 * with Omega_r=9.2364e-5 according to Planck:
 * http://physics.stackexchange.com/questions/94181/where-is-radiation-density-in-the-planck-2013-results
 *==========================================================================*/
double adot(double a)
{
  double omegar0 = 9.2364e-5;
  
  return(sqrt((simu.omega0/a+simu.lambda0*pow2(a)+omegar0/pow2(a))));
}

double dtda3r(double a)
{
  double t1;
  
  t1 = 1./adot(a);
  
  return(pow3(t1));
}

double calc_growthr(double a)
{
  double t1, t2;
  
  t1 = INTEGRATE(dtda3r, 0.0, a, a/10., 1e-8);
  t2 = adot(a)/a;
  
  return(2.5 * simu.omega0 * t1 * t2);
}

/*===========================================================================
* create a timeline for current cosmological model
*===========================================================================*/
void create_timeline(double a_init, double a_fin, tlptr timeline)
{
   int    iloop, iwrite_cosmo;
   double a,t,omega,lambda,rhoc,hubble,super_t,virial,growth,age,dummy,zred;
   
   FILE  *fpcosmo;
   char   cosmo_file[MAXSTRING], string[MAXSTRING];


#ifdef COSMOLOGY
   if(global_io.params->outfile_prefix != NULL)
     {
      strcpy(cosmo_file, global_io.params->outfile_prefix);
      strcat(cosmo_file,"cosmology");
      /*---------------------------------------------
         * check if there is already a .cosmology file
         *---------------------------------------------*/
      if( (fpcosmo = fopen(cosmo_file,"r")) == NULL )
        {
         iwrite_cosmo = TRUE;
        }
      else
        {
         iwrite_cosmo = FALSE;
         fclose(fpcosmo);
        }
     }
   else
      iwrite_cosmo = FALSE;

#else
   /*---------------------------------------------
    * do not write the .cosmology file whatsoever
    *---------------------------------------------*/
   iwrite_cosmo = FALSE;
#endif /* COSMOLOGY */

   /*------------------------------
    * create timeline from scratch
    *------------------------------*/
#ifndef WITH_MPI
   if(io.logfile != NULL)
     {
      fprintf(io.logfile,"\ncreating timeline from a=%g to a=%g (MAXTIME=%d) ... ",a_init,a_fin,MAXTIME);
      fflush(io.logfile);
     }
#endif
      
#ifdef WITH_OPENMP
#pragma omp parallel private(iloop, a, t, super_t, omega, lambda, hubble, rhoc) shared(timeline, simu, a_fin, a_init)
#pragma omp for schedule(static)
#endif
   for(iloop=0; iloop<MAXTIME; iloop++)
     {
      a       = ((double)iloop)/(double)(MAXTIME-1) * (a_fin-a_init) + a_init;
      t       = calc_t(a);
      super_t = calc_super_t(a);
      omega   = calc_omega(a);
      lambda  = calc_lambda(a);
      hubble  = calc_Hubble(a);
      rhoc    = rhoc0 *     (
                             simu.lambda0*(1.-1./pow2(a))           +
                             simu.omega0*(1./pow3(a)-1./pow2(a))    + 1./pow2(a));
      
      timeline->a[iloop]       = a;
      timeline->t[iloop]       = t;
      timeline->super_t[iloop] = super_t;      
      timeline->omega[iloop]   = omega;
      timeline->lambda[iloop]  = lambda;
      timeline->hubble[iloop]  = hubble;
      timeline->rhoc[iloop]    = rhoc;
      timeline->age[iloop]     = t * simu.t_unit*Mpc/Gyr;
      timeline->growth[iloop]  = calc_growth(a);
      
      // do >>not<< call calc_virial() here when using OpenMP as calc_virial() abuses simu.omega0!!!
      
     }
   
   if(iwrite_cosmo)
     {
     fpcosmo = fopen(cosmo_file,"w");
     fprintf(fpcosmo,"#     z(1)           a(2)    t[h^-1 Gyr](3)  Omega(4)    lambda(5)   hubble(6)   RhoCrit(7)     virial(8)      growth(9)        q(10)   super_t(11)\n");
     
     for(iloop=0; iloop<MAXTIME; iloop++)
        {
        /* calc_virial() temporarily changes simu.omega0 */
        timeline->virial[iloop] = calc_virial(timeline->a[iloop]);
          
        fprintf(fpcosmo,"%10.4f %16.10f %12.6f %12.6f %12.6f %12.4f %12.6g %12.6g %12.6g %12.6g  %12.6g\n",
                1.0/timeline->a[iloop]-1,
                timeline->a[iloop],
                timeline->t[iloop] * simu.t_unit*Mpc/Gyr,
                timeline->omega[iloop],
                timeline->lambda[iloop],
                timeline->hubble[iloop],
                timeline->rhoc[iloop],
                timeline->virial[iloop],
                timeline->growth[iloop],
                calc_q(timeline->a[iloop]),
                timeline->super_t[iloop]);
        fflush(fpcosmo);        
        }

      fclose(fpcosmo);
     }
     
#ifndef WITH_MPI
   if(io.logfile != NULL)
     {
      fprintf(io.logfile,"done\n\n");
      fflush(io.logfile);      
     }
#endif
   
#ifdef DEBUG_SUPERCOMOVING
   {
      double a_prev, t_prev, super_t_prev, a_mid;
      double dtc, dTc;
      
      iloop        = 0;
      a_prev       = ((double)iloop)/(double)(MAXTIME-1) * (a_fin-a_init) + a_init;
      t_prev       = calc_t(a_prev);
      super_t_prev = calc_super_t(a_prev);

      for(iloop=1; iloop<MAXTIME; iloop++)
        {
         a       = ((double)iloop)/(double)(MAXTIME-1) * (a_fin-a_init) + a_init;
         t       = calc_t(a);
         super_t = calc_super_t(a);
         
         dtc     = t       - t_prev;
         dTc     = super_t - super_t_prev;
         
         a_mid   = (a+a_prev)/2;
         
         fprintf(stderr,"z=%12.4f  a=%20.10f   t=%20.10f  super_t=%20.10f  dtc=%20.10f   dTc=%20.10f    dTc(dtc)=%20.10f  a=%20.10f     err=%g\n",
                 1./a-1.,a,t,super_t,dtc,dTc,dtc/(pow2(a_mid)),calc_super_a(super_t),fabs(dTc-dtc/(pow2(a_mid))));
         
         a_prev       = a;
         t_prev       = t;
         super_t_prev = super_t;
        }
   }
#endif

#ifdef COSMOLOGY_TERM
   exit(0);
#endif   

}

/* a quick-and-dirty hack for the VDE model which requires H in a table */
double calc_Hubble_VDE(double a)
{
  FILE   *fhubble;
  double *a_array, *h_array;
  double  h, dh;
  int     n_array;
  char    dummyline[MAXSTRING];
  int     jl,ju,jm,j,IORD=2;
  
  fhubble = fopen("hubble_vde.dat","r");
  n_array = -1;
  a_array = calloc(1,sizeof(double));
  h_array = calloc(1,sizeof(double));
  while(fgets(dummyline,MAXSTRING,fhubble) != NULL)
   {
    n_array++;
    a_array = realloc(a_array, (n_array+1)*sizeof(double));
    h_array = realloc(h_array, (n_array+1)*sizeof(double));
    sscanf(dummyline,"%lf %lf", &(a_array[n_array]), &(h_array[n_array]));
    
    /* possible conversions of what has just been read in */
    a_array[n_array]  = 1./(1.+a_array[n_array]);
    h_array[n_array] *= H0;
   }
  
  /* interpolate using IORD points around a_array[j] and h_array[j] */
  jl = 0;
  ju = n_array+1;
  while(ju-jl > 1)
   {
    jm = (ju+jl)/2;
    if(a < a_array[jm]) 
      jl=jm;
    else
      ju=jm;
   }
  j=jl;
  if(j >= n_array-IORD)
    j = n_array-IORD-1;
  polint(&a_array[j],&h_array[j],IORD,a,&h,&dh);
  
  //fprintf(stderr,"%d %g %g\n",j,a,a_array[j],h_array[j]);

  fclose(fhubble);
  
  return(h);
}
