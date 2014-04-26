#include <math.h>
#include <stdio.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"

/*
 * we need to cut-and-paste the create_timeline() routine from ../src/cosmology.c
 * in order to be able to bypass the #ifdef COSMOLOGY statements...
 */
void create_timeline2(double a_init, double a_fin, tlptr timeline);
#define MAXTIME2 10000

int main()
{
   int     no_outputs, i;
   double  zred, z_init, z_final, a_init, a_final, a;
   
   printf("============================\n");
   printf(" calculate age of Universe\n");
   printf("============================\n");
   printf("please give omega0:       ");
   scanf("%lf",&simu.omega0);
   printf("please give lambda0:      ");
   scanf("%lf",&simu.lambda0);
   printf("\n");
   
   z_init  = 500.;
   z_final = 0.;
   a_init  = 1.0/(1.0+z_init);
   a_final = 1.0/(1.0+z_final);
   
   create_timeline2(a_init, a_final, &simu.timeline);
}

/*===========================================================================
* create a timeline for current cosmological model
*===========================================================================*/
void create_timeline2(double a_init, double a_fin, tlptr timeline)
{
   int    iloop;
   double a,t,omega,lambda,rhoc,hubble;
   
   FILE *fpout;
   fpout = fopen("Cosmology.DAT","w");
   fprintf(fpout,"#      z             a       t[h^-1 Gyr]     Omega       lambda      hubble      RhoCrit        virial        growth        q\n");
   
   for(iloop=0; iloop<MAXTIME2; iloop++)
     {
      a      = ((double)iloop+1.0)/(double)MAXTIME2 * (a_fin-a_init) + a_init;
      t      = calc_t(a);
      omega  = calc_omega(a);
      lambda = calc_lambda(a);
      hubble = calc_Hubble(a);
      rhoc   = 3*pow2(hubble)/8./PI/Grav;
      
      timeline->a[iloop]      = a;
      timeline->t[iloop]      = t;
      timeline->omega[iloop]  = omega;
      timeline->lambda[iloop] = lambda;
      timeline->hubble[iloop] = hubble;
      timeline->age[iloop]    = t*Mpc/H0/Gyr;
      timeline->virial[iloop] = calc_virial(a);
      timeline->growth[iloop] = calc_growth(a);
      
      //fprintf(stderr,"%16.8g %16.8g %16.8g   %16.8g %16.8g    %16.8g\n",1.0/a-1,rhoc,3*pow2(hubble)/8./PI/Grav,hubble,calc_Hubble(a),calc_Hubble_VDE(a));
      
      fprintf(fpout,"%12.4f %12.4f %12.6f %12.6f %12.6f %12.4f %12.6g %12.6g %12.6g %12.6g\n",
              1.0/a-1,a,t*Mpc/H0/Gyr,omega,lambda,hubble,rhoc,timeline->virial[iloop],timeline->growth[iloop],calc_q(a));
      fflush(fpout);
     }
}