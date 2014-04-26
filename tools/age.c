#include <math.h>
#include <stdio.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"


/*
 * we need to cut-and-paste the create_timeline() routine from ../src/cosmology.c
 * in order to be able to bypass the #MAXTIME definition...
 */
void create_timeline2(double a_init, double a_fin, tlptr timeline);
#define MAXTIME2 200


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
   printf("please give redshift z:   ");
   scanf("%lf",&zred);
   printf("\n");
   
   z_init  = 100.;
   z_final = 0.;
   a_init  = 1.0/(1.0+z_init);
   a_final = 1.0/(1.0+z_final);
   
   create_timeline2(a_init, a_final, &simu.timeline);
   
   
   a = 1./(1.+zred);
   printf("age = %lf/h Gyr\n",calc_t(a)*Mpc/H0/Gyr);
}

/*===========================================================================
* create a timeline for current cosmological model
*===========================================================================*/
void create_timeline2(double a_init, double a_fin, tlptr timeline)
{
   int    iloop;
   double a,t,omega,lambda,rhoc,hubble;
   
   for(iloop=0; iloop<MAXTIME2; iloop++)
     {
      a      = ((double)iloop+1.0)/(double)MAXTIME2 * (a_fin-a_init) + a_init;
      t      = calc_t(a);
      omega  = calc_omega(a);
      lambda = calc_lambda(a);
      hubble = H0    * sqrt(
                            simu.lambda0 * (1.-1./pow2(a))         +
                            simu.omega0  * (1./pow3(a)-1./pow2(a)) + 1./pow2(a));
      rhoc   = rhoc0 *     (
                            simu.lambda0*(1.-1./pow2(a))           +
                            simu.omega0*(1./pow3(a)-1./pow2(a))    + 1./pow2(a));
      
      timeline->a[iloop]      = a;
      timeline->t[iloop]      = t;
      timeline->omega[iloop]  = omega;
      timeline->lambda[iloop] = lambda;
      timeline->hubble[iloop] = hubble;
      timeline->age[iloop]    = t*Mpc/H0/Gyr;
      timeline->virial[iloop] = calc_virial(a);
     }
}