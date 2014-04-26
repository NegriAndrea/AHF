#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "../src/param.h"

int main()
{
   double npart, omega0, boxsize, mean_dens, inv_mean_dens, pmass;

   printf("=========================================\n");
   printf(" calculate mass of a simulation particle\n");
   printf("=========================================\n");
   printf("please total number of particles:  ");
   scanf("%lf",&npart);
   printf("please give boxsize [Mpc/h]:       ");
   scanf("%lf",&boxsize);
   printf("please give omega0:                ");
   scanf("%lf",&omega0);
   
   if(npart < 32768.)
      npart = npart*npart*npart;

   inv_mean_dens = boxsize*boxsize*boxsize/npart;
   pmass         = omega0*rhoc0*inv_mean_dens;
   
   printf("\npmass = %g [Msun/h] (total mass in box=%g)\n",pmass,omega0*rhoc0*boxsize*boxsize*boxsize);
}
