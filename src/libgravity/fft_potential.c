#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "gravity.h"
#include "../libutility/utility.h"


/*==============================================================================
* CALC_POTENTIAL: 
*
*    rho^ = FFT(rho)
*    Phi^ = G^ * rho^ ;     G^ = -1/k**2 (no adjustment for finite differences)
*    Phi  = FFT(Phi^) 
*
*    uses NR routine fourn.c for FFT's
*==============================================================================*/
void fft_potential(flouble *data, long l1dim)
{
   FILE         *fpout;
   double        Green, sf, sf2;
   double        FFTnorm;
   double        kx, ky, kz, ksquare, kf;
   double        sin2sfkz, sin2sfky;
   long          *k;
   long          k_dim, k_nyq;
   long          i, l, m, n;
   int           forward=1, backward=-1;
   unsigned long nn[NDIM];
   
   /* dimensions of grid */
   nn[0] = l1dim;
   nn[1] = l1dim; 
   nn[2] = l1dim;
   k_dim = l1dim;
   k_nyq = l1dim/2;
   
   /* this will help extracting information from data^ */
   k = (long*) calloc(l1dim,sizeof(long));
   for(i=0; i<l1dim; i++)
     {
      if(i <= k_nyq)
         k[i] = i;
      else
         k[i] = -(l1dim-i);
     }
   
   /* fourn normalisation factor */
   FFTnorm = (double)pow3(l1dim);
   
   /* fundamental wave in box units: kf = 2*pi */
   kf  = TWOPI;
   
   /* adjusting the 7-point-difference Green's function */
   sf  = PI/(double)l1dim/kf;
   sf2 = pow2(sf);
   
   
   /* forward FFT */
   fourn(data-1, nn-1, NDIM, forward);
   
  PkSpectrum.time -= time(NULL);
  PkSpectrum.time += time(NULL);
   
   /* calculation of gravitational potential */
   for(n = 0; n < k_dim; n++) 
     {
      kz       = kf * k[n];
      sin2sfkz = pow2(sin(sf*kz));
      
      for(m = 0; m < k_dim; m++)
        {
         ky       = kf * k[m];
         sin2sfky = pow2(sin(sf*ky));
         
         for(l = 0; l < k_dim; l++)
           {
            kx = kf * k[l];
            
            /* get k^2 => NOT needed for 7-point difference approximation to Green = -1.0 / ksquare */
            /* ksquare = kx * kx + ky * ky + kz * kz; */
            
            /* Green's function */
            if (k[l] == 0 &&  k[m]  == 0 && k[n] == 0)
              {
               Green = 0.0;
              }
            else
              {
               /* 7-point difference approximation to Green = -1.0 / ksquare */
               Green = -sf2/(pow2(sin(sf*kx))+sin2sfky+sin2sfkz);
              }
            
            data[Re(l,m,n,l1dim)] *= Green;
            data[Im(l,m,n,l1dim)] *= Green;
            
           }
        }
     }
   
   
   /* backward FFT */
   fourn(data-1, nn-1, NDIM, backward);
   
   /* normalize FFT output */
   for(i=0; i<2*l1dim*l1dim*l1dim; i++)
      data[i] /= FFTnorm;
   
   /* delete k-array again */
   free(k);
}
