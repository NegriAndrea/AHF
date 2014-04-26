#include <stddef.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

/*==========================================================================================
 * residual:
 *
 * calculate the norm of all residuals on a given grid
 *
 * the options for that norm are:
 *
 *  1. SQNM  ->  sqrt(pow2(cur_residual))
 *  2. ----  ->  |cur_residual|
 *
 *==========================================================================================*/
double residual(gridls *cur_grid)
{
#ifndef AHFlean
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   nptr    tsc_nodes[3][3][3];
   long    x, y, z, no_nodes, ipquad;
   double  FourPiGa;
   double  cur_source, cur_residual, Lphi;
   double  residual;
   

   /* conversion factor to go from density to source term */
   FourPiGa = simu.FourPiG*calc_super_a(cur_grid->timecounter);
   
   /* reset values */
   residual  = 0.0;
   no_nodes  = 0;

   
   /* loop over all nodes */
#ifdef WITH_OPENMP
#pragma omp parallel reduction(+:residual) reduction(+:no_nodes) firstprivate(FourPiGa) private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, tsc_nodes, x, y, z, cur_source, cur_residual, Lphi) shared(cur_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)
     {
       cur_pquad = cur_grid->pquad_array[ipquad];
#else
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
     {
#endif
      for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, z++)
        {
         for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
           {
            for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++)
              {
               for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
                 {
                  for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                     /* convert density to source term */
                     cur_source   = cur_node->dens * FourPiGa;
                                          
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                     
                     /* only use nodes in the interior of the refinement */
                     if(test_tsc(tsc_nodes))
                       {
                        /* count number of nodes used for residual calculation */
                        no_nodes++;
                                                
                        /* obtain residual in cur_node */
                        Lphi         = Laplace_pot(tsc_nodes, cur_grid->spacing2);
                        cur_residual = (double)cur_source - Lphi;
                        
                        /* store residual as well as Lphi temporarily...
                         * note: temp[1] is reserved for R phi already (cf. f2c_pot()! */
                        cur_node->force.temp[0] = cur_residual;
                        cur_node->force.temp[2] = Lphi;
                        
#ifdef SQNM
                        residual += pow2(cur_residual);
#else
                        residual += fabs(cur_residual);
#endif
                       }
                     else
                       {
                        /* reset temporary storage for those nodes that can not be handled properly... */
                        cur_node->force.temp[0] = 0.0;
                        
                        /* remember: "L phi" is supposed to be rho for the correct solution phi */
                        cur_node->force.temp[2] = cur_source;
                       }
                    }
                 }
              }
           }
        }
     }


#ifdef SQNM
   return (sqrt(residual) / (double)no_nodes);
#else
   return (     residual  / (double)no_nodes);
#endif
#else
   fprintf(stderr, "%s should not be called when AHFlean is defined.\n",
           __func__);
   exit(1);
   return 0.0;
#endif /* AHFlean*/
}

