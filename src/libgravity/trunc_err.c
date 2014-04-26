#include <math.h>
#include <stddef.h> 
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

static long unsigned nnodes;
static double        trunc_error;

/*==============================================================================
 *  trunc_err: estimate the truncation error
 *==============================================================================*/
double trunc_err(gridls *fin_grid)
{
#ifndef AHFlean
   gridls       *coa_grid;
   
   long          i, j, k, idim, ipquad;
   
   pqptr         coa_pquad, tcoa_pquad;
   cqptr         coa_cquad, tcoa_cquad;
   nqptr         coa_nquad, tcoa_nquad;
   nptr          coa_node;
   long          coa_x, coa_y, coa_z;
   double        coa_dens;
   nptr          tsc_coanodes[3][3][3];
   
   pqptr         fin_pquad;
   cqptr         fin_cquad, ifin_cquad;
   nqptr         fin_nquad, ifin_nquad;
   nptr          fin_node;
   long          fin_x, fin_y, fin_z;
   nptr          tsc_finnodes[3][3][3];
   
   double        tau, trunc_error;
   
   /* it's your responsibility to make sure that a coarse grid actually exists */
   coa_grid = fin_grid-1;

   
   /* zero force.temp[1] values on coarse grid (also sets >>interier<< coa_pot-values to zero) */
   zero_temp1(coa_grid);

   /*==================================================================================
    * restrict potential to coa_grid                         
    *                                       (fin_node->pot -> coa_node->force.temp[1])
    *                                       (fin_node->pot -> coa_node->pot)
    *==================================================================================*/
   f2c_pot(fin_grid);
   
   /*==================================================================================
    * loop over refined coarse nodes obtaining "Laplace_pot" values
    *                             (coa_node->force.temp[1] -> coa_node->force.temp[2])
    *==================================================================================*/
   cf_Laplace(coa_grid, fin_grid);
   
   /*==================================================================================
    * prolong coa_node "delsquare_pot" values back to fin_nodes
    *                             (coa_node->force.temp[2] -> fin_node->force.temp[2])
    *==================================================================================*/
   c2f_temp2(coa_grid, fin_grid);
   

   /*==================================================================================
    * loop over fin_nodes accumulating truncation error...
    *                                                 (using fin_node->force.temp[2])
    *==================================================================================*/
   trunc_error = 0.0;
   nnodes      = 0;
   
#ifdef WITH_OPENMP
#pragma omp parallel reduction(+:trunc_error) reduction(+:nnodes) private(ipquad, fin_pquad, fin_cquad, ifin_cquad, fin_nquad, ifin_nquad, fin_node, tsc_finnodes, fin_x, fin_y, fin_z, tau) shared(fin_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<fin_grid->no_pquad; ipquad++)
     {
       fin_pquad=fin_grid->pquad_array[ipquad];
#else
   for(fin_pquad=fin_grid->pquad; fin_pquad != NULL; fin_pquad=fin_pquad->next)
     {
#endif
      fin_z = fin_pquad->z;
      
      for(fin_cquad = fin_pquad->loc;
          fin_cquad < fin_pquad->loc + fin_pquad->length; 
          fin_cquad++, fin_z++)  
        {  
         for(ifin_cquad  = fin_cquad; 
             ifin_cquad != NULL; 
             ifin_cquad  = ifin_cquad->next)
           {
            fin_y = ifin_cquad->y;
            
            for(fin_nquad = ifin_cquad->loc;  
                fin_nquad < ifin_cquad->loc + ifin_cquad->length; 
                fin_nquad++, fin_y++) 
              { 
               for(ifin_nquad  = fin_nquad; 
                   ifin_nquad != NULL; 
                   ifin_nquad  = ifin_nquad->next)
                 {
                  fin_x = ifin_nquad->x;
                  
                  for(fin_node = ifin_nquad->loc; 
                      fin_node < ifin_nquad->loc + ifin_nquad->length; 
                      fin_node++, fin_x++)
                    {
                      /* get neighbouring fin_nodes */
                      tsc_finnodes[1][1][1] = fin_node;
                      get_TSCnodes(fin_grid, fin_pquad, ifin_cquad, ifin_nquad, tsc_finnodes, &fin_z, &fin_y, &fin_x);
                      
                      if(test_tsc(tsc_finnodes))
                       {
                        nnodes++;
                        
                        /* compare it to coa_node->dens */
                        tau = (double)fin_node->force.temp[2] - Laplace_pot(tsc_finnodes, fin_grid->spacing2);
                           
#ifdef SQNM
                        trunc_error += pow2(tau);
#else
                        trunc_error += fabs(tau);
                       }
#endif    
                    }   
                 }
              }
           }
        }
     }
   
#ifdef SQNM
   trunc_error = sqrt(trunc_error) /(double) nnodes;
#else
   trunc_error =      trunc_error  /(double) nnodes;
#endif
   
   /* be more concerned with solution on domain grid */
   if(fin_grid >= global.dom_grid)
      trunc_error *= DOMCORRECT;
   
   return (trunc_error);
   
#else
   fprintf(stderr, "%s should not be called when AHFlean is defined.\n",
           __func__);
   exit(1);
   return 0.0;
#endif /* AHFlean*/
  }
