#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

static long unsigned totnodes,innodes;
static double        max_dens;
static double        cur_shift;

void write_residual(gridls *grid, char *prefix)
{
#ifndef AHFlean
   char          filename[100];
   FILE          *resfile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   nptr          tsc_nodes[3][3][3];
   long          x, y, z;
   long          ncount;
   double        r,g,b;
   double        res_mean;
   
   long unsigned no_part;
   double        cur_shift;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)grid->l1dim;

   innodes   = 0;
   totnodes  = 0;
   res_mean  = 0.0;
   
   
   if((resfile = fopen(prefix,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
   
   for(cur_pquad=grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length; 
          cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++, ncount++)
                    {
                     
                     /* count interior nodes */
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes,
                                  &z, &y, &x);
                     if(test_tsc(tsc_nodes) == TRUE)
                        innodes++;
                     
                     totnodes++;
                      res_mean += cur_node->force.temp[0];
                     
                      if(test_tsc(tsc_nodes) == TRUE)
                       {
                        fprintf(resfile,"%ld %ld %ld %g %g %g %g %ld    1\n",
                                x,y,z,
                                ((float)(x)/(float)(grid->l1dim) + cur_shift),
                                ((float)(y)/(float)(grid->l1dim) + cur_shift),
                                ((float)(z)/(float)(grid->l1dim) + cur_shift),
                                cur_node->force.temp[0],
                                x+y*grid->l1dim+z*pow2(grid->l1dim));
                       }
                     else
                       {
                        fprintf(resfile,"%ld %ld %ld %g %g %g %g %ld    0\n",
                                x,y,z,
                                ((float)(x)/(float)(grid->l1dim) + cur_shift),
                                ((float)(y)/(float)(grid->l1dim) + cur_shift),
                                ((float)(z)/(float)(grid->l1dim) + cur_shift),
                                cur_node->force.temp[0],
                                x+y*grid->l1dim+z*pow2(grid->l1dim));
                       }
                     fflush(resfile);
                    }
                 }
              }
           }
        }
     }
   
   
   printf("write_residual:  grid %8ld has %12ld interior nodes\n", grid->l1dim,innodes);
   printf("                 grid %8ld has %12ld nodes in total\n", grid->l1dim,totnodes);
   printf("                 mean residual = %g\n", res_mean/(double)totnodes);

   fclose(resfile);
#else
   fprintf(stderr, "%s should not be called when AHFlean is defined.\n",
           __func__);
   exit(1);
#endif /* AHFlean*/
}
