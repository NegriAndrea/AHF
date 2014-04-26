#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>


/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

//#define STEREO2
//#define BINARY //NOTE: STEREO2 and BINARY are mutually exclusive!

static long unsigned totnodes,innodes;
static double        max_dens;
static double        cur_shift;

void write_density(gridls *grid, char *prefix)
{
   char          filename[100];
   FILE          *densfile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   nptr          tsc_nodes[3][3][3];
   long          x, y, z;
   long          ncount;
   double        r,g,b;
   
   double        max_dens;
   
   long unsigned no_part;
   double        cur_shift;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)grid->l1dim;
   
#ifdef STEREO2
   if(grid->l1dim != global.dom_grid->l1dim)
     {
      r        = pow2(log2((double)(grid->l1dim-global.dom_grid->l1dim))/log2((double)(global.fin_l1dim-global.dom_grid->l1dim)));
      g        = log2((double)(grid->l1dim-global.dom_grid->l1dim))/log2((double)(global.fin_l1dim-global.dom_grid->l1dim));
      b        = 1.0;
     }
   else
     {
      r        = 0.0;
      g        = 0.0;
      b        = 1.0;
     }
#endif
   
   innodes  = 0;
   totnodes = 0;
   
   /* prepare filename */
   write_filename(filename, prefix, grid->l1dim);
   strcat(filename,"-DENS");

#ifdef BINARY
  if((densfile = fopen(filename,"wb")) == NULL)
    {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
    }
#else
   if((densfile = fopen(filename,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
#endif
   
#ifdef STEREO2
   max_dens = -1000;
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
                     if(cur_node->dens > max_dens)
                        max_dens = cur_node->dens;
                    }
                 }
              }
           }
        }
     }
#endif
   
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
                     get_TSCnodes(grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                     
                     if(test_tsc(tsc_nodes) == TRUE)
                       {
                        innodes++;
                        
                        /* interior nodes should not have weird densities! */
                        if(cur_node->dens < -simu.mean_dens)
                           fprintf(stderr,"write_density():  strange density %g < %g (l1dim=%ld x=%g y=%g z=%g)\n",
                                   cur_node->dens, simu.mean_dens, grid->l1dim, (x+0.5)/(double)grid->l1dim, (y+0.5)/(double)grid->l1dim, (z+0.5)/(double)grid->l1dim);
                       }
                     
                     totnodes++;
                     
#ifdef STEREO2
                     r = log10(cur_node->dens+2)/log10(max_dens+2);
                     g = log10(cur_node->dens+2)/log10(max_dens+2);
                     b = 0.1;
                     
                     /* no need to visualize empty nodes... */
                     //if(cur_node->dens > -0.99)
                        fprintf(densfile,"P %f %f %f %g %g %g 6     %g\n",
                                ((float)(x)/(float)(grid->l1dim) + cur_shift)*simu.boxsize,
                                ((float)(y)/(float)(grid->l1dim) + cur_shift)*simu.boxsize,
                                ((float)(z)/(float)(grid->l1dim) + cur_shift)*simu.boxsize,
                                r,g,b,
                                cur_node->dens);
#else
#ifdef BINARY
                     fwrite(&cur_node->dens, sizeof(float), 1, densfile);
#else
                     fprintf(densfile,"%f %f %f %f\n",
                             ((float)(x)/(float)(grid->l1dim) + cur_shift)*simu.boxsize,
                             ((float)(y)/(float)(grid->l1dim) + cur_shift)*simu.boxsize,
                             ((float)(z)/(float)(grid->l1dim) + cur_shift)*simu.boxsize,
                             cur_node->dens);
#endif
#endif
                    }
                 }
              }
           }
        }
     }
   
   
   printf("write_density: grid %8ld has %12ld interior nodes\n", grid->l1dim,innodes);
   printf("               grid %8ld has %12ld nodes in total\n", grid->l1dim,totnodes);

   fflush(densfile);
   fclose(densfile);
}
