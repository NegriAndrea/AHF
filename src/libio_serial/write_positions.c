#include <math.h>

#include <stdio.h>
#include <stdlib.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"

#define STEREO2

void write_positions(gridls *grid)
{
   char          filename[100];
   FILE          *posfile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   long          x, y, z;
   long          ncount;
   double        r,g,b;
   
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
      b        = 1;
     }
#endif

   no_part = 0;
   ncount  = 0;
   
   /* prepare filename */
   write_filename(filename, "Positions.", grid->l1dim);
   
   if((posfile = fopen(filename,"w")) == NULL)
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
                     for(cur_part = cur_node->ll; 
                         cur_part != NULL; 
                         cur_part = cur_part->ll)
                       {
                        no_part++;
#ifdef STEREO2
                        fprintf(posfile,"P %f %f %f %g %g %g  2\n",
                                cur_part->pos[0]*simu.boxsize, 
                                cur_part->pos[1]*simu.boxsize, 
                                cur_part->pos[2]*simu.boxsize,
                                r,g,b);
#else
                        fprintf(posfile,"%ld %f %f %f %f %f %f\n",
                                ncount,
                                (float) x/(float)grid->l1dim + cur_shift,
                                (float) y/(float)grid->l1dim + cur_shift,
                                (float) z/(float)grid->l1dim + cur_shift,
                                cur_part->pos[0], 
                                cur_part->pos[1], 
                                cur_part->pos[2]);
#endif
                       }
                    }
                 }
              }
           }
        }
     }
   
   /*  printf("write_positions: no_part = %ld\n", no_part); */
   
   fflush(posfile);
   fclose(posfile);
}
