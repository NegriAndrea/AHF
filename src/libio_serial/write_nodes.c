#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

#define STEREO2

/*=============================================================================
* simply write the position of all nodes
*=============================================================================*/
void write_nodes(gridls *cur_grid, char *prefix)
{
   char          filename[MAXSTRING];
   FILE          *nodesfile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   long          ipart, ipart_node,x, y, z;
   
   double        cur_shift;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   /* prepare filename */
   write_filename(filename, prefix, cur_grid->l1dim);
   
   strcat(filename,"-NODES");
   
   if((nodesfile = fopen(filename,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
   
   /* loop over all nodes */
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
     {
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
#ifdef STEREO2
                     fprintf(nodesfile,"P %g %g %g 1 0 0 4\n",
                             ((float)(x)/(float)(cur_grid->l1dim) + cur_shift)*simu.boxsize,
                             ((float)(y)/(float)(cur_grid->l1dim) + cur_shift)*simu.boxsize,
                             ((float)(z)/(float)(cur_grid->l1dim) + cur_shift)*simu.boxsize);
#else
                     fprintf(nodesfile,"%ld %ld %ld %f\n",x,y,z,cur_node->dens+simu.mean_dens);
#endif
                    }
                 }
              }
           }
        }
      //fprintf(nodesfile,"# linked pquad...\n");
     }
   fclose(nodesfile);
}
