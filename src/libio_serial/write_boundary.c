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

/*=============================================================================
* access all particles linked to cur_grid
*=============================================================================*/
void write_boundary(gridls *cur_grid)
{
   char          filename[100];
   FILE          *bndfile, *stereofile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   nptr          tsc_nodes[3][3][3];
   long          ipart, ipart_node,x, y, z;
   float         vx_part, vy_part, vz_part;
   float         vx_node, vy_node, vz_node;
   float         vmean_part, vmean_node, direction, dens_node;
   
   /* prepare filename */
   write_filename(filename, "Boundary.", cur_grid->l1dim);
   
   if((bndfile = fopen(filename,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
   
   write_filename(filename, "BoundaryS2.", cur_grid->l1dim);
   
   if((stereofile = fopen(filename,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
   
   /* loop over all nodes -> particles */
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
      
      for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, z++)
         for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
            
            for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++)
               for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
                  
                  for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                     
                    {
                     
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid,cur_pquad,icur_cquad,icur_nquad,
                                  tsc_nodes,&z,&y,&x);
                     
                     if(test_tsc(tsc_nodes) == FALSE)
                       {
                        fprintf(bndfile,   "%ld %ld %ld\n",x,y,z);
                        fprintf(stereofile,"P %ld %ld %ld 1 0 0 10\n",x,y,z);
                       }
                     else
                       {
                        fprintf(stereofile,"P %ld %ld %ld 0 0 1 1\n",x,y,z);
                       }
                    }
                     
}







