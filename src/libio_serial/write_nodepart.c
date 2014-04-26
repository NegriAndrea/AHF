#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"

/*=============================================================================
* access all particles linked to cur_grid
*=============================================================================*/
void write_nodepart(gridls *cur_grid)
{
   char          filename[100];
   FILE          *posfile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   long          ipart, ipart_node,x, y, z;
   
  float cur_shift = 0.5/(double)cur_grid->l1dim;
  
   /* prepare filename */
   write_filename(filename, "NodePart.", cur_grid->l1dim);
   
   if((posfile = fopen(filename,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
   
   ipart = 0;
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
                     ipart_node = 0;
                     
                     for(cur_part = cur_node->ll; 
                         cur_part != NULL; 
                         cur_part = cur_part->ll)
                       {
                        
                        /*-----------------------------------------------------------
                        * from this point on you can do whatever with the particles
                        *-----------------------------------------------------------*/
                         fprintf(posfile,"%f %f %f %f  %f  %f\n",
                                 cur_part->pos[0], 
                                 cur_part->pos[1],
                                 cur_part->pos[2],
                                 cur_part->pos[0]*simu.boxsize, 
                                 cur_part->pos[1]*simu.boxsize,
                                 cur_part->pos[2]*simu.boxsize);
                                 //((float)(cur_part->pos[0])/(float)(cur_grid->l1dim) + cur_shift)*simu.boxsize,
                                 //((float)(cur_part->pos[1])/(float)(cur_grid->l1dim) + cur_shift)*simu.boxsize,
                                 //((float)(cur_part->pos[2])/(float)(cur_grid->l1dim) + cur_shift)*simu.boxsize);
                        ipart++;
                        ipart_node++;
                       }
#ifdef STEREO2
                      //fprintf(posfile,"p %ld %ld %ld 1 0 0\n",x,y,z);
#else
                     /*	if(ipart_node > 1) */
                      // fprintf(posfile,"%ld %ld %ld %ld %f\n",
                      //       x, y, z, ipart_node, cur_node->dens+simu.mean_dens);
#endif
                    }
                     
                     /*  fprintf(stderr,"write_nodepart: %ld = %ld\n",cur_grid->l1dim,ipart); */
                     fclose(posfile);
}
