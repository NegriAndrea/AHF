#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "amr_serial.h"
#include "../libutility/utility.h"

/*==========================================================================
*
* This file contains all routines for relinking particles back
*
* relink         (coa_grid, fin_grid)  -> relink all particles from coa_grid to fin_grid
* relink_back    (coa_grid, fin_grid)  -> relink all particles from fin_grid to coa_grid
* relink_bk_edge (coa_grid, fin_grid)  -> relink only edge particles from fin_grid to coa_grid
*
*==========================================================================*/

static double cur_shift;

/*==============================================================================
* relink: relink particles from coa_grid to fin_grid
*==============================================================================*/
boolean relink(gridls *coa_grid, gridls *fin_grid)
{
   partptr       coa_part, prev_part, next_part;
   long          no_coall, no_finll, no_totll, i;
   
   pqptr         coa_pquad,             tcoa_pquad;
   cqptr         coa_cquad, icoa_cquad, tcoa_cquad;
   nqptr         coa_nquad, icoa_nquad, tcoa_nquad;
   nptr          coa_node;
   long          coa_x, coa_y, coa_z;
   
   pqptr         fin_pquad;
   cqptr         fin_cquad, ifin_cquad;
   nqptr         fin_nquad, ifin_nquad;
   nptr          fin_node;
   long          fin_x, fin_y, fin_z;
   double        x_finnode, y_finnode, z_finnode;
   double        x_part, y_part, z_part;
   double        dx, dy, dz;
   nptr          tsc_nodes[3][3][3];
   
   long          no_finpart, no_finnodes;
   double        ratio;
   long          ipquad;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift   = 0.5/(double)fin_grid->l1dim;
   
   no_finpart  = 0;
   no_finnodes = 0;
   
   /*==============================================================
      * loop over fine grid (and simultanesouly over coarse grid...)
      *==============================================================*/
   for(fin_pquad=fin_grid->pquad; fin_pquad != NULL; fin_pquad=fin_pquad->next)
     {
      fin_z = fin_pquad->z;
      coa_z = fin_z/2;
      
      /* find correct coa_pquad */
      for(tcoa_pquad = coa_grid->pquad; coa_z > tcoa_pquad->z + tcoa_pquad->length; tcoa_pquad = tcoa_pquad->next)
         ;
      
      /* jump to correct cquad */
      coa_cquad = tcoa_pquad->loc + (coa_z - tcoa_pquad->z);
      
      
      for(fin_cquad = fin_pquad->loc;
          fin_cquad < fin_pquad->loc + fin_pquad->length; 
          fin_cquad++, fin_z++)  
        {  
         for(ifin_cquad  = fin_cquad; 
             ifin_cquad != NULL; 
             ifin_cquad  = ifin_cquad->next)
           {
            fin_y = ifin_cquad->y;
            coa_y = fin_y/2;
            
            /* find correct coa_cquad */
            for(tcoa_cquad = coa_cquad; coa_y > tcoa_cquad->y + tcoa_cquad->length; tcoa_cquad = tcoa_cquad->next)
               ;
            
            /* jump to correct nquad */
            coa_nquad = tcoa_cquad->loc + (coa_y - tcoa_cquad->y);
            
            for(fin_nquad = ifin_cquad->loc;  
                fin_nquad < ifin_cquad->loc + ifin_cquad->length; 
                fin_nquad++, fin_y++) 
              { 
               for(ifin_nquad  = fin_nquad; 
                   ifin_nquad != NULL; 
                   ifin_nquad  = ifin_nquad->next)
                 {
                  fin_x = ifin_nquad->x;
                  coa_x = fin_x/2;
                  
                  /* find correct coarse nquad */
                  for(tcoa_nquad = coa_nquad; coa_x > tcoa_nquad->x + tcoa_nquad->length; tcoa_nquad = tcoa_nquad->next)
                     ;
                  
                  /* jump to correct (first!) node */
                  coa_node = tcoa_nquad->loc + (coa_x - tcoa_nquad->x);
                  
                  for(fin_node = ifin_nquad->loc; 
                      fin_node < ifin_nquad->loc + ifin_nquad->length; 
                      fin_node++, fin_x++)
                    {
                     
                     no_finnodes++;
                     
                     /* we now have access to both the fin_node and the mother coa_node */
                     
                     /* are there any particles at this coa_node? */
                     if(coa_node->ll != NULL)
                       {
                        /* transfer only to interior fin_nodes */
                        tsc_nodes[1][1][1] = fin_node;
                        get_TSCnodes(fin_grid, fin_pquad, ifin_cquad, ifin_nquad, tsc_nodes, &fin_z, &fin_y, &fin_x);
                        
                        if(test_tsc(tsc_nodes) == TRUE)
                          {         
                           /* counters for new linked lists (coa_node and fin_node) */
                           no_finll  = 0;
                           no_coall  = 0;
                           
                           /* take first particle attached to coa_node */
                           coa_part = coa_node->ll;
                           
                           /* fin_node coordinates in particle units! */
                           x_finnode = (((double)fin_x) / (double)fin_grid->l1dim) + cur_shift;
                           y_finnode = (((double)fin_y) / (double)fin_grid->l1dim) + cur_shift;
                           z_finnode = (((double)fin_z) / (double)fin_grid->l1dim) + cur_shift;
                           
                           /* calculate distance to fin_node boundaries */
                           x_part = coa_part->pos[X];
                           y_part = coa_part->pos[Y];
                           z_part = coa_part->pos[Z];
                           dx     = fabs(x_part-x_finnode);
                           dy     = fabs(y_part-y_finnode);
                           dz     = fabs(z_part-z_finnode);
                           
                           /* first coa_part lies within fin_node boundaries */
                           while( dx <= cur_shift && dy <= cur_shift && dz <= cur_shift )
                             {
                              /* count no. of particles transferred to fin_node */
                              no_finll++;
                              no_finpart++;
                              
                              /* remove particle from coa_node */
                              coa_node->ll = coa_part->ll;
                              
                              /* insert particle at front (messes up x-sort of linked list!!!) */
                              coa_part->ll = fin_node->ll;
                              fin_node->ll = coa_part;
                              
                              /* move to next particle in linked list */
                              prev_part = coa_part;
                              coa_part  = coa_node->ll;
                              
                              if(coa_part != NULL)
                                {
                                 /* calculate distance to fin_node boundaries again */
                                 x_part = coa_part->pos[X];
                                 y_part = coa_part->pos[Y];
                                 z_part = coa_part->pos[Z];
                                 dx     = fabs(x_part-x_finnode);
                                 dy     = fabs(y_part-y_finnode);
                                 dz     = fabs(z_part-z_finnode);
                                }
                              else
                                {
                                 /* exit while-loop */
                                 dx = 2*cur_shift;
                                }
                             }
                           
                           /* the current coa_part (if existent) lies *not* in fin_node */
                           
                           /* loop over remainder of coa_node particles */
                           while(coa_part != NULL)
                             {
                              /* calculate distance to fin_node boundaries again */
                              x_part = coa_part->pos[X];
                              y_part = coa_part->pos[Y];
                              z_part = coa_part->pos[Z];
                              dx     = fabs(x_part-x_finnode);
                              dy     = fabs(y_part-y_finnode);
                              dz     = fabs(z_part-z_finnode);                              
                              
                              /* particle lies in fin_node */
                              if( dx <= cur_shift && dy <= cur_shift && dz <= cur_shift )
                                {
                                 /* count no. of particles transferred to fin_node */
                                 no_finll++;
                                 no_finpart++;

                                 /* remove particle from coa_node linked-list */
                                 prev_part->ll = coa_part->ll;
                                 /* insert particle at front (messes up x-sort of linked list!!!) */
                                 coa_part->ll = fin_node->ll;
                                 fin_node->ll = coa_part;
                                 
                                 /* prev_part stays the same as coa_part has been removed */
                                 
                                 /* jump to next coa_part in linked-list */
                                 coa_part = prev_part->ll;
                                }
                              else
                                {
                                 /* count no. of particles remaining at coa_node */
                                 no_coall++;
                                 
                                 /* remember coa_part as prev_part */
                                 prev_part = coa_part;
                                 
                                 /* jump to next coa_part in linked-list */
                                 coa_part = coa_part->ll;                                 
                                }
                              
                             } /* while(coa_part!=NULL) */                           
                          } /* test_tsc */
                       } /* coa_part != NULL */
                  
                  /* move to next coa_node */
                  if(is_even(fin_x) == FALSE)
                    {
                     coa_node++;
                     coa_x++;
                    }
                  
                    }
                 }
            
            /* move to next coa_nquad */
            if(is_even(fin_y) == FALSE)
              {
               coa_nquad++;
               coa_y++;
              }
              }
           }
      
      /* move to next coa_cquad */
      if(is_even(fin_z) == FALSE)
        {
         coa_cquad++;
         coa_z++;
        }
        }
     }



   /* now check the ratio of number of particles and number of nodes... */
   fin_grid->size.no_nodes = no_finnodes;
   fin_grid->size.no_part  = no_finpart;


   /* do not allow refinements smaller than MIN_NNODES */
   if(no_finnodes < MIN_NNODES)
      return(FALSE);

   /* check mean ratio of (number of nodes)/(number of particles) on grid */
   if(simu.np_limit == FALSE)
      return(TRUE);
   else
     {
      if(fin_grid->l1dim == 2*simu.NGRID_DOM && simu.Nth_dom < 3.)
         return(TRUE);
      
      ratio = (double)no_finnodes/(double)no_finpart;
      if(ratio > NP_RATIO)
         return(FALSE);
      else
         return(TRUE);
     }

}


/*==============================================================================
* relink_back: relink particles back from fin_grid to coa_grid
*==============================================================================*/
void relink_back(gridls *coa_grid, gridls *fin_grid)
{
   partptr       fin_part, coa_part, for_part;
   part          tran_coa;        /* transfer particle for coarse grid */
   part          tran_fin;        /* transfer particle for fine grid   */
   
   pqptr         coa_pquad,             tcoa_pquad;
   cqptr         coa_cquad, icoa_cquad, tcoa_cquad;
   nqptr         coa_nquad, icoa_nquad, tcoa_nquad;
   nptr          coa_node;
   long          coa_x, coa_y, coa_z;
   
   pqptr         fin_pquad;
   cqptr         fin_cquad, ifin_cquad;
   nqptr         fin_nquad, ifin_nquad;
   nptr          fin_node;
   long          fin_x, fin_y, fin_z;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)fin_grid->l1dim;
   
   /*==============================================================
    * loop over fine grid (and simultanesouly over coarse grid...)
    *==============================================================*/
   for(fin_pquad=fin_grid->pquad; fin_pquad != NULL; fin_pquad=fin_pquad->next)
     {
      fin_z = fin_pquad->z;
      coa_z = fin_z/2;
      
      /* find correct coa_pquad */
      for(tcoa_pquad = coa_grid->pquad; coa_z > tcoa_pquad->z + tcoa_pquad->length; tcoa_pquad = tcoa_pquad->next)
         ;
      
      /* jump to correct cquad */
      coa_cquad = tcoa_pquad->loc + (coa_z - tcoa_pquad->z);
      
      
      for(fin_cquad = fin_pquad->loc;
          fin_cquad < fin_pquad->loc + fin_pquad->length; 
          fin_cquad++, fin_z++)  
        {  
         for(ifin_cquad  = fin_cquad; 
             ifin_cquad != NULL; 
             ifin_cquad  = ifin_cquad->next)
           {
            fin_y = ifin_cquad->y;
            coa_y = fin_y/2;
            
            /* find correct coa_cquad */
            for(tcoa_cquad = coa_cquad; coa_y > tcoa_cquad->y + tcoa_cquad->length; tcoa_cquad = tcoa_cquad->next)
               ;
            
            /* jump to correct nquad */
            coa_nquad = tcoa_cquad->loc + (coa_y - tcoa_cquad->y);
            
            for(fin_nquad = ifin_cquad->loc;  
                fin_nquad < ifin_cquad->loc + ifin_cquad->length; 
                fin_nquad++, fin_y++) 
              { 
               for(ifin_nquad  = fin_nquad; 
                   ifin_nquad != NULL; 
                   ifin_nquad  = ifin_nquad->next)
                 {
                  fin_x = ifin_nquad->x;
                  coa_x = fin_x/2;
                  
                  /* find correct coarse nquad */
                  for(tcoa_nquad = coa_nquad; coa_x > tcoa_nquad->x + tcoa_nquad->length; tcoa_nquad = tcoa_nquad->next)
                     ;
                  
                  /* jump to correct (first!) node */
                  coa_node = tcoa_nquad->loc + (coa_x - tcoa_nquad->x);
                  
                  for(fin_node = ifin_nquad->loc; 
                      fin_node < ifin_nquad->loc + ifin_nquad->length; 
                      fin_node++, fin_x++)
                    {
                     /* we now have access to both the fin_node and the mother coa_node */
                     
                     /* simply transfer all particles back to mother coa_node */
                     /* merge the lists */
                     /* are there actually particles to transfer? */
                     if(fin_node->ll != NULL)
                       {
                        /* are there already particles at coa_node? */
                        if(coa_node->ll != NULL)
                          {
                           /* loop to end of coa_node's linked-list */
                           for(for_part = coa_node->ll; for_part->ll != NULL; for_part=for_part->ll)
                              ;
                           
                           for_part->ll = fin_node->ll;
                          }
                        else
                          {
                           coa_node->ll = fin_node->ll;
                          }
                       }
                     
                     /* there are no more particles at fin_node */
                     fin_node->ll = NULL;
                     
                     /* move to next coa_node */
                     if(is_even(fin_x) == FALSE)
                       {
                        coa_node++;
                        coa_x++;
                       }
                     
                    }
                 }
               
               /* move to next coa_nquad */
               if(is_even(fin_y) == FALSE)
                 {
                  coa_nquad++;
                  coa_y++;
                 }
              }
           }
         
         /* move to next coa_cquad */
         if(is_even(fin_z) == FALSE)
           {
            coa_cquad++;
            coa_z++;
           }
        }
     }
}
/*==============================================================================
* relink_bk_edge: relink edge particles from coa_grid to fin_grid
*==============================================================================*/
void relink_bk_edge(gridls *coa_grid, gridls *fin_grid)
{
   partptr       cur_part, next_part;
   
   pqptr         coa_pquad,             tcoa_pquad;
   cqptr         coa_cquad, icoa_cquad, tcoa_cquad;
   nqptr         coa_nquad, icoa_nquad, tcoa_nquad;
   nptr          coa_node;
   long          coa_x, coa_y, coa_z;
   
   pqptr         fin_pquad;
   cqptr         fin_cquad, ifin_cquad;
   nqptr         fin_nquad, ifin_nquad;
   nptr          fin_node;
   long          fin_x, fin_y, fin_z;
   
   nptr          tsc_nodes[3][3][3];  
   
   cur_shift = 0.5/(double)fin_grid->l1dim;
   
   /*==============================================================
      * loop over fine grid (and simultanesouly over coarse grid...)
      *==============================================================*/
   for(fin_pquad=fin_grid->pquad; fin_pquad != NULL; fin_pquad=fin_pquad->next)
     {
      fin_z = fin_pquad->z;
      coa_z = fin_z/2;
      
      /* find correct coa_pquad */
      for(tcoa_pquad = coa_grid->pquad; coa_z > tcoa_pquad->z + tcoa_pquad->length; tcoa_pquad = tcoa_pquad->next)
         ;
      
      /* jump to correct cquad */
      coa_cquad = tcoa_pquad->loc + (coa_z - tcoa_pquad->z);
      
      
      for(fin_cquad = fin_pquad->loc;
          fin_cquad < fin_pquad->loc + fin_pquad->length; 
          fin_cquad++, fin_z++)  
        {  
         for(ifin_cquad  = fin_cquad; 
             ifin_cquad != NULL; 
             ifin_cquad  = ifin_cquad->next)
           {
            fin_y = ifin_cquad->y;
            coa_y = fin_y/2;
            
            /* find correct coa_cquad */
            for(tcoa_cquad = coa_cquad; coa_y > tcoa_cquad->y + tcoa_cquad->length; tcoa_cquad = tcoa_cquad->next)
               ;
            
            /* jump to correct nquad */
            coa_nquad = tcoa_cquad->loc + (coa_y - tcoa_cquad->y);
            
            for(fin_nquad = ifin_cquad->loc;  
                fin_nquad < ifin_cquad->loc + ifin_cquad->length; 
                fin_nquad++, fin_y++) 
              { 
               for(ifin_nquad  = fin_nquad; 
                   ifin_nquad != NULL; 
                   ifin_nquad  = ifin_nquad->next)
                 {
                  fin_x = ifin_nquad->x;
                  coa_x = fin_x/2;
                  
                  /* find correct coarse nquad */
                  for(tcoa_nquad = coa_nquad; coa_x > tcoa_nquad->x + tcoa_nquad->length; tcoa_nquad = tcoa_nquad->next)
                     ;
                  
                  /* jump to correct (first!) node */
                  coa_node = tcoa_nquad->loc + (coa_x - tcoa_nquad->x);
                  
                  for(fin_node = ifin_nquad->loc; 
                      fin_node < ifin_nquad->loc + ifin_nquad->length; 
                      fin_node++, fin_x++)
                    {
                     /* we now have access to both the fin_node and the mother coa_node */
                     
                     /* are we on the edge? */
                     tsc_nodes[1][1][1] = fin_node;
                     get_TSCnodes(fin_grid, fin_pquad, ifin_cquad, ifin_nquad, tsc_nodes, 
                                  &fin_z, &fin_y, &fin_x);
                     
                     /* current node is on the edge */
                     if(test_tsc(tsc_nodes) == FALSE)
                       {
                        for(cur_part = fin_node->ll; cur_part != NULL; cur_part = next_part)
                          {
                           /* in case store_leaver() calls coagrid_llsearch() */
                           next_part = cur_part->ll;
                          }
                        
                        fin_node->ll = NULL;   /* no particles are linked to fin_node anymore */
                       }
                     
                     /* move to next coa_node */
                     if(is_even(fin_x) == FALSE)
                       {
                        coa_node++;
                        coa_x++;
                       }
                     
                    }
                 }
               
               /* move to next coa_nquad */
               if(is_even(fin_y) == FALSE)
                 {
                  coa_nquad++;
                  coa_y++;
                 }
              }
           }
         
         /* move to next coa_cquad */
         if(is_even(fin_z) == FALSE)
           {
            coa_cquad++;
            coa_z++;
           }
        }
     }
}


