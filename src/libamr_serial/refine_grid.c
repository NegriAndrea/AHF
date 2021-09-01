#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "amr_serial.h"
#include "../libutility/utility.h"


/*================================================================================
*
/*================================================================================
*
*          This file contains everything to refine a given grid
*          ----------------------------------------------------
*
* The visible function is at the very end:
*
* refine_grid(fin_grid, coa_grid)
*
* (a 'fin_grid' structure needs to be available to refine_grid...)
*
*
* This file starts off with declaring all the test functions for
* checking coa_nodes for refinement....
*
* Then follows a hierarchy of functions guiding one through the
* rather complicated process of sweeping through the coarse grid
* and simultaneously creating the fine grid...
*
*
*===============================================================================*/


/*------------------------------------------------------------------------------
* generate a new nquad2-linked-list by simply mirror-ing the nquad1-linked-list
*
*-------------------------------------------------------------------------------*/
void mirror_nquad(nqptr nquad1, nqptr nquad2)
{
   /* we need to mirror every old fin_nquad in the linked list */
   while(nquad1->next != NULL)
     {
      /* copy nquad */
      nquad2->loc    = c_node(nquad1->length);
      nquad2->x      = nquad1->x;
      nquad2->length = nquad1->length;
      nquad2->next   = c_nquad(1);
      
      /* move to next nquad to be copied */
      nquad1 = nquad1->next;
      nquad2 = nquad2->next;
     }
   
   /* copy nquad */
   nquad2->loc    = c_node(nquad1->length);
   nquad2->x      = nquad1->x;
   nquad2->length = nquad1->length;
   nquad2->next   = NULL;
}

/*------------------------------------------------------------------------------
* generate a new cquad2-linked-list by simply mirror-ing the cquad1-linked-list
*
*-------------------------------------------------------------------------------*/
void mirror_cquad(cqptr cquad1, cqptr cquad2)
{
   nqptr i_nquad1, i_nquad2;
   
   /* we need to mirror every old fin_cquad in the linked list */
   while(cquad1->next != NULL)
     {
      /* copy cquad */
      cquad2->loc    = c_nquad(cquad1->length);
      cquad2->y      = cquad1->y;
      cquad2->length = cquad1->length;
      cquad2->next   = c_cquad(1);
      
      /* copy linked nquad's */
      for(i_nquad1 = cquad1->loc, i_nquad2=cquad2->loc;
          i_nquad1 < cquad1->loc+cquad1->length;
          i_nquad1++, i_nquad2++)      
         mirror_nquad(i_nquad1, i_nquad2);
      
      /* move to next cquad to be copied */
      cquad1 = cquad1->next;
      cquad2 = cquad2->next;
     }
   
   /* copy cquad */
   cquad2->loc    = c_nquad(cquad1->length);
   cquad2->y      = cquad1->y;
   cquad2->length = cquad1->length;
   cquad2->next   = NULL;

   /* copy linked nquad's */
   for(i_nquad1 = cquad1->loc, i_nquad2=cquad2->loc;
       i_nquad1 < cquad1->loc+cquad1->length;
       i_nquad1++, i_nquad2++)      
      mirror_nquad(i_nquad1, i_nquad2);
   
}

/*------------------------------------------------------------------------------
 * test_node()
 *              simply test one single node for refinement!
 *------------------------------------------------------------------------------*/
boolean test_node(gridls *cur_grid, pqptr cur_pquad, cqptr cur_cquad, nqptr cur_nquad, nptr cur_node, 
                 long z, long y, long x)
{
   nptr  tsc_nodes[3][3][3];
   int   i,j,k;

#ifdef TEST_CURNODE_ONLY
   if (cur_node->dens >= (cur_grid->critdens-simu.mean_dens))
      return TRUE;
#else
   /* 1. criterion: not an edge node */
   tsc_nodes[1][1][1] = cur_node;
   get_TSCnodes(cur_grid, cur_pquad, cur_cquad, cur_nquad, tsc_nodes, &z, &y, &x);
   if(test_tsc(tsc_nodes) == FALSE)
      return FALSE;
   
   /* 2. criterion: actual node (or any of the surrounding nodes ahead) meets refinement criterion */
   for(k=0; k<NDIM; k++)
      for(j=0; j<NDIM; j++)
         for(i=1; i<NDIM; i++)
            if(tsc_nodes[k][j][i]->dens >= (cur_grid->critdens-simu.mean_dens))
               return TRUE;
#endif
  
   return FALSE;
}
   

/*==============================================================================
*
* now follow all routines for looping over coarse nodes...
* ...and simultaneously creating the fine grid !
*
*==============================================================================*/



/*------------------------------------------------------------------------------
* ref_nquad: nquad refinement function. Add blocks of nodes to supplied nquad.
* Generates additional nquads (linked to "next" pointer) as required and calls
* itself recursively.
*------------------------------------------------------------------------------*/
boolean ref_nquad(gridls *coa_grid, pqptr coa_pquad, cqptr coa_cquad, nqptr coa_nquad,
                  nqptr ifine_nquad, long z, long y)
{
   nptr    icoa_node;            /* pointer to coarse node                     */
   boolean state   = FALSE;      /* on/off refinement state indicator          */
   boolean refined = FALSE;      /* test result variable                       */
   boolean refever = FALSE;      /* have we ever refined along this column     */
   nqptr   icoa_nquad;           /* pointer to coarse cquad                    */
   long    x;                    /* x-coord of current coarse node             */
   long    xlen = 0;             /* length of current refinement (x-direction) */
   long    xoffset;
   int     idummy;
   nptr    tsc_nodes[3][3][3];
   
#ifdef PERIODIC_X
   /* first node is on edge - account for periodicity */
   if(coa_nquad->x == 0)           
     {
      /* find last coarse nquad */
      for(icoa_nquad = coa_nquad; icoa_nquad->next != NULL ; icoa_nquad = icoa_nquad->next)
         ;
      
      /* is the last node in the last coarse nquad on the edge of the grid? */
      if(icoa_nquad->x + icoa_nquad->length == coa_grid->l1dim)
        {
         icoa_node = icoa_nquad->loc + (icoa_nquad->length - 1);
         x         = coa_grid->l1dim - 1;
         
         /* and if test indicates refinement needed... */
         if(test_node(coa_grid,coa_pquad,coa_cquad,icoa_nquad,icoa_node,z,y,x) == TRUE)
           {
            /* ...start a refinement */
            state              = TRUE;
            refever            = TRUE;
            (ifine_nquad)->x   = 0; /* fill in start of refinement */            
           }
         
         /* this is the first node with a periodic counterpart -> it's NOT an edge node */
         xoffset = 0;
        }
      else
        {
         /* this is the first node with no periodic counterpart -> it's an edge node! */
         xoffset = 1;
        }
     }
   else
     {
      /* we are not on the boundary and hence do not try to refine the first node */
      xoffset = 1;
     }
#else /* PERIODIC_X */
   /* do not try to refine the first node as this is an edge node... */
   xoffset = 1;
#endif /* PERIODIC_X */
   
   /* loop over coarse nodes (but the last one!) */
   for(icoa_node = coa_nquad->loc + xoffset, x = coa_nquad->x + xoffset; 
       icoa_node < coa_nquad->loc + (coa_nquad->length - 1); 
       x++, icoa_node++)
     {
      
      /* test if refinement is needed... */
      if(test_node(coa_grid,coa_pquad,coa_cquad,coa_nquad,icoa_node,z,y,x) == TRUE)
        {
         if(state == FALSE)
           {
            state             = TRUE;
            refever           = TRUE;
            (ifine_nquad)->x  = (x * 2); /* fill in start of refinement */
           }
         
         xlen += 2;       /* add 2 fin_grid nodes        */
        }
      
      /* only state indicates refine - we must be concluding a fine_nquad */
      else if(state == TRUE)
        {
         /* no further refining needed */
         state = FALSE;
         
         /* only add ghost cells if not on the edge */
         tsc_nodes[1][1][1] = icoa_node;
         get_TSCnodes(coa_grid,coa_pquad,coa_cquad,coa_nquad,tsc_nodes,&z,&y,&x);
         if(test_tsc(tsc_nodes) == TRUE)
           {
            /* add two (gost) fin_nodes to length */
            xlen += 2;            
           }

         ifine_nquad->length = xlen;            /* fill in length of refinement*/
         ifine_nquad->loc    = c_node(xlen);    /* create the correct no nodes */
         ifine_nquad->next   = c_nquad(1);      /* gen. new linked nquad       */
         ifine_nquad         = (ifine_nquad)->next; /* move to linked nquad    */
         xlen                = 0;               /* reset xlen                  */
        }
     }
      
   /* last node described by coarse nquad */
      
   /* CASE 1. last node, but there is another linked coarse nquad */
   if(coa_nquad->next != NULL)   
     {
      /* conclude current fine_nquad without refining this edge node */
      if(state == TRUE)        
        {
         /* do not increase xlen -> do not refine edge nodes ! */         
         ifine_nquad->length = (xlen);          /* fill in length of refinement*/
         ifine_nquad->loc    = c_node(xlen);        /* generate nodes          */
         ifine_nquad->next   = c_nquad(1);          /* gen. new linked nquad   */
         ifine_nquad         = (ifine_nquad)->next; /* move to linked nquad    */
        }
      
      /* call recursively */
      refined = ref_nquad(coa_grid, coa_pquad, coa_cquad, coa_nquad->next, ifine_nquad, z, y);
     }
   
   /* CASE 2. this last node is at the boundary and refinement is required */
   else if((coa_nquad->x + coa_nquad->length == coa_grid->l1dim) &&
           (test_node(coa_grid, coa_pquad, coa_cquad, coa_nquad, icoa_node, z, y, x) == TRUE))
     {
      /* new refinement */
      if(state == FALSE)  
        {
         refever                  = TRUE;        /* we have refined             */
         (ifine_nquad)->x         = (x * 2);     /* fill in start of refinement */
        }
      
      /* create *two* fin_nodes: cospatial and non-cospatial
       * (OCTSPLIT always requires 2 nodes anyways) */
      xlen += 2;
      
      /* now conclude the fine_nquad */
      ifine_nquad->length   = xlen;
      ifine_nquad->loc      = c_node(xlen);
      ifine_nquad->next     = NULL;
     }

   /* CASE 3. we are not at the boundary but are currently refining -> conclude refinement */
   else if(state == TRUE)     
     {
      /* do not increase xlen -> do not refine edge nodes ! */
      ifine_nquad->length      = xlen;
      ifine_nquad->loc         = c_node(xlen);
      ifine_nquad->next        = NULL;
     }
            
   return (refever || refined);
}

/*--------------------------------------------------------------------------------
* ref_cquad: cquad refinement function. Add blocks of nquads to supplied cquad.
* Generates additional cquads (linked to "next" pointer) as required and calls
* itself recursively.
*--------------------------------------------------------------------------------*/
boolean  ref_cquad(gridls *coa_grid, pqptr coa_pquad, cqptr coa_cquad, 
                   cqptr ifine_cquad, long z)
{
   nqptr   icoa_nquad;           /* pointer to coarse nquad                   */
   boolean refined;              /* refined result indicator                  */
   boolean refever;              /* have we ever refined                      */
   nqptr   i_fine_nq;            /* current fine nquad                        */
   boolean state;                /* on/off state flag                         */
   nqptr   rem_nquad;            /* to remove unused nquads                   */
   long    y;                    /* y-coord of current coarse node            */
   long    ylen;                 /* length of current refinement y-direction) */
   
   cqptr   fst_cquad, lst_cquad; /* do not refine edge nodes and hence calculate "offsets" */
   long    y_low, y_up;

   /* initilise state flags */
   refined = FALSE;
   refever = FALSE;
   state   = FALSE;
   
   /* initialise ylen = [1, ...] */
   ylen = 0;
   
   /*
    * The new cquad requires a nquad be created. Further nquads can be added and
    * blocks of nodess may be created and referenced as the refinement routine
    * requires them.
    */
   ifine_cquad->loc = c_nquad(1);          /* create first nquad                  */
   i_fine_nq        = ifine_cquad->loc;    /* set i_fine_nq to point to new nquad */

   
   /* we should actually avoid looping over the first column,
      but in that case we would need to make sure that the
      first column is *not* periodically wrapped... */ 
   
#ifdef PERIODIC_Y
   /* jump to last cquad */
   fst_cquad = coa_pquad->loc;

   /* jump to last cquad */
   lst_cquad = fst_cquad;
   while(lst_cquad->next != NULL)
      lst_cquad = lst_cquad->next;
   
   /* check whether things are periodically wrapped in y-direction */
   if(fst_cquad->y == 0 && lst_cquad->y+lst_cquad->length == coa_grid->l1dim)
     {
      /* current y-refinement runs through whole box */
      if(coa_cquad->y == 0 && coa_cquad->y+coa_cquad->length == coa_grid->l1dim)
        {
         y_low =  0;
         y_up  =  0;
        }
      /* current y-position is at lower boundary */
      else if(coa_cquad->y == 0)
        {
         y_low =  0;
         y_up  = -1;
        }
      /* current y-position will reach upper boundary */
      else if(coa_cquad->y+coa_cquad->length == coa_grid->l1dim)
        {
         y_low = +1;
         y_up  =  0;
        }
      /* we are somewhere in-between */
      else
        {
         y_low = +1;
         y_up  = -1;
        }
     }
   /* no periodic wrapping whatsoever */
   else
     {
      y_low = +1;
      y_up  = -1;
     }
   
#else /* PERIODIC_Y */
   y_low = +1;
   y_up  = -1;
#endif /* PERIODIC_Y */

   /* loop over all coarse columns (but the last one!) */
   for(icoa_nquad = coa_cquad->loc+y_low, y = coa_cquad->y+y_low; 
       icoa_nquad < coa_cquad->loc + (coa_cquad->length - 1) +  y_up;
       icoa_nquad++, y++)
     {
      /*------------------------
       * test cospatial columns
       *------------------------*/
      
      /* call ref_nquad -> generates a new linked-list of i_fine_nq's (if == TRUE) */
      if(ref_nquad(coa_grid, coa_pquad, coa_cquad, icoa_nquad, i_fine_nq, z, y) == TRUE)
        {
         /* check "next" pointer not NULL */
         if(i_fine_nq->next != NULL)
           {
            
            /* find the last nquad along this column */
            for(rem_nquad = i_fine_nq; (rem_nquad->next)->next != NULL ;
                rem_nquad = rem_nquad->next)
               ;
            
            /*if loc pointer is NULL, destroy*/
            if((rem_nquad->next)->loc == NULL) 
              {
               dest_nquad(&(rem_nquad->next));
               rem_nquad->next = NULL;
              }
           }
         
         /* cospatial column => start refinement */
         if(state == FALSE)
           {
            state          = TRUE;     /* refinement is occurring set state to TRUE */
            refever        = TRUE;     /* refinement has occurred/set refever flag  */
            ifine_cquad->y = (y * 2);  /* fill in start of refinement               */
           }
         
         /* increase length by one */
         ylen++;                     
         
         /* generate another i_fine_nq placeholder */
         ifine_cquad->loc = r_nquad(ifine_cquad->loc, ylen, (ylen + 1));
         
         /* jump to newly created i_fine_nq */
         i_fine_nq        = ifine_cquad->loc + ylen;
         
         /*----------------------------------------------------------
          * simply copy previous i_fine_nq over to newly created one
          *----------------------------------------------------------*/
         mirror_nquad(ifine_cquad->loc+(ylen-1), i_fine_nq);
         
         /* update ylen to account for this newly created fin_nquad-linked-list */
         ylen++;

         /* generate another i_fine_nq placeholder */
         ifine_cquad->loc = r_nquad(ifine_cquad->loc, ylen, (ylen + 1));
         
         /* jump to newly created i_fine_nq */
         i_fine_nq        = ifine_cquad->loc + ylen;
        }
      
      /* we did not refine this cospatial nquad 
       * -> conclude current refinement (if needed, i.e. state==TRUE!) */
      else if(state == TRUE)
        {
         /* set state to FALSE  */
         state = FALSE;

         /* remove the last i_fine_nq placeholder */
         ifine_cquad->loc    = r_nquad(ifine_cquad->loc, (ylen + 1), ylen);

         /* fill in length */
         ifine_cquad->length = ylen;            
         
         /* generate a new ifine_cquad placeholder */
         ifine_cquad->next   = c_cquad(1);         /* add linked fine cquad  */
         ifine_cquad         = ifine_cquad->next;  /* move to linked cquad   */
         
         /* generate another i_fine_nq placeholder */
         ifine_cquad->loc    = c_nquad(1);
         
         /* jump to newly created fin_nquad */
         i_fine_nq           = ifine_cquad->loc;
         
         /* reset ylen */
         ylen                = 0;
        }
     } /* icoa_nquad loop */
   
   /*
    * Deal with last coarse nquad on column (i.e. icoa_nquad = correct exit-value of for-loop!)
    * The following If-Else statement deals
    * with the complicated problem of the last coarse nquad on the column. There
    * are four possibilities:
    * 1. The current cquad is not the last covering the computational domain, and
    *    the refinement flag indicates a refined region is being generated.
    * 2. The current nquad is the last nquad in the computational domain with
    *    ref_nquad producing refinement.
    * 3. Neither of the above are true but the state vector is TRUE and so a
    *    refined region must be ended.
    * 4. All other possibilities - delete unused nquad
    *
    * NOTE:
    *      whatever the situation, the structure of these If-Else statements
    *      entails in any case at least one more call to ref_nquad(...,icoa_nquad,...)
    *      generating a new fin_nquad-linked-list!
    */
   
   /* CASE 1. (there is another linked coa_cquad) 
    * -> conclude current one and perform recursive call to ref_cquad */
   if(coa_cquad->next != NULL)
     {
      /* we are currently refining... */
      if(state == TRUE)  
        {
         /* no finishing call to ref_nquad() -> just conclude! */
                  
         /* conclude refinement */
         ifine_cquad->length = ylen;              /* fill in length of refinement              */         
         ifine_cquad->next   = c_cquad(1);        /* create a linked fine cquad...             */
         ifine_cquad         = ifine_cquad->next; /* ...needed by recursive call to ref_cquad! */
        }
      
      /* CASE 1. -> we have not yet refined (state==FALSE) and won't start now! */
      else
        {
         /* free unused "i_fine_nq==ifine_cquad->loc" placeholder (we already concluded ifine_cquad) */
         dest_nquad(&(ifine_cquad->loc));
        }
      
      /* call ref_cquad recursively for coa_cquad->next 
       * (passing new (or the unused!) ifine_cquad placeholder) */
      refined = ref_cquad(coa_grid, coa_pquad, coa_cquad->next, ifine_cquad, z);
     }

   /* CASE 2. -> we are at the (cospatial) boundary of the box 
    * (whatever the state, do a call to ref_nquad() trying to refine icoa_nquad!) */
   else if((coa_cquad->y + coa_cquad->length == coa_grid->l1dim) &&
           (ref_nquad(coa_grid, coa_pquad, coa_cquad, icoa_nquad, i_fine_nq, z, y) == TRUE))
     {
      /* check "next" pointer not NULL */
      if(i_fine_nq->next != NULL) 
        {
         
         /* find the last nquad along this plane */
         for(rem_nquad = i_fine_nq; (rem_nquad->next)->next != NULL ;
             rem_nquad = rem_nquad->next)
            ;
         
         /* if loc pointer is NULL, destroy */
         if((rem_nquad->next)->loc == NULL) 
           {
            dest_nquad(&(rem_nquad->next));
            rem_nquad->next = NULL;
           }
        }
      
      /* new refinement ? - fill in details */
      if(state == FALSE)    
        {
         /* refined region starts here (possibly because of... */
         /* ...triply-periodic boundary conditions)            */
         refever        = TRUE;
         ifine_cquad->y = (y * 2);
        }
      
      /* update ylen to account for this newly created i_fine_nq-linked-list */
      ylen++;
      
      /* create another i_fine_nq placeholder */
      ifine_cquad->loc = r_nquad(ifine_cquad->loc, ylen, (ylen + 1));
      
      /* jump to newly created i_fine_nq */
      i_fine_nq        = ifine_cquad->loc + ylen;
      
      /* simply copy previous i_fine_nq over to newly created one manually without testing... */
      mirror_nquad(ifine_cquad->loc+(ylen-1), i_fine_nq);
      
      /* update ylen to account for this newly created fin_nquad-linked-list */
      ylen++;
      
      /* fill in the length */
      ifine_cquad->length = ylen;  
     }
   
   /* CASE 3. (same situation as (CASE 1./state==TRUE), but without the recursive call to ref_cquad!) */
   else if(state == TRUE) 
     {
      /* no finishing call to ref_nquad() -> just conclude! */
      
      /* remove the last i_fine_nq placeholder */
      ifine_cquad->loc    = r_nquad(ifine_cquad->loc, (ylen + 1), ylen);
      

      /* fill in length of ifine_cquad */
      ifine_cquad->length = ylen;
     }
   
   /* CASE 4. (same situation as (CASE 1./state==FALSE), but without the recursive call to ref_cquad!) */
   else
     {
      /* free unused "i_fine_nq==ifine_cquad->loc" placeholder (we already concluded ifine_cquad) */
      dest_nquad(&(ifine_cquad->loc));      
     }
   
   /*   printf("ref_cquad:         refever, refined %i %i\n",refever,refined); */
   return (refever || refined);
}

/*-------------------------------------------------------------------------------
* ref_pquad: pquad refinement function. Add blocks of cquads to supplied pquad.
* Generates additional pquads (linked to "next" pointer) as required and calls
* itself recursively.
*-------------------------------------------------------------------------------*/
boolean ref_pquad(gridls *coa_grid, pqptr coa_pquad, pqptr ifine_pquad)
{
   cqptr   icoa_cquad;  /* pointer to coarse cquad, loop pointer                  */
   boolean refined;     /* boolean refinement flag                                */
   boolean refever;     /* boolean refinement flag - recursive flag, have any...  */
                        /* ...previous calls to ref_pquad produced refinement     */
   cqptr   i_fine_cq;   /* current fine cquad                                     */
   boolean state;       /* on/off state flag - indicates whether current fine...  */
                        /* ...cquad is cospatial or not                           */
   cqptr   rem_cquad;   /* pointer to excess cquads - for removal loop            */
   long    z;           /* z-coord of current coarse node                         */
   long    zlen;        /* length of current refinement (z-direction)             */
   
   pqptr   fst_pquad, lst_pquad;  /* do not refine edge nodes and hence calculate "offsets" */
   long    z_low, z_up;
   
   /* initilise state flags */
   refined = FALSE;
   refever = FALSE;
   state   = FALSE;
   
   /* initialise zlen = [1, ...] */
   zlen = 0;
   
   /*
    * The new pquad requires a cquad be created. Further cquads can be added and
    * blocks of nquads may be created and referenced as the refinement routine
    * requires them.
    */
   ifine_pquad->loc = c_cquad(1);      /* create first cquad                       */
   i_fine_cq        = ifine_pquad->loc;/*initialize i_fine_cq to point to new cquad*/
      
   /* we should actually avoid looping over the first plane,
      but in that case we would need to make sure that the
      first column is *not* periodically wrapped... */ 
      
#ifdef PERIODIC_Z
   /* first pquad */
   fst_pquad = coa_grid->pquad;
   
   /* jump to last pquad */
   lst_pquad = fst_pquad;
   while(lst_pquad->next != NULL)
      lst_pquad=lst_pquad->next;

   /* check whether things are periodically wrapped in z-direction */
   if(fst_pquad->z == 0 && lst_pquad->z+lst_pquad->length == coa_grid->l1dim)
     {
      /* current z-refinement runs through whole box */
      if(coa_pquad->z == 0 && coa_pquad->z+coa_pquad->length == coa_grid->l1dim)
        {
         z_low =  0;
         z_up  =  0;
        }
      /* current z-position is at lower boundary */
      else if(coa_pquad->z == 0)
        {
         z_low =  0;
         z_up  = -1;
        }
      /* current z-position will reach upper boundary */
      else if(coa_pquad->z+coa_pquad->length == coa_grid->l1dim)
        {
         z_low = +1;
         z_up  =  0;
        }
      /* we are somewhere in-between */
      else
        {
         z_low = +1;
         z_up  = -1;
        }
     }
   /* no periodic wrapping whatsoever */
   else
     {
      z_low = +1;
      z_up  = -1;
     }
#else /* PERIODIC_Z */
   z_low = +1;
   z_up  = -1;
#endif /* PERIODIC_Z */
   
   /* loop over the coarse planes (but the last one!) */
   for(icoa_cquad = coa_pquad->loc+z_low, z = coa_pquad->z+z_low; 
       icoa_cquad < coa_pquad->loc + (coa_pquad->length - 1) + z_up;
       icoa_cquad++, z++)
     {      
      /*---------------------
      * test cospatial plane
      *----------------------*/
      
      /* call ref_cquad -> generates a new linked-list of i_fine_cq's (if == TRUE)  */
      if(ref_cquad(coa_grid, coa_pquad, icoa_cquad, i_fine_cq, z) == TRUE)
        {
         /* check "next" pointer not NULL */
         if(i_fine_cq->next != NULL) 
           {
            
            /* find the last cquad along this plane */
            for(rem_cquad = i_fine_cq; (rem_cquad->next)->next != NULL ;
                rem_cquad = rem_cquad->next)
               ;
            
            /*if loc pointer is NULL, destroy*/
            if((rem_cquad->next)->loc == NULL)
              {
               dest_cquad(&(rem_cquad->next));
               rem_cquad->next = NULL;
              }
           }
         
         /* cospatial plane => start refinement */
         if(state == FALSE) 
           {
            refever        = TRUE;    /* refinement has occured on this plane     */
            ifine_pquad->z = (z * 2); /* fill in start of refinement              */
            state          = TRUE;    /* refinement state now TRUE - indicates... */
                                      /* ...refinement has already started        */
           }
         
         /* increase length by one */
         zlen++;                     
         
         /* generate another i_fine_cq placeholder */
         ifine_pquad->loc = r_cquad(ifine_pquad->loc, zlen, (zlen + 1));
         
         /* jump to newly created i_fine_cq */
         i_fine_cq      = ifine_pquad->loc + zlen;
         
         /* generate this i_fine_cq-linked-list manually... */
         mirror_cquad(ifine_pquad->loc+(zlen-1), i_fine_cq);
         
         /* update zlen to account for this newly created fin_cquad-linked-list */
         zlen++;         

         /* generate another i_fine_cq placeholder */
         ifine_pquad->loc = r_cquad(ifine_pquad->loc, zlen, (zlen + 1));
         
         /* jump to newly created i_fine_cq */
         i_fine_cq      = ifine_pquad->loc + zlen;
        }
      
      /* we did not refine this cospatial cquad -> conclude current refinement */
      else if(state == TRUE)
        {
         /* set state to FALSE  */
         state = FALSE;
         
         /* remove the last i_fine_cq placeholder */
         ifine_pquad->loc    = r_cquad(ifine_pquad->loc, zlen + 1, zlen);

         /* fill in length */
         ifine_pquad->length = zlen;            
         
         /* generate a new ifine_cquad placeholder */
         ifine_pquad->next   = c_pquad(1);         /* add linked fine cquad  */
         ifine_pquad         = ifine_pquad->next;  /* move to linked cquad   */
         
         /* generate another i_fine_nq placeholder */
         ifine_pquad->loc    = c_cquad(1);
         
         /* jump to newly created fin_nquad */
         i_fine_cq           = ifine_pquad->loc;
         
         /* reset ylen */
         zlen                = 0;
        }
     }
      
   /*
    * Deal with last coarse cquad on plane. (i.e. icoa_cquad = correct exit-value of for-loop!)
    * The following If-Else statement deals
    * with the complicated problem of the last coarse cquad on the plane. There
    * are four possibilities:
    * 1. The current pquad is not the last covering the computational domain, and
    *    the refinement flag indicates a refined region is being generated.
    * 2. The current cquad is the last cquad in the computational domain with
    *    ref_cquad producing refinement.
    * 3. Neither of the above are true but the state vector is TRUE and so a
    *    refined region must be ended.
    * 4. All other possibilities - delete unused cquad
    *
    * NOTE:
    *      whatever the situation, the structure of these If-Else statements
    *      entails in any case at least one more call to ref_cquad(...,icoa_cquad,...)
    *      generating a new fin_cquad-linked-list!
    */

   /* CASE 1. (there is another linked coa_pquad) 
    * -> conclude current one and perform recursive call to ref_pquad */
   if(coa_pquad->next != NULL)
     {
      /* we are currently refining... */
      if(state == TRUE)  
        {
         /* no finishing call to ref_cquad() -> just conclude! */

         /* conclude refinement */
         ifine_pquad->length = zlen;            /* fill in length of refinement */
         
         ifine_pquad->next   = c_pquad(1);      /* create a linked fine pquad   */
         ifine_pquad         = ifine_pquad->next;
        }
      
      /* CASE 1. -> we have not yet refined (state==FALSE) and won't start now! */
      else   
        {
         /* free unused "i_fine_nq==ifine_cquad->loc" placeholder (we already concluded ifine_cquad) */
         dest_cquad(&(ifine_pquad->loc));
        }
      
      /* call ref_pquad recursively for coa_pquad->next */
      refined = ref_pquad(coa_grid, coa_pquad->next, ifine_pquad);
     }

   /* CASE 2. -> we are at the (cospatial) boundary of the box 
    * (whatever the state, do a call to ref_cquad() trying to refine!) */
   else if((coa_pquad->z + coa_pquad->length == coa_grid->l1dim) &&
           (ref_cquad(coa_grid, coa_pquad, icoa_cquad, i_fine_cq, z) == TRUE))
     {
      /* check "next" pointer not NULL */
      if(i_fine_cq->next != NULL) 
        {
         /* find the last cquad along this plane */
         for(rem_cquad = i_fine_cq; (rem_cquad->next)->next != NULL ;
             rem_cquad = rem_cquad->next)
            ;
         
         /* if loc pointer is NULL, destroy */
         if((rem_cquad->next)->loc == NULL)
           {
            dest_cquad(&(rem_cquad->next));
            rem_cquad->next = NULL;
           }
        }
      
      /* new refinement ? - fill in details */
      if(state == FALSE)    
        {
         /* refined region starts here (possibly because of... */
         /* ...triply-periodic boundary conditions)            */
         refever        = TRUE;      /* set refever flag              */
         ifine_pquad->z = (z * 2);   /* enter z coordinate into pquad */
        }
      
      /* update ylen to account for this newly created i_fine_cq-linked-list */
      zlen++;
      
      /* create another i_fine_cq */
      ifine_pquad->loc = r_cquad(ifine_pquad->loc, zlen, (zlen+1));

      /* jump to newly created i_fine_cq */
      i_fine_cq        = ifine_pquad->loc + zlen;
      
      /* simply copy previous i_fine_nq over to newly created one manually without testing... */
      mirror_cquad(ifine_pquad->loc+(zlen-1), i_fine_cq);
      
      /* update ylen to account for this newly created fin_nquad-linked-list */
      zlen++;
         
      /* fill in the length */
      ifine_pquad->length = zlen;   
     }
      
   /* CASE 3. (same situation as (CASE 1./state==TRUE), but without the recursive call to ref_pquad!) */
   else if(state == TRUE)
     {
      /* no finishing call to ref_cquad() -> just conclude! */
      
      /* remove the last i_fine_cq placeholder */
      ifine_pquad->loc    = r_cquad(ifine_pquad->loc, zlen + 1, zlen);
      
      
      /* fill in the length */
      ifine_pquad->length = zlen;   
     }

   /* CASE 4. (same situation as (CASE 1./state==FALSE), but without the recursive call to ref_pquad!) */
   else
     {
      /* destroy unused cquad */
      dest_cquad(&(ifine_pquad->loc));     
     }

   return (refever || refined);
}

/*===============================================================================
* refine_grid: Test coa_grid for refinement, by passing control to ref_pquad. If
* refinement occurs free any excess pquads, else free the first pquad generated
* at start of the function
*===============================================================================*/
boolean refine_grid(gridls *fin_grid, gridls *coa_grid)
{
  pqptr coa_pquad;  /* pquad pointer - old grid                                 */
  pqptr fin_pquad;  /* pquad pointer - new grid                                 */
  pqptr tmp_pquad;  /* pquad pointer - pointer to pquads that will be destroyed */
  boolean refined;  /* refinement flag                                          */
  
  
  /*
   * The new grid requires a pquad be created to other pquads and blocks of
   * cquads may be created and referenced should the refinement routine require
   * them.
   */
  
  fin_grid->pquad    = c_pquad(1);       /* create a pquad                        */
  fin_pquad          = fin_grid->pquad;  /* define the coa_pquad and fin_pquad... */
  coa_pquad          = coa_grid->pquad;  /* ...to point to the first pquads on... */
  
  /* ...the old and new grids respectively */
  
  /*
   * call ref_pquad. If refinement occurs ref_pquad returns the boolean value
   * TRUE or FALSE otherwise.
   */
  refined = ref_pquad(coa_grid, coa_pquad, fin_pquad);
  
  if(refined == TRUE)
    {
      /*
       * Refinement has occured, free potential excess pquad at end of pquad
       * linked list. Loop to end of linked list and test "loc" pointer to test if
       * any cquads are linked to the pquad.
       */
      if(fin_pquad->next != NULL)
        {
          for(tmp_pquad = fin_pquad; (tmp_pquad->next)->next != NULL;
              tmp_pquad = tmp_pquad->next)
            ;
          if((tmp_pquad->next)->loc == NULL)
            {
              dest_pquad(&(tmp_pquad->next));
            }
        }      
    }
  
  /* no refinement hence destroy fin_pquad */
  else 
    {
      dest_pquad(&(fin_grid->pquad));      
    }
    
  
#ifdef WITH_OPENMP
  /* count the number of pquads and make them accessible in a "linear array" */
  /* (this is done in a three-stage process to avoid realloc() calls)        */
  
  /* start from scratch (actually, pquad_array==NULL thanks to free_grid()!? ) */
  if(fin_grid->pquad_array)
    free(fin_grid->pquad_array);
  fin_grid->pquad_array = NULL;
  
  /* 1. count the number of pquad's */
  fin_grid->no_pquad = 0;
  for(tmp_pquad = fin_grid->pquad; tmp_pquad != NULL; tmp_pquad = tmp_pquad->next)
    fin_grid->no_pquad++;
  
  /* 2. allocate linear array to hold all pquad's*/
  fin_grid->pquad_array = (pqptr *) calloc(fin_grid->no_pquad, sizeof(pqptr));
  
  /* 3. copy the pointers into a linear array */
  fin_grid->no_pquad = 0;                    // use this again as a simple counter
  for(tmp_pquad = fin_grid->pquad; tmp_pquad != NULL; tmp_pquad = tmp_pquad->next)
    {
      fin_grid->pquad_array[fin_grid->no_pquad] = tmp_pquad;
      fin_grid->no_pquad++;
    }  
#endif
  
  
  
  /* return refinement flag */
  return (boolean)refined;  
}

