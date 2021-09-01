#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "amr_serial.h"
#include "../libutility/utility.h"

#ifdef EXTRAE_API_USAGE
#include <extrae_user_events.h>
#endif

/*============================================================================
*
* This file contains everything to work with AMIGA's linked-lists
*
*============================================================================*/

/*============================================================================
* build initial linked lists...afterwards this list is only updated
*============================================================================*/
void ll(long unsigned npart, partptr fst_part, gridls *cur_grid) 
{
   partptr       ipart;
   long unsigned jpart;
   unsigned long n[NDIM];     /* particle - node coords        */
   int           i;           /* increment index               */   
   pqptr         cur_pquad;    /* current pquad                 */
   cqptr         cur_cquad;    /* current cquad                 */
   nqptr         cur_nquad;    /* current nquad                 */
   nptr          cur_node;     /* current node                  */

#ifdef EXTRAE_API_USAGE
  Extrae_user_function(1);
#endif

   /* loop over all particles starting with fst_part (first particle) */
   for(jpart = 0; jpart < npart; jpart++)
     {
      ipart = fst_part + jpart;
      
      /*     fprintf(stderr,"ipart = %ld\n",ipart-fst_part); */
#ifdef VERBOSE
      if ( (long unsigned)(ipart-fst_part) > (long unsigned)(fst_part+npart) )
         fprintf(stderr,"%ld %ld %ld\n",ipart-fst_part,npart,jpart);
#endif
      
      /* find cell containing particle */
      for(i = X; i <= Z; i++) 
        {
         n[i] = (unsigned long) ((double)cur_grid->l1dim * (ipart->pos[i]));
         
         /* try to catch 'unperiodic' posititons... */
         if(n[i] > cur_grid->l1dim-1)  n[i] = 0;
         /*if(n[i] < 0)              n[i] = cur_grid->l1dim-1;*/
        }
      
      /*****************************************************************/
      /*****************************************************************/
      
      /* find correct plane */
      cur_pquad  = cur_grid->pquad;       /* pointer to corresponding pquad            */
      cur_cquad  = cur_pquad->loc;     /* pointer to first cquad within that pquad  */
      cur_cquad += n[Z];              /* pointer to plane containing node          */
      
      /* find correct column */
      cur_nquad  = cur_cquad->loc;     /* pointer to first nquad in that cquad      */
      cur_nquad += n[Y];              /* pointer to nquad containing node          */
      
      /* find correct node */
      cur_node   = cur_nquad->loc;     /* pointer to first node in that nquad       */
      cur_node  += n[X];              /* pointer to node containing particle ipart */
      
      /* simply insert at beginning of linked list */
      ipart->ll   = cur_node->ll;
      cur_node->ll = ipart;

     }
#ifdef EXTRAE_API_USAGE
  Extrae_user_function(0);
#endif
}


/*============================================================================
* initialize ll pointers to NULL 
*============================================================================*/
void NULL_ll(gridls *cur_grid)
{
   long    icquad;
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;  
   
   if(cur_grid->l1dim != global.dom_grid->l1dim)
     {
      fprintf(io.logfile,"you are calling NULL_ll() for other than the domain grid!\nEXIT\n");
      exit(0);
     }
   
   for(cur_pquad=cur_grid->pquad;cur_pquad!=NULL;cur_pquad=cur_pquad->next)
      
#ifdef WITH_OPENMP
#pragma omp parallel private(cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node) shared(cur_pquad)
#pragma omp for schedule(static)
      for(icquad=0; icquad<cur_grid->l1dim; icquad++)
        {
         cur_cquad = cur_pquad->loc+icquad;
#else
      for(cur_cquad=cur_pquad->loc;cur_cquad<cur_pquad->loc+cur_pquad->length; 
          cur_cquad++)
        {
#endif
         for(icur_cquad=cur_cquad;icur_cquad!=NULL;icur_cquad=icur_cquad->next)
            
            for(cur_nquad=icur_cquad->loc;cur_nquad<icur_cquad->loc+icur_cquad->length; 
                cur_nquad++)
               for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
                  
                  for(cur_node=icur_nquad->loc;
                      cur_node<icur_nquad->loc+icur_nquad->length;cur_node++)
                     
                     cur_node->ll = NULL;
        }
}

/*============================================================================
* initialize force.new_ll pointers to NULL 
* (otherwise they can't be used !!!)
*============================================================================*/
void NULL_newll(gridls *cur_grid)
{
#ifndef AHFlean
   long    ipquad;
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;  
   
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node) shared(cur_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)
     {
      cur_pquad = cur_grid->pquad_array[ipquad];
#else
   for(cur_pquad=cur_grid->pquad;cur_pquad!=NULL;cur_pquad=cur_pquad->next)
     {
#endif
      for(cur_cquad=cur_pquad->loc;cur_cquad<cur_pquad->loc+cur_pquad->length; 
          cur_cquad++)
         for(icur_cquad=cur_cquad;icur_cquad!=NULL;icur_cquad=icur_cquad->next)
            
            for(cur_nquad=icur_cquad->loc;cur_nquad<icur_cquad->loc+icur_cquad->length; 
                cur_nquad++)
               for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
                  
                  for(cur_node=icur_nquad->loc;
                      cur_node<icur_nquad->loc+icur_nquad->length;cur_node++)
                     
                     cur_node->force.new_ll = NULL;
     }
#else
   fprintf(stderr, "%s should not be called when AHFlean is defined.\n",
           __func__);
   exit(1);
#endif /* AHFlean*/
}

/*===============================================================================
* try to find cur_part to in cur_grid...and maybe next coarser grid
* Note: routine modifies curpart->ll while keeping cur_part unchanged !!!
* (to be used only with move_part() because of cur_node->forces.new_ll !!!)
*===============================================================================*/
void extended_llsearch(partptr cur_part, gridls *cur_grid)
{
   long   n[NDIM];             /* particle - node coords        */
   int    i;                   /* increment index               */
   pqptr  cur_pquad;           /* current pquad                 */
   cqptr  cur_cquad;           /* current cquad                 */
   nqptr  cur_nquad;           /* current nquad                 */
   nptr   cur_node;            /* current node                  */
   int    exp;                 /* exponent for frexp            */
   double edge_shift;          /* coords of edges of first cell */
   
#ifdef VERBOSELOG2
   fprintf(io.logfile,"extended_llsearch: trying to locate particle %ld on grid %ld\n",
           cur_part-global.fst_part,cur_grid->l1dim);
   fflush(io.logfile);
#endif
   
   /* find cell containing particle */
   edge_shift = 0.0;
   
   /* coordinates of node cur_part needs to be linked to */
   for(i = X; i <= Z; i++)
      n[i] = (long) ((double)cur_grid->l1dim * f1mod((double)cur_part->pos[i]-edge_shift+1.0, 1.0));
   
#ifdef VERBOSELOG2
   fprintf(io.logfile,"                   %g %g %g (%g %g %g)  --- %ld %ld %ld\n",
           cur_part->pos[X],cur_part->pos[Y],cur_part->pos[Z],
           cur_part->mom[X],cur_part->mom[Y],cur_part->mom[Z],
           n[X],n[Y],n[Z]);
   fflush(io.logfile);
#endif
   /*-----------------------------------
      * find this node within cur_grid... 
      *-----------------------------------*/
   for(cur_pquad = cur_grid->pquad;
       (n[Z] >= cur_pquad->z + cur_pquad->length) && (cur_pquad->next != NULL); 
       cur_pquad = cur_pquad->next)
      ;
   
   if((n[Z] >= cur_pquad->z) && (n[Z] < cur_pquad->z + cur_pquad->length))
     {
      for(cur_cquad = cur_pquad->loc + (n[Z] - cur_pquad->z);
          (n[Y] >= cur_cquad->y + cur_cquad->length) && (cur_cquad->next != NULL); 
          cur_cquad = cur_cquad->next)
         ;
      
      if((n[Y] >= cur_cquad->y) && (n[Y] < cur_cquad->y + cur_cquad->length))
        {
         for(cur_nquad = cur_cquad->loc + (n[Y] - cur_cquad->y);
             (n[X] >= cur_nquad->x + cur_nquad->length)&&(cur_nquad->next != NULL);
             cur_nquad = cur_nquad->next)
            ;
         
         if((n[X] >= cur_nquad->x) && (n[X] < cur_nquad->x + cur_nquad->length))
           {
            /*------------------------------------
            * eventually reached correct node...
            *------------------------------------*/
            cur_node = cur_nquad->loc + (n[X] - cur_nquad->x);
            
            /* simply insert particle at beginning of linked list */
            cur_part->ll           = cur_node->force.new_ll;
            cur_node->force.new_ll = cur_part;
           }
        }
     }
}


/*==============================================================================
* try to link cur_part to cur_grid...otherwise terminate AMIGA
* Note: routine modifies curpart->ll while keeping cur_part unchanged !!!
* (not to be used with move_part() because of cur_node->forces.new_ll)
*==============================================================================*/
void coagrid_llsearch(partptr cur_part, gridls *cur_grid)
{
   long   n[NDIM];             /* particle - node coords        */
   int    i;                   /* increment index               */   
   pqptr  cur_pquad;           /* current pquad                 */
   cqptr  cur_cquad;           /* current cquad                 */
   nqptr  cur_nquad;           /* current nquad                 */
   nptr   cur_node;            /* current node                  */
   int    exp;                 /* exponent for frexp            */
   double edge_shift;          /* coords of edges of first cell */
   
#ifdef VERBOSELOG2
   fprintf(io.logfile,"coagrid_llsearch: trying to locate particle %ld on %ld grid\n",cur_part-global.fst_part,(cur_grid+1)->l1dim/2);
   fflush(io.logfile);
#endif
   
   /* find cell containing particle */
   edge_shift = 0.0;
   
   /* coordinates of node cur_part needs to be linked to */
   for(i = X; i <= Z; i++)
      n[i] = (long) ((double)cur_grid->l1dim * f1mod((double)cur_part->pos[i]-edge_shift+1.0, 1.0));
   
   /*-----------------------------------
    * find this node within cur_grid... 
    *-----------------------------------*/
   for(cur_pquad = cur_grid->pquad; 
       (n[Z] >= cur_pquad->z + cur_pquad->length) && (cur_pquad->next != NULL); 
       cur_pquad = cur_pquad->next)
      ;
   
   if((n[Z] >= cur_pquad->z) && (n[Z] < cur_pquad->z + cur_pquad->length))
     {
      for(cur_cquad = cur_pquad->loc + (n[Z] - cur_pquad->z); 
          (n[Y] >= cur_cquad->y + cur_cquad->length) && (cur_cquad->next != NULL); 
          cur_cquad = cur_cquad->next)
         ;
      
      if((n[Y] >= cur_cquad->y) && (n[Y] < cur_cquad->y + cur_cquad->length))
        {
         for(cur_nquad = cur_cquad->loc + (n[Y] - cur_cquad->y); 
             (n[X] >= cur_nquad->x + cur_nquad->length)&&(cur_nquad->next != NULL); 
             cur_nquad = cur_nquad->next)
            ;
         
         if((n[X] >= cur_nquad->x) && (n[X] < cur_nquad->x + cur_nquad->length))
           {
            /*------------------------------------
            * eventually reached correct node...
            *------------------------------------*/
            cur_node = cur_nquad->loc + (n[X] - cur_nquad->x);
            
            /* simply insert particle at beginning of linked list */
            cur_part->ll = cur_node->ll;
            cur_node->ll = cur_part;
           }
         else
           {
#ifdef VERBOSELOG
            fprintf(io.logfile,
                    "\ncoagrid_llsearch: particle %ld not in %ld grid (x): %g %g %g\n",
                    cur_part-global.fst_part,
                    cur_grid->l1dim,
                    cur_part->pos[X],cur_part->pos[Y],cur_part->pos[Z]);
            fflush(io.logfile);
#endif
#ifdef TERMINATE
            global.terminate = TRUE;
#endif
            coagrid_llsearch(cur_part, (cur_grid-1));
           }
        }
      else
        {
#ifdef VERBOSELOG
         fprintf(io.logfile,
                 "\ncoagrid_llsearch: particle %ld not in %ld grid (y): %g %g %g\n",
                 cur_part-global.fst_part,
                 cur_grid->l1dim,
                 cur_part->pos[X],cur_part->pos[Y],cur_part->pos[Z]);
         fflush(io.logfile);
#endif
#ifdef TERMINATE
         global.terminate = TRUE;
#endif
         coagrid_llsearch(cur_part, (cur_grid-1));
        }
     }
   else
     {
#ifdef VERBOSELOG
      fprintf(io.logfile,
              "\ncoagrid_llsearch: particle %ld not in %ld grid (z): %g %g %g\n",
              cur_part-global.fst_part,
              cur_grid->l1dim,
              cur_part->pos[X],cur_part->pos[Y],cur_part->pos[Z]);
      fflush(io.logfile);
#endif
#ifdef TERMINATE
      global.terminate = TRUE;
#endif
      coagrid_llsearch(cur_part, (cur_grid-1));
     }
}

/*===============================================================================
* insert cur_part into link list at node new_node
* (returns in cur_part the pointer to the next particle !)
*===============================================================================*/
void add_part_to_ll(nptr new_node, partptr cur_part)
{
   /* simply insert particle at beginning of new linked list */
   cur_part->ll           = new_node->force.new_ll;
   new_node->force.new_ll = cur_part;
   cur_part               = cur_part->ll;
}

