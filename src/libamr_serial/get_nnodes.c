#include <stddef.h>

/* the important definitions have to be included first */
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "amr_serial.h"
#include "../libutility/utility.h"

/*===============================================================================
*
* This file contains the routines to access neighboring nodes
*
* The visible functions are:
*
* get_TSCnodes:     search for all 27  surrounding nodes
* get_MHDnodes:     search for all 125 surrounding nodes
*
*
* Note:   1. x, y,z are the coordinates of that node (in grid units)
*         2. cur_pquad, cur_cquad, cur_nquad are supposed to be the 
*            quad's down the linked list where cur_node lies in...
*
*===============================================================================*/

/* test if all mhd_nodes are present */
boolean test_mhd(nptr mhd_nodes[5][5][5])
{
   int i,j,k;
   
   for(i = 0; i < 5 ; i++)
      for(j = 0; j < 5 ; j++)
         for(k = 0; k < 5 ; k++)
            if(mhd_nodes[i][j][k] == NULL)
               return FALSE;
   return TRUE;
}

/* test if all tsc_nodes are present */
boolean test_tsc(nptr tsc_nodes[3][3][3])
{
   int i,j,k;
   
   for(i = 0; i < 3 ; i++)
      for(j = 0; j < 3 ; j++)
         for(k = 0; k < 3 ; k++)
            if(tsc_nodes[i][j][k] == NULL)
               return FALSE;
   return TRUE;
}

/*===============================================================================
 * count number of actual TSC neighbors
 *===============================================================================*/
int count_tsc(nptr tsc_nodes[3][3][3])
{
   int i,j,k;
   int nn;
   
   nn = 0;
   for(i = 0; i < 3 ; i++)
      for(j = 0; j < 3 ; j++)
         for(k = 0; k < 3 ; k++)
            if(tsc_nodes[i][j][k] != NULL)
               nn++;
   return (nn-1);   /* at least tsc_nodes[1][1][1] != NULL */
}

/*===============================================================================
 * get pointers to six direct neighbors
 *===============================================================================*/
void get_SIXnodes(gridls *cur_grid, pqptr cur_pquad, cqptr cur_cquad, 
                  nqptr cur_nquad, nptr tsc_nodes[3][3][3], int z, int y, int x)
{
   pqptr  s_pquad;      /* pquad to search for nodes */
   cqptr  s_cquad;      /* cquad to search for nodes */
   nqptr  s_nquad;      /* nquad to search for nodes */
   nptr   cur_node;
   double endx;
   double endy;
 
   /* set current node */
   cur_node = tsc_nodes[1][1][1];
  
#ifdef PERIODIC_Z
   /* get nnode[Z][1] */
   /* periodic case */
   
   if(z == cur_grid->l1dim - 1)
     {
      s_pquad = cur_grid->pquad;
      if(s_pquad->z != 0)     /* first plane is not present */
         return;
      
      /* find correct cquad */
      for(s_cquad = s_pquad->loc;
          (y >= (endy = (s_cquad->y+s_cquad->length ))) && (s_cquad->next != NULL); 
          s_cquad = s_cquad->next)
         ;
      
      if((y >= s_cquad->y) && (y < endy))/* is nquad in this cquad */
        {
         for(s_nquad = s_cquad->loc + (y - s_cquad->y);
             (x >= (endx = (s_nquad->x+s_nquad->length ))) && (s_nquad->next!=NULL); 
             s_nquad = s_nquad->next)
            ;
         
         if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
            tsc_nodes[2][1][1] = s_nquad->loc + (x - s_nquad->x);
         else
            return;
        }
      else
         return;
     }
   /* non-periodic case */
   else if(z < cur_pquad->z + cur_pquad->length - 1)
     {
      for(s_cquad = cur_pquad->loc + (z - cur_pquad->z + 1);
          (y >= (endy = (s_cquad->y+s_cquad->length ))) && (s_cquad->next != NULL); 
          s_cquad = s_cquad->next)
         ;
      
      if((y >= s_cquad->y) && (y < endy)) /* is node in this nquad */
        {
         for(s_nquad = s_cquad->loc + (y - s_cquad->y);
             (x >= (endx=(s_nquad->x+s_nquad->length ))) && (s_nquad->next!=NULL); 
             s_nquad = s_nquad->next)
            ;
         
         if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
            tsc_nodes[2][1][1] = s_nquad->loc + (x - s_nquad->x);
         else
            return;
        }
      else
         return;
     }
   else
      return;
   
   /* get nnode[Z][0] */
   /* periodic case */
   if(z == 0)    /* we are at min boundary */
     {
      
      /* search for correct node */
      for(s_pquad = cur_pquad; s_pquad->next != NULL; s_pquad = s_pquad->next)
         ;
      
      if(s_pquad->z + s_pquad->length == cur_grid->l1dim)
        {
         for(s_cquad = s_pquad->loc + s_pquad->length - 1;
             (y >= (endy=(s_cquad->y+s_cquad->length))) && (s_cquad->next!=NULL); 
             s_cquad = s_cquad->next)
            ;
         
         if((y >= s_cquad->y) && (y < endy)) /* is cquad  in this pquad */
           {
            for(s_nquad = s_cquad->loc + (y - s_cquad->y);
                (x >= (endx=(s_nquad->x+s_nquad->length)))&&(s_nquad->next!=NULL);
                s_nquad = s_nquad->next)
               ;
            
            if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
               tsc_nodes[0][1][1] = s_nquad->loc + (x - s_nquad->x);
            else
               return;
           }
         else
            return;
        }
      else
         return;
     }
   /* non-periodic case */
   else if(z > cur_pquad->z)
     {
      for(s_cquad = cur_pquad->loc + (z - cur_pquad->z - 1);
          (y >= (endy = (s_cquad->y + s_cquad->length))) && (s_cquad->next != NULL);
          s_cquad = s_cquad->next)
         ;
      
      if((y >= s_cquad->y) && (y < endy))/* is nquad in this cquad */
        {
         for(s_nquad = s_cquad->loc + (y - s_cquad->y);
             (x >= (endx=(s_nquad->x+s_nquad->length))) && (s_nquad->next != NULL);
             s_nquad = s_nquad->next)
            ;
         
         if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
            tsc_nodes[0][1][1] = s_nquad->loc + (x - s_nquad->x);
         else
            return;
        }
      else
         return;
     }
   else
      return;
   
   /* get nnode[Y][1] */
   /* search for correct node */
   /* move to first cquad along plane */
   
   if(y == cur_grid->l1dim - 1)
     {
      s_cquad = cur_pquad->loc + (z - cur_pquad->z);
      
      if(s_cquad->y == 0)
        {
         for(s_nquad = s_cquad->loc;
             (x >= (endx = (s_nquad->x + s_nquad->length))) && 
             (s_nquad->next != NULL);
             s_nquad = s_nquad->next)
            ;
         
         if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
            tsc_nodes[1][2][1] = s_nquad->loc + (x - s_nquad->x);
         else
            return;
        }
      else
         return;
     }
   else if(y < cur_cquad->y + cur_cquad->length - 1)
     {
      for(s_nquad = cur_cquad->loc + (y - cur_cquad->y + 1);
          (x >= (endx = (s_nquad->x + s_nquad->length))) && (s_nquad->next != NULL);
          s_nquad = s_nquad->next)
         ;
      
      if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
         tsc_nodes[1][2][1] = s_nquad->loc + (x - s_nquad->x);
      else
         return;
     }
   else
      return;
   
   /* get nnode[Y][0] */
   if(y == 0)    /* we are at min boundary */
     {
      /* search for correct node */
      for(s_cquad = cur_pquad->loc + (z - cur_pquad->z);
          s_cquad->next != NULL; s_cquad = s_cquad->next)
         ;
      
      if(s_cquad->y + s_cquad->length == cur_grid->l1dim)
        {
         for(s_nquad = s_cquad->loc + s_cquad->length - 1;
             (x >= (endx = (s_nquad->x+s_nquad->length))) && (s_nquad->next!=NULL);
             s_nquad = s_nquad->next)
            ;
         
         if((x >= s_nquad->x) && (x < endx)) /* is node in this nquad */
           {
            tsc_nodes[1][0][1] = s_nquad->loc + (x - s_nquad->x);
           }
         else
            return;
        }
      else
         return;
     }
   else if(y > cur_cquad->y)
     {
      for(s_nquad = cur_cquad->loc + (y - cur_cquad->y - 1);
          (x >= (endx = (s_nquad->x + s_nquad->length))) && (s_nquad->next != NULL);
          s_nquad = s_nquad->next)
         ;
      
      if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
         tsc_nodes[1][0][1] = s_nquad->loc + (x - s_nquad->x);
      else
         return;
     }
   else
      return;
   
   /* get nnode[X][1] */
   if(x == cur_grid->l1dim - 1)
     {
      s_nquad = cur_cquad->loc + (y - cur_cquad->y);
      
      if(s_nquad->x == 0)
         tsc_nodes[1][1][2] = s_nquad->loc;
      else
         return;
     }
   else if(x < cur_nquad->x + cur_nquad->length - 1)
      tsc_nodes[1][1][2] = cur_node + 1;
   else
      return;
   
   /* get nnode[X][0] */
   if(x == 0)
     {
      for(s_nquad = cur_nquad; s_nquad->next != NULL; s_nquad = s_nquad->next)
         ;
      
      if(s_nquad->x + s_nquad->length == cur_grid->l1dim)
         tsc_nodes[1][1][0] = s_nquad->loc + (s_nquad->length - 1);
      else
         return;
     }
   else if(x > cur_nquad->x)
      tsc_nodes[1][1][0] = cur_node - 1;
   else
      return;
   
#else  /* PERIODIC_Z */
   
   /* get nnode[Z][1] */
   /* periodic case */  
   if(z == cur_grid->l1dim - 1)
      return;
   else if(z < cur_pquad->z + cur_pquad->length - 1)
     {
      for(s_cquad = cur_pquad->loc + (z - cur_pquad->z + 1);
          (y >= (endy = (s_cquad->y+s_cquad->length ))) && (s_cquad->next != NULL); 
          s_cquad = s_cquad->next)
         ;
      
      if((y >= s_cquad->y) && (y < endy)) /* is node in this nquad */
        {
         for(s_nquad = s_cquad->loc + (y - s_cquad->y);
             (x >= (endx=(s_nquad->x+s_nquad->length ))) && (s_nquad->next!=NULL); 
             s_nquad = s_nquad->next)
            ;
         
         if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
            tsc_nodes[2][1][1] = s_nquad->loc + (x - s_nquad->x);
         else
            return;
        }
      else
         return;
     }
   else
      return;
   
   /* get nnode[Z][0] */
   /* periodic case */
   if(z == 0)    /* we are at min boundary */
      return;
   else if(z > cur_pquad->z)
     {
      for(s_cquad = cur_pquad->loc + (z - cur_pquad->z - 1);
          (y >= (endy = (s_cquad->y + s_cquad->length))) && (s_cquad->next != NULL);
          s_cquad = s_cquad->next)
         ;
      
      if((y >= s_cquad->y) && (y < endy))/* is nquad in this cquad */
        {
         for(s_nquad = s_cquad->loc + (y - s_cquad->y);
             (x >= (endx=(s_nquad->x+s_nquad->length))) && (s_nquad->next != NULL);
             s_nquad = s_nquad->next)
            ;
         
         if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
            tsc_nodes[0][1][1] = s_nquad->loc + (x - s_nquad->x);
         else
            return;
        }
      else
         return;
     }
   else
      return;
   
   /* get nnode[Y][1] */
   /* search for correct node */
   /* move to first cquad along plane */
   
   if(y == cur_grid->l1dim - 1)
      return;
   else if(y < cur_cquad->y + cur_cquad->length - 1)
     {
      for(s_nquad = cur_cquad->loc + (y - cur_cquad->y + 1);
          (x >= (endx = (s_nquad->x + s_nquad->length))) && (s_nquad->next != NULL);
          s_nquad = s_nquad->next)
         ;
      
      if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
         tsc_nodes[1][2][1] = s_nquad->loc + (x - s_nquad->x);
      else
         return;
     }
   else
      return;
   
   /* get nnode[Y][0] */
   if(y == 0)    /* we are at min boundary */
      return;
   else if(y > cur_cquad->y)
     {
      for(s_nquad = cur_cquad->loc + (y - cur_cquad->y - 1);
          (x >= (endx = (s_nquad->x + s_nquad->length))) && (s_nquad->next != NULL);
          s_nquad = s_nquad->next)
         ;
      
      if((x >= s_nquad->x) && (x < endx))/* is node in this nquad */
         tsc_nodes[1][0][1] = s_nquad->loc + (x - s_nquad->x);
      else
         return;
     }
   else
      return;
   
   /* get nnode[X][1] */
   if(x == cur_grid->l1dim - 1)
      return;
   else if(x < cur_nquad->x + cur_nquad->length - 1)
      tsc_nodes[1][1][2] = cur_node + 1;
   else
      return;
   
   /* get nnode[X][0] */
   if(x == 0)
      return;
   else if(x > cur_nquad->x)
      tsc_nodes[1][1][0] = cur_node - 1;
   else
      return;
#endif /* PERIODIC_Z */
   
   /* return */
   return;
}

/*--------------------------------------------------------------------------------
* NULL_TSCplane: set all search pointers on plane to NULL 
*--------------------------------------------------------------------------------*/
void NULL_TSCplane(nptr tsc_plane[3][3])
{
   int i,j;         /* loop indices */
   
   for(i = 0; i < 3 ; i++)
      for(j =0; j < 3 ; j++)
         tsc_plane[j][i] = NULL;
}

/*--------------------------------------------------------------------------------
* NULL_TSCcolumnumn: set all search pointers on column to NULL 
*--------------------------------------------------------------------------------*/
void NULL_TSCcolumn(nptr tsc_col[3])
{
   int i;         /* loop index */
   
   for(i = 0; i < 3 ; i++)
      tsc_col[i] = NULL;
}

/*--------------------------------------------------------------------------------
* search_TSCcolumn: search for cur_node--, cur_node++
*--------------------------------------------------------------------------------*/
void search_TSCcolumn(gridls *cur_grid, cqptr cur_cquad, nqptr cur_nquad, 
                   nptr tsc_col[3], long y, long x)
{
   nqptr fstnquad;      /* first nquad on a row */
   
   /* test if cur_node++ is really there */
   if(x < (cur_nquad->x + cur_nquad->length - 1))
     {
      tsc_col[2] = tsc_col[1] + 1;
     }
   
#ifdef PERIODIC_X
   /* deal with periodicity */
   else if(x == (cur_grid->l1dim - 1))
     {
      fstnquad = cur_cquad->loc + (y - cur_cquad->y);
      if(fstnquad->x == 0)
        {
         tsc_col[2] = fstnquad->loc;
        }
      else
         tsc_col[2] = NULL;
     }
#endif /* PERIODIC_X */
   
   else
      tsc_col[2] = NULL;
   
   /* test if cur_node-- is really there */
   if(x > cur_nquad->x)
     {
      tsc_col[0] = tsc_col[1] - 1;
     }
   
#ifdef PERIODIC_X
   /* deal with periodicity */
   else if(x == 0)
     {
      for(fstnquad = cur_nquad; fstnquad->next != NULL; fstnquad = fstnquad->next)
         ;
      
      if((fstnquad->x + fstnquad->length) == cur_grid->l1dim)
         tsc_col[0] = fstnquad->loc + (fstnquad->length - 1);
      else
         tsc_col[0] = NULL;
     }
#endif /* PERIODIC_X */
   
   else
      tsc_col[0] = NULL;
}

/*-------------------------------------------------------------------------------
* search_TSCplane: test the six nodes on same plane - left, right, forward
*-------------------------------------------------------------------------------*/
void search_TSCplane(gridls *cur_grid, pqptr cur_pquad, cqptr cur_cquad, nqptr cur_nquad,
                  nptr tsc_plane[3][3], long z, long y, long x)
{
   cqptr icur_cquad;         /* current coarse cquad    */
   nqptr icur_nquad;         /* current coarse nquad    */
   nptr icur_node;           /* current coarse node     */
   long endx;                /* end coord of icur_nquad */
   long endy;                /* end coord of icur_cquad */
   
   /* 1. test column with cur_node */
   search_TSCcolumn(cur_grid, cur_cquad, cur_nquad, tsc_plane[1], y, x);
   
   
   
   /* 2. test column to left */
   
   /* is the next column there ? */
   if(y < ((cur_cquad->y) + (cur_cquad->length) - 1))
     {
      /* search for correct node */
      for(icur_nquad = cur_cquad->loc + (y - cur_cquad->y) + 1;
          (x >= (endx = (icur_nquad->x + icur_nquad->length))) &&
          (icur_nquad->next != NULL); 
          icur_nquad = icur_nquad->next)
         ;
      
      /* is node in this nquad */
      if((x >=  icur_nquad->x) && (x < endx))
        {
         icur_node       = icur_nquad->loc + (x - icur_nquad->x);
         tsc_plane[2][1] = icur_node;     /* assign middle node */
         search_TSCcolumn(cur_grid, cur_cquad, icur_nquad, tsc_plane[2],
                       modulo(y + 1, cur_grid->l1dim), x);
        }
      else
         NULL_TSCcolumn(tsc_plane[2]);
     }
   
#ifdef PERIODIC_Y
   /* deal with periodic case to left */    
   else if(y == (cur_grid->l1dim - 1))     /* we are at max boundary */
     {
      /* search for correct node */
      
      /* move to first cquad along plane */
      icur_cquad = cur_pquad->loc + (z - cur_pquad->z);
      if(icur_cquad->y == 0)  /* if 1st nquad is on boundary */
        {
         for(icur_nquad = icur_cquad->loc;
             (x >= (endx = (icur_nquad->x + icur_nquad->length ))) &&
             (icur_nquad->next != NULL); 
             icur_nquad = icur_nquad->next)
            ;
         
         /* is node in this nquad */
         if((x >= icur_nquad->x) && (x < endx))
           {
            icur_node       = icur_nquad->loc + (x - icur_nquad->x);
            tsc_plane[2][1] = icur_node;     /* assign middle node */
            search_TSCcolumn(cur_grid, icur_cquad, icur_nquad, tsc_plane[2], 0, x);
           }
         else
            NULL_TSCcolumn(tsc_plane[2]);
        }
      else
         NULL_TSCcolumn(tsc_plane[2]);
     }
#endif /* PERIODIC_Y */
   
   else
      NULL_TSCcolumn(tsc_plane[2]); 
   
   
   
   /* 3. test column to right */
   
   /* is the next column there ? */
   if(y > cur_cquad->y)
     {
      /* search for correct node */
      for(icur_nquad = cur_cquad->loc + (y - cur_cquad->y - 1);
          (x >= (endx = (icur_nquad->x + icur_nquad->length ))) &&
          (icur_nquad->next != NULL); 
          icur_nquad = icur_nquad->next)
         ;
      
      /* is node in this nquad */
      if((x >= icur_nquad->x) && (x < endx))
        {
         icur_node       = icur_nquad->loc + (x - icur_nquad->x);
         tsc_plane[0][1] = icur_node;     /* assign middle node */
         search_TSCcolumn(cur_grid, cur_cquad, icur_nquad, tsc_plane[0],
                       modulo(y - 1, cur_grid->l1dim), x);
        }
      else
         NULL_TSCcolumn(tsc_plane[0]);
      
     }
   
#ifdef PERIODIC_Y
   /* deal with periodic case to right */
   else if(y == 0)    /* we are at min boundary */
     {
      /* search for correct node */
      for(icur_cquad = cur_pquad->loc + (z - cur_pquad->z);
          (( cur_grid->l1dim - 1) >= (endy =(icur_cquad->y + icur_cquad->length))) &&
          (icur_cquad->next != NULL); 
          icur_cquad = icur_cquad->next)
         ;
      
      /* is nquad in this cquad */
      if(((cur_grid->l1dim - 1) >= icur_cquad->y) && ((cur_grid->l1dim-1) < endy))
        {
         for(icur_nquad = icur_cquad->loc + ((cur_grid->l1dim-1) - icur_cquad->y);
             (x >= (endx = (icur_nquad->x + icur_nquad->length ))) &&
             (icur_nquad->next != NULL); 
             icur_nquad = icur_nquad->next)
            ;
         
         /* is node in this nquad */
         if((x >= icur_nquad->x) && (x < endx))
           {
            icur_node       = icur_nquad->loc + (x - icur_nquad->x);
            tsc_plane[0][1] = icur_node;     /* assign middle node */
            search_TSCcolumn(cur_grid, icur_cquad, icur_nquad, tsc_plane[0],
                          (cur_grid->l1dim - 1), x);
           }
         else
            NULL_TSCcolumn(tsc_plane[0]);
        }
      else
         NULL_TSCcolumn(tsc_plane[0]);
     }
#endif /* PERIODIC_Y */
   
   else
      NULL_TSCcolumn(tsc_plane[0]);
}

/*================================================================================
*   get_TSCnodes: search for 26 surrounding nodes
*================================================================================*/
void get_TSCnodes(gridls *cur_grid, pqptr cur_pquad, cqptr cur_cquad, nqptr cur_nquad, 
                  nptr tsc_box[3][3][3], long *z, long *y, long *x)
{
   pqptr icur_pquad;         /* current coarse pquad    */
   cqptr icur_cquad;         /* current coarse cquad    */
   nqptr icur_nquad;         /* current coarse nquad    */
   nptr  icur_node;          /* current coarse node     */
   long  endx;               /* end coord of icur_nquad */
   long  endy;               /* end coord of icur_cquad */
   long  endz;               /* end coord of icur_cquad */
   
   /* 1. search plane with cur_node */
   search_TSCplane(cur_grid, cur_pquad, cur_cquad, cur_nquad, tsc_box[1], *z, *y, *x);
   
   
   /* 2. search plane below cur_node */
   
   /* is plane there ? */
   if(*z > cur_pquad->z)
     {
      /* now look for correct cquad */
      
      /* increment down cquad ll */
      for(icur_cquad = (cur_pquad->loc + (*z - 1 - cur_pquad->z));
          (*y >= (endy = (icur_cquad->y + icur_cquad->length)))
          && (icur_cquad->next != NULL); 
          icur_cquad = icur_cquad->next)
         ;
      
      /* is nquad in this cquad */
      if((*y >=  icur_cquad->y) && (*y < endy))
        {
         /* now look for correct nquad */
         
         /* increment down nquad ll */
         for(icur_nquad = icur_cquad->loc + (*y - icur_cquad->y);
             (*x >= (endx = (icur_nquad->x + icur_nquad->length))) &&
             (icur_nquad->next != NULL); 
             icur_nquad = icur_nquad->next)
            ;
         
         /* is node in this nquad */
         if((*x >=  icur_nquad->x) && (*x < endx))
           {
            icur_node        = icur_nquad->loc + (*x - icur_nquad->x);
            tsc_box[0][1][1] = icur_node;     /* assign middle node */
            search_TSCplane(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_box[0],
                            modulo(*z - 1, cur_grid->l1dim), *y, *x);
           }
         else
            NULL_TSCplane(tsc_box[0]);
        }
      else
        {
         NULL_TSCplane(tsc_box[0]);
        }
     }
   
#ifdef PERIODIC_Z
   /* periodic case */
   else if(*z == 0)
     {
      /* search for correct nquad */
      for(icur_pquad = cur_grid->pquad;
          (( cur_grid->l1dim - 1) >= (endz = (icur_pquad->z + icur_pquad->length)))
          && (icur_pquad->next != NULL); 
          icur_pquad = icur_pquad->next)
         ;
      
      /* is cquad in this pquad */
      if(((cur_grid->l1dim - 1) >= icur_pquad->z) && ((cur_grid->l1dim-1) < endz))
        {
         for(icur_cquad = icur_pquad->loc + ((cur_grid->l1dim-1) - icur_pquad->z);
             (*y >= (endy = (icur_cquad->y + icur_cquad->length ))) &&
             (icur_cquad->next != NULL); 
             icur_cquad = icur_cquad->next)
            ;
         
         /* is nquad in this cquad */
         if((*y >= icur_cquad->y) && (*y < endy))
           {
            for(icur_nquad = icur_cquad->loc + (*y - icur_cquad->y);
                (*x >= (endx = (icur_nquad->x + icur_nquad->length))) &&
                (icur_nquad->next != NULL); 
                icur_nquad = icur_nquad->next)
               ;
            
            /* is node in this nquad */
            if((*x >=  icur_nquad->x) && (*x < endx))
              {
               icur_node        = icur_nquad->loc + (*x - icur_nquad->x);
               tsc_box[0][1][1] = icur_node;     /* assign middle node */
               search_TSCplane(cur_grid, icur_pquad, icur_cquad, icur_nquad, 
                            tsc_box[0], (cur_grid->l1dim - 1), *y, *x);
              }
            else
               NULL_TSCplane(tsc_box[0]);
           }
         else
            NULL_TSCplane(tsc_box[0]);
        }
      else
         NULL_TSCplane(tsc_box[0]); 
     }
#endif /* PERIODIC_Z */
   
   else
      NULL_TSCplane(tsc_box[0]);
   
   
   
   /* 3. search plane above */
   
   if(*z < (cur_pquad->z + cur_pquad->length - 1))
     {
      /* now look for correct cquad */
      
      /* increment down cquad ll */
      for(icur_cquad = (cur_pquad->loc + ((*z + 1) - cur_pquad->z));
          (*y >= (endy = (icur_cquad->y + icur_cquad->length))) &&
          (icur_cquad->next != NULL); 
          icur_cquad = icur_cquad->next)
         ;
      
      /* is nquad in this cquad */
      if((*y >=  icur_cquad->y) && (*y < endy))
        {
         /* now look for correct nquad */
         
         /* increment down nquad ll */
         for(icur_nquad = icur_cquad->loc + (*y - icur_cquad->y);
             (*x >= (endx = (icur_nquad->x + icur_nquad->length))) &&
             (icur_nquad->next != NULL); 
             icur_nquad = icur_nquad->next)
            ;
         
         /* is node in this nquad */
         if((*x >=  icur_nquad->x) && (*x < endx))
           {
            icur_node        = icur_nquad->loc + (*x - icur_nquad->x);
            tsc_box[2][1][1] = icur_node;     /* assign middle node */
            search_TSCplane(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_box[2],
                         modulo((*z + 1), cur_grid->l1dim), *y, *x);
           }
         else
            NULL_TSCplane(tsc_box[2]);
        }
      else
         NULL_TSCplane(tsc_box[2]);
     }
   
#ifdef PERIODIC_Z
   /* deal with periodic case */
   else if(*z == (cur_grid->l1dim - 1))     /* we are at max boundary */
     {
      /* search for correct pquad */
      
      /* move to first pquad  */
      icur_pquad = cur_grid->pquad;
      if(icur_pquad->z == 0)  /* if 1st cquad is on boundary */
        {
         for(icur_cquad = icur_pquad->loc;
             (*y >= (endy = (icur_cquad->y + icur_cquad->length ))) &&
             (icur_cquad->next != NULL); 
             icur_cquad = icur_cquad->next)
            ;
         
         /* is node in this cquad */
         if((*y >= icur_cquad->y) && (*y < endy))
           {
            for(icur_nquad = icur_cquad->loc + (*y - icur_cquad->y);
                (*x >= (endx = (icur_nquad->x + icur_nquad->length))) &&
                (icur_nquad->next != NULL); 
                icur_nquad = icur_nquad->next)
               ;
            
            /* is node in this nquad */
            if((*x >=  icur_nquad->x) && (*x < endx))
              {
               icur_node        = icur_nquad->loc + (*x - icur_nquad->x);
               tsc_box[2][1][1] = icur_node;     /* assign middle node */
               search_TSCplane(cur_grid, icur_pquad, icur_cquad, icur_nquad, 
                            tsc_box[2], 0, *y, *x);
              }
            else
               NULL_TSCplane(tsc_box[2]);
           }
         else
            NULL_TSCplane(tsc_box[2]);
        }
      else
         NULL_TSCplane(tsc_box[2]);
     }
#endif /* PERIODIC_Z */
   
   else
      NULL_TSCplane(tsc_box[2]);
}


/*---------------------------------------------------------------------------------
 * get_MHDcolumn:
 * --------------
 * test for the presence of all 5 nodes in a given column (z, y)
 * (we do not need to check the nodes x+1 and x-1 as they must be there
 *  by construction, i.e. it is impossible to have x-2,x,x+2 without x-1,x+1)
 *
 * return pointer to all relevant 5 nodes in array MHDcolumn[5]
 *
 * NOTE: if the middle node x is missing we do >>not<< search for the others!!!
 *       => this is not a "general purpose routine" 
 *          filling MHDcolumn[5] with the correct pointers
 *---------------------------------------------------------------------------------*/
void get_MHDcolumn(gridls *cur_grid, cqptr cur_cquad, long y, long x,
                  nptr MHDcolumn[5])
{
   nqptr cur_nquad, s_nquad;
      
   /* jump to correct nquad */
   for(cur_nquad = cur_cquad->loc + (y - cur_cquad->y);
       (x >= cur_nquad->x + cur_nquad->length) && (cur_nquad->next != NULL);
       cur_nquad = cur_nquad->next)
      ;
   
   /*--------------------------------------
    * does x actually lie within an nquad?
    *--------------------------------------*/
   if((x >= cur_nquad->x) && (x < cur_nquad->x + cur_nquad->length))
     {
      /*-------------------
       * set middle node x
       *------------------*/
      MHDcolumn[2] = cur_nquad->loc + (x - cur_nquad->x);

      /*-------------------------
       * check for x+1/x+2 nodes
       *-------------------------*/
      if(x == cur_grid->l1dim - 2)
        {
         /* use the chance to set x+1 node ... it >>must<< be there according to refinement strategy! */
         MHDcolumn[3] = cur_nquad->loc + ((x+1) - cur_nquad->x);
         
#ifdef PERIODIC_X
         /* jump back to the beginning of the nquad linked-list */
         s_nquad = cur_cquad->loc + (y - cur_cquad->y);
         
         /* is x+2 node present? */
         if(s_nquad->x == 0)
            MHDcolumn[4] = s_nquad->loc;
         else
            MHDcolumn[4] = NULL;
#else
         MHDcolumn[4] = NULL;
#endif /* PERIODIC_X */
        }
      else if(x == cur_grid->l1dim - 1)
        {
#ifdef PERIODIC_X
         /* jump back to the beginning of the nquad linked-list */
         s_nquad = cur_cquad->loc + (y - cur_cquad->y);
         
         /* is x+2 node present? */
         if(s_nquad->x == 0)
           {
            /* we always refine in packs of 2 and hence both nodes should be there... */
            MHDcolumn[3] = s_nquad->loc;
            MHDcolumn[4] = s_nquad->loc + 1;         
           }
         else
           {
            MHDcolumn[3] = NULL;
            MHDcolumn[4] = NULL;
           }
#else
         MHDcolumn[3] = NULL;
         MHDcolumn[4] = NULL;
#endif /* PERIODIC_X */
        }
      else if(x == cur_nquad->x + cur_nquad->length - 1)
        {
         MHDcolumn[3] = NULL;
         MHDcolumn[4] = NULL;
        }
      else if(x == cur_nquad->x + cur_nquad->length - 2)
        {
         MHDcolumn[3] = cur_nquad->loc + ((x+1)-cur_nquad->x);
         MHDcolumn[4] = NULL;
        }
      else
        {
         MHDcolumn[3] = cur_nquad->loc + ((x+1)-cur_nquad->x);
         MHDcolumn[4] = cur_nquad->loc + ((x+2)-cur_nquad->x);
        }
      
      /*-------------------------
       * check for x-1/x-2 nodes
       *-------------------------*/
      if(x == 1)
        {
         /* use the chance to set x-1 node ... it >>must<< be there according to refinement strategy! */
         MHDcolumn[1] = cur_nquad->loc;
         
#ifdef PERIODIC_X
         for(s_nquad = cur_nquad; s_nquad->next != NULL; s_nquad = s_nquad->next)
            ;
         
         if(s_nquad->x + s_nquad->length == cur_grid->l1dim)
            MHDcolumn[0] = s_nquad->loc + s_nquad->length-1;
         else
            MHDcolumn[0] = NULL;
#else
         MHDcolumn[0] = NULL;
#endif /* PERIODIC_X */
        }
      else if(x == 0)
        {
#ifdef PERIODIC_X
         for(s_nquad = cur_nquad; s_nquad->next != NULL; s_nquad = s_nquad->next)
            ;
         
         if(s_nquad->x + s_nquad->length == cur_grid->l1dim)
           {
            MHDcolumn[1] = s_nquad->loc + s_nquad->length-1;
            MHDcolumn[0] = s_nquad->loc + s_nquad->length-2;
           }
         else
           {
            MHDcolumn[1] = NULL;
            MHDcolumn[0] = NULL;
           }
#else
         MHDcolumn[1] = NULL;
         MHDcolumn[0] = NULL;
#endif /* PERIODIC_X */
        }
      else if(x == cur_nquad->x)
        {
         MHDcolumn[0] = NULL;
         MHDcolumn[1] = NULL;
        }
      else if(x == cur_nquad->x+1)
        {
         MHDcolumn[0] = NULL;
         MHDcolumn[1] = cur_nquad->loc;
        }
      else
        {
         MHDcolumn[0] = cur_nquad->loc + ((x-2)-cur_nquad->x);
         MHDcolumn[1] = cur_nquad->loc + ((x-1)-cur_nquad->x);
        }   
      
     }
   /*-----------------------------------------
    * x lies outside an nquad!
    * (but maybe x-2,x-1,x+1,x+2 lie inside?)
    *-----------------------------------------*/
   else
     {
      MHDcolumn[2] = NULL;
      
      /* in a more general situation we would need to search for all the x-2,x-1,x+1,x+2 nodes,
         but given the nature of MHD equation we >>must<< have access to >>all<<
         x-1, x, x+1 nodes ... if only one is missing we are lost... */
      
      /* it is therefore safe to >>not<< search for the other nodes as we are already missing the x! */
      
      MHDcolumn[0] = NULL;
      MHDcolumn[1] = NULL;
      MHDcolumn[3] = NULL;
      MHDcolumn[4] = NULL;
     }   
}

/*---------------------------------------------------------------------------------
 * NULL_MHDcolumn:
 * --------------
 *
 * simply set all pointers to NULL
 *---------------------------------------------------------------------------------*/
void NULL_MHDcolumn(nptr MHDcolumn[5])
{
   int i;
   
   for(i=0; i<5; i++)
      MHDcolumn[i] = NULL;
}

/*---------------------------------------------------------------------------------
 * NULL_MHDplane:
 * --------------
 *
 * simply set all pointers to NULL
 *---------------------------------------------------------------------------------*/
void NULL_MHDplane(nptr MHDplane[5][5])
{
   int i, j;
   
   for(j=0; j<5; j++)
      for(i=0; i<5; i++)
         MHDplane[j][i] = NULL;
}

/*---------------------------------------------------------------------------------
 * get_MHDplane:
 * -------------
 * return pointer to all relevant 5x5 nodes in array MHDplane[5][5]
 * 
 * NOTE: y and x are always the coordinates of the middle node!
 *---------------------------------------------------------------------------------*/
void get_MHDplane(gridls *cur_grid, pqptr cur_pquad, long z, long y, long x,
                  nptr MHDplane[5][5])
{
   cqptr cur_cquad, s_cquad;
   int   i;

   /* jump to middle cquad */
   for(cur_cquad = cur_pquad->loc + (z - cur_pquad->z);
       (y >= cur_cquad->y + cur_cquad->length) && (cur_cquad->next != NULL);
       cur_cquad = cur_cquad->next)
      ;
   
   /*-------------------------------------
    * does y actually lie within a cquad?
    *-------------------------------------*/
   if((y >= cur_cquad->y) && (y < cur_cquad->y + cur_cquad->length))
     {
      /*---------------------
       * check for y column
       *---------------------*/
      get_MHDcolumn(cur_grid, cur_cquad, y, x, MHDplane[2]);
      
      
      /*---------------------------
       * check for y+1/y+2 columns
       *---------------------------*/
      if (y == cur_grid->l1dim-2)
        {
         /* use the chance to set y+1 column ... it >>must<< be there according to refinement strategy! */
         get_MHDcolumn(cur_grid, cur_cquad, y+1, x, MHDplane[3]);
         
#ifdef PERIODIC_Y
         /* jump back to the beginning of the cquad linked-list */
         s_cquad = cur_pquad->loc + (z - cur_pquad->z);
         if(s_cquad->y == 0)
            get_MHDcolumn(cur_grid, s_cquad, 0, x, MHDplane[4]);
         else
            NULL_MHDcolumn(MHDplane[4]);
#else
         NULL_MHDcolumn(MHDplane[4]);
#endif /* PERIODIC_Y */
        }
      else if (y == cur_grid->l1dim-1)
        {
#ifdef PERIODIC_Y
         /* jump back to the beginning of the cquad linked-list */
         s_cquad = cur_pquad->loc + (z - cur_pquad->z);
         
         if(s_cquad->y == 0)
           {
            get_MHDcolumn(cur_grid, s_cquad, 0, x, MHDplane[3]);
            get_MHDcolumn(cur_grid, s_cquad, 1, x, MHDplane[4]);
           }
         else
           {
            NULL_MHDcolumn(MHDplane[3]);
            NULL_MHDcolumn(MHDplane[4]);
           }
#else
         NULL_MHDcolumn(MHDplane[3]);
         NULL_MHDcolumn(MHDplane[4]);
#endif /* PERIODIC_Y */
        }
      else if (y == cur_cquad->y + cur_cquad->length-1)
        {
         NULL_MHDcolumn(MHDplane[3]);
         NULL_MHDcolumn(MHDplane[4]);
        }
      else if (y == cur_cquad->y + cur_cquad->length-2)
        {
         get_MHDcolumn(cur_grid, cur_cquad, y+1, x, MHDplane[3]);
         NULL_MHDcolumn(MHDplane[4]);
        }
      else
        {
         get_MHDcolumn(cur_grid, cur_cquad, y+1, x, MHDplane[3]);
         get_MHDcolumn(cur_grid, cur_cquad, y+2, x, MHDplane[4]);
        }
      
      
      /*---------------------------
       * check for y-1/y-2 columns
       *---------------------------*/
      if (y == 1)
        {
         /* use the chance to set y-1 column ... it >>must<< be there according to refinement strategy! */
         get_MHDcolumn(cur_grid, cur_cquad, 0, x, MHDplane[1]);

#ifdef PERIODIC_Y
         for(s_cquad = cur_cquad; s_cquad->next != NULL; s_cquad = s_cquad->next)
            ;

         if(s_cquad->y + s_cquad->length == cur_grid->l1dim)
            get_MHDcolumn(cur_grid, s_cquad, cur_grid->l1dim-1, x, MHDplane[0]);
#else
         NULL_MHDcolumn(MHDplane[0]);
#endif /* PERIODIC_Y */
        }
      else if (y == 0)
        {
#ifdef PERIODIC_Y
         for(s_cquad = cur_cquad; s_cquad->next != NULL; s_cquad = s_cquad->next)
            ;
         
         if(s_cquad->y + s_cquad->length == cur_grid->l1dim)
           {
            get_MHDcolumn(cur_grid, s_cquad, cur_grid->l1dim-1, x, MHDplane[1]);
            get_MHDcolumn(cur_grid, s_cquad, cur_grid->l1dim-2, x, MHDplane[0]);
           }
#else
         NULL_MHDcolumn(MHDplane[1]);
         NULL_MHDcolumn(MHDplane[0]);
#endif /* PERIODIC_Y */
        }
      else if(y == cur_cquad->y)
        {
         NULL_MHDcolumn(MHDplane[1]);
         NULL_MHDcolumn(MHDplane[0]);
        }
      else if(y == cur_cquad->y+1)
        {
         get_MHDcolumn(cur_grid, cur_cquad, y-1, x, MHDplane[1]);
         NULL_MHDcolumn(MHDplane[0]);
        }
      else
        {
         get_MHDcolumn(cur_grid, cur_cquad, y-1, x, MHDplane[1]);
         get_MHDcolumn(cur_grid, cur_cquad, y-2, x, MHDplane[0]);
        }
     }
   else
     {
      for(i=0; i<5; i++)
         NULL_MHDcolumn(MHDplane[i]);
     }
}

/*=====================================================================================
 * get_MHDnodes:
 *
 * get pointers to all those 125=5x5x5 neighbouring nodes required for the MHD solver
 *
 * NOTE: this is >>not<< a general purpose routine that fills the array MHDnodes[5][5][5]
 *       with either the pointer or NULL!
 *       -> if one of the relevant nodes is missing you may end up with NULL's in a plane
 *          or column even though if some of the (non-relevant!) nodes are in fact there!
 *
 *=====================================================================================*/
void get_MHDnodes(gridls *cur_grid, pqptr cur_pquad, long z, long y, long x, nptr MHDnodes[5][5][5])
{
   pqptr s_pquad;
   
   /*-------------------------------------
    *    check for z plane
    *-------------------------------------*/
   get_MHDplane(cur_grid, cur_pquad, z, y, x, MHDnodes[2]);
   

   /*-------------------------------------
    *    check for z+1/z+2 planes
    *-------------------------------------*/
   if(z == cur_grid->l1dim-2)
     {
      get_MHDplane(cur_grid, cur_pquad, z+1, y, x, MHDnodes[3]);
#ifdef PERIODIC_Z
      s_pquad = cur_grid->pquad;
      if(s_pquad->z == 0)
         get_MHDplane(cur_grid, s_pquad, 0, y, x, MHDnodes[4]);
      else
         NULL_MHDplane(MHDnodes[4]);
#else
      NULL_MHDplane(MHDnodes[4]);
#endif /* PERIODIC_Z */
     }
   else if(z == cur_grid->l1dim-1)
     {
#ifdef PERIODIC_Z
      s_pquad = cur_grid->pquad;
      if(s_pquad->z == 0)
        {
         get_MHDplane(cur_grid, s_pquad, 0, y, x, MHDnodes[3]);
         get_MHDplane(cur_grid, s_pquad, 1, y, x, MHDnodes[4]);
        }
      else
        {
         NULL_MHDplane(MHDnodes[3]);
         NULL_MHDplane(MHDnodes[4]);
        }
#else
      NULL_MHDplane(MHDnodes[3]);
      NULL_MHDplane(MHDnodes[4]);
#endif /* PERIODIC_Z */
     }
   else if(z == cur_pquad->z + cur_pquad->length-1)
     {
      NULL_MHDplane(MHDnodes[3]);
      NULL_MHDplane(MHDnodes[4]);
     }
   else if(z == cur_pquad->z + cur_pquad->length-2)
     {
      get_MHDplane(cur_grid, cur_pquad, z+1, y, x, MHDnodes[3]);
      NULL_MHDplane(MHDnodes[4]);
     }
   else
     {
      get_MHDplane(cur_grid, cur_pquad, z+1, y, x, MHDnodes[3]);
      get_MHDplane(cur_grid, cur_pquad, z+2, y, x, MHDnodes[4]);
     }
   
   
   /*-------------------------------------
    *    check for z-1/z-2 planes
    *-------------------------------------*/
   if(z == 1)
     {
      get_MHDplane(cur_grid, cur_pquad, z-1, y, x, MHDnodes[1]);
#ifdef PERIODIC_Z
      for(s_pquad = cur_pquad; s_pquad->next != NULL; s_pquad = s_pquad->next)
         ;
      
      if(s_pquad->z + s_pquad->length == cur_grid->l1dim)
         get_MHDplane(cur_grid, s_pquad, cur_grid->l1dim-1, y, x, MHDnodes[0]);
#else
      NULL_MHDplane(MHDnodes[0]);
#endif /* PERIODIC_Z */
     }
   else if(z == 0)
     {
#ifdef PERIODIC_Z
      for(s_pquad = cur_pquad; s_pquad->next != NULL; s_pquad = s_pquad->next)
         ;
      
      if(s_pquad->z + s_pquad->length == cur_grid->l1dim)
        {
         get_MHDplane(cur_grid, s_pquad, cur_grid->l1dim-1, y, x, MHDnodes[1]);
         get_MHDplane(cur_grid, s_pquad, cur_grid->l1dim-2, y, x, MHDnodes[0]);
        }
      else
        {
         NULL_MHDplane(MHDnodes[1]);
         NULL_MHDplane(MHDnodes[0]);
        }
#else
      NULL_MHDplane(MHDnodes[1]);
      NULL_MHDplane(MHDnodes[0]);
#endif /* PERIODIC_Z */
     }
   else if(z == cur_pquad->z)
     {
      NULL_MHDplane(MHDnodes[1]);
      NULL_MHDplane(MHDnodes[0]);
     }
   else if(z == cur_pquad->z+1)
     {
      get_MHDplane(cur_grid, cur_pquad, z-1, y, x, MHDnodes[1]);
      NULL_MHDplane(MHDnodes[0]);
     }
   else
     {
      get_MHDplane(cur_grid, cur_pquad, z-1, y, x, MHDnodes[1]);
      get_MHDplane(cur_grid, cur_pquad, z-2, y, x, MHDnodes[0]);
     }
   
}  


/*==============================================================================================
 * check if the node with coordinates (x_node,y_node,z_node) corresponds to (nptr node)!
 *==============================================================================================*/
void check_node(gridls *cur_grid, nptr node, long z_node, long y_node, long x_node)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   long    x, y, z;
   int     ifound = FALSE;
   
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
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
                      cur_node++, x++)
                    {
                     if(x == x_node && y == y_node && z == z_node)
                       {
                        ifound = TRUE;
                        
                        if(cur_node != node)
                          {
                           fprintf(stderr,"check_node:     l1dim=%ld     x=%ld y=%ld z=%ld   x_node=%ld y_node=%ld z_node=%ld    cur_node=%ld   check_node=%ld\n",
                                   cur_grid->l1dim,x,y,z,x_node,y_node,z_node,(long)cur_node,(long)check_node);
                          }
                        goto exit_routine;
                       }
                    }
                 }
              }
           }
        }
     }
   
   if(node != NULL && ifound == FALSE)
     {
      fprintf(stderr,"check_node:     NOT FOUND   (x=%ld y=%ld z=%ld)\n",x_node,y_node,z_node);
     }
   
exit_routine:
      ;
}


#ifdef WITH_MPI

extern nptr
get_node_from_key(gridls *curgrid,
                  sfc_key_t key,
                  sfc_curve_t ctype,
                  uint32_t bits)
{
	uint32_t pos[3];
	pqptr pq;
	cqptr cq;
	nqptr nq;

	sfc_curve_calcPos(ctype, key, bits, pos);

	return get_node_from_pos(curgrid, pos);
}

extern nptr
get_node_from_pos(gridls *curgrid,
                  uint32_t *pos)
{
	pqptr pq;
	cqptr cq;
	nqptr nq;

	pq = curgrid->pquad;
	if (pq->z > pos[2])
		return NULL;
	while (pq->z+pq->length-1 < pos[2]) {
		pq = pq->next;
		if ( (pq == NULL) || (pq->z > pos[2]) )
			return NULL;
	}
	cq = pq->loc + (pos[2] - pq->z);

	if (cq->y > pos[1])
		return NULL;
	while (cq->y+cq->length-1 < pos[1]) {
		cq = cq->next;
		if ( (cq == NULL) || (cq->y > pos[1]) )
			return NULL;
	}
	nq = cq->loc + (pos[1] - cq->y);

	if (nq->x > pos[0])
		return NULL;
	while (nq->x+nq->length-1 < pos[0]) {
		nq = nq->next;
		if ( (nq == NULL) || (nq->x > pos[0]) )
			return NULL;
	}

	return nq->loc + (pos[0] - nq->x);
}

#endif /* WITH_MPI */
