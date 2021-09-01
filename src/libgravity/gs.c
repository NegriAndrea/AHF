#include <stddef.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

static double FourPiGa;

/*-------------------------------------------------------------------------------
 * gs_node: do gauss-siedel differencing on node 
 *-------------------------------------------------------------------------------*/
void gs_node(gridls *cur_grid, pqptr cur_pquad, cqptr cur_cquad, nqptr cur_nquad, 
             nptr cur_node, long z, long y, long x)
{
#ifndef AHFlean
  double cur_source;
  nptr   tsc_nodes[3][3][3];
  
  /* convert density to source term */
  cur_source = cur_node->dens * FourPiGa;
  
  tsc_nodes[1][1][1] = cur_node;
  get_TSCnodes(cur_grid, cur_pquad, cur_cquad, cur_nquad, tsc_nodes, &z, &y, &x);
  
  if(test_tsc(tsc_nodes))
    {
#ifdef NO_SOR
      cur_node->pot = 0.1666666666667 *
      ((double)tsc_nodes[2][1][1]->pot + (double)tsc_nodes[0][1][1]->pot +
       (double)tsc_nodes[1][2][1]->pot + (double)tsc_nodes[1][0][1]->pot +
       (double)tsc_nodes[1][1][2]->pot + (double)tsc_nodes[1][1][0]->pot -
       (cur_source * cur_grid->spacing2));
#else
      cur_node->pot = (W_SOR * (0.1666666666667 *
                                ((double)tsc_nodes[2][1][1]->pot + (double)tsc_nodes[0][1][1]->pot +
                                 (double)tsc_nodes[1][2][1]->pot + (double)tsc_nodes[1][0][1]->pot +
                                 (double)tsc_nodes[1][1][2]->pot + (double)tsc_nodes[1][1][0]->pot -
                                 (cur_source * cur_grid->spacing2)))) +
      (1.0-W_SOR) * (double)cur_node->pot;
#endif
    }   
#else
   fprintf(stderr, "%s should not be called when AHFlean is defined.\n",
           __func__);
   exit(1);
#endif /* AHFlean*/
}

/*------------------------------------------------------------------------------- 
 * nquad control routine 
 *-------------------------------------------------------------------------------*/
void gs_nquad(gridls *cur_grid, pqptr cur_pquad, cqptr cur_cquad,
              nqptr fst_nquad, boolean sweep_type, long z, long y)
{
  nqptr cur_nquad;      /* current nquad                       */
  nqptr old_nquad;      /* polonger to old nquad               */
  nptr  cur_node;       /* current nquad                       */
  boolean o_e;          /* odd or even? - result of nquad test */
  long x=0;             /* x-coord                             */
  
  /* deal with periodic case */
#ifdef PERIODIC_X
  if(fst_nquad->x == 0)
    {
      if(sweep_type == TRUE)
        {
          /* test if last nquad exists */
          for(cur_nquad = fst_nquad; cur_nquad->next != NULL;
              cur_nquad = cur_nquad->next)
            ;
          
          if((cur_grid->l1dim == (cur_nquad->x + cur_nquad->length)) &&
             (fst_nquad->length > 1))
            {
              gs_node(cur_grid, cur_pquad, cur_cquad, fst_nquad, fst_nquad->loc, 
                      z, y, 0);
            }
        }
    }
#endif /* PERIODIC_X */
  
  /* loop over nquads on this column */
  for(cur_nquad = fst_nquad, old_nquad = NULL, cur_node = NULL; cur_nquad != NULL; 
      old_nquad = cur_nquad, cur_nquad = cur_nquad->next)
    {
      /* test first node in nquad for odd/even */
      o_e = is_even(cur_nquad->x);
      if(o_e == sweep_type)
        {
          cur_node = cur_nquad->loc + 2;
          x        = cur_nquad->x   + 2;
        }
      else
        {
          cur_node = cur_nquad->loc + 1;
          x        = cur_nquad->x   + 1;
        }
      
      /* loop over inner nodes */
      while(cur_node < (cur_nquad->loc + cur_nquad->length - 1))
        {
          gs_node(cur_grid, cur_pquad, cur_cquad, cur_nquad, cur_node, z, y, x);
          cur_node += 2;
          x        += 2;
        }
    }
  
  /* deal with periodic case */
#ifdef PERIODIC_X
  if(sweep_type == FALSE)
    {
      if(old_nquad->x + old_nquad->length == cur_grid->l1dim)
        {
          /* test if first node exists */
          if(fst_nquad->x == 0)
            gs_node(cur_grid, cur_pquad, cur_cquad, old_nquad, cur_node, z, y, x);
        }
    }
#endif /* PERIODIC_X */
}


/*------------------------------------------------------------------------------- 
 * cquad control routine 
 *-------------------------------------------------------------------------------*/
void gs_cquad(gridls *cur_grid, pqptr cur_pquad, cqptr fst_cquad,
              boolean sweep_type, long z)
{
  cqptr cur_cquad;      /* current cquad                       */
  cqptr old_cquad;      /* old cquad                           */
  nqptr cur_nquad;      /* current nquad                       */
  boolean cur_type=0;   /* current boolean value of sweep_type */
  long y=0;             /* y-coord of cur_nquad                */
  
  /* deal with periodic case */ 
#ifdef PERIODIC_Y
  if(fst_cquad->y == 0)
    {
      /* test if last nquad exists */
      for(cur_cquad = fst_cquad; cur_cquad->next != NULL; cur_cquad=cur_cquad->next)
        ;
      if((cur_grid->l1dim == (cur_cquad->y + cur_cquad->length)) &&
         (fst_cquad->length > 1))
        {
          gs_nquad(cur_grid,cur_pquad, fst_cquad, fst_cquad->loc, sweep_type, z, 0);
        }
    }
#endif /* PERIODIC_Y */
  
  /* loop over cquads on this plane */
  for(cur_cquad = fst_cquad, old_cquad = NULL; cur_cquad != NULL; 
      old_cquad = cur_cquad, cur_cquad = cur_cquad->next)
    {
      /* test first nquad in cquad for odd/even */  
      if((is_even(cur_cquad->y)) == TRUE)
        cur_type = set_opposite(sweep_type);
      else
        cur_type = sweep_type;
      
      /* loop over inner nquads */  
      for(cur_nquad = cur_cquad->loc + 1, y = cur_cquad->y + 1;
          cur_nquad < (cur_cquad->loc + cur_cquad->length - 1);
          cur_nquad++, y++)
        {
          gs_nquad(cur_grid, cur_pquad, cur_cquad, cur_nquad, cur_type, z, y);
          cur_type = set_opposite(cur_type);
        }
    }
  
  /* deal with periodic case */
#ifdef PERIODIC_Y
  if(old_cquad->y + old_cquad->length == cur_grid->l1dim)
    {
      /* test if first nquad exists */
      if(fst_cquad->y == 0)
        {
          cur_nquad = old_cquad->loc + (old_cquad->length - 1);
          gs_nquad(cur_grid, cur_pquad, old_cquad, cur_nquad, cur_type, z, y);
        }
    }
#endif /* PERIODIC_Y */
}

/*-------------------------------------------------------------------------------
 * pquad control routine
 *-------------------------------------------------------------------------------*/
void gs_pquad(gridls *cur_grid, pqptr fst_pquad, boolean sweep_type)
{
  long    ipquad;
  pqptr   cur_pquad;        /* current pquad                       */
  pqptr   old_pquad;        /* old pquad                           */
  cqptr   cur_cquad;        /* current cquad                       */
  boolean cur_type=0;       /* current boolean value of sweep_type */
  long    z=0;              /* z-coord of cur_cquad                */
    
  /* deal with z==0 case */  
#ifdef PERIODIC_Z
  if(fst_pquad->z == 0)
    {
      /* test if last cquad exists */
      for(cur_pquad = fst_pquad; cur_pquad->next != NULL; cur_pquad=cur_pquad->next)
        ;
      
      if((cur_grid->l1dim == (cur_pquad->z + cur_pquad->length)) && (fst_pquad->length > 1))
          gs_cquad(cur_grid, fst_pquad, fst_pquad->loc, sweep_type, 0);
    }
#endif /* PERIODIC_Z */
  
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(sweep_type) private(cur_type, cur_pquad, cur_cquad, z, old_pquad) shared(fst_pquad, cur_grid)
#pragma omp for schedule(static)
  /* loop over pquads */
  for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)
    {
      cur_pquad = cur_grid->pquad_array[ipquad];
#else
  /* loop over pquads */
  for(cur_pquad = fst_pquad; cur_pquad != NULL; cur_pquad = cur_pquad->next)
    {
#endif
      /* test first cquad in pquad for odd/even */
      if((is_even(cur_pquad->z)) == TRUE)
        cur_type = set_opposite(sweep_type);
      else
        cur_type = sweep_type;
      
      /* loop over inner cquads [z+1, z+(length-1)]*/
      for(cur_cquad = cur_pquad->loc + 1, z = cur_pquad->z + 1;
          cur_cquad < (cur_pquad->loc + cur_pquad->length - 1);
          cur_cquad++, z++)
        {
          gs_cquad(cur_grid, cur_pquad, cur_cquad, cur_type, z);
          cur_type = set_opposite(cur_type);
        }
      
      /* remember last pquad */
      old_pquad = cur_pquad;
    }
      
  /* deal with z==l1dim */ 
#ifdef PERIODIC_Z
      
#ifdef WITH_OPENMP
  /* we cannot trust the exit values of the OpenMP loop and hence need to determine these numbers again:
   * old_pquad, cur_cquad, cur_type, z
   */
  for(cur_pquad = fst_pquad; cur_pquad != NULL; cur_pquad = cur_pquad->next)
    {
      /* test first cquad in pquad for odd/even */
      if((is_even(cur_pquad->z)) == TRUE)
        cur_type = set_opposite(sweep_type);
      else
        cur_type = sweep_type;
      
      /* loop over inner cquads [z+1, z+(length-1)]*/
      for(cur_cquad = cur_pquad->loc + 1, z = cur_pquad->z + 1;
          cur_cquad < (cur_pquad->loc + cur_pquad->length - 1);
          cur_cquad++, z++)
        cur_type = set_opposite(cur_type);
      
      /* remember last pquad */
      old_pquad = cur_pquad;
    }
#endif
      
  if(old_pquad->z + old_pquad->length == cur_grid->l1dim)
    {
      /* test if first cquad exists */
      if(fst_pquad->z == 0)
        {
          cur_cquad = old_pquad->loc + (old_pquad->length - 1);
          gs_cquad(cur_grid, old_pquad, cur_cquad, cur_type, z);
        }
    }
#endif /* PERIODIC_Z */
}

/*=============================================================================
 * perform Gauss-Seidel sweeps in chessboard manner
 *=============================================================================*/
void gs(gridls *cur_grid)
{
  
  /* conversion factor to go from density to source term */
  FourPiGa = simu.FourPiG*calc_super_a(cur_grid->timecounter);
  
  /* Red */
  gs_pquad(cur_grid, cur_grid->pquad, TRUE);

  /* Black */
  gs_pquad(cur_grid, cur_grid->pquad, FALSE);
}



