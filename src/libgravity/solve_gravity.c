#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "gravity.h"
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

/*
 * there are two ways to stop the GS iteration procedure:
 *
 * 1. the residuals are smaller than the truncation error
 *
 * 2. the residuals are smaller than LIMIT given below
 *
 */
#define LIMIT     0.00001

#define ONE_THIRD 0.33333333333333333333

#include "gravity.h"

/*==============================================================================
 * solve on coarsest grid
 *==============================================================================*/
void solve_cg(gridls *cur_grid)
{
#ifndef AHFlean   
  flouble *dens_array;         /* density array pointer */
  pqptr    cur_pquad;          /* current pquad         */
  cqptr    cur_cquad;          /* current cquad         */
  nqptr    cur_nquad;          /* current nquad         */
  nptr     cur_node;           /* current node          */
  long     i, j, k, l1dim, FFTarray_length;
  double   FourPiGa;
  
  /* conversion factor to go from densito to source term */
  FourPiGa = simu.FourPiG*calc_super_a(cur_grid->timecounter);
  
  /* array dimension */
  l1dim           = cur_grid->l1dim;
  FFTarray_length = 2*l1dim*l1dim*l1dim;
  
  /* generate complex (!) density array for FFT */
  if((dens_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
   {
    fprintf(io.logfile,"solve_cg: could not allocate density array for FFT\n");
    fflush(io.logfile);
    fclose(io.logfile);
    exit(1);
   }
  
  /* fill density array for FFT ... no need for quad-ll's !! */
  cur_pquad = cur_grid->pquad;
  for(k = 0, cur_cquad = cur_pquad->loc; k < l1dim; k++, cur_cquad++)
    for(j = 0, cur_nquad = cur_cquad->loc; j < l1dim; j++, cur_nquad++)
      for(i = 0, cur_node = cur_nquad->loc; i < l1dim; i++, cur_node++)
       {
        dens_array[Re(i,j,k,l1dim)] = cur_node->dens * FourPiGa;  /* real part      */
        dens_array[Im(i,j,k,l1dim)] = 0.0;                         /* imaginary part */
       }
  
  /* solve by FFT */
  fft_potential(dens_array, l1dim);
  
  /* fill node potential values */
  for(k = 0, cur_cquad = cur_pquad->loc; k < cur_grid->l1dim; k++, cur_cquad++)
    for(j = 0, cur_nquad = cur_cquad->loc; j < cur_grid->l1dim; j++, cur_nquad++)
      for(i = 0, cur_node = cur_nquad->loc; i < cur_grid->l1dim; i++,cur_node++)
        cur_node->pot = dens_array[Re(i,j,k,l1dim)];
  
  /* destroy memory assigned to dens_array */
  free(dens_array);
  
#else
  fprintf(stderr, "%s should not be called when AHFlean is defined.\n",
          __func__);
  exit(1);
#endif /* AHFlean*/
}


/*============================================================================
 * test convergence of current grid
 *============================================================================*/
boolean converged(gridls *cur_grid)
{
  double trunc_error;
  double residual;
  
  residual     =             cur_grid->cur_resid;
  trunc_error  = ONE_THIRD * cur_grid->trunc_err;
  
  if     (residual <= LIMIT)
    return TRUE;
  else if(residual <= (trunc_error * CONVCRIT)) 
    return TRUE;
  else
    return FALSE;
}

/*=============================================================================
 * test for slow convergence
 *=============================================================================*/
boolean slow_conv(gridls *cur_grid)
{
#ifdef NO_MULTIGRID
  return FALSE;
#endif
  if(cur_grid->cur_resid > (ETA * cur_grid->old_resid))
    return TRUE;
  else
    return FALSE;
}

/*=============================================================================
 * solve for potential on domain grids
 *=============================================================================*/
void solve_dom_gravity(gridls *grid_list)
{
  int     grid_no;             /* the current grid number     */
  int     lstgrid_no;          /* no_grids minus one          */
  int     i;                   /* index for GS sweep loop     */
  gridls *cur_grid;            /* pointer to the current grid */
  boolean cg_conv;             /* has current grid converged? */
  int     gs_sweeps;
  
  gs_sweeps = DOMSWEEPS;
  
  solve_cg(global.dom_grid);
  return;
  
}


/*==============================================================================
 * solve for potential on refinements
 *==============================================================================*/
void solve_ref_gravity(gridls *grid_list, int grid_no)
{
  int     i;                   /* index for GS sweep loop     */
  gridls *cur_grid;            /* pointer to the current grid */
  boolean cg_conv;             /* has current grid converged? */
  int     gs_sweeps;
  
  /* set cur_grid */
  cur_grid = grid_list+grid_no;
  
  /* get boundary values from next coarser level */
  if(cur_grid->multistep == 3)
    go_down(grid_list+grid_no-1);
  
  cur_grid->no_sweeps  = 0;
  cur_grid->cur_resid  = 0.;
  cg_conv              = FALSE;
  gs_sweeps            = REFSWEEPS;
  
  /* do GS sweeps until converged */
  while(cg_conv == FALSE)
   {
    
    /* do GS relaxation sweep */
    for(i = 0; i < gs_sweeps; i++)
      gs(cur_grid);
    
    cur_grid->no_sweeps += REFSWEEPS;
    
    /* calculate current residual (stored for each node in force.temp[0]) */
    cur_grid->cur_resid  = residual(cur_grid);
    
    
#ifdef FIXED_SWEEPS
    cg_conv = TRUE;
#else
    /* calculate current truncation error (uses force.temp[1] and force.temp[2])  */
    cur_grid->trunc_err = trunc_err(cur_grid);
    
    /* check for convergence, i.e. compare residual against truncation error */
    cg_conv             = converged(cur_grid);
    
    
#ifdef SWEEP_TEST
    fprintf(stderr,"grid=%10ld   resid=%16.8g  trunc_err=%16.8g   -> cg_conv=%10d\n",
            cur_grid->l1dim, cur_grid->cur_resid, CONVCRIT*ONE_THIRD*cur_grid->trunc_err, cg_conv);
#endif
    
#endif
    
    if(fabs(cur_grid->old_resid - cur_grid->cur_resid) < ZERO)
      gs_sweeps += REFSWEEPS;   /* just in case we are somehow stuck */
    
   }   /* cg_conv == FALSE */
  
  /* do GS sweeps one more time */
  for(i = 0; i < REFSWEEPS; i++)
    gs(cur_grid);
  
  cur_grid->cur_resid = residual(cur_grid);
  
}


/*==============================================================================
 * control routine for multi grid solver
 *==============================================================================*/
void solve_gravity(gridls *grid_list, int curgrid_no)
{
  if(curgrid_no == global.domgrid_no)
    solve_dom_gravity(grid_list);
  else
    solve_ref_gravity(grid_list, curgrid_no);
}


