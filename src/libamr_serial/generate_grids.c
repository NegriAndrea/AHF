#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

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

/*===============================================================================
 * gen_domgrids: generate the domain grids
 *===============================================================================*/
gridls *gen_domgrids(int *no_grids)
{
  gridls *grid_list;       /* pointer to array of grid entries          */
  gridls *cur_grid;        /* pointer to current grid                   */
  double  mass2dens;        /* conversion factor to get mass density     */
  double  mass2partdens;    /* conversion factor to get particle density */
  double  frac, dl1dim, dno_part;
  int     exp;
  int     l1dim;

#ifdef EXTRAE_API_USAGE
  Extrae_user_function(1);
#endif

  /* calculate the total number of domain grids */
  frac        = (double) simu.NGRID_DOM / (double) simu.NGRID_MIN;
  exp         = (float) (log(frac)/log(2.));
  *no_grids   = exp + 1;
  
#ifdef VERBOSE
  fprintf(stderr,"\n\ngen_domgrids: total number of domain grids = %d (%g, %d)\n\n",*no_grids,frac,exp);
#endif
  
  /*======================================
   * there are '*no_grids' to be created:
   *        a) simu.NGRID_MIN^3 grid
   *            ...
   *        x) simu.NGRID_DOM^3 grid
   *=======================================*/
  
  /* allocate memory for all grids */
  grid_list = (gridls *) calloc(*no_grids, sizeof(gridls));
  
  /*
   * 'grid_list' is the pointer to the first gridls-structure
   * out of '*no_grids' of such structures
   */
  
  /* generate the all grids and initialize the structure with the corresponding data */
  for(cur_grid = grid_list, l1dim = simu.NGRID_MIN; l1dim <= simu.NGRID_DOM; l1dim *= 2, cur_grid++)
    {     
      /* mass2dens <=> 1/rho_bar [code units] */
      dl1dim        = (double)l1dim;
      dno_part      = (double)simu.no_part;
      mass2dens     = pow3(dl1dim)/simu.no_vpart;
      mass2partdens = pow3(dl1dim)/dno_part;
      
      /*===========================
       * fill in grid_list details
       *===========================*/
      cur_grid->timecounter    = global.super_t;
      cur_grid->l1dim          = l1dim;
      cur_grid->spacing        = (double)1.0 / (double) cur_grid->l1dim;
      cur_grid->spacing2       = pow2(cur_grid->spacing);

      cur_grid->masstodens     = mass2dens;
      cur_grid->masstopartdens = mass2partdens;
#ifdef MULTIMASS
      cur_grid->critdens       = simu.Nth_dom * mass2partdens;
#else
      cur_grid->critdens       = simu.Nth_dom * mass2dens;
#endif
      cur_grid->old_resid      = 0.0;
      cur_grid->cur_resid      = 0.0;
      cur_grid->no_sweeps      = 0;
      
      cur_grid->time.potential = 0;
      cur_grid->time.density   = 0;
      cur_grid->time.DK        = 0;
      cur_grid->time.grid      = 0;
      cur_grid->time.hydro     = 0;
      
      cur_grid->multistep             = 0;
      cur_grid->next                  = FALSE;
      
      /*========================================================
       * generate the actual pquad, cquad, and nquad structures
       *========================================================*/
      alloc_quads(cur_grid, l1dim);
    }
    
  /*===========================
   * keep track of domain grid
   *===========================*/
  global.dom_grid   = grid_list + (*no_grids - 1);
  global.domgrid_no = (*no_grids - 1);

#ifdef EXTRAE_API_USAGE
  Extrae_user_function(0);
#endif

  return(grid_list);
}


/*==============================================================================
 * gen_refgrid: 
 *--------------------
 * 1. test current grid for refinement, and generate adaptive grid (if required) 
 * 2. re-link particles 
 * 3. ajdust densities as neccessary.
 *==============================================================================*/
boolean gen_refgrid(gridls **grid_list, int *no_grids)
{
  gridls *coa_grid, *fin_grid;     /* old/new grid pointer                    */
  boolean refined;                 /* refined YES/NO flag                     */
  boolean mg;                      /* does refinement hold enough particles   */
  int     idim;
  
#ifdef REF_TEST
  fprintf(stderr,"\ngen_refgrid:         current coarse grid         = %ld\n",
          (*grid_list + (*no_grids-1))->l1dim);
#endif
  
  /* set coarse grid pointer */
  coa_grid = *grid_list + *no_grids - 1;
  fin_grid = *grid_list + *no_grids;
  
  /* timing... */
  coa_grid->time.grid -= time(NULL);
  
  if(coa_grid->next == FALSE)  /* fin_grid does not exist ? */
    {
      /* reallocate the grid list array (including space for one additional grid) */
      (*grid_list) = (gridls *) realloc(*grid_list, (*no_grids + 1) * sizeof(gridls));
      
      if((*grid_list) == NULL)
        {
          fprintf(io.logfile,"gen_refgrid: error reallocating grid list\n");
          fflush(io.logfile);
          fclose(io.logfile);
          exit(1);
        }
      
      /* grid_list block might have moved in memory */
      global.dom_grid = *grid_list + global.domgrid_no;
      coa_grid        = *grid_list + *no_grids - 1;
      fin_grid        = *grid_list + *no_grids;
      
      /* fill in new grid's entries */
      fin_grid->masstodens     = coa_grid->masstodens     * CRITMULTI;
      fin_grid->masstopartdens = coa_grid->masstopartdens * CRITMULTI;
#ifdef MULTIMASS
      fin_grid->critdens       = simu.Nth_ref*fin_grid->masstopartdens;
#else
      fin_grid->critdens       = simu.Nth_ref*fin_grid->masstodens;
#endif
      fin_grid->timecounter    = coa_grid->timecounter;
      fin_grid->timestep       = coa_grid->timestep/2;
      fin_grid->l1dim          = coa_grid->l1dim * 2;
      fin_grid->spacing        = (double)1.0 / (double)fin_grid->l1dim;
      fin_grid->spacing2       = fin_grid->spacing * fin_grid->spacing;
      
      fin_grid->no_pquad       = 0;
      fin_grid->pquad_array    = NULL;

      /* we measure the time throughout the two DKD cycles of this fin_grid */
      fin_grid->time.potential = 0;
      fin_grid->time.density   = 0;
      fin_grid->time.DK        = 0;
      fin_grid->time.grid      = 0;
      fin_grid->time.hydro     = 0;
       
      fin_grid->no_sweeps      = 0;
      fin_grid->cur_resid      = 0.;
      
      /* remember that fin_grid exists from now on */
      fin_grid->next = FALSE;
      coa_grid->next = TRUE;
    }
  
  
  /* (re-)set values that are changing within the two DKD cycles the grid lives through */
  fin_grid->old_resid             = 0.0;
  fin_grid->cur_resid             = 0.0;
  
  fin_grid->multistep             = 0;
  
  
  /*------------------------------
   * now call refinement rountine 
   *------------------------------*/
  refined = refine_grid(fin_grid, coa_grid);
  
  if(refined)
    {
#ifdef REF_TEST
      fprintf(stderr,"                     placed new refinement       = %ld\n", fin_grid->l1dim);
#endif
      
      /* update number of grids */
      *no_grids += 1;
      
      /* now relink particles to new grid */
      mg = relink(coa_grid, fin_grid);        /* performs: unassign_part */
      
      /* now assign mass to new grid: THIS IS NEEDED AS IT RETURNS "mg" ! */
      zero_dens        (fin_grid);
      mg = assign_npart(fin_grid);
      
#ifdef REF_TEST
      fprintf(stderr,"             grid %6ld:  nnodes=%ld npart=%ld => %g (<! %g)\n",
              fin_grid->l1dim,fin_grid->size.no_nodes,fin_grid->size.no_part,
              (double)fin_grid->size.no_nodes/(double)fin_grid->size.no_part, CRITMULTI);
#endif
      
      /* only use grid if it holds enough particles */
      if(mg == FALSE)
        {
          refined = FALSE;
          relink_back(coa_grid, fin_grid);
          fin_grid->size.no_nodes         = 0;
          fin_grid->size.no_part          = 0;
          
          free_grid(fin_grid, no_grids);
          
#ifdef REF_TEST
          fprintf(stderr,"                     destroyed new refinement    = %ld (%g)\n",
                  fin_grid->l1dim, fin_grid->critdens);
#endif
        }
      else
        {
          /* get potential onto fine grid (important for boundaries !) */
          go_down(coa_grid);          
        }
    }
  
  else
    {
#ifdef REF_TEST
      fprintf(stderr,"                     NO new refinement grid !\n\n");
#endif
      /*
       * no need to call free_grid() 
       * -> temporary pquads, cquads, nquads already destroyed within ref_grid() 
       *
       * no need to realloc grid_list
       * -> grid will be kept in memory from now on...
       */
    }
  
  /* grid_list block might have changed position in memory */
  global.dom_grid = *grid_list + global.domgrid_no;
  
  /* ...timing */
  coa_grid->time.grid += time(NULL);
  
  return(refined);
}

/*================================================================================
 * recursively generate the AMR hierarchy
 *================================================================================*/
boolean gen_AMRhierarchy(gridls **grid_list, int *no_grids)
{
  gridls *cur_grid;                   /* current grid pointer                  */
  gridls *fin_grid;                   /* fine grid pointer                     */
  gridls *for_grid;                   /* sometimes needed for debugging...     */
  int     fingrid_no, curgrid_no;     /* index of finest and current grid      */
  int     grid_no;
  
  boolean runAHF;  /* only run AHF when the finest grid level has been reached */
  boolean refined; /* refinement flag                                          */
  
#ifdef EXTRAE_API_USAGE
  Extrae_user_function(1);
#endif

  /* up to now we plan to run AHF */
  runAHF = TRUE;
  
  
  /*=========================================================================
   * try to refine finest grid:
   *
   * 1. get (number!) density field on currently finest grid right
   * 2. try to refine grid using 'number of particles per node' criterion
   *=========================================================================*/
  
  /* currently the finest grid... */
  fin_grid = *grid_list + *no_grids - 1;
  
  if(!(((global.fst_cycle == FALSE) && (fin_grid->l1dim == global.fin_l1dim)) || 
       ((fin_grid->l1dim) == simu.NGRID_MAX)))
   {
    timing.genrefgrids  -= time(NULL);
    
    /* get number density: "npart/node" */
    fin_grid->time.grid -= time(NULL);
    zero_dens   (fin_grid);
    assign_npart(fin_grid);
    fin_grid->time.grid += time(NULL);
    
    /* try to refine grid under "npart/node" criterion */
    refined = gen_refgrid(grid_list, no_grids);
    
    timing.genrefgrids  += time(NULL);
   }
  else
   {
    refined = FALSE;
   } 
  
  /*--------------------------------------------------
   * when refined we make a recursive call to step()
   *--------------------------------------------------*/
  if(refined == TRUE) 
   {
    /* recursive call to gen_AMRhierarchy() */
    runAHF = gen_AMRhierarchy(grid_list, no_grids);
    
    curgrid_no      = *no_grids - 1;
    fingrid_no      = *no_grids - 1;
    global.dom_grid = *grid_list + global.domgrid_no;
    cur_grid        = *grid_list + curgrid_no;
   }
  
  /*---------------------------------------------
   * otherwise stay with the actual finest grid
   *---------------------------------------------*/
  else /* refined == FALSE */
   {
    /* no new refinement => stay with current finest grid */
    curgrid_no      = *no_grids - 1;
    fingrid_no      = *no_grids - 1;
    global.dom_grid = *grid_list + global.domgrid_no;
    cur_grid        = *grid_list + curgrid_no;
    
    /* store force resolution in header info */
    io.header.cur_reflevel = (float)(curgrid_no-global.domgrid_no);
    io.header.cur_frcres   = (float)(simu.boxsize/(double)cur_grid->l1dim);
    io.header.cur_frcres  *= 3000.;
    
    
    energy.time -= time(NULL);
    energy.time += time(NULL);
    
    if(global.fst_cycle == TRUE)
     {
      global.fin_l1dim = cur_grid->l1dim;
      global.fst_cycle = FALSE;
     }
   } /* if(refined) */
  
  
  
  /*=========================================================================
   *
   *         WE NOW HAVE ACCESS TO THE FULL BLOWN GRID HIERACHY
   *
   *=========================================================================*/
  
  
  /*=========================================================================
   *          here are those parts dealing with AHFpotcentre,
   *              the actual AHF is happening in main.c
   *=========================================================================*/
  ahf.time -= time(NULL);
  /* we only run AHF when having reached the end of the grid hierarchy */
  if(runAHF)
   {
    
    /* get the density absolutely right on all amr levels */
    /*====================================================*/
    timing.densrecovery -= time(NULL);
#ifdef AHFdensrecovery
    /* re-assign density right on all grids */
    for(for_grid = cur_grid; for_grid >= global.dom_grid; for_grid--)
      restore_dens(for_grid);      
    if(curgrid_no != fingrid_no) 
      refill_dens(cur_grid);
    if(cur_grid != global.dom_grid)
      stack_dens(cur_grid, global.dom_grid+1);
#endif
    timing.densrecovery += time(NULL);
    
    
    /* obtain gravitational potential */
    /*================================*/
    timing.potcentre -= time(NULL);
#ifdef AHFpotcentre    
    /* solve for potential on domain grids */
    fprintf(stderr,"AHFpotcentres:   solving on domain grid %d ... ",global.dom_grid->l1dim);
    solve_dom_gravity(*grid_list);
    fprintf(stderr,"done\n");
    
    /* solve for potential on refinement grids */
    fprintf(stderr,"                 solving on refinement grids ... ");
    if(curgrid_no > global.domgrid_no)
      for(grid_no = global.domgrid_no+1; grid_no <= curgrid_no; grid_no++)
       {
        fprintf(stderr," %d  ",(global.dom_grid+grid_no)->l1dim);
        solve_ref_gravity(*grid_list, grid_no);
       }
    fprintf(stderr,"done\n");    
#endif /* AHFpotcentre */
    timing.potcentre += time(NULL);
    
   }
  ahf.time += time(NULL);

#ifdef EXTRAE_API_USAGE
  Extrae_user_function(0);
#endif

  /* return FALSE indicating that AHF has been run */
  return(FALSE);
} 

