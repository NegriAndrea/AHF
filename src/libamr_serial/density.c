#include <stddef.h>
#include <stdio.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "amr_serial.h"

#ifdef EXTRAE_API_USAGE
#include <extrae_user_events.h>
#endif

/*===============================================================================
*
* This file contains all routines to assign/unassign particles to grids
* as well as to re-correct density fields before solving for the potential
*
*
* The visible functions are:
*
* zero_dens(cur_grid)     -> zero density field on given grid
* assign_dens(cur_grid)   -> assign all particles to a given grid giving density
* assign_npart(cur_grid)  -> assign all particles to a given grid
* assign_part()           -> assign a single particle to a grid
* unassign_part()         -> unassign a single particle from a grid
*
* (the assignment scheme (either NGP, CIC, or TSC) is controlled via
   *  #define NGP or #define CIC or #define TSC (cf. Makefile and 00README))
*
*
* restore_dens(cur_grid)  -> restore P1/P2 density field on coarser grids
* refill_dens(coa_grid)   -> there is a fin_grid which causes a hole in coa_grid
* reflect_dens(cur_grid)  -> reflect dens. via f2c interpolation to coarser grids
* overlap_dens(cur_grid)  -> add P1/P2 contribution from coarse to fine grid
*
*===============================================================================*/

/*===========================================================================
* assign "mass density" to grid
*===========================================================================*/
boolean assign_dens(gridls *cur_grid)
{
   pqptr cur_pquad;
   cqptr cur_cquad, icur_cquad;
   nqptr cur_nquad, icur_nquad;
   nptr  cur_node;
   nptr  tsc_nodes[3][3][3];
   long  x, y, z;
   
   partptr cur_part;           /* current particle                   */
   int     idim;               /* coord changing index               */
   int     i,j,k;              /* indices for 3D arrays              */
   double  pnarg_a;            /* pn sep temp arg                    */
   double  pnarg_b;            /* pn sep temp arg                    */
   double  tpnarg;             /* temp to calc pn args               */
   dvect   xyz_coords;         /* actual float coords of node        */
   dvect   pn_sep;             /* particle-node separation / spacing */
   dvect   weights[3];         /* weights in each dimension          */
   dvect   temp_coords;        /* float vector - calc un_mod coords  */   
   double  mass2dens, ratio, cur_shift;

   long    no_part, no_nodes;
   long    ipquad;
   
   /* reset counters for no_parts and no_nodes on cur_grid */
   no_part  = 0;
   no_nodes = 0;

   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;

   /* loop over whole grid */
#ifdef WITH_OPENMP2
  /* NOTE: this parallelisation is not 100% functional!
   *
   *       Actually it should be safe as we assign different pquad's to different threads and hence we should be able
   *       to safely modify tsc_nodes[][][]->dens without any collisions!
   *       But with the current refinement scheme (i.e. adding a boundary layer of 2 cells) it can happen that
   *       two linked quad's are actually spatially connected which is >>not<< reflected in the quad-structure!
   *
   *       We could use "#pragma omp atomic" to ensure proper "+=" for the density, but tests have shown that this
   *       dramatically slows down the runtime!
   *
   *       In that regards, feel free to switch this OpenMP parallelisation on at your own responsibility...
   */
#pragma omp parallel reduction(+:no_part, no_nodes) private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, x, y, z, tsc_nodes, cur_part, idim, i, j, k, pnarg_a, pnarg_b, tpnarg, xyz_coords, pn_sep, weights, temp_coords, mass2dens) shared(cur_shift, cur_grid, simu)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)
     {
      cur_pquad=cur_grid->pquad_array[ipquad];
#else
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
     {
#endif
      z = cur_pquad->z;
      
      for(cur_cquad=cur_pquad->loc;
          cur_cquad<cur_pquad->loc+cur_pquad->length;
          cur_cquad++, z++){
         for(icur_cquad=cur_cquad;icur_cquad!=NULL;icur_cquad=icur_cquad->next){
            y = icur_cquad->y;
            
            for(cur_nquad=icur_cquad->loc;
                cur_nquad<icur_cquad->loc+icur_cquad->length;
                cur_nquad++, y++){
               for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next){
                  x = icur_nquad->x;
                  
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                     /* increment number of nodes associated with current grid */
                      no_nodes++;
                     
                     /* calculate realspace coords of cur_node */
                     temp_coords[X] = (((double)x) / (double)cur_grid->l1dim) + cur_shift;
                     temp_coords[Y] = (((double)y) / (double)cur_grid->l1dim) + cur_shift;
                     temp_coords[Z] = (((double)z) / (double)cur_grid->l1dim) + cur_shift;
                     
                     xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
                     xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
                     xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);
                     
                     /* obtain 26 neighbours for TSC mass assignment */
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                     
                     /* loop over all particles attached to cur_node */
                     for(cur_part = cur_node->ll; cur_part != NULL; cur_part = cur_part->ll)
                       {
                        /* increment number of particles associated with current grid */
                        no_part++;
                       
                        /* calc fraction to be assigned to each tsc node */
                        for(idim = 0; idim < 3; idim++)
                          {
                           pn_sep[idim] = ((double)cur_part->pos[idim] - xyz_coords[idim]) * (double)cur_grid->l1dim;
                           
                           /* deal with periodic boundary conditions */
                           if(fabs(pn_sep[idim]) > 0.5*(double)cur_grid->l1dim)
                             {
                              tpnarg       = (double)cur_part->pos[idim] + 0.5;
                              pnarg_a      = f1mod(tpnarg+1.0, 1.0);
                              tpnarg       = xyz_coords[idim] + 0.5;
                              pnarg_b      = f1mod(tpnarg+1.0, 1.0);
                              pn_sep[idim] = (pnarg_a - pnarg_b) * (double)cur_grid->l1dim;
                             }
                           
#ifdef TSC
                           weights[1][idim] = 0.75 - pow2(pn_sep[idim]);
                           weights[0][idim] = pow2((0.5 - pn_sep[idim]))/2;
                           weights[2][idim] = pow2((0.5 + pn_sep[idim]))/2;
#endif
#ifdef CIC
                           if(pn_sep[idim] > 0.)
                             {
                              weights[0][idim] = 0.0;
                              weights[1][idim] = 1.0 - pn_sep[idim];
                              weights[2][idim] =       pn_sep[idim];
                             }
                           else
                             {
                              weights[0][idim] =      -pn_sep[idim];
                              weights[1][idim] = 1.0 + pn_sep[idim];
                              weights[2][idim] = 0.0;
                             }
#endif
#ifdef NGP
                           weights[0][idim] = 0.0;
                           weights[1][idim] = 1.0;
                           weights[2][idim] = 0.0;
#endif
                          }
                        
                        /* calculate correct mass to density conversion factor according to AMIGA setup */
                        mass2dens = 
#ifdef MULTIMASS
                           (double)cur_part->weight * 
#endif
                           cur_grid->masstodens;
                        
                        /* assign mass */
                        for(k = 0; k < 3; k++){
                           for(j = 0; j < 3; j++){
                              for(i = 0; i < 3; i++){
                                 if(tsc_nodes[k][j][i] != NULL)
                                   {
                                    tsc_nodes[k][j][i]->dens += mass2dens * weights[k][Z]*weights[j][Y]*weights[i][X];
                                   } 
                              } } } 
                        
                       } /* cur_part */
                    } /* cur_node */
                     
               }
            }
         }
      }
   }
                     
        
   cur_grid->size.no_part  = no_part;
   cur_grid->size.no_nodes = no_nodes;

   /* do not allow refinements smaller than MIN_NNODES */
   if(cur_grid->size.no_nodes < MIN_NNODES)
      return(FALSE);
   
   /* check mean ratio of (number of nodes)/(number of particles) on grid */
   if(simu.np_limit == FALSE)
      return(TRUE);
   else
     {
      if(cur_grid->l1dim == 2*simu.NGRID_DOM && simu.Nth_dom < 3.)
         return(TRUE);
      
      ratio = (double)cur_grid->size.no_nodes/(double)cur_grid->size.no_part;
      if(ratio > NP_RATIO)
         return(FALSE);
      else
         return(TRUE);
     }
   
}

/*===========================================================================
* assign "number of particle density" to grid
*===========================================================================*/
boolean assign_npart(gridls *cur_grid)
{
   pqptr cur_pquad;
   cqptr cur_cquad, icur_cquad;
   nqptr cur_nquad, icur_nquad;
   nptr  cur_node;
   nptr  tsc_nodes[3][3][3];
   long  x, y, z;
   
   partptr cur_part;           /* current particle                   */
   int     idim;               /* coord changing index               */
   int     i,j,k;              /* indices for 3D arrays              */
   double  pnarg_a;            /* pn sep temp arg                    */
   double  pnarg_b;            /* pn sep temp arg                    */
   double  tpnarg;             /* temp to calc pn args               */
   dvect   xyz_coords;         /* actual float coords of node        */
   dvect   pn_sep;             /* particle-node separation / spacing */
   dvect   weights[3];         /* weights in each dimension          */
   dvect   temp_coords;        /* float vector - calc un_mod coords  */   
   double  mass2dens, ratio, cur_shift;
   
   long    no_part, no_nodes;
   long    ipquad;
   
#ifdef EXTRAE_API_USAGE
  Extrae_user_function(1);
#endif

   /* reset counters for no_parts and no_nodes on cur_grid */
   no_part  = 0;
   no_nodes = 0;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   /* loop over whole grid */
#ifdef WITH_OPENMP2
  /* NOTE: this parallelisation is not 100% functional!
   *
   *       Actually it should be safe as we assign different pquad's to different threads and hence we should be able
   *       to safely modify tsc_nodes[][][]->dens without any collisions!
   *       But with the current refinement scheme (i.e. adding a boundary layer of 2 cells) it can happen that
   *       two linked quad's are actually spatially connected which is >>not<< reflected in the quad-structure!
   *
   *       We could use "#pragma omp atomic" to ensure proper "+=" for the density, but tests have shown that this
   *       dramatically slows down the runtime!
   *
   *       In that regards, feel free to switch this OpenMP parallelisation on at your own responsibility...
   */
#pragma omp parallel reduction(+:no_part, no_nodes) firstprivate(cur_shift) private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, x, y, z, tsc_nodes, cur_part, idim, i, j, k, pnarg_a, pnarg_b, tpnarg, xyz_coords, pn_sep, weights, temp_coords, mass2dens) shared(cur_grid)
#pragma omp for schedule(static)
      for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)
        {
         cur_pquad=cur_grid->pquad_array[ipquad];
#else
      for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
        {
#endif
            z = cur_pquad->z;
            
            for(cur_cquad=cur_pquad->loc;
                cur_cquad<cur_pquad->loc+cur_pquad->length;
                cur_cquad++, z++){
               for(icur_cquad=cur_cquad;icur_cquad!=NULL;icur_cquad=icur_cquad->next){
                  y = icur_cquad->y;
                  
                  for(cur_nquad=icur_cquad->loc;
                      cur_nquad<icur_cquad->loc+icur_cquad->length;
                      cur_nquad++, y++){
                    for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next){
                      x = icur_nquad->x;
                      
                      for(cur_node = icur_nquad->loc;
                          cur_node < icur_nquad->loc + icur_nquad->length;
                          cur_node++, x++)
                      {
                        /* increment number of nodes associated with current grid */
                        no_nodes++;
                        
                        /* calculate realspace coords of cur_node */
                        temp_coords[X] = (((double)x) / (double)cur_grid->l1dim) + cur_shift;
                        temp_coords[Y] = (((double)y) / (double)cur_grid->l1dim) + cur_shift;
                        temp_coords[Z] = (((double)z) / (double)cur_grid->l1dim) + cur_shift;
                        
                        xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
                        xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
                        xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);
                        
                        /* obtain 26 neighbours for TSC mass assignment */
                        tsc_nodes[1][1][1] = cur_node;
                        get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                        
                        /* loop over all particles attached to cur_node */
                        for(cur_part = cur_node->ll; cur_part != NULL; cur_part = cur_part->ll)
                        {
#if (defined GAS_PARTICLES && defined AHFdmonlypeaks)
                          if(fabs(cur_part->u-PDM) < ZERO)
#endif
                          {
                            
                            /* increment number of particles associated with current grid */
                            no_part++;
                            
                            /* calc fraction to be assigned to each tsc node */
                            for(idim = 0; idim < 3; idim++)
                            {
                              pn_sep[idim] = ((double)cur_part->pos[idim] - xyz_coords[idim]) * (double)cur_grid->l1dim;
                              
                              /* deal with periodic boundary conditions */
                              if(fabs(pn_sep[idim]) > 0.5*(double)cur_grid->l1dim)
                              {
                                tpnarg       = (double)cur_part->pos[idim] + 0.5;
                                pnarg_a      = f1mod(tpnarg+1.0, 1.0);
                                tpnarg       = xyz_coords[idim] + 0.5;
                                pnarg_b      = f1mod(tpnarg+1.0, 1.0);
                                pn_sep[idim] = (pnarg_a - pnarg_b) * (double)cur_grid->l1dim;
                              }
                              
#ifdef TSC
                              weights[1][idim] = 0.75 - pow2(pn_sep[idim]);
                              weights[0][idim] = pow2((0.5 - pn_sep[idim]))/2;
                              weights[2][idim] = pow2((0.5 + pn_sep[idim]))/2;
#endif
#ifdef CIC
                              if(pn_sep[idim] > 0.)
                              {
                                weights[0][idim] = 0.0;
                                weights[1][idim] = 1.0 - pn_sep[idim];
                                weights[2][idim] =       pn_sep[idim];
                              }
                              else
                              {
                                weights[0][idim] =      -pn_sep[idim];
                                weights[1][idim] = 1.0 + pn_sep[idim];
                                weights[2][idim] = 0.0;
                              }
#endif
#ifdef NGP
                              weights[0][idim] = 0.0;
                              weights[1][idim] = 1.0;
                              weights[2][idim] = 0.0;
#endif
                            }
                            
                            /* calculate correct particle to density conversion factor according to AMIGA setup */
#if (defined REFINE_BARYONIC_MASS && defined GAS_PARTICLES)
                            if(cur_part->u >= PGAS || cur_part->u == PSTAR)
                            mass2dens = (double)cur_part->weight * cur_grid->masstodens;
                            else
                            mass2dens = cur_grid->masstopartdens;
#else /* REFINE_BARYONIC_MASS && GAS_PARTICLES */
                            mass2dens = cur_grid->masstopartdens;
#endif /* REFINE_BARYONIC_MASS && GAS_PARTICLES */
                            
                            /* assign mass */
                            for(k = 0; k < 3; k++){
                              for(j = 0; j < 3; j++){
                                for(i = 0; i < 3; i++){
                                  if(tsc_nodes[k][j][i] != NULL)
                                  {
                                    tsc_nodes[k][j][i]->dens += mass2dens * weights[k][Z]*weights[j][Y]*weights[i][X];
                                  }
                                } } }
                          } /* if(DM particle) */
                        } /* cur_part */
                      } /* cur_node */
                      
                    }
                  }
               }
            }
           }
          
        
   cur_grid->size.no_part  = no_part;
   cur_grid->size.no_nodes = no_nodes;

#ifdef EXTRAE_API_USAGE
  Extrae_user_function(0);
#endif

   /* do not allow refinements smaller than MIN_NNODES */
   if(cur_grid->size.no_nodes < MIN_NNODES)
      return(FALSE);
   
   /* check mean ratio of (number of nodes)/(number of particles) on grid */
   if(simu.np_limit == FALSE)
      return(TRUE);
   else
     {
      if(cur_grid->l1dim == 2*simu.NGRID_DOM && simu.Nth_dom < 3.)
         return(TRUE);
      
      ratio = (double)cur_grid->size.no_nodes/(double)cur_grid->size.no_part;
      if(ratio > NP_RATIO)
         return(FALSE);
      else
         return(TRUE);
     }
   
}

/*===========================================================================
 * zero density field on given grid
 *===========================================================================*/
void zero_dens(gridls *cur_grid)
{
   pqptr cur_pquad;
   cqptr cur_cquad, icur_cquad;
   nqptr cur_nquad, icur_nquad;
   nptr  cur_node;
   long  ipquad;
   
#ifdef EXTRAE_API_USAGE
  Extrae_user_function(1);
#endif

#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node) shared(cur_grid)
#pragma omp for schedule(static)
  for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)     
    {
     cur_pquad=cur_grid->pquad_array[ipquad];
#else
  for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
    {
#endif
         
      for(cur_cquad=cur_pquad->loc;
          cur_cquad<cur_pquad->loc+cur_pquad->length;
          cur_cquad++)
         for(icur_cquad=cur_cquad;icur_cquad!=NULL;icur_cquad=icur_cquad->next)
            
            for(cur_nquad=icur_cquad->loc;
                cur_nquad<icur_cquad->loc+icur_cquad->length;
                cur_nquad++)
               for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
                  
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++)
                    {
                     cur_node->dens          = -simu.mean_dens;
#ifndef AHFlean
                     cur_node->force.temp[0] = 0.0;
                     cur_node->force.temp[1] = 0.0;
#else
                     cur_node->force.tempdens = 0.0;
#endif
                    }
    }
#ifdef EXTRAE_API_USAGE
  Extrae_user_function(0);
#endif
}


/*===============================================================================
* assign mass of a single particle to given grid
*===============================================================================*/
void assign_part(gridls *cur_grid, pqptr fst_pquad, cqptr fst_cquad, 
                 nqptr fst_nquad, nptr cur_node, long z, long y, long x, 
                 partptr cur_part)
{
   nptr   tsc_nodes[3][3][3]; /* nodes to assign to                 */
   dvect  xyz_coords;         /* actual float coords of node        */
   dvect  pn_sep;             /* particle-node separation / spacing */
   int    idim;               /* coord changing index               */
   int    i,j,k;              /* indices for 3D arrays              */
   dvect  weights[3];         /* weights in each dimension          */
   double pnarg_a;            /* pn sep temp arg                    */
   double pnarg_b;            /* pn sep temp arg                    */
   double tpnarg;             /* temp to calc pn args               */
   dvect  temp_coords;        /* float vector - calc un_mod coords  */
   pqptr  cur_pquad;          /* current pquad                      */
   cqptr  cur_cquad;          /* current cquad                      */
   nqptr  cur_nquad;          /* current nquad                      */
   double mass2dens;
   double cur_shift;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   /* current node in centre */
   tsc_nodes[1][1][1] = cur_node;
   
   /* calculate realspace coords of cur_node */
   temp_coords[X] = (((double)x) / (double)cur_grid->l1dim) + cur_shift;
   temp_coords[Y] = (((double)y) / (double)cur_grid->l1dim) + cur_shift;
   temp_coords[Z] = (((double)z) / (double)cur_grid->l1dim) + cur_shift;
   
   xyz_coords[X] = f1mod(temp_coords[X]+1.0, 1.0);
   xyz_coords[Y] = f1mod(temp_coords[Y]+1.0, 1.0);
   xyz_coords[Z] = f1mod(temp_coords[Z]+1.0, 1.0);
   
   /* get correct pointers to current quads */
   for(cur_pquad = fst_pquad; z > cur_pquad->z + cur_pquad->length; 
       cur_pquad	= cur_pquad->next)
      ;
   for(cur_cquad = fst_cquad; y > cur_cquad->y + cur_cquad->length; 
       cur_cquad	= cur_cquad->next)
      ;
   for(cur_nquad = fst_nquad; x > cur_nquad->x + cur_nquad->length; 
       cur_nquad	= cur_nquad->next)
      ;
   
   /* get pointers to tsc nodes */
   get_TSCnodes(cur_grid, cur_pquad, cur_cquad, cur_nquad, tsc_nodes, &z, &y, &x);
   
   /* calc fraction to be assigned to each tsc node */
   
   for(idim = 0; idim < NDIM; idim++)
     {
      pn_sep[idim] = ((double)cur_grid->l1dim) * 
      ((double)cur_part->pos[idim]-xyz_coords[idim]);
      
      if(fabs(pn_sep[idim]) > 0.5*(double)cur_grid->l1dim)
        {
         tpnarg          = (double)cur_part->pos[idim] + 0.5;
         pnarg_a         = f1mod(tpnarg+1.0, 1.0);
         tpnarg          = xyz_coords[idim] + 0.5;
         pnarg_b         = f1mod(tpnarg+1.0, 1.0);
         pn_sep[idim] = (pnarg_a - pnarg_b) * (double)cur_grid->l1dim;
        }
#ifdef TSC
      weights[1][idim] = 0.75 - pow2(pn_sep[idim]);
      weights[0][idim] = pow2((0.5 - pn_sep[idim]))/2;
      weights[2][idim] = pow2((0.5 + pn_sep[idim]))/2;
#endif
#ifdef CIC
      if(pn_sep[idim] > 0.)
        {
         weights[0][idim] = 0.0;
         weights[1][idim] = 1.0 - pn_sep[idim];
         weights[2][idim] =       pn_sep[idim];
        }
      else
        {
         weights[0][idim] =      -pn_sep[idim];
         weights[1][idim] = 1.0 + pn_sep[idim];
         weights[2][idim] = 0.0;
        }
#endif
#ifdef NGP
      weights[0][idim] = 0.0;
      weights[1][idim] = 1.0;
      weights[2][idim] = 0.0;
#endif
     }
   
   /* calculate correct mass to density conversion factor according to AMIGA setup */
   mass2dens = 
#ifdef MULTIMASS
      (double)cur_part->weight * 
#endif
      cur_grid->masstodens;
   
   /* assign mass */  
   for(k = 0; k < 3; k++){
      for(j = 0; j < 3; j++){
         for(i = 0; i < 3; i++){
            if(tsc_nodes[k][j][i] != NULL)
              {
               tsc_nodes[k][j][i]->dens += mass2dens * weights[k][Z]*weights[j][Y]*weights[i][X];
              } 
         } } }
}

/*===============================================================================
* un-assign mass of a single particle from given grid
*===============================================================================*/
void unassign_part(gridls *cur_grid, pqptr fst_pquad, cqptr fst_cquad, 
                   nqptr fst_nquad, nptr cur_node, long z, long y, long x, 
                   partptr cur_part)
{
   nptr   tsc_nodes[3][3][3];  /* nodes to assign to                 */
   int    idim;                /* coord changing index               */
   int    i,j,k;               /* indices for 3D arrays              */
   pqptr  cur_pquad;           /* current pquad                      */
   cqptr  cur_cquad;           /* current cquad                      */
   nqptr  cur_nquad;           /* current nquad                      */
   double pnarg_a;             /* pn sep temp arg                    */
   double pnarg_b;             /* pn sep temp arg                    */
   double tpnarg;              /* temp to calc pn args               */
   dvect  temp_coords;         /* float vector - calc un_mod coords  */
   dvect  xyz_coords;          /* actual float coords of node        */
   dvect  pn_sep;              /* particle-node separation / spacing */
   dvect  weights[3];          /* weights in each dimension          */
   double mass2dens;
   double cur_shift;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   /* current node in centre */
   tsc_nodes[1][1][1] = cur_node;
   
   /* calculate realspace coords of cur_node */
   temp_coords[X] = (((double)x) / (double)cur_grid->l1dim) + cur_shift;
   temp_coords[Y] = (((double)y) / (double)cur_grid->l1dim) + cur_shift;
   temp_coords[Z] = (((double)z) / (double)cur_grid->l1dim) + cur_shift;
   
   xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
   xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
   xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);
   
   /* get correct pointers to current quads */
   for(cur_pquad = fst_pquad; z > cur_pquad->z+cur_pquad->length; cur_pquad = cur_pquad->next)
      ;
   for(cur_cquad = fst_cquad; y > cur_cquad->y+cur_cquad->length; cur_cquad = cur_cquad->next)
      ;
   for(cur_nquad = fst_nquad; x > cur_nquad->x+cur_nquad->length; cur_nquad = cur_nquad->next)
      ;
   
   /* get pointers to tsc nodes */
   get_TSCnodes(cur_grid, cur_pquad, cur_cquad, cur_nquad, tsc_nodes, &z, &y, &x);
   
   /* calc fraction to be un-assigned from each tsc node */
   for(idim = 0; idim < NDIM; idim++)
     {
      pn_sep[idim] = ((double)cur_grid->l1dim) * (cur_part->pos[idim] - xyz_coords[idim]);
      
      if(fabs(pn_sep[idim]) > 0.5*(double)cur_grid->l1dim)
        {
         tpnarg          = (double)cur_part->pos[idim] + 0.5;
         pnarg_a         = f1mod(tpnarg+1.0, 1.0);
         tpnarg          = xyz_coords[idim] + 0.5;
         pnarg_b         = f1mod(tpnarg+1.0, 1.0);
         pn_sep[idim]    = (pnarg_a - pnarg_b) * (double)cur_grid->l1dim;
        }
#ifdef TSC
      weights[1][idim] = 0.75 - pow2(pn_sep[idim]);
      weights[0][idim] = pow2((0.5 - pn_sep[idim]))/2;
      weights[2][idim] = pow2((0.5 + pn_sep[idim]))/2;
#endif
#ifdef CIC
      if(pn_sep[idim] > 0.)
        {
         weights[0][idim] = 0.0;
         weights[1][idim] = 1.0 - pn_sep[idim];
         weights[2][idim] =       pn_sep[idim];
        }
      else
        {
         weights[0][idim] =     - pn_sep[idim];
         weights[1][idim] = 1.0 + pn_sep[idim];
         weights[2][idim] = 0.0;
        }
#endif
#ifdef NGP
      weights[0][idim] = 0.0;
      weights[1][idim] = 1.0;
      weights[2][idim] = 0.0;
#endif
     }
   
   /* calculate correct mass to density conversion factor according to AMIGA setup */
   mass2dens = 
#ifdef MULTIMASS
      (double)cur_part->weight * 
#endif
      cur_grid->masstodens;
   
   /* un-assign mass */
   for(k = 0; k < 3; k++){
      for(j = 0; j < 3; j++){
         for(i = 0; i < 3; i++){
            if(tsc_nodes[k][j][i] != NULL)
              {
               tsc_nodes[k][j][i]->dens -= mass2dens * weights[k][Z]*weights[j][Y]*weights[i][X];
              }
         } } }
}


/*==========================================================================
* add P1/P2 contribution to already existing P3 density
*==========================================================================*/
void overlap_dens(gridls *cur_grid)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    tsc_nodes[3][3][3];
   nptr    cur_node;
   long    x, y, z, ipquad;
   
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, tsc_nodes, x, y, z) shared(cur_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)   
     {
      cur_pquad=cur_grid->pquad_array[ipquad];
#else
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
     {
#endif
       
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
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                     if(test_tsc(tsc_nodes) == TRUE)
                       {
                        /* eventually add contribution from coa_grid to fin_grid */
#ifndef AHFlean
                        cur_node->dens += cur_node->force.temp[0];
#else
                        cur_node->dens += cur_node->force.tempdens;
#endif
                        
                        /* reset temporary storage */
#ifndef AHFlean
                        cur_node->force.temp[0] = 0.0;
#else
                        cur_node->force.tempdens = 0.0;
#endif
                       }
                    }
     }
}


/*============================================================================
* restore P1/P2 density on coa_grid by assign_dens()
*============================================================================*/
void restore_dens(gridls *cur_grid)
{
   /* re-assign density on next coarser level */
   zero_dens(cur_grid);
   assign_dens(cur_grid);
}

/*=============================================================================
* refill the density whole within cur_grid caused be an existing finer grid
*=============================================================================*/
void refill_dens(gridls *coa_grid)
{
   gridls *fin_grid;
   
   fin_grid = coa_grid + 1;
   
   /* re-calculate P3 density field on fine grid
      (fin_node->dens) */
   zero_dens(fin_grid);
   assign_dens(fin_grid);
   
   /* store (P1/P2) density spilling from coarse to fine grid 
      (coa_node->dens -> fin_node->force.temp[0]) */
   fc_overlap(fin_grid);
   
   /* interpolate P3 density from fine to coarse grid
      (fin_node->dens -> coa_node->dens) */
   f2c_dens(fin_grid);
   
   /* add stored P1/P2 contribution to already existing P3 density
      (fin_node->dens += fin_node->force.temp[0]) */
   overlap_dens(fin_grid);
}

/*==============================================================================
* reflect density on cur_grid onto *all* coarser grids
*==============================================================================*/
void reflect_dens(gridls *cur_grid, boolean reflect)
{
   gridls *for_grid;
   gridls *first_grid, *last_grid;
   
   
   /* no re-assignment of cur_grid density field needed...  */
   /* ...we just re-calculated 'P3 density' within step() ! */
   
   /* injection to next coarser levels needed */
   first_grid = cur_grid;
   if(cur_grid == global.dom_grid+1) /* fc_overlap etc. is accessing for_grid-1 */
      last_grid = cur_grid;
   else
      last_grid = cur_grid-1;
   
   for(for_grid = first_grid; for_grid >= last_grid; for_grid--)
     {
      /* store (P1/P2) density spilling from next coarser grid to for_grid 
      ([for-1]_node->dens -> for_node->force.temp[0]) */
      fc_overlap(for_grid);
      
      if(reflect)
         /* interpolate density from for_grid to next coarser grid 
         ([for-1]_node->dens = for_node->dens) */
         f2c_dens(for_grid);
      
      /* add stored P1/P2 contribution to already existing P3 density 
         (for_node->dens += for_node->force.temp[0]) */
      overlap_dens(for_grid);
     }
}

/*==============================================================================
* get density correct on all grids from first_grid to last_grid
*==============================================================================*/
void stack_dens(gridls *first_grid, gridls *last_grid)
{
   gridls *for_grid;
   
   last_grid = MAX(global.dom_grid+1,last_grid);
   
   for(for_grid = first_grid; for_grid >= last_grid; for_grid--)
     {
      /* store (P1/P2) density spilling from next coarser grid to for_grid */
      fc_overlap(for_grid);
      
      /* interpolate density from for_grid to next coarser grid */
      f2c_dens(for_grid);
      
      /* add stored P1/P2 contribution to already existing P3 density */
      overlap_dens(for_grid);
     }
}

#ifndef AHFlean
/*===========================================================================
 * zero potential field on given grid
 *===========================================================================*/
void zero_temp1(gridls *cur_grid)
 {
  pqptr cur_pquad;
  cqptr cur_cquad, icur_cquad;
  nqptr cur_nquad, icur_nquad;
  nptr  cur_node;
  
  nptr tsc_nodes[3][3][3];
  long x, y, z;
  
  long ipquad;
  
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, tsc_nodes, x, y, z) shared(cur_grid)
#pragma omp for schedule(static)
  for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)     
   {
    cur_pquad=cur_grid->pquad_array[ipquad];
#else
  for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
   {
#endif
    z = cur_pquad->z;
    
    for(cur_cquad=cur_pquad->loc;
        cur_cquad<cur_pquad->loc+cur_pquad->length;
        cur_cquad++, z++)
      for(icur_cquad=cur_cquad;icur_cquad!=NULL;icur_cquad=icur_cquad->next)
       {
        y = icur_cquad->y;
        
        for(cur_nquad=icur_cquad->loc;
            cur_nquad<icur_cquad->loc+icur_cquad->length;
            cur_nquad++, y++)
          for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
           {
            x = icur_nquad->x;
            
            for(cur_node = icur_nquad->loc; 
                cur_node < icur_nquad->loc + icur_nquad->length; 
                cur_node++, x++)
             {
              cur_node->force.temp[1] = 0.0;
              
#ifdef F2C_POT
              tsc_nodes[1][1][1] = cur_node;
              get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
              
              /* do not tinker with boundary values! */
              if(test_tsc(tsc_nodes) == TRUE)
                cur_node->pot = 0.0;
#endif
             }
           }
       }
   }
}
#endif /* AHFlean */

