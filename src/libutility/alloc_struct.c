#include <stdio.h>
#include <stdlib.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "utility.h"

/*-----------------------------------------------------------------
* dest_node: destroy node pointed at by *ptr and set *ptr to NULL
*-----------------------------------------------------------------*/
void dest_node(nptr *ptr)
{
   free(*ptr);
   *ptr = NULL;
}

/*-------------------------------------------------------------------
* dest_nquad: destroy nquad pointed at by *ptr and set *ptr to NULL
*-------------------------------------------------------------------*/
void dest_nquad(nqptr *ptr)
{
   free(*ptr);
   *ptr = NULL;
}

/*-------------------------------------------------------------------
* dest_cquad: destroy cquad pointed at by *ptr and set *ptr to NULL
*-------------------------------------------------------------------*/
void dest_cquad(cqptr *ptr)
{
   free(*ptr);
   *ptr = NULL;
}

/*-------------------------------------------------------------------
* dest_pquad: destroy pquad pointed at by *ptr and set *ptr to NULL
*-------------------------------------------------------------------*/
void dest_pquad(pqptr *ptr)
{
   free(*ptr);
   *ptr = NULL;
}

/*----------------------------------------------------------------------------
* c_part:
* Create block_size particles and initialize structure members to zero. This
* function is essentially a cover function for calloc(block_size,
                                                      * sizeof(part)). The position and momentum members of all particle structures
* are set equal to 0.0 and the linked list pointer is set to NULL.
*----------------------------------------------------------------------------*/
partptr c_part(long block_size)
{
   partptr   pptr;     /* pointer to new particles                     */
   partptr   i_pptr;   /* pointer to increment down block of particles */
   int       i;        /* coordinate index (X,Y,Z)                     */
   
   /* generate block_size particles */
   if((pptr = (partptr) calloc(block_size, sizeof(part))) == NULL)
     {
      fprintf(stderr,"c_part: error callocing %ld particles, sizeof(part)=%d\n",block_size,sizeof(part));
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"c_part: error callocing %ld particles, sizeof(part)=%d\n",block_size,sizeof(part));
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /*
    * loop over nodes filling pot and dens with 0.0 and setting linked-list
    * pointer to NULL.
    */
   
//   for(i_pptr = pptr; i_pptr < pptr + block_size; i_pptr++)
//     {
//      for(i = X; i <= Z; i++)
//        {
//         i_pptr->pos[i] = 0.0;
//         i_pptr->mom[i] = 0.0;
//        }
//      i_pptr->ll = NULL;
//     }

  return pptr;
}

/*----------------------------------------------------------------------------
 * c_gas:
 * Create block_size gas particles and initialize structure members to zero.
 *----------------------------------------------------------------------------*/
gasptr c_gas(long block_size)
{
   gasptr   pptr;     /* pointer to new particles                     */
   gasptr   i_pptr;   /* pointer to increment down block of particles */
   int      i;        /* coordinate index (X,Y,Z)                     */
   
   /* generate block_size particles */
   if((pptr = (gasptr) calloc(block_size, sizeof(gas))) == NULL)
     {
      fprintf(stderr,"c_gas: error callocing gas particles\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"c_gas: error callocing gas particles\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /*
    * loop over nodes filling pot and dens with 0.0 and setting linked-list
    * pointer to NULL.
    */
//   for(i_pptr = pptr; i_pptr < pptr + block_size; i_pptr++)
//     {
//      for(i = X; i <= Z; i++)
//        {
//         i_pptr->u = 0.0;
//        }
//     }

  return pptr;
}

/*-----------------------------------------------------------------------------
* c_node:
* Create block_size nodes and initialize structure members to zero. This
* function is essentially a cover function for calloc(block_size,
                                                      * sizeof(node)). The position and momentum members of all particle structures
* are set equal to 0.0 and the linked list pointer is set to NULL.
*-----------------------------------------------------------------------------*/
nptr c_node(long block_size)
{
   nptr node_ptr;     /* pointer to new nodes */
   nptr i_node_ptr;   /* pointer to increment down block */
   int i;             /* coordinate index (X,Y,Z) */
   
   /* generate nodes */
   if((node_ptr = (nptr) calloc(block_size, sizeof(node))) == NULL)
     {
      fprintf(stderr,"c_node: error callocing nodes\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"c_node: error callocing nodes\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /*
    * loop over nodes filling pot, dens and force.forces vector components with
    * 0.0. The head-of-chain pointer (ll) is set to NULL.
    */
   
//   for(i_node_ptr = node_ptr; i_node_ptr < node_ptr + block_size; i_node_ptr++)
//     {
//#ifndef AHFlean
//      i_node_ptr->pot     = 0.0;
//#endif
//      i_node_ptr->dens    = 0.0;
//#ifndef AHFlean
//      for(i = X; i <=Z; i++)
//         i_node_ptr->force.forces[i] = 0.0;
//#endif
//      i_node_ptr->ll      = NULL;
//     }

  return node_ptr;
}

/*------------------------------------------------------------------------------
* r_node:
* Change size of block of nodes to new_size and initialize structure members
* to zero. This function is a cover function for realloc(ptr, new_size *
                                                         * sizeof(node)). The pot, dens and force.forces vector of ALL node
* structures are set equal to 0.0 and the linked list pointer is set to NULL.
*------------------------------------------------------------------------------*/
nptr r_node(nptr node_ptr, long new_size)
{
   int i;              /* loop index */
   nptr i_node_ptr;    /* pointer to point at new nodes */
   
   /* generate new nodes */
   if((node_ptr = (nptr) realloc(node_ptr, new_size * sizeof(node))) == NULL)
     {
      fprintf(stderr,"r_node: error reallocing nodes\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"r_node: error reallocing nodes\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /* loop over all nodes setting structure members to zero */
   if(new_size > 0)
     {
      for(i_node_ptr = node_ptr; i_node_ptr < node_ptr + new_size; i_node_ptr++)
        {
#ifndef AHFlean
         i_node_ptr->pot     = 0.0;
#endif
         i_node_ptr->dens    = 0.0;
#ifndef AHFlean
         for(i = X; i <=Z; i++)
            i_node_ptr->force.forces[i] = 0.0;
#endif
         i_node_ptr->ll      = NULL;
        }
     }
   return node_ptr;
}

/*-------------------------------------------------------------------------------
* c_nquad:
* Create block_size nquads and initialize structure members. This function is
* a cover function for calloc(block_size, sizeof(nquad)). The loc and next
* pointers are set equal to NULL, while the x coordinate and length members of
* the nquad structure are set equal to 0.
*-------------------------------------------------------------------------------*/
nqptr c_nquad(long block_size)
{
   nqptr nquad_ptr;     /* pointer to new nquads           */
   nqptr i_nquad_ptr;   /* pointer to increment down block */
   
   /* generate nquads */
   if((nquad_ptr = (nqptr) calloc(block_size, sizeof(nquad))) == NULL)
     {
      fprintf(stderr,"c_nquad: error callocing nquads\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"c_nquad: error callocing nquads\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /* loop over nquads initializing members of structure */
//   for(i_nquad_ptr = nquad_ptr; i_nquad_ptr < nquad_ptr + block_size; i_nquad_ptr++)
//     {
//      i_nquad_ptr->loc    = NULL;
//      i_nquad_ptr->x      = 0;
//      i_nquad_ptr->length = 0;
//      i_nquad_ptr->next   = NULL;
//     }

  return nquad_ptr;
}

/*-------------------------------------------------------------------------------
* r_nquad:
* Change size of block of nquads to new_size and initialize structure members.
* This function is a cover function for realloc(ptr, new_size * sizeof(nquad)).
* NEW nquads are initialized. The loc and next pointers are set equal to NULL,
* while the x coordinate and length members of the nquad structure are set
* equal to 0.
*-------------------------------------------------------------------------------*/
nqptr r_nquad(nqptr nquad_ptr, long old_size, long new_size)
{
   nqptr new_nquad_ptr, old_nquad_ptr;
   
   /* destroy unused pointers */
   if(new_size < old_size)
     {
      for(old_nquad_ptr = nquad_ptr+new_size; old_nquad_ptr < nquad_ptr+old_size;
          old_nquad_ptr++)
         free_nquad(old_nquad_ptr);
     }
   
   /* re-size block of nquads to new_size */
   if((nquad_ptr = (nqptr) realloc(nquad_ptr, new_size * sizeof(nquad))) == NULL)
     {
      fprintf(stderr,"r_nquad: error reallocing nquad\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"r_nquad: error reallocing nquad\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /* If there are new nquads initialize them */
   if(new_size > old_size)
     {
      for(new_nquad_ptr = nquad_ptr + old_size; new_nquad_ptr < nquad_ptr + new_size;
          new_nquad_ptr++)
        {
         new_nquad_ptr->loc    = NULL;
         new_nquad_ptr->x      = 0;
         new_nquad_ptr->length = 0;
         new_nquad_ptr->next   = NULL;
        }
     }
   return nquad_ptr;
}

/*-------------------------------------------------------------------------------
* c_cquad:
* Create block_size cquads and initialize structure members. This function is
* a cover function for calloc(block_size, sizeof(cquad)). The loc and next
* pointers are set equal to NULL, while the y coordinate and length members of
* the cquad structure are set equal to 0.
*-------------------------------------------------------------------------------*/
cqptr c_cquad(long block_size)
{
   cqptr cquad_ptr;     /* pointer to new cquads */
   cqptr i_cquad_ptr;   /* pointer to increment down block */
   
   /* generate block_size cquads */
   
   if((cquad_ptr = (cqptr) calloc(block_size, sizeof(cquad))) == NULL)
     {
      fprintf(stderr,"c_cquad: error callocing cquads\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"c_cquad: error callocing cquads\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /* loop over cquads initializing */
//   for(i_cquad_ptr = cquad_ptr; i_cquad_ptr < cquad_ptr + block_size; i_cquad_ptr++)
//     {
//      i_cquad_ptr->loc = NULL;
//      i_cquad_ptr->y = 0;
//      i_cquad_ptr->length = 0;
//      i_cquad_ptr->next = NULL;
//     }

  return cquad_ptr;
}

/*-------------------------------------------------------------------------------
* r_cquad:
* Change size of block of cquads to new_size and initialize structure members.
* This function is a cover function for realloc(ptr, new_size * sizeof(cquad)).
* NEW cquads are initialized. The loc and next pointers are set equal to NULL,
* while the y coordinate and length members of the nquad structure are set
* equal to 0.
*-------------------------------------------------------------------------------*/
cqptr r_cquad(cqptr cquad_ptr, long old_size, long new_size)
{
   cqptr new_cquad_ptr, old_cquad_ptr;
   
   /* destroy unused cquads */
   if(new_size < old_size)
     {
      for(old_cquad_ptr = cquad_ptr+new_size; old_cquad_ptr < cquad_ptr+old_size;
          old_cquad_ptr++)
         free_cquad(old_cquad_ptr);
     }
   
   
   /* generate new cquads */
   if((cquad_ptr = (cqptr) realloc(cquad_ptr, new_size*sizeof(cquad)))  == NULL)
     {
      fprintf(stderr,"r_cquad: error reallocing cquad\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"r_cquad: error reallocing cquad\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /* If there are new cquads initialize them */
   if(new_size > old_size)
     {
      for(new_cquad_ptr = cquad_ptr + old_size; 
          new_cquad_ptr < cquad_ptr + new_size; 
          new_cquad_ptr++)
        {
         new_cquad_ptr->loc    = NULL;
         new_cquad_ptr->y      = 0;
         new_cquad_ptr->length = 0;
         new_cquad_ptr->next   = NULL;
        }
     }
   return cquad_ptr;
}

/*-------------------------------------------------------------------------------
* c_pquad:
* Create block_size pquads and initialize structure members. This function is
* a cover function for calloc(block_size, sizeof(pquad)). The loc and next
* pointers are set equal to NULL, while the y coordinate and length members of
* the pquad structure are set equal to 0.
*-------------------------------------------------------------------------------*/
pqptr c_pquad(long block_size)
{
   pqptr pquad_ptr;     /* pointer to new pquads */
   pqptr i_pquad_ptr;   /* pointer to increment down block */
   
   /* generate new pquads */
   
   if((pquad_ptr = (pqptr) calloc(block_size, sizeof(pquad))) == NULL)
     {
      fprintf(stderr,"c_pquad: error callocing pquads\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"c_pquad: error callocing pquads\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /* loop over pquads initializing */
//   for(i_pquad_ptr = pquad_ptr; i_pquad_ptr < pquad_ptr + block_size; i_pquad_ptr++)
//     {
//      i_pquad_ptr->loc = NULL;
//      i_pquad_ptr->z = 0;
//      i_pquad_ptr->length = 0;
//      i_pquad_ptr->next = NULL;
//     }

  return pquad_ptr;
}


/*-------------------------------------------------------------------------------
* r_pquad:
* Change size of block of pquads to new_size and initialize structure members.
* This function is a cover function for realloc(ptr, new_size * sizeof(pquad)).
* NEW pquads are initialized. The loc and next pointers are set equal to NULL,
* while the z coordinate and length members of the cquad structure are set
* equal to 0.
*-------------------------------------------------------------------------------*/
pqptr r_pquad(pqptr pquad_ptr, long old_size, long new_size)
{
   pqptr new_pquad_ptr, old_pquad_ptr;
   
   /* destroy unused pquads */
   if(new_size < old_size)
     {
      for(old_pquad_ptr = pquad_ptr+new_size; old_pquad_ptr < pquad_ptr+old_size;
          old_pquad_ptr++)
         free_pquad(old_pquad_ptr);
     }
   
   
   /* generate new pquads */
   if((pquad_ptr = (pqptr) realloc(pquad_ptr, new_size*sizeof(pquad)))  == NULL)
     {
      fprintf(stderr,"r_pquad: error reallocing pquad\n");
      fflush(stderr);
      fclose(stderr);
      fprintf(io.logfile,"r_pquad: error reallocing pquad\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
   
   /* If there are new pquads initialize them */
   if(new_size > old_size)
     {
      for(new_pquad_ptr = pquad_ptr + old_size; 
          new_pquad_ptr < pquad_ptr + new_size; 
          new_pquad_ptr++)
        {
         new_pquad_ptr->loc    = NULL;
         new_pquad_ptr->z      = 0;
         new_pquad_ptr->length = 0;
         new_pquad_ptr->next   = NULL;
        }
     }
   return pquad_ptr;
}

/*=============================================================================
 * alloc_quads: generate the full quad-hierarchy for a domain grid of size l1dim
 *=============================================================================*/
void alloc_quads(gridls *cur_grid, long l1dim)
{
  pqptr cur_pquad;          /* domain-grid pquad pointer */
  cqptr cur_cquad;          /* domain-grid cquad pointer */
  nqptr cur_nquad;          /* domain-grid nquad pointer */
  int   idim;
  
  /*======================================================================
   * create the grid, i.e. one pquad and all subsequent cquads and nquads
   *======================================================================*/
  cur_grid->pquad      = c_pquad(1);           /* generate pquad      */
  cur_pquad            = cur_grid->pquad;      /* set cur_pquad       */
  cur_pquad->length    = l1dim;                /* set length to l1dim */
  cur_pquad->loc       = c_cquad(l1dim);       /* generate cquads     */
  
  /* loop over cquads */
  for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + l1dim; cur_cquad++)
    {
      cur_cquad->length = l1dim;          /* set length to l1dim    */
      cur_cquad->loc    = c_nquad(l1dim); /* generate nquads        */
      
      /* loop over nquads */
      for(cur_nquad = cur_cquad->loc; cur_nquad < cur_cquad->loc + l1dim; cur_nquad++)
        {
          cur_nquad->length = l1dim;         /* set length to l1dim       */
          cur_nquad->loc    = c_node(l1dim); /* generate nodes            */
          /* -> l1dim^3 nodes in total */
        }
    }
  
  /* make pquad available in a linear array (actually only used with WITH_OPENMP) */
  cur_grid->no_pquad       = 1;
  cur_grid->pquad_array    = (pqptr *) calloc(1, sizeof(pqptr));
  cur_grid->pquad_array[0] = cur_grid->pquad;  
}



/*==============================================================================
 * free_nquad: free nodes attached to cur_nquad
 *==============================================================================*/
void free_nquad(nqptr cur_nquad)
{
  nptr fst_node;
  
  if(cur_nquad->next != NULL)
    {
      free_nquad(cur_nquad->next);
      dest_nquad(&(cur_nquad->next));
    }
  
  fst_node = cur_nquad->loc;
  dest_node(&(fst_node));
}

/*==============================================================================
 * free_cquad: free nquads attached to cur_cquad
 *==============================================================================*/
void free_cquad(cqptr cur_cquad)
{
  nqptr cur_nquad;       /* current nquad */
  
  if(cur_cquad->next != NULL)
    {
      free_cquad(cur_cquad->next);
      dest_cquad(&(cur_cquad->next));
    }
  
  /* loop over the nquads */
  for(cur_nquad = cur_cquad->loc; 
      cur_nquad < cur_cquad->loc + cur_cquad->length; 
      cur_nquad++)
    free_nquad(cur_nquad);
  
  cur_nquad = cur_cquad->loc;
  dest_nquad(&(cur_nquad));   /* destroy the block of nquads */
}


/*==============================================================================
 * free_pquad: free cquads attached to cur_pquad
 *==============================================================================*/
void free_pquad(pqptr cur_pquad)
{
  cqptr cur_cquad;       /* current nquad */
  
  if(cur_pquad != NULL)
    {
      if(cur_pquad->next != NULL)
        {
          free_pquad(cur_pquad->next);
          dest_pquad(&(cur_pquad->next));
        }
      
      /* loop over the cquads */
      for(cur_cquad = cur_pquad->loc; 
          cur_cquad<cur_pquad->loc + cur_pquad->length; 
          cur_cquad++)
        free_cquad(cur_cquad);
      
      cur_cquad = cur_pquad->loc;
      dest_cquad(&(cur_cquad));   /* destroy the block of cquads */
    }
}

/*==============================================================================
 * free_grid: free memory associated with supplied grid
 *==============================================================================*/
void free_grid(gridls *cur_grid, int *no_grids)
{
  /* pass cur_pquad to free_pquad */
  free_pquad(cur_grid->pquad);
  dest_pquad(&(cur_grid->pquad));
  
  /* also erase that linear pquad array */
  cur_grid->no_pquad = 0;
  if(cur_grid->pquad_array)
    free(cur_grid->pquad_array);
  cur_grid->pquad_array = NULL;
  
  /* reset all sorts of things to zero */
  cur_grid->multistep = 0;
    
  /*======================================================================
   * do NOT reset time and size of cur_grid...it's needed for the logfile
   *======================================================================*/
  
  /* reduce the number of grids */
  (*no_grids) -= 1;
}



#if (defined AHF || defined AHF2)
/*----------------------------------------------------------------------------
* c_profile:
* allocate memory for halo profile
*----------------------------------------------------------------------------*/
void c_profile(HALO *cur_halo, int nbins)
{
   cur_halo->prof.nbins   = nbins;
   cur_halo->prof.r       = (double *)       calloc(nbins, sizeof(double));
   cur_halo->prof.npart   = (long unsigned*) calloc(nbins, sizeof(long unsigned));
   cur_halo->prof.nvpart  = (double *)       calloc(nbins, sizeof(double));
   cur_halo->prof.ovdens  = (double *)       calloc(nbins, sizeof(double));
   cur_halo->prof.dens    = (double *)       calloc(nbins, sizeof(double));
   cur_halo->prof.v2_circ = (double *)       calloc(nbins, sizeof(double));
   cur_halo->prof.v_esc2  = (double *)       calloc(nbins, sizeof(double));
   cur_halo->prof.sig_v   = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.Ekin    = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.Epot    = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.Lx      = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.Ly      = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.Lz      = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.axis1   = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E1x     = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E1y     = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E1z     = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.axis2   = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E2x     = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E2y     = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E2z     = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.axis3   = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E3x     = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E3y     = (double *)       calloc(nbins, sizeof(double)); 
   cur_halo->prof.E3z     = (double *)       calloc(nbins, sizeof(double)); 
#ifdef AHFphspdens
#ifdef AHFmeanvelocities
   cur_halo->prof.sigma2_vx_sh = malloc(nbins*18*sizeof(double));
#else
   cur_halo->prof.sigma2_vx_sh = malloc(nbins*6*sizeof(double));
#endif
   cur_halo->prof.sigma2_vy_sh = cur_halo->prof.sigma2_vx_sh + 1 *nbins;
   cur_halo->prof.sigma2_vz_sh = cur_halo->prof.sigma2_vx_sh + 2 *nbins;
   cur_halo->prof.sigma2_vr_sh = cur_halo->prof.sigma2_vx_sh + 3 *nbins;
   cur_halo->prof.sigma2_vtheta_sh = cur_halo->prof.sigma2_vx_sh + 4 *nbins;
   cur_halo->prof.sigma2_vphi_sh = cur_halo->prof.sigma2_vx_sh + 5 *nbins;
#ifdef AHFmeanvelocities
   cur_halo->prof.mean_vx_sh = cur_halo->prof.sigma2_vx_sh + 6 *nbins;
   cur_halo->prof.mean_vy_sh = cur_halo->prof.sigma2_vx_sh + 7 *nbins;
   cur_halo->prof.mean_vz_sh = cur_halo->prof.sigma2_vx_sh + 8 *nbins;
   cur_halo->prof.mean_vr_sh = cur_halo->prof.sigma2_vx_sh + 9 *nbins;
   cur_halo->prof.mean_vtheta_sh = cur_halo->prof.sigma2_vx_sh + 10*nbins;
   cur_halo->prof.mean_vphi_sh = cur_halo->prof.sigma2_vx_sh + 11*nbins;
   cur_halo->prof.mean_vx_sp = cur_halo->prof.sigma2_vx_sh + 12*nbins;
   cur_halo->prof.mean_vy_sp = cur_halo->prof.sigma2_vx_sh + 13*nbins;
   cur_halo->prof.mean_vz_sp = cur_halo->prof.sigma2_vx_sh + 14*nbins;
   cur_halo->prof.mean_vr_sp = cur_halo->prof.sigma2_vx_sh + 15*nbins;
   cur_halo->prof.mean_vtheta_sp = cur_halo->prof.sigma2_vx_sh + 16*nbins;
   cur_halo->prof.mean_vphi_sp = cur_halo->prof.sigma2_vx_sh + 17*nbins;
#endif
#endif
#ifdef GAS_PARTICLES
  cur_halo->prof.M_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.M_star  = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.u_gas   = (double*) calloc(nbins, sizeof(double));
#ifdef AHFdisks
  cur_halo->prof.k_gas     = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.Ekin_gas  = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.Lx_gas    = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.Ly_gas    = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.Lz_gas    = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.axis1_gas = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E1x_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E1y_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E1z_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.axis2_gas = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E2x_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E2y_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E2z_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.axis3_gas = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E3x_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E3y_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E3z_gas   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.k_star     = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.Ekin_star  = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.Lx_star    = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.Ly_star    = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.Lz_star    = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.axis1_star = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E1x_star   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E1y_star   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E1z_star   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.axis2_star = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E2x_star   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E2y_star   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E2z_star   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.axis3_star = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E3x_star   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E3y_star   = (double*) calloc(nbins, sizeof(double));
  cur_halo->prof.E3z_star   = (double*) calloc(nbins, sizeof(double));
#endif /* AHFdisks */
#endif /* GAS_PARTICLES */
#ifdef METALHACK
  cur_halo->prof.z_gas   = malloc(nbins*2*sizeof(double));
  cur_halo->prof.z_star  = cur_halo->prof.z_gas + nbins;
#endif
}

/*----------------------------------------------------------------------------
* dest_profile:
* free all pointers associated with halo profile
*----------------------------------------------------------------------------*/
void dest_profile(HALO *cur_halo)
{
  if(cur_halo->prof.r)
   {
    free(cur_halo->prof.r);
    free(cur_halo->prof.npart);
    free(cur_halo->prof.nvpart);
    free(cur_halo->prof.ovdens);
    free(cur_halo->prof.dens);
    free(cur_halo->prof.v2_circ);
    free(cur_halo->prof.v_esc2);
    free(cur_halo->prof.sig_v);
    free(cur_halo->prof.Ekin);
    free(cur_halo->prof.Epot);
    free(cur_halo->prof.Lx);
    free(cur_halo->prof.Ly);
    free(cur_halo->prof.Lz);
    free(cur_halo->prof.axis1);
    free(cur_halo->prof.E1x);
    free(cur_halo->prof.E1y);
    free(cur_halo->prof.E1z);
    free(cur_halo->prof.axis2);
    free(cur_halo->prof.E2x);
    free(cur_halo->prof.E2y);
    free(cur_halo->prof.E2z);
    free(cur_halo->prof.axis3);
    free(cur_halo->prof.E3x);
    free(cur_halo->prof.E3y);
    free(cur_halo->prof.E3z);
#ifdef AHFphspdens
    free(cur_halo->prof.sigma2_vx_sh);
    free(cur_halo->prof.sigma2_vy_sh);
    free(cur_halo->prof.sigma2_vz_sh);
    free(cur_halo->prof.sigma2_vr_sh);
    free(cur_halo->prof.sigma2_vtheta_sh);
    free(cur_halo->prof.sigma2_vphi_sh);
#endif
#ifdef AHFmeanvelocities
    free(cur_halo->prof.mean_vx_sh);
    free(cur_halo->prof.mean_vy_sh);
    free(cur_halo->prof.mean_vz_sh);
    free(cur_halo->prof.mean_vr_sh);
    free(cur_halo->prof.mean_vtheta_sh);
    free(cur_halo->prof.mean_vphi_sh);
    free(cur_halo->prof.mean_vx_sp);
    free(cur_halo->prof.mean_vy_sp);
    free(cur_halo->prof.mean_vz_sp);
    free(cur_halo->prof.mean_vr_sp);
    free(cur_halo->prof.mean_vtheta_sp);
    free(cur_halo->prof.mean_vphi_sp);
#endif
#ifdef GAS_PARTICLES
    free(cur_halo->prof.M_gas);
    free(cur_halo->prof.M_star);
    free(cur_halo->prof.u_gas);
#ifdef AHFdisks
    free(cur_halo->prof.Ekin_gas);
    free(cur_halo->prof.Lx_gas);
    free(cur_halo->prof.Ly_gas);
    free(cur_halo->prof.Lz_gas);
    free(cur_halo->prof.axis1_gas);
    free(cur_halo->prof.E1x_gas);
    free(cur_halo->prof.E1y_gas);
    free(cur_halo->prof.E1z_gas);
    free(cur_halo->prof.axis2_gas);
    free(cur_halo->prof.E2x_gas);
    free(cur_halo->prof.E2y_gas);
    free(cur_halo->prof.E2z_gas);
    free(cur_halo->prof.axis3_gas);
    free(cur_halo->prof.E3x_gas);
    free(cur_halo->prof.E3y_gas);
    free(cur_halo->prof.E3z_gas);
    free(cur_halo->prof.Ekin_star);
    free(cur_halo->prof.Lx_star);
    free(cur_halo->prof.Ly_star);
    free(cur_halo->prof.Lz_star);
    free(cur_halo->prof.axis1_star);
    free(cur_halo->prof.E1x_star);
    free(cur_halo->prof.E1y_star);
    free(cur_halo->prof.E1z_star);
    free(cur_halo->prof.axis2_star);
    free(cur_halo->prof.E2x_star);
    free(cur_halo->prof.E2y_star);
    free(cur_halo->prof.E2z_star);
    free(cur_halo->prof.axis3_star);
    free(cur_halo->prof.E3x_star);
    free(cur_halo->prof.E3y_star);
    free(cur_halo->prof.E3z_star);
#endif
#endif
#ifdef METALHACK
    free(cur_halo->prof.z_gas);
    //free(cur_halo->prof.z_star); DO NOT FREE THIS AS THIS IS NOT MALLOC'ED! see c_profile() above
#endif
   }
   
   cur_halo->prof.nbins   = 0;
   cur_halo->prof.r       = NULL;
   cur_halo->prof.npart   = NULL;
   cur_halo->prof.nvpart  = NULL;
   cur_halo->prof.ovdens  = NULL;
   cur_halo->prof.dens    = NULL;
   cur_halo->prof.v2_circ = NULL;
   cur_halo->prof.v_esc2  = NULL;
   cur_halo->prof.sig_v   = NULL;
   cur_halo->prof.Ekin    = NULL;
   cur_halo->prof.Epot    = NULL;
   cur_halo->prof.Lx      = NULL;
   cur_halo->prof.Ly      = NULL;
   cur_halo->prof.Lz      = NULL;
   cur_halo->prof.axis1   = NULL;
   cur_halo->prof.E1x     = NULL;
   cur_halo->prof.E1y     = NULL;
   cur_halo->prof.E1z     = NULL;
   cur_halo->prof.axis2   = NULL;
   cur_halo->prof.E2x     = NULL;
   cur_halo->prof.E2y     = NULL;
   cur_halo->prof.E2z     = NULL;
   cur_halo->prof.axis3   = NULL;
   cur_halo->prof.E3x     = NULL;
   cur_halo->prof.E3y     = NULL;
   cur_halo->prof.E3z     = NULL;
#ifdef AHFphspdens
   cur_halo->prof.sigma2_vx_sh = NULL;
   cur_halo->prof.sigma2_vy_sh = NULL;
   cur_halo->prof.sigma2_vz_sh = NULL;
   cur_halo->prof.sigma2_vr_sh = NULL;
   cur_halo->prof.sigma2_vtheta_sh = NULL;
   cur_halo->prof.sigma2_vphi_sh = NULL;
#ifdef AHFmeanvelocities
   cur_halo->prof.mean_vx_sh = NULL;
   cur_halo->prof.mean_vy_sh = NULL;
   cur_halo->prof.mean_vz_sh = NULL;
   cur_halo->prof.mean_vr_sh = NULL;
   cur_halo->prof.mean_vtheta_sh = NULL;
   cur_halo->prof.mean_vphi_sh = NULL;
   cur_halo->prof.mean_vx_sp = NULL;
   cur_halo->prof.mean_vy_sp = NULL;
   cur_halo->prof.mean_vz_sp = NULL;
   cur_halo->prof.mean_vr_sp = NULL;
   cur_halo->prof.mean_vtheta_sp = NULL;
   cur_halo->prof.mean_vphi_sp = NULL;
#endif
#endif
#ifdef GAS_PARTICLES
  cur_halo->prof.M_gas   = NULL;
  cur_halo->prof.M_star  = NULL; 
  cur_halo->prof.u_gas  = NULL; 
#endif
#ifdef METALHACK
  cur_halo->prof.z_gas  = NULL;
  cur_halo->prof.z_star = NULL;
#endif
   
}
#endif
