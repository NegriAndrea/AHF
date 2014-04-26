#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

#ifdef DEBUG_GRIDS
#include "../libio_serial/io_serial.h"
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

static long dummy_variable;
void dummy_ahf_gridinfo(void)
{
}

#ifdef AHF

#include "ahf.h"

SRINDEX *spatialRefIndex;
int     *numIsoRef;          /* Number of spatially isolated refinements */
int      totnumIsoRef;       /* Number of spatially isolated refinements */


/*********************************************************************
 *********************************************************************
 * Identify the spatially connected refinements colouring the nodes that are connected
 * Dump the node - positions - density - spatial tag
 */


void ahf_gridinfo(gridls *grid_list, int curgrid_no)
{

  gridls	  *cur_grid;
  pqptr     cur_pquad;
  cqptr     cur_cquad, icur_cquad;
  nqptr     cur_nquad, icur_nquad;
  nptr      cur_node;

  int 			refinecounter;
	
  long			x,y,z;

  /* The colouring tools */
  int				colCounter=0;
  int				colREPLACE;
  int				currentCol;

  nptr   		tsc_nodes[3][3][3]; /* nodes to assign to */

  int              firstCOLOUR;

  SPATIALREF      *spatialRefHead;
  SPATIALREF		*tmpSpatialRef;
  SPATIALREF 		*current;

#ifdef AHFgridinfofile
  FILE          *gridinfofile;
  char          filename1[100];
#endif

  int i,tmp;
  int boundCount;
  int numWrongNodes, colRight;

#ifdef AHFDEBUG
  char          filename2[100];
  FILE          *fout;
#endif

  intXYZ tmpPeriodic;
  int		colour;
  SPATIALREF 		*previous;

  int iterate, startAHFgrid;
  double a, a3, omega, ovlim, rho_vir;
  double fl1dim, refine_len, refine_vol, refine_ovdens, refine_mass;

  double  cur_shift;
   
#ifdef VERBOSE
  /***************************************************************************/
  fprintf(stderr,"################## ahf_gridinfo ###################\n");
  fprintf(stderr,"Number of grids           = %d\n", curgrid_no+1);
#endif
  fprintf(io.logfile,"################## ahf_gridinfo ###################\n");
  fprintf(io.logfile,"Number of grids           = %d\n", curgrid_no+1);
  fflush(io.logfile);
  
  /* total number of grids to consider for AHF */
  /*  ahf.no_grids = curgrid_no - global.domgrid_no + 1;*/
  ahf.no_grids = curgrid_no - global.domgrid_no;
  
#ifdef VERBOSE
  fprintf(stderr,"Number of refinements     = %d\n", ahf.no_grids);
  fprintf(stderr,"global.domgrid_no         = %d\n", global.domgrid_no);
#endif
  fprintf(io.logfile,"Number of refinements     = %d\n", ahf.no_grids);
  fprintf(io.logfile,"global.domgrid_no         = %d\n", global.domgrid_no);
  
  
  
#ifdef DEBUG_GRIDS
  for(cur_grid=global.dom_grid+1; cur_grid<=global.dom_grid+ahf.no_grids; cur_grid++)
   {
    fprintf(stderr,"writing l1dim = %ld ... ",cur_grid->l1dim);
    //output_grid(cur_grid,0);
    write_density(cur_grid,"cur_grid-");
    fprintf(stderr,"done\n");
   }
  
  exit(0);
#endif

  
  
	/* cosmology related stuff */
	a        = global.a;
	a3       = pow3(a);
	omega    = calc_omega(a);
	ovlim    = calc_virial(a);
  rho_vir  = a3 * calc_rho_vir(a); /* comoving(!) density used to normalize densities */
  
  /* get the number of the grid satisfying virial overdensity criterion */
  refine_mass  = simu.pmass*simu.med_weight;
  cur_grid     = global.dom_grid;
  startAHFgrid = 1;
  for ( i=0; i<=ahf.no_grids; i++) 
    {
					
     fl1dim = ((double)(cur_grid->l1dim));
     cur_grid++;
     
     refine_len    = (double)(simu.boxsize/((double)(fl1dim)));
     refine_vol    = pow3(refine_len);
     refine_ovdens = (simu.Nth_ref*refine_mass/refine_vol) / rho_vir;
#ifdef VERBOSE
     fprintf(stderr,"l1dim = %16.0f refine_ovdens = %16.4g ovlim = %16.4g\n", fl1dim,refine_ovdens,ovlim);
#endif
     fprintf(io.logfile,"l1dim = %16.0f refine_ovdens = %16.4g ovlim = %16.4g\n", fl1dim,refine_ovdens,ovlim);
     if ( refine_ovdens <  ovlim ) 
        startAHFgrid = i;     
    }
  
  /* store maximum overdensity possible with current refinement hierarchy (used by ahf_halos()!) */
  global.max_ovdens = refine_ovdens;
#ifdef VERBOSE
  fprintf(stderr,"max_ovdens = %g\n",global.max_ovdens);
#endif
  fprintf(io.logfile,"max_ovdens = %g\n",global.max_ovdens);
  
  /* the first refinement grid to consider */
  ahf.min_ref = startAHFgrid+AHF_MIN_REF_OFFSET;
  
  /* total number of grids to be used with AHF */
  ahf.no_grids = ahf.no_grids - ahf.min_ref + 1; 
  
#ifdef VERBOSE
  fprintf(stderr,"min_ref = %d    (ahf_nogrids = %d)\n", ahf.min_ref,ahf.no_grids);
#endif
  fprintf(io.logfile,"min_ref = %d    (ahf_nogrids = %d)\n", ahf.min_ref,ahf.no_grids);
  fflush(io.logfile);
  
  /* Initialising important variables */
#ifdef VERBOSE
  fprintf(stderr,"  initialising force.color on all grids...\n");
#endif
  spatialRefHead = NULL; 	/* Setting the HEAD of the link list of SPATIALREF's */
  spatialRefTail = NULL;
  colREPLACE     = FALSE;
  firstCOLOUR    = 1;

#ifdef WITH_OPENMP
  //
  // this loop over all grids and all refinements generates a linked-list starting
  // with spatialRefHead that contains the spatially connected refinements independent
  // of the actual grid (that information is contained in spatialRefHead->refLevel)
  //
  // to parallelize one could do this generation independently on each grid,
  // but would need to merge all linked-lists afterwards
  // 
  //#pragma omp parallel private(cur_grid, ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, firstCOLOUR, spatialRefHead, spatialRefTail, cur_node) shared(refinecounter)
  //#pragma omp for schedule(static)
#endif
  for( refinecounter=0; refinecounter<ahf.no_grids; refinecounter++ ) 
    {
      cur_grid = global.dom_grid+ahf.min_ref+refinecounter;
      
#ifdef VERBOSE
     /* fprintf() to double-check that we set min_ref and ahf.no_grids correctly! */
     fprintf(stderr,"%12d || l1dim = %12ld\n", refinecounter, cur_grid->l1dim);
#endif
					
     /************************************************************
     * Initialising all the colour tages in preperation for the 
     * Spatial refinement identification.
     * Therefore, we are zeroing all the colour tags on the current refinement
     ************************************************************/
     for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
      {
        for(cur_cquad = cur_pquad->loc;
            cur_cquad < cur_pquad->loc + cur_pquad->length; 
            cur_cquad++)  
          {  
           for(icur_cquad  = cur_cquad; 
               icur_cquad != NULL; 
               icur_cquad  = icur_cquad->next)
             {
              for(cur_nquad = icur_cquad->loc;  
                  cur_nquad < icur_cquad->loc + icur_cquad->length; 
                  cur_nquad++) 
                { 
                 for(icur_nquad  = cur_nquad; 
                     icur_nquad != NULL; 
                     icur_nquad  = icur_nquad->next)
                   {
                    for(cur_node = icur_nquad->loc; 
                        cur_node < icur_nquad->loc + icur_nquad->length; 
                        cur_node++)
                      {
                       /* zeroing the colour */
                       cur_node->force.colour = 0;
                      }
                   }
                }
             }
          }
       } 
     
     /************************************************************
     * Start looping through the refinement identifying the spatial refinements
     ************************************************************/
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
                       
                       
                       /* The current colour of the node */
                       currentCol = cur_node->force.colour;
                       
                       /* We have finished updating the colours */
                       if ( (colREPLACE==TRUE) && (currentCol==0) )
                          colREPLACE=FALSE;
                       
                       /* Replacing the colours */
                       if ( colREPLACE ) 
                         { 		
                          
                          if ( col.numReplace == 2 ) 
                            {
                             
                             if (currentCol==col.holder[2])
                                cur_node->force.colour=col.holder[1];
                             
                            } 
                          else if ( col.numReplace == 3 ) 
                            {
                             
                             if (currentCol==col.holder[2])
                                cur_node->force.colour=col.holder[0];
                             
                            } 
                          else if ( col.numReplace == 4 ) 
                            {
                             
                             if (currentCol==col.holder[2])
                                cur_node->force.colour=col.holder[0];
                             
                             if (currentCol==col.holder[1])
                                cur_node->force.colour=col.holder[0];
                             
                            } 
                          else 
                            {
                             fprintf(stderr,"Something is wrong with colREPLACE");
                            }
                          
                         } 
                       /* New Node to check */
                       else 
                         { 					 
                          
                          /* Getting the relevant nodes and filling colHolder 	*/
                          /* search for pointers if any are null return NULL 	*/
                          tsc_nodes[1][1][1] = cur_node;
                          get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                          
                          /* Gathering the colour information for this node 		*/
                          colourInfo(tsc_nodes);
                          
                          /************************************************/
                          /* all nodes are zero, numReplace=0 */
                          if ( col.numReplace == 0 )
                            {
                             /* Start a new colour and keep going */
                             if ( (tmpSpatialRef = malloc(sizeof(SPATIALREF))) == NULL ) 
                               {
                                fprintf(stderr,"No memory for SPATIALREF\n");
                                exit(1);
                               }
                             
                             /* Colour the node with a unique name */
                             colCounter++;
                             
                             cur_node->force.colour = colCounter;
                             
                             /* Point SPATIALREF to the right information */
                             tmpSpatialRef->name 			  = colCounter;
                             tmpSpatialRef->periodic.x	= 0;
                             tmpSpatialRef->periodic.y	= 0;
                             tmpSpatialRef->periodic.z	= 0;
                             tmpSpatialRef->refLevel		= refinecounter;
                             tmpSpatialRef->cur_pquad	  = cur_pquad;
                             tmpSpatialRef->cur_cquad	  = cur_cquad;
                             tmpSpatialRef->icur_cquad	= icur_cquad;
                             tmpSpatialRef->cur_nquad	  = cur_nquad;
                             tmpSpatialRef->icur_nquad	= icur_nquad;
                             tmpSpatialRef->cur_node		= cur_node;
                             tmpSpatialRef->x				    = x;
                             tmpSpatialRef->y				    = y;
                             tmpSpatialRef->z				    = z;
                             
                             /* Join the new  SPATIALREF to the link list */
                             spatialRefTail = insertColour(tmpSpatialRef, firstCOLOUR);
                             
                             /* Now just for the first time we need to make head also point to this newColour */
                             if ( firstCOLOUR==1 ) {
                                spatialRefHead = spatialRefTail; 
                                firstCOLOUR=0;
                             }
                             
                            }
                          
                          /************************************************/
                          /* One node is non-zero, numReplace=1 */
                          if ( col.numReplace == 1 ) 
                            {
                             
                             /* set using the current colour and keep going */
                             cur_node->force.colour = col.holder[2];
                             
                            }
                          /************************************************/
                          /* Two nodes are non-zero and are not equal, numReplace=2 */
                          if ( col.numReplace == 2 ) 
                            {
                             
                             /* Replace the previous colours with the lowest colour */
                             /* set colour using the lowest value : in this case 1 */
                             cur_node->force.colour = col.holder[1];
                             
                             if (!colREPLACE) 
                               { 
                                /* Find and destroy the redundant SPATIALREF */
                                tmpSpatialRef=deleteColour(spatialRefHead, col.holder[2]);
                                
                                /* STU :: This is where you start looking for the deleted colours */
                                cur_pquad	   = tmpSpatialRef->cur_pquad;
                                cur_cquad	   = tmpSpatialRef->cur_cquad;
                                icur_cquad	 = tmpSpatialRef->icur_cquad;
                                cur_nquad	   = tmpSpatialRef->cur_nquad;
                                icur_nquad	 = tmpSpatialRef->icur_nquad;
                                cur_node 	   = tmpSpatialRef->cur_node;
                                x            = tmpSpatialRef-> x;
                                y 				   = tmpSpatialRef-> y;
                                z 				   = tmpSpatialRef-> z;
                                
                                free(tmpSpatialRef);
                                tmpSpatialRef = NULL;
                                 
                                colREPLACE=TRUE;
                                
                                /* When going back to replace the colour the for loop progress
                                   the node by one - which we don't want to happen 
                                   Hence, we have to do an inital sweep on the first node with  
                                   the colour that is going to be deleted*/
                                /* Getting the relevant nodes and filling colHolder 	*/
                                /* search for pointers if any are null return NULL 	*/
                                currentCol         = cur_node->force.colour;
                                tsc_nodes[1][1][1] = cur_node;
                                get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                                /* When going back to replace the colour the for loop progress
                                   the node by one - which we don't want to happen */
                                if ( col.numReplace == 2 ) 
                                  {
                                   
                                   if (currentCol==col.holder[2])
                                      cur_node->force.colour=col.holder[1];
                                   
                                  } 
                                else if ( col.numReplace == 3 ) 
                                  {
                                   
                                   if (currentCol==col.holder[2])
                                      cur_node->force.colour=col.holder[0];
                                   
                                  } 
                                else if ( col.numReplace == 4 ) 
                                  {
                                   
                                   if (currentCol==col.holder[2])
                                      cur_node->force.colour=col.holder[0];
                                   
                                   if (currentCol==col.holder[1])
                                      cur_node->force.colour=col.holder[0];
                                   
                                  } 
                                else 
                                  {
                                   fprintf(stderr,"Something is wrong with colREPLACE");
                                  }
                                /*****************************************************************/
                                /*****************************************************************/
                                
                               }
                             
                            }
                          
                          /************************************************/
                          /* Two nodes are non-zero and are not equal, numReplace=2 */
                          /* With the last two colour elements equal */
                          if ( col.numReplace == 3 ) 
                            {
                             
                             /* Replace the previous colours with the lowest colour */
                             /* set colour using the lowest value : in this case 1 */
                             cur_node->force.colour = col.holder[0];
                             
                             if (!colREPLACE) 
                               { 
                                /* Find and destroy the redendent SPATIALREF */
                                tmpSpatialRef=deleteColour(spatialRefHead, col.holder[2]);
                                
                                cur_pquad	   = tmpSpatialRef->cur_pquad;
                                cur_cquad	   = tmpSpatialRef->cur_cquad;
                                icur_cquad	 = tmpSpatialRef->icur_cquad;
                                cur_nquad	   = tmpSpatialRef->cur_nquad;
                                icur_nquad	 = tmpSpatialRef->icur_nquad;
                                cur_node 	   = tmpSpatialRef->cur_node;
                                x 				   = tmpSpatialRef-> x;
                                y 				   = tmpSpatialRef-> y;
                                z 				   = tmpSpatialRef-> z;
                                
                                free(tmpSpatialRef);
                                tmpSpatialRef = NULL;
                                colREPLACE=TRUE;
                                
                                
                                /*****************************************************************/
                                /* Getting the relevant nodes and filling colHolder 	*/
                                /* search for pointers if any are null return NULL 	*/
                                currentCol         = cur_node->force.colour;
                                tsc_nodes[1][1][1] = cur_node;
                                get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                                /* When going back to replace the colour the for loop progress
                                   the node by one - which we don't want to happen */
                                if ( col.numReplace == 2 ) 
                                  {
                                   
                                   if (currentCol==col.holder[2])
                                      cur_node->force.colour=col.holder[1];
                                   
                                  } 
                                else if ( col.numReplace == 3 ) 
                                  {
                                   
                                   if (currentCol==col.holder[2])
                                      cur_node->force.colour=col.holder[0];
                                   
                                  } 
                                else if ( col.numReplace == 4 ) 
                                  {
                                   
                                   if (currentCol==col.holder[2])
                                      cur_node->force.colour=col.holder[0];
                                   
                                   if (currentCol==col.holder[1])
                                      cur_node->force.colour=col.holder[0];
                                   
                                  } 
                                else 
                                  {
                                   fprintf(stderr,"Something is wrong with colREPLACE");
                                  }
                                /*****************************************************************/
                                /*****************************************************************/
                                
                               }
                             
                            }
                          
                          /************************************************/
                          /* Three node are non-zero and are not equal, numReplace=3*/
                          if ( col.numReplace == 4 ) 
                            {
                             
                             /* Loop back and replace the dupulacates */
                             /* set colour using the lowest value : in this case 0 */
                             cur_node->force.colour = col.holder[0];
                             
                             if (!colREPLACE) 
                               {
                                /* Find and destroy the redendent SPATIALREF's */
                                tmpSpatialRef=deleteColour(spatialRefHead, col.holder[2]);
                                free(tmpSpatialRef);
                                tmpSpatialRef=deleteColour(spatialRefHead, col.holder[1]);
                                
                                cur_pquad	   = tmpSpatialRef->cur_pquad;
                                cur_cquad	   = tmpSpatialRef->cur_cquad;
                                icur_cquad	 = tmpSpatialRef->icur_cquad;
                                cur_nquad	   = tmpSpatialRef->cur_nquad;
                                icur_nquad	 = tmpSpatialRef->icur_nquad;
                                cur_node 	   = tmpSpatialRef->cur_node;
                                x 				   = tmpSpatialRef-> x;
                                y 				   = tmpSpatialRef-> y;
                                z 				   = tmpSpatialRef-> z;
                                
                                free(tmpSpatialRef);
                                colREPLACE=TRUE;
                                
                                
                                /*****************************************************************/
                                /*****************************************************************/
                                /*****************************************************************/
                                /*****************************************************************/
                                /* Getting the relevant nodes and filling colHolder 	*/
                                /* search for pointers if any are null return NULL 	*/
                                currentCol         = cur_node->force.colour;
                                tsc_nodes[1][1][1] = cur_node;
                                get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                                /* When going back to replace the colour the for loop progress
                                   the node by one - which we don't want to happen */
                                if ( col.numReplace == 2 )
                                  {
                                   
                                   if (currentCol==col.holder[2])
                                      cur_node->force.colour=col.holder[1];
                                   
                                  } 
                                else if ( col.numReplace == 4 ) 
                                  {
                                   
                                   if (currentCol==col.holder[2])
                                      cur_node->force.colour=col.holder[0];
                                   
                                   if (currentCol==col.holder[1])
                                      cur_node->force.colour=col.holder[0];
                                   
                                  } 
                                else 
                                  {
                                   fprintf(stderr,"Something is wrong with colREPLACE");
                                  }
                                /*****************************************************************/
                                /*****************************************************************/
                                
                               }	
                             
                            }
                          
                          /* ERROR */
                          if ( col.numReplace > 4 )
                             fprintf(stderr,"ERROR with counting col.numReplace");
                          
                         }
                       
                      }
                   }
                }
             }
          }
       } 
     
     
    }
	
	
  /*******************************************************************************************************/
  /* Checking that the grid is coloured correctly */
  /*******************************************************************************************************/
#ifdef AHFDEBUG2
  /* Printing the PVIEW file */	
  sprintf(filename2,"%shalo_pview3",global_io.params->outfile_prefix);
  fprintf(stderr,"%s\n",filename2);
  if((fout = fopen(filename2,"w")) == NULL) 
    {
     fprintf(stderr,"could not open %s\n", filename2);
     exit(1);
    }
#endif


  refinecounter=0;
  for( 	iterate=0, cur_grid=global.dom_grid+ahf.min_ref;  iterate<ahf.no_grids; iterate++, cur_grid++) 
    {
     /* shift of cell centre as compared to edge of box [grid units] */
     cur_shift = 0.5/(double)cur_grid->l1dim;

     /************************************************************
     * Dealing with the extents of the periodic refinements
     ************************************************************/
     numWrongNodes = 0;
     boundCount    = 0;
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
                      { /* current node */
                       
                       colour = cur_node->force.colour;
                       
                       tsc_nodes[1][1][1] = cur_node;
                       get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                       
                       colRight = 0;
                       colRight = checkColourInfo(tsc_nodes);
                       numWrongNodes += colRight;
                       
#ifdef AHFDEBUG2
                       if ( colRight == 1 )
                          fprintf(fout,"p %g %g %g 0.0 1.0 0.5\n",
                                  ((float) x/(float)cur_grid->l1dim + cur_shift)*simu.boxsize,
                                  ((float) y/(float)cur_grid->l1dim + cur_shift)*simu.boxsize,
                                  ((float) z/(float)cur_grid->l1dim + cur_shift)*simu.boxsize);
#endif
                       
                       /********************************************************************************/	
                       /* PERIODIC BOUNDARY CONDITIONS */
                       /* We are on the 'bottom' boundary */
                        if ( (x == 0) || (y == 0) || (z == 0) ) {
                          
                          /* Gathering periodic information */ 
                          tmpPeriodic = testBound(tsc_nodes,x,y,z);
                          
                          /* This isolated refinement is also periodic */
                          if ( (tmpPeriodic.x == 1) || (tmpPeriodic.y == 1) || (tmpPeriodic.z == 1) ) {
                             
                             /* finding the corresponding isolated spatial refinement in the spatalREF linklist */
                             current=spatialRefHead;
                             while ( current != NULL ) {
                                
                                previous = current;
                                
                                if ( current->name == colour )
                                   break;
                                
                                current = current->next;
                             } 
                             
                             /* Setting the periodic elements */
                             if ( tmpPeriodic.x == 1 )
                                previous->periodic.x = 1;
                             
                             if ( tmpPeriodic.y == 1 )
                                previous->periodic.y = 1;
                             
                             if ( tmpPeriodic.z == 1 )
                                previous->periodic.z = 1;
                             
                          }
                          
                          boundCount++;
                       }
                       
                      }/* current node */
                   }
                }
             }
          }
       }
#ifdef AHFDEBUG
     fprintf(stderr,"Num Boundary (%d) = %d\n",refinecounter,boundCount);
#endif
     refinecounter++;
    }
#ifdef AHFDEBUG2
  fclose(fout);
#endif

								

  /****************************************************************
   ****************************************************************
   * Creating the spatial Refinement Index
   *
   * And boundary condition information!!!
   * 
   */

  spatialRefIndex = NULL;
  if ((spatialRefIndex = calloc(colCounter+1,sizeof(SRINDEX)))==NULL) 
    {
    fprintf(stderr,"Error in allocating the memory for spatialRefIndex array\n");
    exit(0);
  }	
#ifdef VERBOSE
  fprintf(stderr,"colCounter(%d)\n",colCounter);
#endif

  /* Creating a tempary counting array */
  numIsoRef=NULL;
  if ((numIsoRef = calloc(ahf.no_grids,sizeof(int)))==NULL) 
    {
     fprintf(stderr,"Error in allocating the memory for numIsoRef array\n");
     exit(0);
    }
  for ( i=0; i<ahf.no_grids; i++)
    numIsoRef[i]=0;
 
  /* Setting the spatialRefIndex to -1 */
  for ( i=0; i<colCounter; i++) 
    {
     spatialRefIndex[i].refLevel = -1;
     spatialRefIndex[i].isoRefIndex = -1;
    }
 

  /* Loop through the spatial refinements */
  /* Record their refinement level and isolated refinement index */
  totnumIsoRef=0;
  current=spatialRefHead;
  while ( current!=NULL ) 
    {
     
     spatialRefIndex[current->name].refLevel    = current->refLevel;
     spatialRefIndex[current->name].isoRefIndex = numIsoRef[current->refLevel];
     spatialRefIndex[current->name].periodic.x  = current->periodic.x;
     spatialRefIndex[current->name].periodic.y  = current->periodic.y;
     spatialRefIndex[current->name].periodic.z  = current->periodic.z;
     
     numIsoRef[current->refLevel]++;
     totnumIsoRef++;
     
     current 	= current->next;
    }

	
  /****************************************************************
   ****************************************************************
   * Looping throught the Grid structure and dumping the gridInfo file
   * positions - density - spatial tag
   */


#ifdef AHFgridinfofile
  fprintf(stderr,"### PRINTING GRIDINFO FILES\n");
  /* Note:: do this for each refinement level starting with the first refinment */
  for( 	iterate=0, cur_grid=global.dom_grid+ahf.min_ref;  iterate<ahf.no_grids; iterate++, cur_grid++) 
    {
     
     /* prepare filename */
     sprintf(filename1,"%sz%.3f.grid.%06d",global_io.params->outfile_prefix,global.z,cur_grid->l1dim);
     fprintf(stderr,"%s\n",filename1);
     
     if((gridinfofile = fopen(filename1,"w")) == NULL) 
       {
        fprintf(stderr,"could not open %s\n", filename1);
        exit(1);
       }
     
     /* shift of cell centre as compared to edge of box [grid units] */
     cur_shift = 0.5/(double)cur_grid->l1dim;

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
                       fprintf(gridinfofile,"p %lf %lf %lf 0.0 0.0 1.0 %lf %d\n",  
                               simu.boxsize*1000*((double) x/(double)cur_grid->l1dim + cur_shift),
                               simu.boxsize*1000*((double) y/(double)cur_grid->l1dim + cur_shift),
                               simu.boxsize*1000*((double) z/(double)cur_grid->l1dim + cur_shift),
                               cur_node->dens, cur_node->force.colour);
                      }	
                   }
                }
             }
          }
       }
     
     /* Flusing and closing this grids file */
     fflush(gridinfofile);
     fclose(gridinfofile);
     
    }
#endif

  /*	fprintf(stderr,"I'm finished here, my work is done.... time to exit :)\n");
   *	exit(1); */
  fprintf(io.logfile,"################## ahf_gridinfo finished ##################\n\n");
  fflush(io.logfile);
  
  /* free() all spatialRefHead memory */
  current = spatialRefHead;
  while(current != NULL)
   {
    previous = current->next;
    free(current);
    current = previous;
   }
}



/*
************************************************************
************************************************************
*  Delete colour / spatial refinment from the link list
*/
SPATIALREF* deleteColour(SPATIALREF* spatialRefHead, int name)
{

  SPATIALREF *current;
  SPATIALREF *previous;
  SPATIALREF *tmpSpatialRef;
	
  /***********************************************/ 
  /* Looking for the SPATIALREF with the same name 			*/
  /* Now reset the current pointers to the lowest colour 	*/
  /* And destroy the redendent SPATIALREF 						*/
  /***********************************************/ 
  if ( (tmpSpatialRef = malloc(sizeof(SPATIALREF))) == NULL ) 
    {
     fprintf(stderr,"No memory for SPATIALREF\n");
     exit(1);
    }

  /* ERROR :: head and tail pointing to the same thing */
  if ( (spatialRefHead->name == name) && (spatialRefTail->name == name) ) 
    {
     
     fprintf(stderr,"ERROR :: head and tail never point to the same thing");
     exit(1);
     
    } /* ERROR :: Deleting the first element of the link list */ 
  else if ( spatialRefHead->name == name ) 
    {

     fprintf(stderr,"ERROR :: You never delete the first element of the link list");
     exit(1);
     
    } 
  /* Deleting the last element of the link list */
  else if ( spatialRefTail->name == name ) 
    {
     current  = spatialRefHead;
     previous = current;
     while ( current!=spatialRefTail ) 
       {
        previous = current;
        current  = current->next;
       }
     /*	fprintf(stderr,"HELP\n"); */
     previous->next = NULL;
     spatialRefTail = previous;
     
     /* Now reset the current pointers to the lowest colour */
     tmpSpatialRef->cur_pquad  = current->cur_pquad;
     tmpSpatialRef->cur_cquad  = current->cur_cquad;
     tmpSpatialRef->icur_cquad = current->icur_cquad;
     tmpSpatialRef->cur_nquad  = current->cur_nquad;
     tmpSpatialRef->icur_nquad = current->icur_nquad;
     tmpSpatialRef->cur_node   = current->cur_node;
     tmpSpatialRef->x          = current->x;
     tmpSpatialRef->y          = current->y;
     tmpSpatialRef->z          = current->z;
     
     /* Removing the pointer */
     free(current);
    }
  
  /* The colour in safely inside the link list */
  else 
    {
     /* start at the beginning */
     current = spatialRefHead;
     while ( current->name != name ) 
       { 
        previous = current;
        current  = current->next;
       }
     previous->next = current->next;
     
     /* Now reset the current pointers to the lowest colour */
     tmpSpatialRef->cur_pquad   = current->cur_pquad;
     tmpSpatialRef->cur_cquad   = current->cur_cquad;
     tmpSpatialRef->icur_cquad  = current->icur_cquad;
     tmpSpatialRef->cur_nquad   = current->cur_nquad;
     tmpSpatialRef->icur_nquad  = current->icur_nquad;
     tmpSpatialRef->cur_node    = current->cur_node;
     tmpSpatialRef->x           = current->x;
     tmpSpatialRef->y           = current->y;
     tmpSpatialRef->z           = current->z;
     
     /* Removing the pointer */
     free(current);
    }

  /* Returning the current pointers to the lowest colour to start the replace search */
  return(tmpSpatialRef);
}



/*
************************************************************
************************************************************
*  Inserting colour/spatial refinment into the link list
*/
SPATIALREF* insertColour( SPATIALREF* newColour, int firstCOLOUR )
{

  /* Point to NULL since it is the end element */
  newColour->next = NULL;


  /* This makes what tail points to, point to newColour 	*/
  /* Only do it when it is not the first time it is done 	*/
  if ( firstCOLOUR==0 )
    spatialRefTail->next = newColour;

  return(newColour);

}

/*
************************************************************
************************************************************
*  Checking that the nodes are coloured correctly 
*/
int	checkColourInfo(nptr tsc_nodes[3][3][3])
{

  int	i;
  int	r,g,b;
  int 	numzero;

  int	colRight=0;

  /* Looking back to fill the colour Holder	*/

  if ( tsc_nodes[1][0][1] != NULL ) {
    r = tsc_nodes[1][0][1]->force.colour;
  } else {
    r=0; /* If there is nothing there just set it to 0 */
  }

  if ( tsc_nodes[1][1][0] != NULL ) {
    g = tsc_nodes[1][1][0]->force.colour;
  } else {
    g=0;
  }

  if ( tsc_nodes[0][1][1] != NULL ) {
    b = tsc_nodes[0][1][1]->force.colour;
  } else {
    b=0;
  }

  colRight=0;
	
  col.holder[0] = r;
  col.holder[1] = g;
  col.holder[2] = b;

  /*	qsort(col.holder,NDIM,sizeof(int),(int (* )(void *, void *))colourCompare);*/
  qsort(col.holder,NDIM,sizeof(int),colourCompare);

  r = col.holder[0];
  g = col.holder[1];
  b = col.holder[2];

  numzero=0;
  for( i=0; i<NDIM; i++ ) {
    if ( col.holder[i] == 0 ) {
      numzero++;
    }
  }

  if ( numzero == 3 ) {
    col.numReplace = 0;
    colRight = 0;

  } else if ( numzero == 2 ) {

    col.numReplace = 1;
    colRight = 0;

  } else if ( numzero == 1 ) {

    if ( g==b ) {
      col.numReplace = 1;
      colRight = 0;
    } else {
      col.numReplace = 2;
      colRight = 1;
    }

  } else {

    if ( (r == g) && (r == b) ) { 					
      col.numReplace = 1;
      colRight = 0;
    } else if ( (r==g) && (r!=b) ) { 				
      col.numReplace = 2;
      colRight = 1;
    } else if ( (r==b) && (r!=g) ) {
      col.numReplace = 2;
      colRight = 1;
    } else if ( (g==b) && (g!=r) ) {
      col.numReplace = 3;
      colRight = 1;
    } else if ( (r!=g) && (r!=b) ) { 				
      col.numReplace = 4;
      colRight = 1;
    } else {
      fprintf(stderr,"# There should not be any other else!\n");
      colRight = 1;
    }

  }

	
  return(colRight);
	

}
	
/*
************************************************************
************************************************************
*  Looking in all 6 directions to see if there are different nodes
*/
intXYZ	testBound(nptr tsc_nodes[3][3][3], int x, int y, int z)
{
  
  intXYZ periodic;
  
  periodic.x=0;
  periodic.y=0;
  periodic.z=0;
	
  /* YYYYYYYYYYYYYYY */
  if ( y == 0 ) {
    if ( tsc_nodes[1][0][1] != NULL ) /* Looking back - and finding something there => periodic */
      periodic.y = 1;
  }
	
  /* XXXXXXXXXXXXXXXX */
  if ( x == 0 ) {
    if ( tsc_nodes[1][1][0] != NULL )
      periodic.x = 1;
  }
  
  /* ZZZZZZZZZZZZZZZ */
  if ( z == 0 ) {
    if ( tsc_nodes[0][1][1] != NULL )
      periodic.z = 1;
  }
  
  return(periodic);
}


/*
************************************************************
************************************************************
*  Gathering all the information we need to know about the colour of the node
*/
void	colourInfo(nptr tsc_nodes[3][3][3])
{

  int	i,j;
  int	r,g,b,t,h,n;
  int numzero;
  int tmpCol[6];
  /*	int tmpCol[3]; */
  int	NEW;

  /* Cleaning the colour holder */
  col.holder[0] = 0; col.holder[1] = 0; col.holder[2] = 0;
  col.numReplace = 0;
  col.numCol = 0;
//  free(col.info);
  col.info = NULL;
	

  /* Looking back to fill the colour Holder	*/

  /*********************************************************
   * Y :: node below me colour */
  if ( tsc_nodes[1][0][1] != NULL ) {
    r = tsc_nodes[1][0][1]->force.colour;
  } else {
    r=0; /* If there is nothing there just set it to 0 */
  }
  /*********************************************************
   * Y :: node above me colour */
  if ( tsc_nodes[1][2][1] != NULL ) {
    t = tsc_nodes[1][2][1]->force.colour;
  } else {
    t=0; /* If there is nothing there just set it to 0 */
  }


  /*********************************************************
   *  X :: right node colour */	
  if ( tsc_nodes[1][1][0] != NULL ) {
    g = tsc_nodes[1][1][0]->force.colour;
  } else {
    g=0;
  }
  /*********************************************************
   *  X :: left node colour */	
  if ( tsc_nodes[1][1][2] != NULL ) {
    h = tsc_nodes[1][1][2]->force.colour;
  } else {
    h=0;
  }

	
  /*********************************************************
   *  Z :: node to my back colour */
  if ( tsc_nodes[0][1][1] != NULL ) {
    b = tsc_nodes[0][1][1]->force.colour;
  } else {
    b=0;
  }
  /*********************************************************
   *  Z :: node to my front colour */
  if ( tsc_nodes[2][1][1] != NULL ) {
    n = tsc_nodes[2][1][1]->force.colour;
  } else {
    n=0;
  }

  /*********************************************************/
  /* list of unique colours */

  tmpCol[0] = r; tmpCol[1] = g; tmpCol[2] = b;
  tmpCol[3] = t; tmpCol[4] = h; tmpCol[5] = n;

  numzero=0;
  for (i=0;i<6;i++) {
    /*	for (i=0;i<3;i++) { */

    if ( tmpCol[i] != 0 ) { /* entre if not a zero */

      /* Is this the first colour we have come accross? */
      if ( col.numCol == 0 ) {

	if ((col.info = malloc(sizeof(COLGATH)))==NULL) {
	  fprintf(stderr,"malloc failed in allocating the col.info\n");
	  exit(-1);
	}
	col.info[0].name = tmpCol[i];
	col.info[0].num  = 1;
	col.numCol++;

      } else {

	/* have we seen this colour before? */
	NEW = 0; /* yes - we want a new colour */
	for (j=0;j<col.numCol;j++) {
								
	  if ( tmpCol[i] == col.info[j].name ) { /* It is an old colour - iterate the counters */
	    NEW = 1; /* We have seen this colour before */
	    col.info[j].num++;
	  }
	}

	/* It is a new colour - create a space for it*/
	if ( NEW == 0 ) {
								
	  col.numCol++;
	  if ((col.info = realloc(col.info,(col.numCol)*sizeof(COLGATH)))==NULL) {
	    fprintf(stderr,"realloc failed in allocating the col.info\n");
	    exit(-1);
	  }
	  col.info[col.numCol-1].name = tmpCol[i];
	  col.info[col.numCol-1].num  = 1;
	}
      }
    } else {
      numzero++;
    }
  }


  /*********************************************************/

#ifdef VERBOSE					
  if ( col.numCol > 3 ) {
    fprintf(stderr,"col.numCol\n");
    for ( i=0;i<col.numCol;i++ ) 
      fprintf(stderr,"%d ",col.info[i].name);		
    fprintf(stderr,"\n");
  }
#endif
	
  for ( i=0;i<col.numCol;i++ ) 
    col.holder[i] = col.info[i].name;

  /* Sorting the color holder from smallest to greatest */
  qsort(col.holder,NDIM,sizeof(int),colourCompare);

  r = col.holder[0];
  g = col.holder[1];
  b = col.holder[2];

  numzero=0;
  for( i=0; i<NDIM; i++ ) {
    if ( col.holder[i] == 0 ) {
      numzero++;
    }
  }
	
  /*********************************************************/


  if ( numzero == 3 ) {
    /* new colour */
    col.numReplace = 0;
  } else if ( numzero == 2 ) {
    /* set colour col.holder[2] */
    col.numReplace = 1;
  } else if ( numzero == 1 ) {

    if ( g==b ) {
      /* set colour col.holder[2] */
      col.numReplace = 1;
    } else {
      /* replace colour col.holder[2] with col.holder[1] */
      col.numReplace = 2;
    }

  } else {

    if ( (r == g) && (r == b) ) { 						/* All the same */
      /* set colour col.holder[2] */
      col.numReplace = 1;
    } else if ( (r==g) && (r!=b) ) { /* If two are the same */
      /* replace colour col.holder[2] with col.holder[1] */
      col.numReplace = 2;
    } else if ( (r==b) && (r!=g) ) {
      /* replace colour col.holder[2] with col.holder[1] */
      col.numReplace = 2;
    } else if ( (g==b) && (g!=r) ) {
      /* replace colour col.holder[2] with col.holder[0] */
      col.numReplace = 3;
    } else if ( (r!=g) && (r!=b) ) { 		/* If they are all diferent */ 
      /* replace colour col.holder[2] and col.holder[1] with col.holder[0] */
      col.numReplace = 4;
    } else {
      fprintf(stderr,"# There should not be any other else!\n");
    }

  }
  
  free(col.info);
}


/*
************************************************************
************************************************************
*  qsort stuff
*/
int colourCompare(const void *p1, const void *p2) {
		
  if ( *(int *)p1 < *(int *)p2 )
    return(-1);
  else if ( *(int *)p1 > *(int *)p2 )
    return(1);
  else
    return(0);
}

#endif /* AHF */

