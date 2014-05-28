#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "common.h"
#ifdef WITH_MPI
#	include <mpi.h>
#	include "libutility/loadbalance.h"
#	include "comm.h"
#endif

#ifdef AHF
#include "libahf/ahf.h"
#endif

#ifdef AHF2
#include "libahf2/ahf.h"
#endif

#include "libsfc/sfc.h"
#include "startrun.h"
#include "libutility/utility.h"
#include "libgravity/gravity.h"
//#include "libio_serial/io_serial.h"

#ifdef NEWAMR
#include "libtree/tree.h"
#else
#include "libamr_serial/amr_serial.h"
#endif


#ifdef WITH_MPI
static void
local_communicate_all(void);
#endif

#ifdef AHFrfocus
#include <assert.h>
static void
local_focusSphere(void);
#endif

/*==============================================================================
 * MAIN: where everything starts ....
 *==============================================================================*/
int main(int argc, char **argv)
{
  gridls  *grid_list;        /* pointer to list of grids            */
  
  int     no_grids;          /* total number of grids               */
  int     no_timestep;       /* number of coarse grid timesteps     */
  int     no_first_timestep; /* number of initial timestep          */
  
  double  timecounter;       /* time variable                       */
  double  timestep;          /* timestep size                       */
  double  timecounter_final; /* for all sorts of tests...           */
  
  char     AMIGA_input[MAXSTRING];
  
#ifdef WITH_MPI
  uint64_t newparts;
#endif

#ifdef NEWAMR
  ahf2_patches_t *patches=NULL;
#endif
  
  /*============================================================ 
   * we always read the relevant parameters from an input file!
   *===========================================================*/
  if(argc<2)
   {
    fprintf(stderr,"usage:    %s AMIGA.input\n", argv[0]);
    fprintf(stderr,"       or %s --parameterfile\n", argv[0]);
    exit(1);
   }
  
  /*============================================================ 
   * maybe the user only wants the parameterfile?
   *===========================================================*/
  if(strcmp(argv[1],"--parameterfile") == 0)
   {
    global_io.params                 = (io_parameter_t) calloc(1,sizeof(io_parameter_struct_t));
    global_io.params->outfile_prefix = (char *) calloc(MAXSTRING,sizeof(char));
    global.a                         = 1;
    strcpy(global_io.params->outfile_prefix,"AHF");
    write_parameterfile();
    exit(0);
   }
  else
   {
    strcpy(AMIGA_input, argv[1]);
   }
  
  
  
  
  /* check for some DEFINEFLAGS mistakes */
#if (defined NCPUREADING_EQ_NFILES && defined BCASTHEADER)
  fprintf(stderr,"you cannot define NCPUREADING_EQ_NFILES and BCASTHEADER at the same time\nABORTING\n");
  exit(1);
#endif
  
  
  
  
#if (!defined WITH_MPI)
  WRITEAHFLOGO(stderr);
#endif
  
  /* how much memory per node and particle for this particular run */
  global.bytes_node = sizeof(struct node);
  global.bytes_part = sizeof(struct particle);
  
#	ifdef WITH_MPI
	/* Initialize the MPI environment */
	common_initmpi(&argc, &argv);
#		ifdef MPI_TIMING
	global_mpi.start = MPI_Wtime();
#		endif
#	endif
  
  /*======================================================== 
   * startrun:    input the initial data from infile 
   *========================================================*/
  timing.io       -= time(NULL);
  
  timing.startrun -= time(NULL);
	startrun((argc > 1) ? argv[1] : NULL, &timecounter, &timestep, &no_first_timestep);
  timing.startrun += time(NULL);
  
  
#ifdef DEBUG_STARTRUN
  /*===========================================================
   * DEBUG_STARTRUN:
   * we simply check if the particles have been read correctly
   *===========================================================*/
 {
  FILE *fpout;
  char outname[MAXSTRING];
  partptr cur_part;
  
#ifdef WITH_MPI
  sprintf(outname,"test-%d.ascii",global_mpi.rank);
#else
  sprintf(outname,"test.ascii");
#endif
  
  fpout = fopen(outname,"w");
  
  for(cur_part=global_info.fst_part; cur_part<(global_info.fst_part+global_info.no_part); cur_part++)
    fprintf(fpout,"%e %e %e\n",cur_part->pos[X]*simu.boxsize,cur_part->pos[Y]*simu.boxsize,cur_part->pos[Z]*simu.boxsize);
  
  fclose(fpout);
#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  //exit(0);
 }
#endif /* DEBUG_STARTRUN */
  
  
  /*==========================================================================================
   * AHFptfocus:
   *
   * only use a certain type of particles ("pt") and focus ("focus") the AHF analysis on them
   *
   *==========================================================================================*/
#if (defined AHFptfocus && defined MULTIMASS && defined GAS_PARTICLES)
  /* global_info.no_part
   * global_info.fst_part
   *                       => the no. of particles and relevant pointer for this CPU */
  timing.ptfocus -= time(NULL);
 {
#ifdef DEBUG_AHFptfocus
  FILE *fp;
  fp = fopen("ascii.dat","w");
#endif
  long unsigned no_part;
  partptr       fst_part, cur_part, new_part;
  int           ikeep;
  
  fprintf(stderr,"\n==================================================================\n");
  fprintf(stderr,"                          AHFptfocus\n");
  fprintf(stderr,"               ? ARE YOU SURE ABOUT THIS FLAG ?\n");
  fprintf(stderr,"==================================================================\n");
  fprintf(stderr,"AHF will now remove all particles whose type is not %d\n",AHFptfocus);
  fprintf(stderr,"starting with %ld particles -> ",global_info.no_part);
  
  /* 1. count number of particles to keep */
  no_part  = 0;
  for(cur_part=global_info.fst_part; cur_part<(global_info.fst_part+global_info.no_part); cur_part++)
   {
    /* we only want ot keep those particles with type AHFptfocus */
    if(AHFptfocus == 0)
     {
      if(cur_part->u >= AHFptfocus)
        no_part++;
     }
    else
     {
      if(fabs(cur_part->u+AHFptfocus) < ZERO)
        no_part++;
     }
    
    /* only keep the high-resolution particles */
    // if(cur_part->u >= 0 || cur_part->u == PDM || cur_part->u == PSTAR)
    //   no_part++;
    
   }
  
  /* allocate memory for new particles */
  fst_part = c_part(no_part);
  
  /* 2. remove all other particles */
  new_part = fst_part;
  for(cur_part=global_info.fst_part; cur_part<(global_info.fst_part+global_info.no_part); cur_part++)
   {
    ikeep = 0;
    
    /* we only want ot keep those particles with type AHFptfocus */
    if(AHFptfocus == 0)
     {
      if(cur_part->u >= AHFptfocus)
        ikeep = 1;
     }
    else
     {
      if(fabs(cur_part->u+AHFptfocus) < ZERO)
        ikeep = 1;
     }
    
    /* only keep the high-resolution particles */
    // if(cur_part->u >= 0 || cur_part->u == PDM || cur_part->u == PSTAR)
    //   ikeep = 1;
    
    if(ikeep)
     {
      new_part->pos[X] = cur_part->pos[X];
      new_part->pos[Y] = cur_part->pos[Y];
      new_part->pos[Z] = cur_part->pos[Z];
      new_part->mom[X] = cur_part->mom[X];
      new_part->mom[Y] = cur_part->mom[Y];
      new_part->mom[Z] = cur_part->mom[Z];
      new_part->weight = cur_part->weight;
      new_part->u      = cur_part->u;
#if (!(defined AHF_NO_PARTICLES && defined AHFlean))
      new_part->id     = cur_part->id;
#endif
      
#ifdef DEBUG_AHFptfocus
      fprintf(fp,"%f %f %f %f %f %f %f %f\n",
              new_part->pos[X],
              new_part->pos[Y],
              new_part->pos[Z],
              new_part->mom[X],
              new_part->mom[Y],
              new_part->mom[Z],
              new_part->weight,
              new_part->u);
#endif
      
      new_part++;
     }
   }
  
  /* erase old particle list and store new one */
  free(global_info.fst_part);
  global_info.fst_part = fst_part;
  
  /* update global.no_part parameter */
  global_info.no_part  = no_part;
  fprintf(stderr,"ended with %ld particles\n\n",global_info.no_part);
#ifdef DEBUG_AHFptfocus
  fclose(fp);
#endif
 }
  timing.ptfocus += time(NULL); 
#endif /* AHFptfocus */
  
  
  
  
#ifdef AHFrfocus
  /*====================================================================
   * This is for focussing on a Sphere defined in param.h
   * This assumes that periodicity can be neglected for deciding
   * whether a particle is inside the selected sphere or not.
   *====================================================================*/
  timing.rfocus -= time(NULL);
	local_focusSphere();
  timing.rfocus += time(NULL);
#endif 
  
  
#		if (defined WITH_MPI && defined MPI_TIMING)
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1), "Startrun done in %fs", global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#		endif
  
#		ifdef WITH_MPI
  timing.loadbalance -= time(NULL);
	/* Sort the particles in a particle block structure */
	io_logging_section(global_io.log, "Initial Load-Balancing and Particle Distribution");
  io_logging_subsection(global_io.log, "Loadbalancing");
  loadbalance_update(global_io.log, global_info.loadbal, global_info.fst_part, global_info.no_part);
#			ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1), "Loadbalance done in %fs", global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#			else
  io_logging_msg(global_io.log, INT32_C(1), "Loadbalance done.");
#			endif
  loadbalance_log(global_io.log, global_info.loadbal);
  timing.loadbalance += time(NULL);
#		else
	/* Generate the SFC keys for all particles */
  timing.sfckey -= time(NULL);
	for (uint64_t i=0; i<global_info.no_part; i++) {
		partptr part=global_info.fst_part+i;
		part->sfckey = sfc_curve_calcKey(global_info.ctype,
		                                 (double)(part->pos[0]),
		                                 (double)(part->pos[1]),
		                                 (double)(part->pos[2]),
		                                 BITS_PER_DIMENSION);
	}
	/* Sorting all particles to have fast access later on */
	qsort(global_info.fst_part,
	      global_info.no_part,
	      sizeof(part),
	      &cmp_sfckey_part);
  timing.sfckey += time(NULL);
#		endif /* WITH_MPI*/
  
#		ifdef WITH_MPI
  timing.distribution -= time(NULL);
  
  /* Do a first sort of the particles, required for distributing */
  io_logging_subsection(global_io.log, "Sorting particles");
  qsort(global_info.fst_part,
        global_info.no_part,
        sizeof(part),
        &cmp_sfckey_part);
#			ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1), "Sorting done in %fs", global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#			else
  io_logging_msg(global_io.log, INT32_C(1), "Sorting done.");
#			endif
  
  /* Distribute the particles */
  io_logging_subsection(global_io.log, "Distributing particles");
  io_logging_msg(global_io.log, INT32_C(0), "Currently having %"PRIu64" particles.", global_info.no_part);
  comm_dist_part(global_io.log,
                 &(global_info.fst_part),
                 &(global_info.no_part),
                 global_info.loadbal);
#			ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1), "Distributing done in %fs", global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#			else
  io_logging_msg(global_io.log, INT32_C(1), "Distributing done.");
#			endif
  io_logging_msg(global_io.log, INT32_C(0), "Having %"PRIu64" particles!", global_info.no_part);
  
	/* Do the AHF distribution*/
	io_logging_subsection(global_io.log, "AHF distribution (duplicating)");
	newparts = comm_dist_part_ahf(global_io.log,
	                              &(global_info.fst_part),
	                              &(global_info.no_part),
	                              global_info.loadbal);
	io_logging_msg(global_io.log, INT32_C(0), "Received %"PRIu64" new particles.", newparts);
	/* We need to sort the particles again */
	qsort(global_info.fst_part, global_info.no_part, sizeof(part), &cmp_sfckey_part);
#				ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1), "AHF distribution done in %fs", global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#				else
	io_logging_msg(global_io.log, INT32_C(1), "AHF distribution done.");
#				endif
  
  timing.distribution += time(NULL);
#		endif /* WITH_MPI */
  

#ifdef AHFsplit_only
  /*====================================================================
   * we only split the data using the SFC and
   * dump the data into multilpe files
   *====================================================================*/
 {
  io_file_t dumpf;
  io_file_strg_struct_t strg;
  char *fname;
  
  /* Start tge section */
  io_logging_section(global_io.log, "Dumping AHF chunk to file");
  
  /* First generate the filename */
  fname = (char *)malloc( sizeof(char) *( strlen(global_io.params->outfile_prefix)+30));
  if (fname == NULL) {
    io_logging_memfatal(global_io.log, "filename string");
    common_terminate(EXIT_FAILURE);
  }
  sprintf(fname, "%s.chunk.%04i.dump", global_io.params->outfile_prefix, global_mpi.rank);
  io_logging_msg(global_io.log, UINT32_C(0), "Used filename: %s", fname);
  
  fflush(NULL);
  MPI_Barrier(MPI_COMM_WORLD);
  
  /* Assign particles to structure */
  strg.posx.val = (void *)(global_info.fst_part->pos);
  strg.posx.stride =   (char *)((global_info.fst_part+1)->pos  ) - (char *)(global_info.fst_part->pos);
  strg.posy.val = (void *)(global_info.fst_part->pos+1);
  strg.posy.stride =   (char *)((global_info.fst_part+1)->pos+1) - (char *)(global_info.fst_part->pos+1);
  strg.posz.val = (void *)(global_info.fst_part->pos+2);
  strg.posz.stride =   (char *)((global_info.fst_part+1)->pos+2) - (char *)(global_info.fst_part->pos+2);
  strg.momx.val = (void *)(global_info.fst_part->mom);
  strg.momx.stride =   (char *)((global_info.fst_part+1)->mom  ) - (char *)(global_info.fst_part->mom);
  strg.momy.val = (void *)(global_info.fst_part->mom+1);
  strg.momy.stride =   (char *)((global_info.fst_part+1)->mom+1) - (char *)(global_info.fst_part->mom+1);
  strg.momz.val = (void *)(global_info.fst_part->mom+2);
  strg.momz.stride =   (char *)((global_info.fst_part+1)->mom+2) - (char *)(global_info.fst_part->mom+2);
#	ifdef MULTIMASS
  strg.weight.val = (void *)&(global_info.fst_part->weight);
  strg.weight.stride =   (char *)&((global_info.fst_part+1)->weight) - (char *)&(global_info.fst_part->weight);
#	else
  strg.weight.val = NULL;
  strg.weight.stride = (ptrdiff_t)0;
#	endif /* MULTIMASS */
#	ifdef GAS_PARTICLES
  strg.u.val = (void *)&(global_info.fst_part->u);
  strg.u.stride =   (char *)&((global_info.fst_part+1)->u) - (char *)&(global_info.fst_part->u);
#	else
  strg.u.val = NULL;
  strg.u.stride = (ptrdiff_t)0;
#	endif /* GAS_PARTICLE */
#	if (defined AHFlean && defined AHF_NO_PARTICLES)
  strg.id.val = NULL;
  strg.id.stride = (ptrdiff_t)0;
#	else
  strg.id.val = &(global_info.fst_part->id);
  strg.id.stride =   (char *)&((global_info.fst_part+1)->id) - (char *)&(global_info.fst_part->id);
#	endif
  strg.bytes_float = sizeof(global_info.fst_part->pos[0]);
#	if (defined AHFlean && defined AHF_NO_PARTICLES)
  strg.bytes_int = 0;
#	else
  strg.bytes_int = sizeof(global_info.fst_part->id);
#	endif
  
  /* Open the dump file now */
  dumpf = io_file_open(global_io.log, fname, IO_FILE_ARES, IO_FILE_UNKOWN_SWAPPING, IO_FILE_WRITE, 0);
  
  /* Write the particles */
  io_file_writepart(global_io.log, dumpf, 0, global_info.no_part, strg);
  
  /* Set the header values */
  ((io_ares_t)dumpf)->header->no_part = (uint64_t)simu.no_part;
  ((io_ares_t)dumpf)->header->no_species = UINT64_C(0);
  ((io_ares_t)dumpf)->header->no_vpart = simu.no_vpart;
  ((io_ares_t)dumpf)->header->boxsize = simu.boxsize;
  ((io_ares_t)dumpf)->header->omega0 = simu.omega0;
  ((io_ares_t)dumpf)->header->lambda0 = simu.lambda0;
  ((io_ares_t)dumpf)->header->pmass = simu.pmass;
  ((io_ares_t)dumpf)->header->minweight = simu.min_weight;
  ((io_ares_t)dumpf)->header->maxweight = simu.max_weight;
  ((io_ares_t)dumpf)->header->a_initial = simu.a_initial;
  ((io_ares_t)dumpf)->header->a_current = global.a;
  ((io_ares_t)dumpf)->header->timestep = timestep;
  ((io_ares_t)dumpf)->header->minkey = global_info.loadbal->fstkey[global_mpi.rank];
  ((io_ares_t)dumpf)->header->maxkey = global_info.loadbal->lstkey[global_mpi.rank];
  ((io_ares_t)dumpf)->header->lb_level = global_info.loadbal->level;
  ((io_ares_t)dumpf)->header->rank = global_mpi.rank;
  ((io_ares_t)dumpf)->header->size = global_mpi.size;
  
  /* Log the file */
  io_file_log(global_io.log, dumpf);
  
  /* Close the file and clean up*/		
  io_file_close(global_io.log, &dumpf);
  free(fname);
 }
	common_terminate(EXIT_SUCCESS);
#endif /*  AHFsplit_only */
  
  
  
  
#ifdef AHF_DUMP_AFTER_READ_TO_ASCII
  /*====================================================================
   * write an ASCII file of the data just read
   *====================================================================*/
 {
  FILE *dumpf;
  char *fname;
  
  /* First generate the filename */
  fname = (char *)malloc( sizeof(char) *( strlen(global_io.params->outfile_prefix)+35));
  if (fname == NULL) {
    io_logging_memfatal(global_io.log, "filename string");
    common_terminate(EXIT_FAILURE);
  }
#ifdef WITH_MPI
  sprintf(fname, "%s.chunk.%04i.ascii", global_io.params->outfile_prefix, global_mpi.rank);
#else
  sprintf(fname, "%s.DUMP.ascii", global_io.params->outfile_prefix);
#endif
  io_logging_msg(global_io.log, UINT32_C(0), "Used filename: %s", fname);
  fflush(NULL);
#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  dumpf = fopen(fname, "w");
  fprintf(dumpf, "# x y z  vx vy vz  ID\n");
  for (uint64_t i=0L; i<global_info.no_part; i++) {
    fprintf(dumpf, "%15e %15e %15e   %15e %15e %15e   %lu\n",
            global_info.fst_part[i].pos[0],
            global_info.fst_part[i].pos[1],
            global_info.fst_part[i].pos[2],
            global_info.fst_part[i].mom[0],
            global_info.fst_part[i].mom[1],
            global_info.fst_part[i].mom[2],
            (unsigned long)global_info.fst_part[i].id);
  }
  fclose(dumpf);
  common_terminate(EXIT_SUCCESS);
 }
#endif

  
#ifdef WITH_MPI
	loadbalance_minimalMemory(global_io.log, global_info.loadbal);
#endif
	io_logging_msg(global_io.log, INT32_C(5), "amiga_main:  running with %" PRIu64 " particles", global_info.no_part);
	io_logging_part(global_io.log, "Handing over logging to AMIGA");
  
  timing.io       += time(NULL);
  
  
  
  
  /*===================================================================== 
   * at this point we completely read in the data file
   * and are ready to proceed with generating the
   * grid hierarchy or what else we plan to do...
   *=====================================================================*/
  
  
  
  
  
  
  /*====================================================================
   *  GENERATE THE FULL BLOWN AMR HIERARCHY AND ORGANIZE IT INTO A TREE
   *====================================================================*/
  
#ifdef NEWAMR
{
  io_logging_msg(global_io.log, INT32_C(0), "### Executing with NEWAMR ###\n");
  
  /* 1. organize the particles into a tree */
  fprintf(stderr,"[main] Calling generate_tree...\n");
#ifndef AHF2_read_spatialRef
  patches=generate_tree(global_info.no_part, global_info.fst_part, simu.Nth_dom, simu.AHF_MINPART);
#else
  patches = (ahf2_patches_t *) calloc(1, sizeof(ahf2_patches_t));
#endif
  fprintf(stderr,"[main] Exit from generate_tree\n");

	//Generate patchtree.out
  //patch_connection_review(patches);
  
  /* 2. moving things around */
	global.fst_part     = global_info.fst_part;
  global.no_part      = global_info.no_part;
  global.total_time   = 0.;
  global.output_count = 0;
  global.no_timestep  = no_first_timestep;

  /* 3. perform halo analysis */
  fprintf(stderr,"\nmain: calling ahf_halos():\n");
  ahf.time -= time(NULL);
  timing.ahf_halos -= time(NULL);
  ahf_halos(patches);
  timing.ahf_halos += time(NULL);
  ahf.time += time(NULL);  
}
#else /* NEWAMR */
  
  /*===================================================================== 
   * generate the domain grids: simu.NGRID_MIN^3, ...., simu.NGRID_DOM^3 
   *=====================================================================*/
  timing.gendomgrids -= time(NULL);
  grid_list = gen_domgrids(&no_grids);   
  timing.gendomgrids += time(NULL);
  
  /*===================================================================== 
   * build initial linked list 
   *=====================================================================*/
  timing.ll -= time(NULL);
	ll(global_info.no_part, global_info.fst_part, global.dom_grid);
	global.fst_part = global_info.fst_part;
  global.no_part  = global_info.no_part;
  timing.ll += time(NULL);
  
  /*================================================================
   * assign particles to the domain grid with simu.NGRID_DOM^3 nodes 
   *================================================================*/
  zero_dens(global.dom_grid);
  assign_npart(global.dom_grid);
  
  /*================================================================
   * initialize some counters 
   *================================================================*/
  no_timestep         = no_first_timestep+1;  /* count total number of integration steps */
  global.total_time   = 0.;                   /* cumulative total time for simulation    */
  global.output_count = 0;                    /* count the number of outputs             */
  
  /* make *current* time step available to AHF/etc. routines */
  global.no_timestep = no_first_timestep;
  
  /*========================================================================================= 
   * recursively call gen_AMRhierarchy() to generate the AMR hierarchy...
   *=========================================================================================*/
  global.fst_cycle = TRUE;
  gen_AMRhierarchy(&grid_list, &no_grids);
  
  /*========================================================================================= 
   * eventually perform AHF analysis of AMR hierarchy
   *=========================================================================================*/
  ahf.time -= time(NULL);
  
  /* get spatially connected refinement patches */
  timing.ahf_gridinfo -= time(NULL);
  ahf_gridinfo(grid_list, no_grids-1);
  timing.ahf_gridinfo += time(NULL);
  
  
  /* get AHF halos */
  timing.ahf_halos -= time(NULL);
  ahf_halos(grid_list);
  timing.ahf_halos += time(NULL);    
  
  ahf.time += time(NULL);
 
  /*========================================================================================= 
   * update logfile and say bye-bye
   *=========================================================================================*/
  write_logfile(timecounter, timestep, no_timestep);
  
  /* free all allocated memory... */
  free(grid_list);
#endif /* NEWAMR */  
  
  
  /*============================================================================
   *                                   BYE BYE!
   *============================================================================*/  
  free(io.icfile_name);
  free(io.dumpfile_name);
  free(io.logfile_name);
  free(io.outfile_prefix);
  free(global.termfile_name);
  
  free(global.fst_part);
  if(global.fst_gas)
    free(global.fst_gas);
  if(global.fst_star)
    free(global.fst_star);
  
  
  fprintf(io.logfile, "==========================================================\n");
  fprintf(io.logfile, "                       FINISHED (v%3.1f/%03d)\n",VERSION,BUILD);
  fprintf(io.logfile, "==========================================================\n");
  fclose(io.logfile);
  
#	ifdef WITH_MPI
	/* Gracefully terminate MPI */
	MPI_Finalize();
#	endif
  
  return EXIT_SUCCESS;
}

#ifdef WITH_MPI
/*====================================================================
 *  MPI communications
 *====================================================================*/
static void
local_communicate_all(void)
{
	if (global_mpi.rank == 0) {
		/*********************\
		 *      SENDING      *
     \*********************/
		/** This is where to send to */
		int target;
    
		for (target = 1; target < global_mpi.size; target++) {
			MPI_Send(&global, sizeof(struct info_global),
			         MPI_BYTE,
			         target, 0, MPI_COMM_WORLD);
		}
		/*********************\
		 *    DONE SENDING   *
     \*********************/
	} else {
		/*********************\
		 *     RECIEVING     *
     \*********************/
		/** The status */
		MPI_Status stat;
    
		MPI_Recv(&global, sizeof(struct info_global),
		         MPI_BYTE,
		         0, 0, MPI_COMM_WORLD, &stat);
    
		/*********************\
		 *  DONE  RECIEVING  *
     \*********************/
	}
  
	return;
}
#endif

#ifdef AHFrfocus
/*====================================================================
 *  remove all particles outside a spherical region about centre[]
 *
 *  (though called AHFrfocus it is a multi-purpose routine)
 *====================================================================*/
static void
local_focusSphere(void)
{
	uint64_t oldNumPart = global_info.no_part;
	uint64_t newNumPart = UINT64_C(0);
	double center[3];
	double radSqr;
	int8_t   *tags;
	uint64_t i, j;
	partptr newParts;
  
	/* The sphere to focus on, give in Mpc/h (cf. param.h) */
	center[X] = AHFrfocusX;
	center[Y] = AHFrfocusY;
	center[Z] = AHFrfocusZ;
	radSqr    = AHFrfocusR;
  
	fprintf(stderr, "ATTENTION!  This is AHFrfocus calling:\n");
	fprintf(stderr, "This will remove all particles outside this sphere:\n");
	fprintf(stderr, "Center: (%lf %lf %lf) Mpc/h",
	        center[X], center[Y], center[Z]);
	fprintf(stderr, "Radius: %lf Mpc/h\n", radSqr);
	fprintf(stderr,"starting with %ld particles -> ",global_info.no_part);
  
	/* Convert center and radius to AHF units */
	center[X] /= simu.boxsize;
	center[Y] /= simu.boxsize;
	center[Z] /= simu.boxsize;
	radSqr /= simu.boxsize;
	radSqr *= radSqr;
  
	/* Keeps a record whether the particle is kept or not. */
	tags = malloc(sizeof(int8_t) * oldNumPart);
  
	/* Check which particles to keep (ignoring periodicity) */
	for (i=UINT64_C(0); i<oldNumPart; i++) {
		double dpos[3];
		double rSqr;
    
		dpos[X] = global_info.fst_part[i].pos[X] - center[X];
		dpos[Y] = global_info.fst_part[i].pos[Y] - center[Y];
		dpos[Z] = global_info.fst_part[i].pos[Z] - center[Z];
		rSqr  = dpos[X] * dpos[X] + dpos[Y] * dpos[Y]
    + dpos[Z] * dpos[Z];
		if (rSqr <= radSqr) {
			tags[i] = 1;
			newNumPart++;
		} else {
			tags[i] = 0;
		}
	}
  
	/* Get new particles and copy the old ones over. */
	newParts = c_part(newNumPart);
	for (j = i = UINT64_C(0); i<oldNumPart; i++) {
		if (tags[i] == 1) {
			memcpy(newParts + j, global_info.fst_part + i,
			       sizeof(struct particle));
			j++;
		}
	}
	assert(j == newNumPart);
  
	/* Clean */
	free(tags);
	free(global_info.fst_part);
  
	/* Activate the new particles. */
	global_info.fst_part = newParts;
	global_info.no_part  = newNumPart;
	fprintf(stderr,"ended with %ld particles\n\n",global_info.no_part);
}
#endif


