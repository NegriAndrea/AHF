#include "define.h"


#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "libutility/utility.h"
#include "libio/io.h"
#include "libio/io_file.h"

#ifdef WITH_OPENMP
#include <omp.h>
#endif

/**
 * \brief Helper function for reading the parameters.
 */
static void
local_startrunParams(char *paramfile);

/**
 * \brief Helper function for setting up the logging.
 */
static void
local_startrunLog(void);

/**
 * \brief Helper function for opening and initializing the IC file.
 */
static void
local_startrunFopen(void);

/**
 * \brief Helper function for reading the IC file.
 */
static void
local_startrunRead(void);

/**
 * \brief Helper function for setting the simulation parameters.
 */
static void
local_startrunSimparams();

/**
 * \brief Helper function for setting the return values.
 */
static void
local_startrunRetset(double *timecounter,
                     double *timestep,
                     int32_t *no_first_timestep);


/*====================================================================================================
 * startrun()
 *====================================================================================================*/
extern void startrun(char *paramfile, double *timecounter, double *timestep, int32_t *no_first_timestep)
{
	/* Read the parameters */
	local_startrunParams(paramfile);


  /* Now set up the logging */
	local_startrunLog();
	io_logging_part(global_io.log, "Setting up the run");
  
	/* FIXME Do some glueing FIXME */
	io.logfile = global_io.log->logfile;
  
	/* Open the file */
	local_startrunFopen();
  
	/* Set global_info */
	global_info.fst_part = NULL;
	global_info.no_part = UINT64_C(0);
#	ifdef WITH_MPI
	global_info.loadbal = loadbalance_new(global_io.log,
	                                      LOADBALANCE_EQUALPART,
	                                      SFC_CURVE_HILBERT,
                                        //	                                      LOADBALANCE_DOMAIN_LEVEL,
	                                      global_io.params->lb_level,
	                                      global_mpi.size);
#	endif
  
	/* Now read the initial conditions (sets the rest of global_info) */
	local_startrunRead();
  
	/* Load the dark energy tables for interpolation if DE is defined */
#ifdef DARK_ENERGY
	read_dark_energy_table(global_io.params->defile_name);
#endif

	/* Set the simulation parameters */
	local_startrunSimparams();
  
	/* Now set the returned simulation counter stuff thingies */
	local_startrunRetset(timecounter, timestep, no_first_timestep);
  
	/* Set global time counter */
	global.super_t = *timecounter;                      
	global.a       = calc_super_a(global.super_t);   
	global.t       = calc_t(global.a);            
	global.z       = 1./global.a - 1.;
  
#ifdef NO_EXPANSION
  global.a       = 1.0;
  global.t       = *timecounter;
  global.z       = 0.0;
  global.super_t = -1.;
#endif
  
	/* And now lets be so nice and close the file... */
	io_logging_section(global_io.log, "Tidying");
	io_file_close(global_io.log, &(global_io.file));


  
  write_parameterfile();
  
	return;
}

extern void
stoprun(void)
{	
	io_parameter_del(&(global_io.params));
	io_logging_stop(&(global_io.log));
  
	return;
}

static void
local_startrunParams(char *paramfile)
{
	global_io.params = io_parameter_get(paramfile);
	if (global_io.params == NULL) {
		common_terminate(EXIT_FAILURE);
	}
  
	return;
}

static void
local_startrunLog(void)
{
#ifdef WITH_OPENMP
  int nthreads, ncpus, tid;
#endif
  
	global_io.log = io_logging_start(global_io.params->outfile_prefix, INT32_C(VERBOSITY), IO_LOGGING_FLAGS_DUPLICATE_CRITICAL);
	if (global_io.log == NULL) {
		io_parameter_del(&(global_io.params));
		common_terminate(EXIT_FAILURE);
	}
  
	/* And since we can log now, log a bit. */
	io_logging_hello(global_io.log, VERSION, BUILD);
	io_logging_msg(global_io.log, INT32_C(2), " User Input:");
	io_logging_msg(global_io.log, INT32_C(2), "=============");
	io_logging_msg(global_io.log, INT32_C(2), "ic_filename       = %s (%s, %i)", global_io.params->icfile_name, io_file_typestr(global_io.params->ic_filetype),(int)(global_io.params->ic_filetype));
	io_logging_msg(global_io.log, INT32_C(2), "outfile_prefix    = %s", global_io.params->outfile_prefix);
	io_logging_msg(global_io.log, INT32_C(2), "LgridDomain       = %i", global_io.params->NGRID_DOM);
	io_logging_msg(global_io.log, INT32_C(2), "LgridMax          = %i", global_io.params->NGRID_MAX);
	io_logging_msg(global_io.log, INT32_C(2), "NminPerHalo       = %i", global_io.params->AHF_MINPART);
	io_logging_msg(global_io.log, INT32_C(2), "VescTune          = %g", global_io.params->AHF_VTUNE);
	io_logging_msg(global_io.log, INT32_C(2), "NperDomCell       = %g", global_io.params->Nth_dom);
	io_logging_msg(global_io.log, INT32_C(2), "NperRefCell       = %g", global_io.params->Nth_ref);
	io_logging_msg(global_io.log, INT32_C(2), "RhoVir            = %i", global_io.params->UseRhoBack);
	io_logging_msg(global_io.log, INT32_C(2), "Dvir              = %g", global_io.params->UserDvir);
	io_logging_msg(global_io.log, INT32_C(2), "MaxGatherRad      = %g Mpc/h", global_io.params->MaxGatherRad);
#ifdef WITH_MPI
	io_logging_msg(global_io.log, INT32_C(2), "LevelDomainDecomp = %i", global_io.params->lb_level);
	io_logging_msg(global_io.log, INT32_C(2), "NcpuReading       = %" PRIu32 " ", global_io.params->reader);
#endif
#ifdef AHF_LRSI
	io_logging_msg(global_io.log, INT32_C(2), "lrsi_beta         = %g", global_io.params->lrsi_beta);
	io_logging_msg(global_io.log, INT32_C(2), "lrsi_r_s          = %g", global_io.params->lrsi_r_s);
#endif
#ifdef DARK_ENERGY
	io_logging_msg(global_io.log, INT32_C(2), "DarkEnergyFile    = %s", global_io.params->defile_name);
#endif
  
#ifdef WITH_OPENMP
	io_logging_msg(global_io.log, INT32_C(2), "\n\n OpenMP Info:");
	io_logging_msg(global_io.log, INT32_C(2), "==============");
#pragma omp parallel private(nthreads, tid)
  {
   /* Obtain thread number */
   tid = omp_get_thread_num();
   
   /* Only master thread does this */
   if (tid == 0)
    {
     ncpus    = omp_get_num_procs();
     nthreads = omp_get_num_threads();
     io_logging_msg(global_io.log, INT32_C(2), "Number of available processors = %d", ncpus);
     io_logging_msg(global_io.log, INT32_C(2), "Number of threads in use       = %d", nthreads);
    }
  }
#endif /* WITH_OPENMP */
  
	return;
}

static void
local_startrunFopen(void)
{
	io_logging_section(global_io.log, "Opening the data file");
  
	global_io.file = io_file_open(global_io.log,
	                              global_io.params->icfile_name,
	                              global_io.params->ic_filetype,
#	ifdef BYTESWAP
	                              IO_FILE_IS_SWAPPED,
#	else
	                              IO_FILE_UNKOWN_SWAPPING,
#	endif
                                
	                              IO_FILE_READ,
	                              global_io.params->reader);
  
	if (global_io.file == NULL)
		common_terminate(EXIT_FAILURE);
  
	/* Set the scaling
   * NOTE, no reading/writing to any file will be performed; only the file structure will be accessed and updated */
	if (global_io.file->ftype == IO_FILE_GADGET)
		io_gadget_resetscale(global_io.log,
		                     (io_gadget_t)global_io.file,
		                     global_io.params->GADGET_l2Mpch,
		                     global_io.params->GADGET_m2Msunh);
	if (global_io.file->ftype == IO_FILE_MGADGET)
		io_mgadget_resetscale(global_io.log,
		                      (io_mgadget_t)global_io.file,
		                      global_io.params->GADGET_l2Mpch,
		                      global_io.params->GADGET_m2Msunh);
  

	/* Init the file */
	io_file_init(global_io.log, global_io.file);
  

	/* Now dump the file information */
	io_file_log(global_io.log, global_io.file);
  
	return;
}

static void
local_startrunRead(void)
{
	io_file_strg_struct_t strg;
	uint64_t pskip = UINT64_C(0);
	uint64_t pread = UINT64_MAX;
	partptr tmppart;
  
	io_logging_section(global_io.log, "Reading data from file");
  
	/* See if we are supposed to read anything at all */
	if (global_io.file->ftype == IO_FILE_EMPTY) {
		global_info.no_part = 0;
		global_info.fst_part = NULL;
		return;
	}
    
  /* First create particle storage */
	io_logging_subsection(global_io.log, "Creating Storage");
	global_info.no_part = io_file_get_numpart(global_io.log,
	                                          global_io.file,
	                                          &pskip, &pread);
	global_info.fst_part = c_part((long)global_info.no_part);
 
	/* Create the description of the storage */
	strg.posx.val = (void *)(global_info.fst_part->pos);
	strg.posx.stride =   (char *)((global_info.fst_part+1)->pos)
  - (char *)(global_info.fst_part->pos);
	strg.posy.val = (void *)(global_info.fst_part->pos+1);
	strg.posy.stride =   (char *)((global_info.fst_part+1)->pos+1)
  - (char *)(global_info.fst_part->pos+1);
	strg.posz.val = (void *)(global_info.fst_part->pos+2);
	strg.posz.stride =   (char *)((global_info.fst_part+1)->pos+2)
  - (char *)(global_info.fst_part->pos+2);
	strg.momx.val = (void *)(global_info.fst_part->mom);
	strg.momx.stride =   (char *)((global_info.fst_part+1)->mom)
  - (char *)(global_info.fst_part->mom);
	strg.momy.val = (void *)(global_info.fst_part->mom+1);
	strg.momy.stride =   (char *)((global_info.fst_part+1)->mom+1)
  - (char *)(global_info.fst_part->mom+1);
	strg.momz.val = (void *)(global_info.fst_part->mom+2);
	strg.momz.stride =   (char *)((global_info.fst_part+1)->mom+2)
  - (char *)(global_info.fst_part->mom+2);
#	ifdef MULTIMASS
	strg.weight.val = (void *)&(global_info.fst_part->weight);
	strg.weight.stride =   (char *)&((global_info.fst_part+1)->weight)
  - (char *)&(global_info.fst_part->weight);
#	else
	strg.weight.val = NULL;
	strg.weight.stride = (ptrdiff_t)0;
#	endif
#	if (!(defined AHFlean && defined AHF_NO_PARTICLES))
	strg.id.val = &(global_info.fst_part->id);
	strg.id.stride =   (char *)&((global_info.fst_part+1)->id)
  - (char *)&(global_info.fst_part->id);
#	else
	strg.id.val = NULL;
	strg.id.stride = (ptrdiff_t)0;
#	endif
#	ifdef GAS_PARTICLES
	strg.u.val = &(global_info.fst_part->u);
	strg.u.stride =   (char *)&((global_info.fst_part+1)->u)
  - (char *)&(global_info.fst_part->u);
#	else
	strg.u.val = NULL;
	strg.u.stride = (ptrdiff_t)0;
#	endif
	strg.bytes_float = sizeof(global_info.fst_part->pos[0]);
#	if (!(defined AHFlean && defined AHF_NO_PARTICLES))
	strg.bytes_int = sizeof(global_info.fst_part->id);
#	else
	strg.bytes_int = 0;
#	endif
  
#	ifdef METALHACK
	strg.z.val = &(global_info.fst_part->z);
	strg.z.stride =   (char *)&((global_info.fst_part+1)->z)
  - (char *)&(global_info.fst_part->z);
	strg.age.val = &(global_info.fst_part->age);
	strg.age.stride =   (char *)&((global_info.fst_part+1)->age)
  - (char *)&(global_info.fst_part->age);
#	endif
  
#	ifdef VERBOSE
	/* Print the description */
	io_file_strg_log(global_io.log, strg);
#	endif
  
  /* Now read the particles */
	io_logging_subsection(global_io.log, "Reading");
	if (io_file_readpart(global_io.log, global_io.file,
	                     pskip, pread, strg) == UINT64_C(0) ) {
		/* We read 0 particles from the file, this is an error */
		common_terminate(EXIT_FAILURE);
	}
  
	/* Print the first two particles to the logfile */
	io_logging_subsection(global_io.log, "Short sanity check");
	tmppart = global_info.fst_part;
	io_logging_msg(global_io.log, INT32_C(5),
	               "First particle:");
	io_logging_msg(global_io.log, INT32_C(5),
	               "    positions (x,y,z):      %g  %g  %g",
	               tmppart->pos[0],
	               tmppart->pos[1],
	               tmppart->pos[2]);
	io_logging_msg(global_io.log, INT32_C(5),
	               "    velocities (vx,vy,vz):  %g  %g  %g",
	               tmppart->mom[0],
	               tmppart->mom[1],
	               tmppart->mom[2]);
#	ifdef MULTIMASS
	io_logging_msg(global_io.log, INT32_C(5),
	               "    weight:                 %g",
	               tmppart->weight);
#	endif
#	if	(!(defined AHFlean && defined AHF_NO_PARTICLES))
	io_logging_msg(global_io.log, INT32_C(5),
	               "    ID:                     %" PRIpartid,
	               tmppart->id);
#	endif
#	ifdef GAS_PARTICLES
	io_logging_msg(global_io.log, INT32_C(5),
	               "    energy:                 %g",
	               tmppart->u);
#	endif
#	ifdef METALHACK
	io_logging_msg(global_io.log, INT32_C(5),
	               "    metallicity:            %g",
	               tmppart->z);
	io_logging_msg(global_io.log, INT32_C(5),
	               "    age:                    %g",
	               tmppart->age);
#	endif
	tmppart = global_info.fst_part+global_info.no_part-1;
	io_logging_msg(global_io.log, INT32_C(5),
	               "Last particle:");
	io_logging_msg(global_io.log, INT32_C(5),
	               "    positions (x,y,z):      %g  %g  %g",
	               (tmppart)->pos[0],
	               (tmppart)->pos[1],
	               (tmppart)->pos[2]);
	io_logging_msg(global_io.log, INT32_C(5),
	               "    velocities (vx,vy,vz):  %g  %g  %g",
	               (tmppart)->mom[0],
	               (tmppart)->mom[1],
	               (tmppart)->mom[2]);
#	ifdef MULTIMASS
	io_logging_msg(global_io.log, INT32_C(5),
	               "    weight:                 %g",
	               (tmppart)->weight);
#	endif
#	if	(!(defined AHFlean && defined AHF_NO_PARTICLES))
	io_logging_msg(global_io.log, INT32_C(5),
	               "    ID:                     %" PRIpartid,
	               (tmppart)->id);
#	endif
#	ifdef GAS_PARTICLES
	io_logging_msg(global_io.log, INT32_C(5),
	               "    energy:                 %g",
	               (tmppart)->u);
#	endif
#	ifdef METALHACK
	io_logging_msg(global_io.log, INT32_C(5),
	               "    metallicity:            %g",
	               tmppart->z);
	io_logging_msg(global_io.log, INT32_C(5),
	               "    age:                    %g",
	               tmppart->age);
#	endif
  
#	ifdef MPI_DEBUG
	io_logging_subsection(global_io.log, "Longer sanity check");
	io_logging_msg(global_io.log, INT32_C(0),
	               "Fileobject after reading particles (supposedly "
	               "with correct multimass information now).");
	io_file_log(global_io.log, global_io.file);
#	endif
  
  return;
}

static void
local_startrunSimparams()
{
	io_logging_section(global_io.log, "Setting simulation parameter");
  
	io_logging_subsection(global_io.log, "Information from file");
#	ifdef WITH_MPI
	if (global_mpi.rank != 0) {
		io_logging_msg(global_io.log, INT32_C(4), "Not setting up myself, will receive.");
    fflush(NULL);
	} else {
#	else
   {
#	endif
		int32_t no_timestep;
    
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_BOXSIZE, (void *)&(simu.boxsize));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_OMEGA0, (void *)&(simu.omega0));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_OMEGAL, (void *)&(simu.lambda0));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_PMASS, (void *)&(simu.pmass));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_NOPART, (void *)&(simu.no_part));
#		ifdef MULTIMASS
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_NOVPART, (void *)&(simu.no_vpart));
#		else
		simu.no_vpart = (double)(simu.no_part);
#		endif
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_NOSPECIES, (void *)&(simu.no_species));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_AINITIAL, (void *)&(simu.a_initial));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_DOUBLE, (void *)&(simu.double_precision));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_MMASS, (void *)&(simu.multi_mass));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_MINWEIGHT, (void *)&(simu.min_weight));
		io_file_get(global_io.log, global_io.file, IO_FILE_GET_MAXWEIGHT, (void *)&(simu.max_weight));
    
		/* Copy over the information contained in the parameter file */
		simu.NGRID_DOM     = global_io.params->NGRID_DOM;
		simu.NGRID_MIN     = simu.NGRID_DOM;
		simu.Nth_dom       = global_io.params->Nth_dom;
		simu.Nth_ref       = global_io.params->Nth_ref;
    simu.lb_level      = global_io.params->lb_level;
    
    //fprintf(stderr,"simu.lb_level=%d global_io.params->lb_level=%d\n",simu.lb_level,global_io.params->lb_level);
    
    simu.MaxGatherRad  = global_io.params->MaxGatherRad;
    simu.UserDvir      = global_io.params->UserDvir;
    simu.UseRhoBack    = global_io.params->UseRhoBack;
		simu.NGRID_MAX     = global_io.params->NGRID_MAX;
		simu.AHF_MINPART   = global_io.params->AHF_MINPART;
		simu.AHF_VTUNE     = global_io.params->AHF_VTUNE;
    
    simu.GADGET_m2Msunh= global_io.params->GADGET_m2Msunh;
    simu.GADGET_l2Mpch = global_io.params->GADGET_l2Mpch;

#		ifdef AHF_LRSI
		simu.lrsi_beta     = global_io.params->lrsi_beta;
		simu.lrsi_r_s      = global_io.params->lrsi_r_s;
#		endif

#if (defined AHFmixHaloIDandSnapID || defined SUSSING2013)
    simu.isnap         = global_io.params->isnap;
#endif
    
		/* Set quantities given by constants */
#		ifdef NP_LIMIT
		simu.np_limit = TRUE;
#		else
		simu.np_limit = FALSE;
#		endif
		simu.mean_dens = (double) 1.0;
    
#ifdef MULTIMASS
    simu.multi_mass = 1;
#endif
    
		simu.mmfocus  = 0;
		simu.hydro    = 0;
		simu.magneto  = 0;
    
		/* Set the time unit */
		simu.t_unit = 1/H0; // we assume that we only ever deal with
		                    // cosmological simulations...
    
		/* Set derived quantities */
		simu.SHIFT     = ((double)0.5000000/(double) simu.NGRID_DOM);
		simu.z_initial = (double)1.0/simu.a_initial - (double)1.0;
		simu.a_final   = (double)1.0/((double)1.0 + simu.z_final);
		simu.FourPiG   = 1.5*simu.omega0;
    
		/* Do some sanity checks */
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_NOTSTEP, (void *)&(no_timestep));
    
		if ( isless(fabs(simu.a_initial-simu.a_final), ZERO) ) {
			io_logging_warn(global_io.log, INT32_C(3),
			                "Since a_initial = %g is equal to "
			                "a_final = %g, create_timeline will not "
			                "function correctly, setting "
			                "a_initial = .1 * a_final = %g",
			                simu.a_initial, simu.a_final,
			                simu.a_final / 10.0);
			simu.a_initial = simu.a_final / 10.0;
		}
		if ( simu.a_initial > simu.a_final ) {
			io_logging_warn(global_io.log, INT32_C(3),
			                "Since a_initial = %g is greater than "
			                "a_final = %g, create_timeline will not "
			                "function correctly, setting "
			                "a_initial = 0.001",
			                simu.a_initial, simu.a_final);
			simu.a_initial = 0.001;
      simu.z_initial = 1./simu.a_initial - 1.;
		}
   } /* End of stuff done solely by process 0 */
    
    io_logging_subsection(global_io.log, "Gathering from reading processes");
#	ifdef WITH_MPI
    io_logging_msg(global_io.log, INT32_C(4), "Broadcast of simulation parameters!");
    MPI_Bcast(&simu, sizeof(struct param_simu),
              MPI_BYTE,
              0,
              MPI_COMM_WORLD);
    io_logging_msg(global_io.log, INT32_C(4), "Broadcast done.");
#	endif
    
    
    /* Create timeline */
    io_logging_subsection(global_io.log, "Local setup");
    io_logging_msg(global_io.log, INT32_C(2), "Creating timeline from a = %g to a = %g",
                   simu.a_initial/10., simu.a_final);
    create_timeline(simu.a_initial/10., simu.a_final, &simu.timeline);
    io_logging_msg(global_io.log, INT32_C(2), "Timeline created");
    
    /* Set the SFC information */
    io_logging_msg(global_io.log, INT32_C(2), "Setting volume boundaries");
#	ifdef AHFrestart
    if(global_io.params->ic_filetype != IO_FILE_ARES)
     {
      fprintf(stderr,"AHFrestart only works together with ic_filetype=5\nPlease make sure that your input file is of this type and adjust AHF.input\nExiting now!\n");
      exit(0);
     }
    global_info.minkey = (sfc_key_t)( ((io_ares_t)(global_io.file))->header->minkey);
    global_info.maxkey = (sfc_key_t)( ((io_ares_t)(global_io.file))->header->maxkey);
    global_info.level = ((io_ares_t)(global_io.file))->header->lb_level;
#	else
    //    global_info.level = LOADBALANCE_DOMAIN_LEVEL;
    global_info.level = global_io.params->lb_level;
    global_info.minkey = (sfc_key_t)0;
    global_info.maxkey = (sfc_key_t)((1<<(3*global_info.level))-1);
#	endif
    global_info.ctype = SFC_CURVE_HILBERT;
    io_logging_msg(global_io.log, INT32_C(2),  "  minkey: %"SFC_PRIkey, global_info.minkey);
    io_logging_msg(global_io.log, INT32_C(2),  "  maxkey: %"SFC_PRIkey, global_info.maxkey);
    io_logging_msg(global_io.log, INT32_C(2),  "  level : %i", global_info.level);
    io_logging_msg(global_io.log, INT32_C(2),  "  ctype : %s", sfc_curve_typestr(global_info.ctype));
    
    /* Now that we have the timeline, set the time variables */
    simu.super_t_initial = calc_super_t(simu.a_initial);
    simu.super_t_final   = calc_super_t(simu.a_final);
    simu.t_initial       = calc_t(simu.a_initial);
    simu.t_final         = calc_t(simu.a_final);
    
    
    /* FIXME
     * Not set or not properly set simu-structure members:
     *
     * Hydro variables:
     * gamma
     * omegab
     * omegaDM
     * f_b
     * H_frac
     * T_init
     * e_init
     * med_weight
     * l_unit
     * m_unit
     *
     * AHF variable:
     * no_halos
     *
     * Unkown:
     * ifdef GAS_PARTICLES: no_gas
     * ifdef GADGET: no_stars
     *
     */
    /* Will use dummy values */
    simu.gamma    = 0.0;
    simu.omegab   = 0.0;
    simu.omegaDM  = simu.omega0;
    simu.f_b      = 0.0;
    simu.H_frac   = 0.0;
    simu.T_init   = 0.0;
    simu.B_init   = 0.0;
    simu.e_init   = 0.0;
    simu.no_halos = 0;
    simu.med_weight = simu.max_weight; // TODO: this is very conservative yet leads to more credible halos
    simu.l_unit = 0.0;
    simu.m_unit = 0.0;
    simu.no_gas   = 0;
    simu.no_stars = 0;
    
    //#	ifdef VERBOSE
    /* Be so kind and write everything to the logfile */
    io_logging_subsection(global_io.log, "Used simulation parameters");
    io_logging_msg(global_io.log, INT32_C(5), "simu.omega0          :  %g", simu.omega0);
    io_logging_msg(global_io.log, INT32_C(5), "simu.lambda0         :  %g", simu.lambda0);
    io_logging_msg(global_io.log, INT32_C(5), "simu.boxsize         :  %g", simu.boxsize);
    io_logging_msg(global_io.log, INT32_C(5), "simu.a_initial       :  %g", simu.a_initial);
    io_logging_msg(global_io.log, INT32_C(5), "simu.a_final         :  %g", simu.a_final);
    io_logging_msg(global_io.log, INT32_C(5), "simu.z_initial       :  %g", simu.z_initial);
    io_logging_msg(global_io.log, INT32_C(5), "simu.z_final         :  %g", simu.z_final);
    io_logging_msg(global_io.log, INT32_C(5), "simu.t_initial       :  %g", simu.t_initial);
    io_logging_msg(global_io.log, INT32_C(5), "simu.t_final         :  %g", simu.t_final);
    io_logging_msg(global_io.log, INT32_C(5), "simu.super_t_initial :  %g", simu.super_t_initial);
    io_logging_msg(global_io.log, INT32_C(5), "simu.super_t_final   :  %g", simu.super_t_final);
    io_logging_msg(global_io.log, INT32_C(5), "simu.mean_dens       :  %g", simu.mean_dens);
    io_logging_msg(global_io.log, INT32_C(5), "simu.FourPiG         :  %g", simu.FourPiG);
    io_logging_msg(global_io.log, INT32_C(5), "simu.pmass           :  %g", simu.pmass);
    io_logging_msg(global_io.log, INT32_C(5), "simu.t_unit          :  %g", simu.t_unit);
    io_logging_msg(global_io.log, INT32_C(5), "simu.gamma           :  %g", simu.gamma);
    io_logging_msg(global_io.log, INT32_C(5), "simu.timeline (ptr)  :  %p", (void*)&(simu.timeline));
    io_logging_msg(global_io.log, INT32_C(5), "simu.no_part         :  %lu", simu.no_part);
    io_logging_msg(global_io.log, INT32_C(5), "simu.no_vpart        :  %g", simu.no_vpart);
    io_logging_msg(global_io.log, INT32_C(5), "simu.no_species      :  %i", simu.no_species);
    io_logging_msg(global_io.log, INT32_C(5), "simu.no_halos        :  %lu", simu.no_halos);
    io_logging_msg(global_io.log, INT32_C(5), "simu.NGRID_DOM       :  %i", simu.NGRID_DOM);
    io_logging_msg(global_io.log, INT32_C(5), "simu.NGRID_MIN       :  %i", simu.NGRID_MIN);
    io_logging_msg(global_io.log, INT32_C(5), "simu.NGRID_MAX       :  %i", simu.NGRID_MAX);
    io_logging_msg(global_io.log, INT32_C(5), "simu.Nth_dom         :  %g", simu.Nth_dom);
    io_logging_msg(global_io.log, INT32_C(5), "simu.Nth_ref         :  %g", simu.Nth_ref);
    io_logging_msg(global_io.log, INT32_C(5), "simu.MaxGatherRad    :  %g", simu.MaxGatherRad);
    io_logging_msg(global_io.log, INT32_C(5), "simu.lb_level        :  %d", simu.lb_level);
    io_logging_msg(global_io.log, INT32_C(5), "simu.min_weight      :  %g", simu.min_weight);
    io_logging_msg(global_io.log, INT32_C(5), "simu.max_weight      :  %g", simu.max_weight);
    io_logging_msg(global_io.log, INT32_C(5), "simu.np_limit        :  %i", simu.np_limit);
    io_logging_msg(global_io.log, INT32_C(5), "simu.mmfocus         :  %i", simu.mmfocus);
    io_logging_msg(global_io.log, INT32_C(5), "simu.multi_mass      :  %i", simu.multi_mass);
    io_logging_msg(global_io.log, INT32_C(5), "simu.double_precision:  %i", simu.double_precision);
    //#	endif /* VERBOSE */
    
    
    //fprintf(stderr,"simu.lb_level=%d global_io.params->lb_level=%d\n",simu.lb_level,global_io.params->lb_level);

    
    return;
  }
  
  static void
  local_startrunRetset(double *timecounter,
                       double *timestep,
                       int32_t *no_first_timestep)
 {
	io_logging_subsection(global_io.log, "Setting time counter");
  
#	ifdef WITH_MPI
	if (global_mpi.rank == 0) {
#	else
   {
#	endif
		double a_current;
    
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_NOTSTEP, (void *)no_first_timestep);
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_TSTEP, (void *)timestep);
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_A, (void *)&a_current);
		*timecounter = calc_super_t(a_current);
   }
    
#	ifdef WITH_MPI
    MPI_Bcast(timecounter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(timestep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(no_first_timestep, 1, MPI_INT, 0, MPI_COMM_WORLD);
#	endif
    
    return;
  }
  
