#ifndef DEFINE_INCLUDED
#define DEFINE_INCLUDED

/*=============================================================================
 * this is written into the logfile just for information
 *=============================================================================*/
#define VERSION 1.0
#define BUILD   116

/*=============================================================================
 * here we switch on/off various features of AMIGA (i.e. DEFINEFLAGS)
 * (#define statements without actually defining a value...)
 *=============================================================================*/

/*--------------------------------------------------- 
 *                  MISC
 *--------------------------------------------------*/

#define PERIODIC                     // use periodic boundary conditions 
#define VERBOSE                      // let the user know what's going on
//#define VERBOSE2                     // dump as much runtime information as possible

//#define MULTIMASS                    // you MUST switch this on when the simulation features particles of different masses
//#define BYTESWAP                     // forces a byteswap of the input file
 
//#define GAS_PARTICLES                // you MUST switch this on when the simulation contains gas and/or star particles
                                       /* a few more words about this switch:
                                        *   - for historical reasons it is called GAS_PARTICLES even though it actually
                                        *     deals with baryons in the simulation
                                        *   - if you do not switch */
//#define WITH_MPI                     // switch on MPI domain decomposition
//#define WITH_OPENMP                  // switch on OpenMP parallisation of for-loops
//#define REFINE_BARYONIC_MASS         // use mass as refinement criterion for baryons (but number density for dark matter!)
//#define CHECK_RLIMIT_NOFILE          // uses system functions to increase file descriptor limitation if needed
#define FOPENCLOSE                   // open/close files, rather than opening multiple simulation files at the same time
#define BCASTHEADER                  // only one MPI task will read all the relevant header information and then broadcast
//#define NCPUREADING_EQ_NFILES        // this should speed up I/O of multiple snapshot files, but only works for this condition

//#define EXTRAE_API_USAGE              // activate the Extrae API calls to trace user functions with Extrae_user_function()


/*--------------------------------------------------
 *                     AHF
 *--------------------------------------------------*/

/* removes all variables related to the potential solver */
#define AHFlean                      /* will be deactivated for AHFpotcentre                                       */

/* decide whether you want to re-assign the density on all levels after generating all refinements */
//#define AHFdensrecovery              /* will be deactivated below for AHFcomcentre and activated for AHFpotcentre! */

/* whether or not to dump the *.AHF_substructure file    */
#define AHFsubstructure            /* dump substructure information (based upon halo-tree!) to file    */

/* whether or not to dump the *.AHF_disks file (in ASCII format at the moment) */
//#define AHFdisks

/* write catalogues in binary format (AHF_substructure and AHF_particlesSTARDUST not implemented yet)  */
/* (use AHFbinary2ascii.c to convert and merge MPI files!)                                             */
//#define AHFbinary

// only use DM particles to define density peaks
//#define AHFdmonlypeaks

/* miscellaneous flags */
#define AHFdmonly_Rmax_r2          /* base determination of Rmax and r2 upon DM only                 */
#define AHFparticle_Rmax_r2        /* base determination of Rmax and r2 upon sorted particle list    */
//#define AHFsorthalosbymass        /* this will sort the halos by mass when writing the output files */
//#define AHFaddDMonlyproperties    /* additionally calculates halo properties based upon DM only     */
//#define AHFnoremunbound           /* do not perform an unbinding procedure                          */
//#define AHFnoHubbleDrag           /* will not consider the Hubble term when unbinding               */
//#define AHFignore_ugas            /* ignores thermal energy of gas particles                        */
//#define AHFreducedinertiatensor   /* use reduced inertia tensor for shape deterimations             */
//#define AHFprofilerise            /* checks for rising profile in rem_outsideRvir() and chops halo  */
//#define AHFsplinefit              /* this will use a spline-interpolation to get Rmax, Vmax, and r2 (use with caution as it may not work for the low-mass haloes!) */
//#define AHFphspdens               /* add additional information about phase-space density to output */
//#define AHFvmbp                   /* interested in the velocity of the most bound particle?         */
//#define AHFundoPositionShiftAndScale /* dumps halo positions in coordinates found in input file     */
//#define AHFshellshape             /* base moment of inertia tensor on particles in shells           */
//#define AHFsplit_only             /* splits input simulation data into multiple files only          */
//#define AHFrestart                /* activate this when analysing output from AHFsplit_only         */
//#define AHFexciseSubhaloStars     /* this writes an additional _particlesSTARDUSTexcised file       */
//#define AHFnewHaloIDs             /* assigns a unique ID to each halo even across MPI tasks         */
//#define AHFcR1                    /* adds concentration calculation ala Wang et al. (arxiv:2310.00200) to _halos       */

/* restrict analysis to certain particles only */
//#define AHFhiresfocus             /* only keep particles of type [0], [1], and [4]                  */
//#define AHFptfocus  1             /* only keep particles of type 1                                  */
//#define AHFrfocus                 /* restrict analysis to a spherical region (specified in param.h) */

/* the following flags decide how to make Parent-Daughter Assignment in analyseRef() */
/* (you MUST specify one of them...)                                                 */
//#define PARDAU_DISTANCE            /* use sub-grid with closest distance to follow host  */
//#define PARDAU_NODES               /* use sub-grid with most nodes to follow host        */
#define PARDAU_PARTS               /* use sub-grid with most particles to follow host    */

/* the following centre flags determine how to calculate the halo centre... (you MUST specify one and only one of them...)                                        */
//#define AHFmaxdenscentre          /* use cell with maximum density as halo centre                       */
//#define AHFgeomcentre             /* use geometrical centre as halo centre                              */
//#define AHFpotcentre              /* use potential weighted centre as halo centre                       */
#define AHFcomcentre              /* use centre-of-mass of particles on finest refinement as halo centre */

/* the following flags are mainly used for debugging purposes */
//#define AHFgeom                   /* writes .AHF_halos.geom files                      */
//#define AHFgridinfofile           /* writes grid information into a new file           */
//#define AHFcentrefile             /* writes a file containing substructure info        */
//#define AHFcentrefileBASIC        /* writes a basic file containing substructure info  */
//#define AHFgridtreefile           /* file containing information about grid-tree       */






















/*=================================================================================
 *
 *
 *  the standard user of AHF does not need to touch any of the features below...
 *
 *
 *=================================================================================*/


/*--------------------------------------------------
 *               LEAVE AS IS, PLEASE...
 *--------------------------------------------------*/
/* resolve some sensible and important dependencies */
#if (defined AHFcomcentre || defined AHFgeomcentre)     // safe to switch off the density recovery in step() prior to AHF
#undef AHFdensrecovery
#endif

#ifdef AHFmaxdenscentre  // we need to have absolutely correct density values
#define AHFdensrecovery
#undef AHFlean
#endif

#ifdef AHFpotcentre // undefine all other options to determine the centre
#undef AHFmaxdenscentre
#undef AHFcomcentre
#undef AHFgeomcentre
#undef AHFlean
#define AHFdensrecovery
#endif

#ifdef AHFbinary
#undef AHFsubstructure // no binary format supported yet
#undef AHFdisks        // no binary format supported yet
#endif

/*--------------------------------------------------
 *            -DAHFfast
 *--------------------------------------------------*/
#ifdef AHFfast                                          // careful, these parallelisations may lead to inaccuracies
#define WITH_OPENMP2
#define WITH_OPENMP3
#undef AHFparticle_Rmax_r2
#endif

/*--------------------------------------------------
 *            -DAHFexciseSubhaloStars
 *--------------------------------------------------*/
#ifdef AHFexciseSubhaloStars
#ifndef METALHACK
#define METALHACK
#endif
#endif

/*--------------------------------------------------
 *           -DWITH_MPI or -DAHFrestart
 *--------------------------------------------------*/
#if (defined WITH_MPI || defined AHFrestart)
#  undef  REF_TEST
#  undef  VERBOSE
#  define AHFnewHaloIDs
#endif


/*--------------------------------------------------
 *                    -DAHFmixHaloIDandSnapID
 *--------------------------------------------------*/
#ifdef AHFmixHaloIDandSnapID
#ifndef AHFnewHaloIDs
#define AHFnewHaloIDs
#endif
#endif

/*--------------------------------------------------
 *                    -DSUSSING2013
 *--------------------------------------------------*/
// this flag simply leads to haloids of the type
//       haloid = snapid*1e12+ihalo
#ifdef SUSSING2013
#undef AHFnewHaloIDs  // SUSSING2013 overwrites AHFnewHaloIDs
#endif

// using SUSSING2013_particles will write a different _particles file
//  now also including information about the total energy and distance to the host
#ifdef SUSSING2013_particles
#ifndef SUSSING2013
#define SUSSING2013
#endif
#endif

/*--------------------------------------------------
 *                    -DGADGET
 *--------------------------------------------------*/
#ifdef GADGET2  // keep this flag for legacy reasons in case someone still counts on it...
#define GADGET
#endif
#ifdef GADGET
 #ifndef MULTIMASS
  #define MULTIMASS	/* avoid unnecessary compiler warnings */
 #endif
#define GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                    -DTIPSY
 *--------------------------------------------------*/
#ifdef TIPSY
#define  MULTIMASS
#define  GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                    -DDEVA
 *--------------------------------------------------*/
#ifdef DEVA2
#define DEVA
//#define DEVA2_QHULL_FILE "../snapshots/cube00.qhull"
#endif

#ifdef DEVA
#define MULTIMASS
#define GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                -DMARE_NOSTRUM
 *--------------------------------------------------*/
#ifdef MARE_NOSTRUM
#define MULTIMASS
#define GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                -DMETALHACK
 *--------------------------------------------------*/
#ifdef METALHACK
#	ifndef MULTIMASS
#		define MULTIMASS
#	endif
#	ifndef GAS_PARTICLES
#		define GAS_PARTICLES;
#	endif
#endif

/*--------------------------------------------------
 *                -DGAS_PARTICLES
 *--------------------------------------------------*/
#ifdef GAS_PARTICLES
#ifndef MULTIMASS
#define MULTIMASS
#endif
#endif

/*--------------------------------------------------
 *                     misc.
 *--------------------------------------------------*/
#ifdef NO_GAS
#undef GAS_PARTICLES
#endif

#ifdef VERBOSE2
#define VERBOSE
#endif

#ifdef VERBOSE
#define REF_TEST
#endif

#ifdef NCPUREADING_EQ_NFILES
#undef BCASTHEADER
#endif

#ifdef BCASTHEADER
#ifndef FOPENCLOSE
#define FOPENCLOSE
#endif
#endif

#ifdef AHFcentrefileBASIC
#ifndef AHFcentrefile
#define AHFcentrefile
#endif
#endif

/*--------------------------------------------
 * more transparent to read in source-code...
 *--------------------------------------------*/
#ifndef TSC
#ifndef CIC   /* forgotten to define mass assignemnt scheme ? => use TSC then... */
#ifndef NGP
#define TSC
#endif
#endif
#endif

#ifndef CONTINUE
#define TERMINATE
#define TERMINATE2  /* used in leavers.c */
#define VERBOSELOG
#endif /* CONTINUE */

#ifdef PERIODIC
#define PERIODIC_X
#define PERIODIC_Y
#define PERIODIC_Z
#endif /* PERIODIC */

#endif

