#ifndef PARAM_INCLUDED
#define PARAM_INCLUDED

#include <float.h>

/*================================================================================
 * switch on/off various features of AMIGA not passed via DEFINEFLAGS in Makefile
 *================================================================================*/
#include "define.h"

/*=============================================================================
* AHF related parameters
*=============================================================================*/
#define AHF_MINPART_GAS     10     /* min. number of gas for spin and shape calculation                                          */
#define AHF_MINPART_STARS   10     /* min. number of stars for spin and shape calculation                                        */
#define AHF_MINPART_SHELL   10     /* min. number of particles in a profile shell for using AHFshellshape                        */
#define AHF_NBIN_MULTIPLIER 1      /* (integer value) increases the number of radial bins by this factor from the standard value */
#define AHF_HIRES_DM_WEIGHT 1.0
#define AHF_HOSTHALOLEVEL   1      /* first level to be considered as credible to spawn subbaloes                                */
#define AHF_HOSTSUBOVERLAP  0.5    /* how far should the subhalo have entered into the host                                      */
#define AHF_MIN_REF_OFFSET  0      /* offset for first refinement to be used by AHF                                              */
#define AHF_RISE            1.00   /* Rho > AHF_RISE*Rho_prev -> rising density                                                  */
#define AHF_SLOPE           0.99   /* outer halo profile at least like r^-AHF_SLOPE                                              */
#define AHF_MAXNRISE        2      /* try to catch variations in density                                                         */
#define AHF_Rmax_r2_NIGNORE 5      /* how many central particle to ignore with AHFparticle_RMax_r2 feature                       */
#define PGAS                0.0   /* identifier for gas particles; has to be exactly 0.0!!!!                                     */
#define PDM                -1.0   /* identifier for dm particles; whatever negative value                                        */
#define PSTAR              -4.0   /* identifier for star particles; whatever negative value                                      */
#define PDMbndry           -5.0   /* identifier for bndry particles; whatever negative value                                     */

/* spherical region used with -DAHFrfocus (cf. define.h), units in Mpc/h */
#define AHFrfocusX 57.06044914
#define AHFrfocusY 52.61864923
#define AHFrfocusZ 48.70489744
#define AHFrfocusR 0.5


/*=============================================================================
 * some physical constants...
 *=============================================================================*/
#define Gyr       3.1558e16         /* [sec]                   */
#define Mpc       3.08567782e19     /* [km]                    */
#define H0        100.              /* [h*km]/[sec*Mpc]        */
#define rhoc0     2.7755397e11      /* [h^2*Msun]/[Mpc^3]      */
#define Grav      4.3006485e-9      /* [Mpc*km^2]/[Msun*sec^2] */
#define cH0	      2998.0		        /* c/H0 (in h^-1 Mpc)      */
#define kB_per_mp 0.825481286614E-2 /* [(km/sec)^2/K]          */
#define kBoltzman 6.9416792         /* [(km/sec)^2 Msun/K]     */
#define Msun      1.9891e30         /* [kg]                    */

/*============================================================================
 * grid parameters
 *============================================================================*/
#define MIN_NNODES 125            /* smallest grid-block 5x5x5     */

/*============================================================================
 * handling output file names etc.
 *============================================================================*/
#define MAXSTRING    2048   /* used for char statement, i.e. filenames etc.     */
#define AMIGAHEADER  2048   /* maximum size (in bytes) for output file header   */
#define HEADERSTRING 256    /* no. of characters for header string in outfiles  */
#define HEADERSIZE  (HEADERSTRING*sizeof(char)+2*sizeof(long)+6*sizeof(int)+46*sizeof(double))
#define FILLHEADER  (AMIGAHEADER-HEADERSIZE)

/*=============================================================================
 * general parameters
 *=============================================================================*/
#define TERMINATE_AMIGA {"terminateAMIGA"}
#define MAXTIME      10000   /* interpolation points for timeline (cf. tdef.h)   */
#define bytes2GB  9.313225746154785e-10

/*============================================================================
 * potential solver
 *============================================================================*/
#define DOMSWEEPS  10     /* number of domain and refinement sweeps...   10   */
#define REFSWEEPS  10     /* ...before checking for convergence          10   */
#define W_SOR      1.34   /* successive over-relaxation parameter        1.34 */
#define ETA        0.625  /* parameter for slow convergence              0.625*/  
#define CONVCRIT   0.1    /* convergence criterion                       0.1  */
#define DOMCORRECT 0.25   /* constrained criterion on domain grid        0.25 */
#define SPEEDFRAC  0.01   /* fraction of particles allowed to speed      0.01 */

/*=============================================================================
 * MPI parameters
 *=============================================================================*/
#ifdef MPI_DEBUG
#undef MAXTIME
#define MAXTIME 100
#endif

/** Number of bits used per dimension for the calculation of the Hilbert key */
#define BITS_PER_DIMENSION 21

/** The maximum amount of particles that can be send in on flush */
#define MAX_SEND_PARTICLES 1000000

/** Defines the verbosity level for the io_logging_* functions if
 * NEWSTARTUN is used. Depending on this value some messages might not
 * appear and hence this can be used to reduce the chatter produced by
 * the io_logging_* function. During the starup this will be passed to
 * the logging module (in startrun.c). The lower the number, the less
 * output will be produced. */
#define VERBOSITY 6

/*============================================================================
 * time stepping
 *============================================================================*/
#define NSTEPS       1000    /* number of (initial) steps for time stepping     */
#define CA_CRIT      0.15    /* restricts timestep due to da/a criterion    0.05*/ 
#define CELLFRAC_MAX 0.2     /* how far are particles allowed to move       0.2*/
#define CELLFRAC_MIN 0.05    /* how far should particles move at least      0.05*/
#define CF_MEAN      ((CELLFRAC_MAX+CELLFRAC_MIN)/2)

/*=============================================================================
 * finally some convenient abreviations...
 *=============================================================================*/
#define NDIM      3          /* DO NOT EVER TOUCH THIS NUMBER!  */
#define CRITMULTI 8.0
#define NP_RATIO  7.5
#ifdef DOUBLE 
#define ZERO         (1E-12)
#define MACHINE_ZERO (5E-16)
#else
#define ZERO         (1e-6)
#define MACHINE_ZERO (5e-16)
#endif
#define MZERO        (1e-10)  /* used when reading GADGET files */


/*=============================================================================
 * GADGET related parameters (only relevant for libio_serial.a !!!)
 *=============================================================================*/
#define GADGET_MUNIT 1.0e10    /* GADGET mass unit in Msol/h                   */
#ifdef GADGET_LUNIT_KPC
#define GADGET_LUNIT 1.0e-3
#else                          /* GADGET length unit in Mpc/h                  */
#define GADGET_LUNIT 1.0
#endif

/*============================================================================= 
 * no code is complete without defining these numbers ;-)
 *=============================================================================*/
#define PI    3.14159265358979323846264338
#define TWOPI 6.28318530717958647692528677
#define SQRT2 1.41421356237309504880168872

/*=============================================================================
 * The numbers 0, 1 and 2 are used to represent the X, Y and Z coordinates of
 * a three dimensional object (eg vector). To aid readability these numbers are
 * replaced by the names X, Y and Z when used individually.
 *=============================================================================*/
#define X 0     /* x-coord symbol */
#define Y 1     /* y-coord symbol */
#define Z 2     /* z-coord symbol */

/*=============================================================================
 * Boolean parameters
 *=============================================================================*/
#define YES   1
#define NO    0

#define ON    1
#define OFF   0

#define TRUE  1
#define FALSE 0

#endif

