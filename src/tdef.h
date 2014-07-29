#include <stdio.h>
#include <time.h>
#include "define.h"
#include "param.h"


#ifndef TDEF_INCLUDED
#define TDEF_INCLUDED

#	include <stdint.h>
#	include "libsfc/sfc.h"


/*--------------
 * Useful Types
 *--------------*/

typedef float  fvect[NDIM];      /* "vector" of floats  */
typedef double dvect[NDIM];      /* "vector" of doubles */
typedef char   boolean;          /* C has no boolean data type */
#ifdef DOUBLE
typedef double flouble;
#else
typedef float  flouble;
#endif

/*=========================================================================
 *  particle structure type
 *=========================================================================*/

/** The ID type */
typedef uint64_t part_id_t;
#	define PRIpartid PRIu64

typedef struct particle *partptr;
typedef struct particle
{
  /* linked-list coord  
   *    :- NOTE -: 
   * for the linkedlist sort to work this element *must* be first in the struct
   */
  partptr ll;
  
  
  flouble   pos[NDIM];                	/* position vector             */
  flouble   mom[NDIM];                	/* momentum vector             */
#ifdef MULTIMASS
  flouble   weight;                    /* weight of particle          */
#endif
    
  sfc_key_t sfckey;

#if (!(defined AHF_NO_PARTICLES && defined AHFlean))
  part_id_t id;
#endif
#if (defined GAS_PARTICLES)
  flouble    u;
#endif 
  
#if (defined METALHACK)
  flouble    z;
  flouble    age;
#endif
  
#ifdef DEVA
  long itype;
#endif
  
#ifdef SUSSING2013
  flouble E;
#endif
  
#ifdef STORE_MORE
  flouble rho;
  flouble eps;
#endif
}part;

/*=========================================================================
 *  gas particle structure type
 *=========================================================================*/
typedef struct gas *gasptr;
typedef struct gas
{
  /* positions and masses are stored along with "normal" particles... */
  flouble u;
}gas;


/*=========================================================================
 *  star particle structure type
 *=========================================================================*/
typedef struct star *starptr;
typedef struct star
{
  /* positions and masses are stored along with "normal" particles... */
  flouble dummy;
}star;



/*========================================================================
 * Grid structures and associated pointers
 *========================================================================*/

/*---------------------------------------------------
 * The node or gridpoint structure for refinements.
 *---------------------------------------------------*/
typedef struct node *nptr;
typedef struct node
{
#ifndef AHFlean
  flouble  pot;              	    /* gravitational potential              */
#endif
  flouble  dens;                  /* mass density                         */
  partptr  ll;              	    /* head of particle linked-list         */
  
  /*
   * The forces at the node are only required at the end of a timestep.
   * Hence the space allocated to the forces is used to store temporary
   * data at other times.
   */
  union forcestorage
 {
#ifndef AHFlean
  flouble  forces[NDIM];       /* grid force components                   */
  flouble  temp[NDIM];         /* old potential, source term, whatever... */
#endif
#ifdef AHFlean
  flouble  tempdens;
#endif
  partptr  new_ll;             /* pointer for rebuilding linked list      */
  int      colour;             /* used with -DAHF                         */
 }force;
  
}node;


/*-------------------------
 * QUAD logical structures
 *-------------------------*/

typedef struct nquad *nqptr;
typedef struct nquad
{
  nptr loc;                       /* pointer to first node in block*/
#ifndef AHFlean
  long x;                         /* x-coord of first node in block */
  long length;                    /* length of block of nodes */
#else
  int x;
  int length;
#endif
  nqptr next;               	  /* pointer to next nquad (or null) */
}nquad;

typedef struct cquad *cqptr;
typedef struct cquad
{
  nqptr loc;                      /* pointer to first nquad in block */
#ifndef AHFlean
  long y;                         /* y-coord of first nquad in block */
  long length;                    /* length of block of nquads */
#else
  int y;
  int length;
#endif
  cqptr next;               	  /* pointer to next cquad (or null) */
}cquad;

typedef struct pquad *pqptr;
typedef struct pquad
{
  cqptr loc;                	  /* pointer to first cquad in block */
#ifndef AHFlean
  long z;                         /* z-coord of first cquad in block */
  long length;                    /* length of block of cquads */
#else
  int z;
  int length;
#endif
  pqptr next;                     /* pointer to next pquad (or null) */
}pquad;

/*---------------------------
 * the grid structure itself
 *---------------------------*/

typedef struct gridlist
{
  pqptr   pquad;                 /* pointer to first pquad on grid           */
  long    no_pquad;              /* how many pquad's are there?              */
  pqptr  *pquad_array;           /* align pquads in memory when using OpenMP */
  
  long unsigned    l1dim;        /* virtual grid length                      */
  double  spacing;               /* grid spacing                             */
  double  spacing2;              /* grid spacing squared                     */
  double  critdens;              /* critical density on this grid            */
  double  masstodens;            /* mass to density conversion factor        */
  double  masstopartdens;        /* mass to partdens converion factor        */
  
  double  old_resid;             /* old residual                             */
  double  cur_resid;             /* current residual                         */
  double  trunc_err;             /* the truncation error                     */
  long    no_sweeps;             /* number of sweeps done on grid            */
  
  double  timecounter;           /* integration variable                     */
  double  timestep;              /* current timestep for integration variable*/
  int     multistep;             /* keep track of current multistep-phase    */
  
  struct  grid_size
 {
  long unsigned   no_part;     /* number of particles linked to grid       */
  long unsigned   no_nodes;    /* number of nodes actually present         */
 }size;
  
  struct  grid_time
 {
  time_t potential;             /* deriving the potenital by GS sweeps      */
  time_t density;               /* deriving the proper densities            */
  time_t DK;                    /* drifting and kicking particles           */
  time_t grid;                  /* everything related to grid hierarchy     */
  time_t hydro;                 /* time spent for the hydro-solver          */
 }time;
    
  boolean  next;                  /* the grid structure is never destroyed    */
  
#ifdef WITH_MPI
	/** Number of boundary nodes*/
	uint32_t num_bound;
	/** The boundary nodes */
	nptr bound;
#endif
  
}gridls;

/*=============================================================================
 * main:   structure carrying useful common information
 *=============================================================================*/
typedef struct info_global
{
  /* access to the domain grid */
  gridls *dom_grid;      /* actual pointer to domain grid strcuture             */
  int     domgrid_no;    /* locate domain grid within grid_list                 */
  
  /* access to all particles currently in use (NOTE: global.no_part is not necessarily equal to simu.no_part!) */
  partptr       fst_part;
  long unsigned no_part;       /* here we store the actual number of particles        */
  /* simu.no_part gives the total number of particles in the simulation while global.no_part counts the actual particles currently used */
  gasptr        fst_gas;
  long unsigned no_gas;        /* how many gas particles               */
  long          offset_gas;    /* for those input files that contain gas particles         */
  
  starptr       fst_star;
  long unsigned no_stars;      /* how many star particls               */
  long          offset_stars;  /* for those input files that contain gas particles         */
  
  
  /* access to current timecounter */
  double  a;
  double  t;
  double  z;
  double  super_t;
  int     no_timestep;
  
  long    fin_l1dim;     /* finest refinement level reached per domain cycle    */
  boolean fst_cycle;
      
  int     output_count;
  boolean restart;
  
  double  total_time;
  
  double  bytes_node;
  double  bytes_part;
  
  boolean       terminate;    /* flag that terminates code after current step     */
  char         *termfile_name;
  
  int           architecture; /* little or big endian machine                     */
  
  long unsigned ndummy;
  double        fdummy;
  
  
  double  ovlim;         /* virial overdensity parameter                             */
  double  rho_b;         /* stores comoving(!) background density                    */
  double  rho_vir;       /* stores comoving(!)   virial   density                    */
  double  max_ovdens;    /* maximum overdensity possible with currrent AMR hierarchy */
  
  int     ioflag;        /* flag that indicates that we wrote an output file         */
  
} info_global;

/*=============================================================================
 * tline:  structure carrying the functions a(t), omega(t), and lambda(t)
 *=============================================================================*/
typedef struct tline
{
  double age[MAXTIME];
  double hubble[MAXTIME];
  double omega[MAXTIME];
  double lambda[MAXTIME];
  double virial[MAXTIME];
  double growth[MAXTIME];
  double a[MAXTIME];
  double t[MAXTIME];
  double super_t[MAXTIME];
  double rhoc[MAXTIME];
}tline, *tlptr;

/*=============================================================================
 * simu:   structure carrying all sorts of information on the simulation
 *=============================================================================*/
typedef struct param_simu
{
  /* user supplied values */
  int           NGRID_DOM;
  int           NGRID_MAX;
  double        Nth_dom;
  double        Nth_ref;
  int           UseRhoBack;     /* flag indicating to use rho_b for halo edge defition      */
  double        UserDvir;       /* if >0 this will be the ovlim value used for halo edge    */
  double        MaxGatherRad;
  int           lb_level;
  double        AHF_VTUNE;
  int           AHF_MINPART;
  
  double        GADGET_m2Msunh;
  double        GADGET_l2Mpch;  // the GADGET units

  
  /* cosmological information */
  double        omega0;
  double        lambda0;
  tline         timeline;   
  
  
  /* timestepping information */
  double        a_initial;
  double        a_final;
  double        z_initial;
  double        z_final;
  double        t_initial;
  double        t_final;
  double        super_t_initial;
  double        super_t_final;
  
  
  /* unit stuff */
  double        mean_dens;
  double        FourPiG;       /* 4piG factor when using -DISOLATED    */
  double        boxsize;       /* the size of the cosmological box     */
  double        pmass;         /* mass of a single particle            */
  double        t_unit;        /* the time unit                        */
  
  double        l_unit;        /* currently boxsize and pmass are the units...          */
  double        m_unit;        /* ...but at some stage this should be clearly separated */
  
  double        pos_shift[3];  /* keep track of any shifts and scaling */
  double        pos_scale;     /* applied to the positions             */
  
  
  /* HYDRO information */
  double        gamma;        /* adiabatic index for equation-of-state */
  double        omegab;       /* baryonic matter content               */
  double        omegaDM;      /* dark matter content                   */
  double        f_b;          /* baryon fraction omegab/omega0         */
  double        H_frac;       /* mass fraction of molecular hydrogen   */
  double        T_init;       /* initial gas temperature               */
  double        B_init;
  double        e_init;       /* initial internal energy density       */
  
  
  /* numerical details */
  double        no_vpart;      /* virtual number of particles          */
  long unsigned no_part;       /* real (physical) number of particles  */
  long unsigned no_species;    /* number of different particle species */
  long unsigned no_gas;        /* how many gas particles               */
  long unsigned no_stars;      /* how many star particls               */
  
  
  long unsigned no_halos;      /* number of halos (used by AHF)  */
  
  int           NGRID_MIN;
  double        SHIFT;
  double        min_weight;    /* minimum particle mass in internal units     */ 
  double        max_weight;    /* maximum particle mass in internal units     */
  double        med_weight;
  int           np_limit;
  
  int           mmfocus;
  int           multi_mass;
  int           double_precision;
  int           hydro;
  int           magneto;
    
#ifdef AHF_LRSI
  double lrsi_beta;
  double lrsi_r_s;
#endif
  
#if (defined AHFmixHaloIDandSnapID || SUSSING2013)
  uint64_t isnap;
#endif
  
} param_simu;

/*=============================================================================
 * files:   structure carrying all information on output/input files
 *=============================================================================*/
typedef struct info_io
{
  /*------------------------------------------------------------------------
   * this first block is filled by copying the appropriate user_data!
   *------------------------------------------------------------------------*/
  
  
  /*------------------------------------------------------------------------
   * information for managing input and output files
   *   - filenames
   *   - pointer to logfile for global access
   *   - no. of output files and their respective expansion factors
   *------------------------------------------------------------------------*/
  FILE         *logfile;
  
  char         *icfile_name;           /* name of file with IC's             */
  char         *outfile_prefix;        /* prefix for output files            */
  char         *dumpfile_name;         /* name of dump file                  */
  char         *logfile_name;          /* name for log file                  */
  
  
  /*-------------------------------------------------------------------------------------------------
   * here we store the particles that are being read in from the input file
   *
   * fst_part is an array of length no_part that holds >>all<< particles
   * (NOTE: a particle structure only carries pos, mom, and weight!)
   *
   * fst_gas is an array of length no_gas that holds >>additional<< information for gas particles
   * (NOTE: the mass/weight is stored in fst_part[offset_gas -> offset_gas+no_gas-1]
   *
   * fst_star is an array of length no_stars that holds >>additional<< information for stars
   * (NOTE: the mass/weight is stored in fst_part[offset_gas+offset_stars -> offset_gas+offset_gas+no_stars-1]
   *------------------------------------------------------------------------------------------------*/
  
  
  /* here we store the particles read in from file (no_part is stored in io.header, too!) */
  partptr       fst_part;      // fst_part gives access to >>all<< particles
  long unsigned no_part;       // = no_DMpart + no_gas + no_stars
  
  gasptr        fst_gas;       // additional information other than mass for gas particles
  long unsigned no_gas;
  long          offset_gas;    // offset in fst_part[] to access gas particles
  
  starptr       fst_star;      // additional information other than mass for star particles
  long unsigned no_stars;
  long          offset_stars;  // offset in fst_part[] to access gas particles
  
  
  /*----------------------------------------------------------------------------
   *                            AMIGA file header
   *
   * here we store information that is relevant for 
   *  - restarting a simulation
   *  - obtaining information about the status of the simulation
   *    at the time the output file was written
   *
   *
   * this header gets written to file and is being constantly extended
   *
   * therefore, the ordering of variables has no physical meaning and
   * is more a reflection of the development of AMIGA than anything else
   *
   * if you require downwards compatibility, please do not change the ordering
   *----------------------------------------------------------------------------*/
  struct io_header_block
 {
  char          header[HEADERSTRING];    // a string at your disposal ;-)
  
  int           multi_mass;
  int           double_precision;        // remember previous settings
  
  long unsigned no_part;                 // this is the total number of particles (=DM+gas+stars)
  long unsigned no_species;
  double        no_vpart;                // information about particles (cf. min/max_weight below!)
  
  double        timestep;
  int           no_timestep;             // stage of simulation and last timestep used
  
  double        boxsize;                 // [length unit]
  double        omega0;
  double        lambda0;                 // information about the cosmological model
  
  double        pmass;                   // [mass unit]
  
  double        cur_reflevel;
  double        cur_frcres;              // information about currently finest refinement level
  
  double        a_initial;               // for historical reasons we store the expansion factor a
  double        a_current;               // rather than supercomoving time in the input/output file
  
  double        K_initial;
  double        K_current;
  double        U_initial;
  double        U_current;
  double        Eintegral;
  double        Econst;                  // needed for layzer_irvine() energy check
  
  double        paramNSTEPS;
  double        paramNGRID_DOM;
  double        paramNth_dom;
  double        paramNth_ref;
  double        paramE_UPDATE;
  double        paramCELLFRAC_MAX;
  double        paramCELLFRAC_MIN;
  double        paramCA_CRIT;
  double        paramMAX_L1DIM;
  double	      paramDOMSWEEPS;
  double        paramREFSWEEPS;          // technical details
  
  double        paramAHF_MINPART;
  double		    paramAHF_VTUNE;
  double        paramAHF_RISE;
  double        paramAHF_SLOPE;
  double        paramAHF_MAXNRISE;       // technical details about AHF
  
  double        min_weight;
  double        max_weight;
  double        t_unit;
  double        B_init;                  // relevant for MHD solver
  double        param_dummy5;
  double        param_dummy6;
  double        param_dummy7;
  double        param_dummy8;            // empty space that once was used (downwards compatibility!)
  
  double        version;
  int           build;                   // the code is under constant development ;-)
  
  double        omegab;
  double        gamma;
  double        H_frac;
  double        T_init;                  // relevant for the HYDRO solver
  
  int           hydro;
  int           magneto;
  
  double        med_weight;
  
  char          dummy[FILLHEADER];       // empty space for future additions (e.g. MHD, ...)
 } header;
  
} info_io;



/*======================================
 * definitions used by AHF
 *======================================*/
typedef struct {
  double x,y,z;
} XYZ;

typedef struct {
  double min,max;
} MINMAX;

/* this structure is used for DM, gas, and star properties only */
typedef struct {
  long unsigned npart;
  XYZ           pos_com;
  XYZ           pos_mbp;
  XYZ           vel;
  double        Mass;
  double        lambda;
  double        lambdaE;
  XYZ           AngMom;
  XYZ           axis;
  XYZ           E1;
  XYZ           E2;
  XYZ           E3;
  double        Ekin;
  double        Epot;
} SPECIESPROP;


/* the HALOPROFILE structure */
typedef struct {
  int             nbins;
  double         *r;
  long unsigned  *npart;
  double         *nvpart;
  double         *ovdens;
  double         *dens;
  double         *v2_circ;
  double         *v_esc2;
  double         *sig_v;
  double         *Ekin;
  double         *Epot;
  double         *Lx;
  double         *Ly;
  double         *Lz;
  double         *axis1;
  double         *E1x;
  double         *E1y;
  double         *E1z;
  double         *axis2;
  double         *E2x;
  double         *E2y;
  double         *E2z;
  double         *axis3;
  double         *E3x;
  double         *E3y;
  double         *E3z;
  
#ifdef GAS_PARTICLES
  double *M_gas;
  double *M_star;
  double *u_gas;
#ifdef AHFdisks
  double *k_gas;
  double *Ekin_gas;
  double *Lx_gas;
  double *Ly_gas;
  double *Lz_gas;
  double *axis1_gas;
  double *E1x_gas;
  double *E1y_gas;
  double *E1z_gas;
  double *axis2_gas;
  double *E2x_gas;
  double *E2y_gas;
  double *E2z_gas;
  double *axis3_gas;
  double *E3x_gas;
  double *E3y_gas;
  double *E3z_gas;
  double *k_star;
  double *Ekin_star;
  double *Lx_star;
  double *Ly_star;
  double *Lz_star;
  double *axis1_star;
  double *E1x_star;
  double *E1y_star;
  double *E1z_star;
  double *axis2_star;
  double *E2x_star;
  double *E2y_star;
  double *E2z_star;
  double *axis3_star;
  double *E3x_star;
  double *E3y_star;
  double *E3z_star;
#endif /* AHFdisks */
#endif /* GAS_PARTICLES */
  
#ifdef METALHACK
  double *z_gas;
  double *z_star;
#endif

#ifdef AHFphspdens
  double *sigma2_vx_sh;
  double *sigma2_vy_sh;
  double *sigma2_vz_sh;
  double *sigma2_vr_sh;
  double *sigma2_vtheta_sh;
  double *sigma2_vphi_sh;
#ifdef AHFmeanvelocities
  double *mean_vx_sh;
  double *mean_vy_sh;
  double *mean_vz_sh;
  double *mean_vr_sh;
  double *mean_vtheta_sh;
  double *mean_vphi_sh;
  double *mean_vx_sp;
  double *mean_vy_sp;
  double *mean_vz_sp;
  double *mean_vr_sp;
  double *mean_vtheta_sp;
  double *mean_vphi_sp;
#endif /* AHFmeanvelocities */
#endif /* AHFphspdens */
} HALOPROFILE;


/* the HALO structure */
typedef struct {   
  /*--------------------------------
   * access to the halo's particles
   *--------------------------------*/
  long unsigned  npart;       /* total number of particles                        */
  long unsigned *ipart;       /* indizes of halo particles                        */
  long unsigned  nll;         /* temporary storage of no. of parts attached to ll */
  partptr        ll;          /* spatialRef2halos() still requires linked-lists   */

#ifdef AHFexciseSubhaloStars
  long unsigned  npart_uniquestars;
  long unsigned *ipart_uniquestars; /* the stars only belonging to and only to the halo */
  double         Mstar_excised;
  double         mean_z_star_excised;
#endif
    
  /*----------------------------
   * integral properties
   *----------------------------*/
  XYZ     pos;           /* position                                                     */
  XYZ     vel;           /* bulk velocity                                                */
  double  M_vir;         /* mass                                                         */
  double  R_vir;         /* radius                                                       */
  double  sigV;          /* velocity dispersion                                          */
  double  v_esc2;        /* escape velocity at R_vir                                     */
  double  V2_max;        /* peak of rotation curve                                       */
  double  R_max;         /* position of rotation curve peak                              */
  double  r2;            /* position in density profile where rho*r^2 peaks              */
  double  ovdens;        /* over-density at R_vir                                        */
  double  lambda;        /* spin paramerer ala Bullock et al. (2001)                     */
  double  lambdaE;       /* energy based spin parameter ala Peebles                      */
  double  Ekin;          /* total kinetic energy                                         */
  double  Epot;          /* total potential energy                                       */
  double  SurfP;         /* surface pressure term ala Shaw et al. (2006)                 */
  double  Phi0;          /* potential normalisation constant                             */
  XYZ     pos_com;
  double  com_offset;    /* offset of centre-of-mass     to halo centre                  */
  XYZ     pos_mbp;       /* position of most bound particle                              */
  XYZ     vel_mbp;       /* velocity of most bound particle                              */
  double  mbp_offset;    /* offset of most-bound-particle to halo centre                 */
  XYZ     AngMom;        /* angular momentum vector                                      */
  XYZ     axis;
  XYZ     E1;            /* all the moment of inertia stuff                              */
  XYZ     E2;
  XYZ     E3;
  double  fMhires;       /* fraction of hi-res to all DM particles                       */ 
  double  cNFW;          /* NFW concentration according to Eq.(9) in Prada et al. (2012) */
  double  R_edge;        /* ignore this quantity!                                        */

#ifdef METALHACK
  double mean_z_gas;     /* mean metallicity of the gas particles                        */
  double mean_z_star;    /* mean metallicity of the star particles                       */
#endif    

#ifdef GAS_PARTICLES
#ifndef AHFaddDMonlyproperties
  SPECIESPROP DM_only;
#endif
  SPECIESPROP gas_only;  /* some integral quantities based upon gas particles only       */
  SPECIESPROP stars_only;/* some integral quantities based upon star particles only      */
#endif
  
  
  /*----------------------------
   * the halo profile
   *----------------------------*/
  HALOPROFILE prof;
  
  
  /*----------------------------
   * halo-tree information
   *----------------------------*/  
  int hostHalo;          /* Index to the halo that hosts this one based on the grids          */
  int hostHaloLevel;     /* keeps track on which level the host halo was assigned             */
  int numSubStruct;      /* number of subhalos                                                */
  int	*subStruct;	  	   /* indizes of the substructure halos                                 */
  
  double gatherRad;      /* radius used to gather initial set of particles                    */
  double spaRes;         /* The spatial resolution of the halo                                */			
  int    refLev;         /* The deepest refinement level:  0 - is the what we start AHF on    */			
  
  int    numNodes;       /* number of nodes in this halo                                      */
  
#	if (defined WITH_MPI || defined AHFrestart)
  boolean ignoreme;
#	endif
  
#ifdef AHFnewHaloIDs
  /*--------------------------------------------
   * unique IDs to be used for the output files
   *--------------------------------------------*/
  uint64_t haloID;
  uint64_t hostHaloID;
#endif
  
} HALO;

typedef struct {
  HALO    *halos;
  uint64_t nhalos;
} halo_s;

typedef struct {
  int x,y,z;
} intXYZ;

/* An index array to save computation */
typedef struct {
  int    refLevel;
  int    isoRefIndex;
  intXYZ periodic;
} SRINDEX;


/*=============================================================================
 * in ahf we are going to store some global information about AHF
 *=============================================================================*/
typedef struct info_ahf
{
  int       no_grids;    /* Number of grids used in AHF    */
  int       min_ref;     /* coarsest grid to be considered */
  
  /* timing information */
  time_t    time;  
} info_ahf;

/*=============================================================================
 * in timing we are going to store some global information about timings
 *=============================================================================*/
typedef struct info_timing
{
  /* timing information */
  time_t    io;
  time_t      startrun;
  time_t      ptfocus;
  time_t      rfocus;
  time_t      loadbalance;
  time_t      sfckey;
  time_t      distribution;
  
  // AHF1 only timings
  time_t    gendomgrids;
  time_t    ll;
  time_t    genrefgrids;
  time_t    densrecovery;
  time_t    potcentre;
  time_t    ahf_gridinfo;
  time_t      RefCentre;
  time_t      analyseRef;
  time_t      spatialRef2halos;

  // AHF2 only timings
  time_t      patchtree2halos;
  
  // AHF1 and AHF2
  time_t      ahf_halos;
  time_t      ahf_halos_sfc_constructHalo;
  time_t      ahf_io;
  time_t      generate_tree;
#ifdef AHF2
  double     generate_tree_v2;
#endif
  
} info_timing;



/*==============================================================*/
/* the periodic cubes list for the lightcone                    */
/*              (used with -DLIGHTCONE)                         */
/*==============================================================*/
typedef struct cubes *cubeptr;
typedef struct cubes {
  int I,J,K;
  int used;
} cube;


/*=============================================================================
 * energy:   structure carrying all information on energy conservation
 *=============================================================================*/
typedef struct info_layzer_irvine
{
  double  K_initial, K_current, Kold;
  double  U_initial, U_current, Uold;
  double  aold;
  double  integral;
  double  econst;
  double  echeck;
  
  time_t  time;
} info_energy;

/*=============================================================================
 * powerspectrum:   structure carrying all information about on-the-fly P(k)
 *=============================================================================*/
typedef struct info_Pk
{
  double    *rk;
  double    *Pk;
  
  double    Pk_ini;
  double    Pk_now;
  
  double    Dgrowth_ini;
  double    Dgrowth_now;
  
  int       dump_Pk;
  
  time_t    time;
} info_Pk;

/*=============================================================================
 * info_gadget:   structure carrying all information about GADGET files
 *=============================================================================*/
#define SIZEOFGADGETHEADER 256    /* Size of GADGET header in bytes */
typedef struct info_gadget
{
  int        no_gadget_files;
  int        i_gadget_file;
  int      *(np[6]);
  int        nall1;
  partptr    lst_part;
  gasptr     lst_gas;
  
#ifdef LGADGET
  long long  IDmin;
  long long  IDmax;
#else
  int        IDmin;
  int        IDmax;
#endif
  double     mmin;
  
  gasptr     fst_gas;
  
  starptr    fst_star;
  
  struct io_gadget_header
 {
  int      np[6];
  double   massarr[6];
  double   expansion;
  double   redshift;
  int      flagsfr;
  int      flagfeedback;
  int      nall[6];
  int      flagcooling;
  int      NumFiles;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     unused[SIZEOFGADGETHEADER- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
 } header; 
  
} info_gadget;

/*==============================================================*/
/* Variables for reading in the TIPSY format */
/*===============================================================*/
#ifdef TIPSY
double tipsy_boxsize,tipsy_omega0,tipsy_lambda0,tipsy_initalz,tipsy_currentimeno;
#define MAXDIM 3
#define forever for(;;)

typedef float Real;

struct gas_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real rho;
  Real temp;
  Real hsmooth;
  Real metals ;
  Real phi ;
} ;

struct gas_particle *gas_particles;

struct dark_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real eps;
  Real phi ;
} ;

struct dark_particle *dark_particles;

struct star_particle {
  Real mass;
  Real pos[MAXDIM];
  Real vel[MAXDIM];
  Real metals ;
  Real tform ;
  Real eps;
  Real phi ;
} ;

struct star_particle *star_particles;

struct tipsy_dump {
  double time;
  int    nbodies;
  int    ndim;
  int    nsph;
  int    ndark;
  int    nstar;
  
  /* Jeremy Balin addition */
  char align[ (32 - sizeof(double) - 5 * sizeof(int)) / sizeof(char) ];
  /* total size should be 32 bytes to make alignment okay with
   * 64-bit architectures (ie. alphas) */
  
} ;
#endif /* TIPSY */


#endif

