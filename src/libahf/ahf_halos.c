/*==========================================================================
 *
 * This is a re-write from the bottom up of AHF!!
 * In this file I have merged _centres with _halos
 *
 * 08.02.2009: small adaptations for working with HaloTracker as well
 *
 *
 *==========================================================================*/

#ifdef AHF

/***************************************************************************
 *   Includes                                                              *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "ahf_halos.h"
#include "ahf_io.h"
#include "../libutility/utility.h"
#include "../libutility/timer.h"
#include "../libamr_serial/amr_serial.h"

/* Use the new particle finding methods */
#  include "ahf_halos_sfc.h"

#ifdef AHF_SQL
#  include "ahf_io_sql.h"
#endif

#ifdef EXTRAE_API_USAGE
#include <extrae_user_events.h>
#endif

/***************************************************************************
 *   Macros and assorted definitions                                       *
 ***************************************************************************/

/* This puts coordinate B in the restframe of A */
#define RESTFRAME(A, B) \
    ((A)-(B) > 0.5 ? (A)-(B)-1. : ((A)-(B) < -0.5 ? (A)-(B)+1. : (A)-(B)))


/***************************************************************************
 *   Structure and type definitions                                        *
 ***************************************************************************/
/* General */
typedef struct {
	double x, y, z, norm;
} XYZnorm;

typedef struct {
	int refLevel, isoRefIndex;
} INDEX;

typedef struct {
	double eigenVectors[3][3];
	double eigenValues[3];
	double e1, e2;
	double T;
} INERTIA;

typedef struct {
	XYZ    L;
	XYZ    pos;
	double dens;
	int    colour;
} GRIDATACOL;

/* Spatial Refinement */
struct spatREF {
	int    colour;                 /* The colour of the isolated refinement                 */
	int    numNodes;               /* The number of nodes in the refinement                 */
	double volume;                 /* Volume of the refinement                              */
	double maxDens;                /* Maximum density of node on the refinement             */

	XYZ           boundRefDiv;     /* The divider to calculate the min-max                  */

	XYZ           centre;          /* The actual centre to be used as halo centre           */
	XYZnorm       centreCMpart;    /* The centre of mass of particles on spatialRef[][]     */
	XYZnorm       centreGEOM;      /* The geometric centre of spatialRef[][]                */
	XYZnorm       centreDens;      /* The density weighted centre of spatialRef[][]         */
	XYZnorm       centrePot;       /* The potential weighted centre of spatialRef[][]       */
	MINMAX        x, y, z;         /* the extent of the refinement	                      	*/

	long unsigned numParts;        /* Number of particles in the refinement                 */
	partptr       ll;              /* head of particle linked-list                          */

	int           numSubStruct;    /* Number of sub halos                                   */
	INDEX         *subStruct;      /* Indexes of the substructure halos                     */

	int           numParDom;       /* Number of partent domains                             */
	INDEX         *parDom;         /* Indexes of the parent domains                         */

	INDEX         daughter;        /* Indexes of the daughter refinemnt                     */

	int           haloIndex;       /* Which halo do I belong to ?                           */

	double        closeRefDist;    /* Distance to the closest refinement                    */

	/*  GRIDATA	*nodedata;	 The node data ------ not sure if will include?
	 *  INERTIA	inertia;	 unity inertia tensor
	 *  INERTIA	winertia;	 The Weighted inertia tensor
	 *  double	dx,dy,dz;	 Length of the REf in the x,y,z dir
	 **/
};
typedef struct spatREF SPATIALREF;


/***************************************************************************
 *   Prototypes of local functions                                         *
 ***************************************************************************/
/* General functions */
int nodeCompare(const void *, const void *);
int refCompare(const void *, const void *);
int haloCompare(const void *, const void *);

/* Specific functions */
int     RefCentre                   (gridls *, int, SPATIALREF **);
int     analyseRef                  (int, SPATIALREF **);
int     spatialRef2halos            (int, SPATIALREF **);
int     WriteGridtreefile           (const char *, int, SPATIALREF **);
double  get_haloedge                (double *, double *, int, int);
int     compare                     (struct particle *, struct particle *, XYZ *);
boolean check_subhalo               (HALO *host, HALO *sub);
void    exciseSubhaloStars          (HALO *halos, long ihost);

/***************************************************************************
 *   Global variables                                                      *
 ***************************************************************************/
SRINDEX        *spatialRefIndex;
int            totnumIsoRef; /* Number of spatially isolated refinements */
int            *numIsoRef;   /* Number of spatially isolated refinements */

static int     numHalos;
static HALO    *halos;

static int     removecount = 0;
static int     movecount   = 0;
static int     numDensZero = 0;
static int     numPartZero = 0;
static SRINDEX *densZero;

/* Arrays */
static double *gridl1dim;

/* Inertia tensor */
static double matrix[NDIM][NDIM];        /* The input inertia tensor    */
static double eigenv[NDIM][NDIM];        /* The calculates eigenvectors	*/

/* conversion factors etc. */
double r_fac, x_fac, v_fac, m_fac, rho_fac, phi_fac, u_fac, Hubble;


/***************************************************************************
 *   Implementation of exported functions                                  *
 ***************************************************************************/
void
ahf_halos(gridls *grid_list)
{
	long unsigned *idxtmp, *idx, ipart, ifrac, nhalos;
	double        *fsort;
	int           i, j, k, l, slen, kstar;
	FILE          *fout;
	char          filename[MAXSTRING], fprefix[MAXSTRING],
	              file_no[MAXSTRING];
	partptr       current, cur_part;
	double        r, g, b;
	double        rad, t_relax, age;
	int           ibin, r_conv_i;
	double        a, a2, a3, rho_crit, omega, lambda, ovlim, rho_b, rho_vir;
	int           numGoodHalos;
	double        Xhost, Yhost, Zhost, Xsat, Ysat, Zsat, dX, dY, dZ, Rhost2,
	              Dsat2;
	long          nsat, *isat;
	SPATIALREF    **spatialRef;
  int           numSubStruct, *SubStruct, ihost, isub;

#ifdef EXTRAE_API_USAGE
  Extrae_user_function(1);
#endif

	/* pointless trying to find halos when there are no useable grids */
	if (ahf.no_grids <= 0)
		return;

	/* conversion factors to get physical units */
	r_fac   = simu.boxsize * global.a;
	x_fac   = simu.boxsize;
	v_fac   = simu.boxsize / simu.t_unit / global.a;   //NOTE: AHF stores a^2 \dot{x} as velocity and hence we need to divide by a !!!!
	m_fac   = simu.pmass;
  u_fac   = pow2(simu.boxsize / simu.t_unit);
	rho_fac = simu.pmass / pow3(simu.boxsize);
	phi_fac = Grav * simu.pmass / (simu.boxsize * global.a);

	/* cosmology related stuff */
	a        = global.a;
	a3       = pow3(a);
	omega    = calc_omega(a);
	lambda   = calc_lambda(a);
	ovlim    = calc_virial(a);
	Hubble   = calc_Hubble(a);              /* in km/sec/Mpc */
	rho_crit = a3* calc_rho_crit(a);        /* comoving(!) critical density   */
	rho_b    = omega * rho_crit;            /* comoving(!) background density */
  rho_vir  = a3*calc_rho_vir(a);          /* comoving(!) density used to normalize densities */
  
	/* Properly store those values in the global structure */
	global.ovlim   = ovlim;
	global.rho_b   = rho_b;
	global.rho_vir = rho_vir;

#ifdef VERBOSE
	fprintf(stderr, "#################### ahf_halos ####################\n");
	fprintf(stderr, "AHF_MINPART  = %d\n", simu.AHF_MINPART);
	fprintf(stderr, "AHF_VTUNE    = %g\n", simu.AHF_VTUNE);
	fprintf(stderr, "z            = %g\n", 1./a-1.);
	fprintf(stderr, "Omega(z)     = %g\n", omega);
	fprintf(stderr, "OmegaL(z)    = %g\n", lambda);
	fprintf(stderr, "rho_crit(z)  = %g\n", rho_crit);
	fprintf(stderr, "rho_back(z)  = %g\n", rho_b);
	fprintf(stderr, "rho_vir(z)   = %g (actual normalisation density)\n", rho_vir);
	fprintf(stderr, "Delta_vir(z) = %g\n", ovlim);
	fprintf(stderr, "Hubble(z)    = %g\n", Hubble);
#endif
	fprintf(io.logfile,
	        "#################### ahf_halos ####################\n");
	fprintf(io.logfile, "AHF_MINPART  = %d\n", simu.AHF_MINPART);
	fprintf(io.logfile, "AHF_VTUNE    = %g\n", simu.AHF_VTUNE);
	fprintf(io.logfile, "z            = %g\n", 1./a-1.);
	fprintf(io.logfile, "Omega(z)     = %g\n", omega);
	fprintf(io.logfile, "OmegaL(z)    = %g\n", lambda);
	fprintf(io.logfile, "rho_crit(z)  = %g\n", rho_crit);
	fprintf(io.logfile, "rho_back(z)  = %g\n", rho_b);
	fprintf(io.logfile, "rho_vir(z)   = %g (actual normalisation density)\n", rho_vir);
	fprintf(io.logfile, "Delta_vir(z) = %g\n", ovlim);
	fprintf(io.logfile, "Hubble(z)    = %g\n", Hubble);
	fflush(io.logfile);
  

	/* prepare output filenames */
#  ifdef WITH_MPI
	snprintf(fprefix, MAXSTRING, "%s.%04d.", global_io.params->outfile_prefix, global_mpi.rank);
#  else
	snprintf(fprefix, MAXSTRING, "%s.", global_io.params->outfile_prefix);
#  endif
#  ifdef AHFDEBUG
	io_logging_msg(global_io.log, INT32_C(0), "Using %s as the file prefix for AHF outputs.", fprefix);
#  endif

	/* OUTPUT CONVENTION: 3 digits*/
	sprintf(file_no, "z%.3f", fabs(global.z));
	strcat(fprefix, file_no);


	/**************************************************************************
	 * Allocating the memory for the isolated refinements array
	 * This information is only a temporary holder and can be set free
	 * once the refinements have been joined to form halos
	 ***************************************************************************/
	spatialRef = NULL;
	if ((spatialRef = calloc(ahf.no_grids, sizeof(SPATIALREF *))) == NULL) {
		fprintf(stderr, "Error in allocating the memory for spatialRef array\n");
		exit(0);
	}
	for (i = 0; i < ahf.no_grids; i++) {
		spatialRef[i] = NULL;
		if ((spatialRef[i] = calloc(numIsoRef[i], sizeof(SPATIALREF))) == NULL) {
			fprintf(stderr, "Error in allocating the memory for spatialRef array\n");
			exit(0);
		}
	}

	/* Initialising the spatialRef array */
	for (i = 0; i < ahf.no_grids; i++) {
		for (j = 0; j < numIsoRef[i]; j++) {
			spatialRef[i][j].numNodes             = 0;
			spatialRef[i][j].maxDens              = -1.0;

			spatialRef[i][j].centreDens.x         = 0.0;
			spatialRef[i][j].centreDens.y         = 0.0;
			spatialRef[i][j].centreDens.z         = 0.0;
			spatialRef[i][j].centreDens.norm      = 0.0;

			spatialRef[i][j].centrePot.x          = 0.0;
			spatialRef[i][j].centrePot.y          = 0.0;
			spatialRef[i][j].centrePot.z          = 0.0;
			spatialRef[i][j].centrePot.norm       = 0.0;

			spatialRef[i][j].centreGEOM.x         = 0.0;
			spatialRef[i][j].centreGEOM.y         = 0.0;
			spatialRef[i][j].centreGEOM.z         = 0.0;
			spatialRef[i][j].centreGEOM.norm      = 0.0;

			spatialRef[i][j].centreCMpart.x       = 0.0;
			spatialRef[i][j].centreCMpart.y       = 0.0;
			spatialRef[i][j].centreCMpart.z       = 0.0;
			spatialRef[i][j].centreCMpart.norm    = 0.0;

			spatialRef[i][j].x.min                = 100000.0;
			spatialRef[i][j].x.max                = -100000.0;
			spatialRef[i][j].y.min                = 100000.0;
			spatialRef[i][j].y.max                = -100000.0;
			spatialRef[i][j].z.min                = 100000.0;
			spatialRef[i][j].z.max                = -100000.0;

			spatialRef[i][j].haloIndex            = -1;
			spatialRef[i][j].daughter.refLevel    = -1;
			spatialRef[i][j].daughter.isoRefIndex = -1;

			spatialRef[i][j].numParts             = 0;
			spatialRef[i][j].ll                   = NULL;
			spatialRef[i][j].numSubStruct         = 0;
			spatialRef[i][j].subStruct            = NULL;
			spatialRef[i][j].numParDom            = 0;
			spatialRef[i][j].parDom               = NULL;
			spatialRef[i][j].closeRefDist         = -1.0;

			spatialRef[i][j].boundRefDiv.x        = -1.0;
			spatialRef[i][j].boundRefDiv.y        = -1.0;
			spatialRef[i][j].boundRefDiv.z        = -1.0;
		}
	}

	/**************************************************************************
	 * Collect the isolated refinements */
#ifdef VERBOSE
	fprintf(stderr, "\nCollecting the isolated refinements ... ");
	fprintf(io.logfile, "\nCollecting the isolated refinements ... ");
	fflush(io.logfile);
#endif
  timing.RefCentre -= time(NULL);
	if (RefCentre(grid_list, ahf.no_grids, spatialRef) == 0) {
		fprintf(stderr, "Stuffed up collecting the grid data\n");
		exit(-1);
	}
  timing.RefCentre += time(NULL);
#ifdef VERBOSE
	fprintf(stderr, "done\n");
	fprintf(io.logfile, "done\n");
	fflush(io.logfile);
#endif

	/**************************************************************************
	 * Analyse the isolated refinements */
#ifdef VERBOSE
	fprintf(stderr, "\nAnalysing the isolated refinements:\n");
	fprintf(io.logfile, "\nAnalysing the isolated refinements:\n");
	fflush(io.logfile);
#endif
  timing.analyseRef -= time(NULL);
	if (analyseRef(ahf.no_grids, spatialRef) == 0) {
		fprintf(stderr, "Stuffed up analysing the isolated refinements\n");
		exit(-1);
	}
  timing.analyseRef += time(NULL);
  
  timing.generate_tree += time(NULL);

#ifdef AHFgridtreefile
  /**************************************************************************
   * Writing the spatialRef[][] tree to a file */
  fprintf(stderr, "\nWriting the spatialRef[][] tree to a file\n");
  fprintf(io.logfile, "\nWriting the spatialRef[][] tree to a file\n");
  fflush(io.logfile);
  
  WriteGridtreefile(fprefix, ahf.no_grids, spatialRef);

  //exit(0);
#endif

	/**************************************************************************
	 * Connect the isolated refinements */
#ifdef VERBOSE
	fprintf(stderr, "\nConverting the spatialRef[][] tree to a halos[] array:\n");
	fprintf(io.logfile, "\nConverting the spatialRef[][] tree to a halos[] array:\n");
	fflush(io.logfile);
#endif
  timing.spatialRef2halos -= time(NULL);
	if (spatialRef2halos(ahf.no_grids, spatialRef) == 0) {
		fprintf(stderr, "Stuffed up converting spatialRef[][] to halos[] \n");
		exit(-1);
	}
  timing.spatialRef2halos += time(NULL);

  
#ifdef AHFwritePreliminaryHalos
 {
  FILE *fp;
  char outfile[MAXSTRING];
  long unsigned *idx, *idxtmp;
  double        *fsort;
  
  // sort haloes by npart
  idx    = (long unsigned *)calloc(numHalos, sizeof(long unsigned));
  idxtmp = (long unsigned *)calloc(numHalos, sizeof(long unsigned));
  fsort  = (double *)       calloc(numHalos, sizeof(double));
  for (i = 0; i < numHalos; i++)
    fsort[i] = (double)halos[i].npart;
  indexx(numHalos, fsort-1, idxtmp-1);
  
  /* indexx sorts ascending and gives indizes starting at 1 */
  for (i = 0; i < numHalos; i++)
    idx[numHalos - i - 1] = idxtmp[i] - 1;
  free(idxtmp);
  free(fsort);
  
  sprintf(outfile,"%s.AHF_preliminaryhalos",fprefix);
  fp = fopen(outfile,"w");
  
  fprintf(fp,"%d\n",numHalos);
  for(k=0; k<numHalos; k++) {
    i = idx[k];
    fprintf(fp,"%lu %18.14lf %18.14lf %18.14lf %18.14lf %18.14lf %18.14lf %d %d %d %d %d\n",
            halos[i].npart,
            halos[i].pos.x,
            halos[i].pos.y,
            halos[i].pos.z,
            halos[i].gatherRad,
            halos[i].R_vir,
            halos[i].spaRes,
            halos[i].refLev,
            halos[i].numNodes,
            halos[i].hostHalo,
            halos[i].hostHaloLevel,
            halos[i].numSubStruct);
    
#ifdef HaloTree
    if(halos[i].numSubStruct > 0) {
      for(k=0; k<halos[i].numSubStruct; k++) {
        fprintf(fp,"  %d\n",halos[i].subStruct[k]);
      }
    }
#endif
  }
  free(idx);
  
  fclose(fp);
  //exit(0);
 }
#endif

	/**************************************************************************
	 * Freeing the spatialRef array ALEX */
	for (i = 0; i < ahf.no_grids; i++) {
    for (j = 0; j < numIsoRef[i]; j++) {
      if(spatialRef[i][j].subStruct) free(spatialRef[i][j].subStruct);
      if(spatialRef[i][j].parDom)    free(spatialRef[i][j].parDom);
    }
		free(spatialRef[i]);
   }
	free(spatialRef);


#ifndef AHFphi_infty
	{
		/* the grids are no longer needed => free memory */
		gridls *for_grid;
		int    no_grids;

    no_grids = ahf.no_grids    + ahf.min_ref - 1;
		for_grid = global.dom_grid + no_grids;
#    ifdef VERBOSE
		fprintf(stderr, "\nFree'ing all grid structures:\n");
		fprintf(io.logfile, "\nFree'ing all grid structures:\n");
#    endif
		while (no_grids >= 0) {
#    ifdef VERBOSE
			fprintf(stderr, "  free'ing grid %ld\n", for_grid->l1dim);
			fprintf(io.logfile, "  free'ing grid %ld\n", for_grid->l1dim);
			fflush(io.logfile);
#    endif
			free_grid(for_grid, &no_grids);
			for_grid--;
		}
	}
#endif   /* AHFphi_infty */



  
  /**************************************************************************/
  /*                CALCULATE ALL RELEVANT HALO PROPERTIES                  */
	/**************************************************************************/
  timing.ahf_halos_sfc_constructHalo -= time(NULL);
#ifdef VERBOSE
	fprintf(stderr, "\nConstructing Halos (%d)\n", numHalos);
	fprintf(stderr,   "===================================\n\n");
	fprintf(io.logfile, "\nConstructing Halos (%d)\n", numHalos);
	fprintf(io.logfile,   "===================================\n\n");
	fflush(io.logfile);
#endif

#ifdef WITH_OPENMP
#  pragma omp parallel for schedule (dynamic)	shared(halos, numHalos) private(i)
#endif
  /* this construct integral as well as profile properties of each halo */
	for (i = 0; i < numHalos; i++) {
		ahf_halos_sfc_constructHalo(halos + i);
	}  
  
#if (defined WITH_MPI || defined AHFrestart)
  /**************************************************************************/
  /*            FLAG/REMOVE HALOES NOT IN ACTUAL MPI DOMAIN                 */
	/**************************************************************************/
  timing.ahf_halos_sfc_constructHalo -= time(NULL);
#ifdef VERBOSE
  //	fprintf(stderr, "Removing haloes not in MPI domain\n");
  //	fprintf(stderr, "=================================\n\n");
  //	fprintf(io.logfile, "Removing haloes not in MPI domain\n");
  //	fprintf(io.logfile, "=================================\n\n");
  //	fflush(io.logfile);
	fprintf(stderr, "Flagging haloes not in MPI domain\n");
	fprintf(stderr, "=================================\n\n");
	fprintf(io.logfile, "Flagging haloes not in MPI domain\n");
	fprintf(io.logfile, "=================================\n\n");
	fflush(io.logfile);
#endif
	/* Removing the haloes in the boundary zone */
	//rem_boundary_haloes();
  
	/* Flagging the haloes in the boundary zone */
	flag_boundary_haloes();
  
  timing.ahf_halos_sfc_constructHalo += time(NULL);
#endif

	/**************************************************************************/
  /*                     RE-HASH SUBHALO INFORMATION                        */
  /*          (does not properly function for MPI version yet)              */
	/**************************************************************************/
#ifdef VERBOSE
	fprintf(stderr, "Re-hashing substructure information using final radii\n");
	fprintf(stderr, "=====================================================\n\n");
	fprintf(io.logfile, "Re-hashing substructure information using final radii\n");
	fprintf(io.logfile, "=====================================================\n\n");
	fflush(io.logfile);
#endif
#ifdef WITH_OPENMP
#  pragma omp parallel for schedule (dynamic)	shared(halos,numHalos,stderr) private(i,k,numSubStruct,SubStruct,ihost,isub) default(none)
#endif
	for (i = 0; i < numHalos; i++) {
    
    /* deal with substructure IDs and lists */
    numSubStruct = 0;
    SubStruct    = NULL;
    for(k=0; k<halos[i].numSubStruct; k++)
     {
      /* quick-and-easy access to host and subhalo in halos[] */
      isub  = halos[i].subStruct[k];
      ihost = halos[isub].hostHalo; // this should be identical to i !?
      
//      if(ihost != i) fprintf(stderr,"ihost=%ld i=%ld\n",ihost,i);
//      if(isub < 0) fprintf(stderr,"ihost=%ld i=%ld isub=%ld\n",ihost,i,isub);
      
      /* if subhalo spawned from AHF_HOSTHALOLEVEL (or below) we check for distance (and mass!) */
      if(halos[isub].hostHaloLevel >= AHF_HOSTHALOLEVEL)
       {
        if(check_subhalo(halos+ihost, halos+isub) == TRUE)
         {
          numSubStruct++;
          SubStruct = (int *) realloc(SubStruct, numSubStruct*sizeof(int));
          SubStruct[numSubStruct-1] = isub;

#ifdef MPI_SUBHALO_FIX
          // if the host halo is not written, also do not write the subhalo
          if(halos[i].ignoreme == TRUE)
            halos[isub].ignoreme = TRUE;
          
          // but if the host halo is written, also write the subhalo
          else
            halos[isub].ignoreme = FALSE;
#endif
         }
        else
         {
          /* if it lies outside mark it as field halo */
          halos[isub].hostHalo = -1;
         }
       }
            
      /* if it spawned from above AHF_HOSTHALOLEVEL mark it as field halo */
      else
       {
        halos[isub].hostHalo = -1;
       }
     }
    
    /* copy the new substructure list over to host halo structure */
    halos[i].numSubStruct = numSubStruct;
    if(numSubStruct>0)
     {
      // remove old subStruct[] array from halos[].
      if(halos[i].subStruct) free(halos[i].subStruct);
      
      // put the new memory for halos[].subStruct[] into place
      halos[i].subStruct = SubStruct;

      // old version that copies the data over
//      for (k = 0; k < halos[i].numSubStruct; k++)
//        halos[i].subStruct[k] = SubStruct[k];
//      free(SubStruct);
      
     }
  }
  timing.ahf_halos_sfc_constructHalo += time(NULL);
    
  
#ifdef AHFnewHaloIDs
  /**************************************************************************/
  /*                   ASIGN A UNIQUE ID TO EACH HALO                       */
  /*   this is constructing a unique halo ID using the number of particles  */
  /*     and position encoding it into a 64-bit unsigned integer number     */
  /*     this has to be done before the MPI boundary haloes are removed!    */
	/**************************************************************************/
  timing.ahf_halos_sfc_constructHalo -= time(NULL);
  
  /* first get unique ID for each halo */
	for (i = 0; i < numHalos; i++) {
    halos[i].haloID = getHaloID(halos,i);
  }

  /* second update hostHalo to reflect the new IDs */
	for (i = 0; i < numHalos; i++) {
    if(halos[i].hostHalo > -1) {
      halos[i].hostHaloID = halos[halos[i].hostHalo].haloID;
    }
  }
  
  timing.ahf_halos_sfc_constructHalo += time(NULL);
#endif
  
#ifdef AHFexciseSubhaloStars // note: this flag also switched on METALHACK and GAS_PARTICLES!!!
  /**************************************************************************/
  /*              REMOVE SUBHALO STAR PARTICLES FROM HOST HALO              */
  /*     this is a hack to allow for the calculation of the luminosity      */
  /*        of central galaxies as the halo of the central normally         */
  /*          contains all star particles from the subhaloes, too           */
	/**************************************************************************/
  timing.ahf_halos_sfc_constructHalo -= time(NULL);
#ifdef VERBOSE
	fprintf(stderr, "Excising star particles in subhaloes from their hosts\n");
	fprintf(stderr, "=====================================================\n\n");
	fprintf(io.logfile, "Excising star particles in subhaloes from their hosts\n");
	fprintf(io.logfile, "=====================================================\n\n");
	fflush(io.logfile);
#endif
  
#ifdef WITH_OPENMP
#  pragma omp parallel for schedule (dynamic) shared(halos, numHalos, simu) private(i)
#endif
	for (i = 0; i < numHalos; i++) {
    /* it only makes sense to excise for haloes that
     *
     * a) contain stars
     * b) have subhaloes
     *    a) will be checked here, but b) will be checked inside exciseSubhaloStars()
     *
     *  Note, excising has to happen *before* flagging/removing outside MPI domain
     *
     */
    if(halos[i].npart >= simu.AHF_MINPART && halos[i].stars_only.npart > 0) {
      exciseSubhaloStars(halos,i);
    }
  }
  timing.ahf_halos_sfc_constructHalo += time(NULL);
#endif /* AHFexciseSubhaloStars */
  
  
  
  
  timing.ahf_io -= time(NULL);
#ifdef AHFcentrefile
  /**************************************************************************/
  /*                        JUST DUMP ALL CENTRES                           */
  /*      we want to write that file here after rem_boundary_haloes()       */
  /*        and after the calculation of all properties, actually!          */
  /* (Note, the only difference to the standard writing of _halos is that   */
  /*  here we do not check for AHF_MINPART and MPI boundaries)              */
	/**************************************************************************/
	ahf_io_WriteCenterfile(fprefix, halos, numHalos);
#endif   /* AHFcentrefile */
  
  
	/**************************************************************************/
  /*        ORDER HALOS WITH RESPECT TO NUMBER OF PARTICLES OR MASS         */
	/**************************************************************************/
	if (numHalos > 0) {
		idx    = (long unsigned *)calloc(numHalos, sizeof(long unsigned));
		idxtmp = (long unsigned *)calloc(numHalos + 1, sizeof(long unsigned));
		fsort  = (double *)calloc(numHalos + 1, sizeof(double));
		for (i = 0; i < numHalos; i++) {
#ifdef AHFsorthalosbymass
			fsort[i + 1] = (double)halos[i].M_vir;
#else
      fsort[i + 1] = (double)halos[i].npart;
#endif
    }
		indexx(numHalos, fsort, idxtmp);

		/* indexx sorts ascending, but we want descending */
		for (i = 0; i < numHalos; i++)
			idx[numHalos - i - 1] = idxtmp[i + 1] - 1;

		free(idxtmp);
		free(fsort);
	} else {
		/* If there are no halos, have a proper empty idx array */
		idx = NULL;
	}

	/**************************************************************************/
  /*                           OUTPUT FILES                                 */
	/**************************************************************************/
#ifdef VERBOSE
  fprintf(stderr, "Writing all AHF output files\n");
  fprintf(stderr, "=============================\n");
#endif
#ifdef AHF_SQL
	{
		ahf_io_sql_t sql;
		uint64_t     idOffset = 0;
		double timingTotal, timingSpecific;

#ifdef VERBOSE
   fprintf(stderr, "### GENERATING AHF SQLITE3 DATABASE\n");
#endif
   
   timingTotal = timingSpecific = -timer_getTime();
#  ifdef AHF_SQL_ONEDB_PER_TABLE
		sql = ahf_io_sql_new(fprefix, AHF_IO_SQL_MODE_MULTIPLE, idOffset);
#  else
		sql = ahf_io_sql_new(fprefix, AHF_IO_SQL_MODE_SINGLE, idOffset);
#  endif
		timingSpecific += timer_getTime();
		io_logging_msg(global_io.log, INT32_C(0),
		               "  Used %12.5gs to create SQL object",
		               timingSpecific);
		timingSpecific = - timer_getTime();
		ahf_io_sql_createTables(sql);
		timingSpecific += timer_getTime();
		io_logging_msg(global_io.log, INT32_C(0),
		               "  Used %12.5gs to create tables", timingSpecific);
		timingSpecific = - timer_getTime();
		ahf_io_sql_writeHalos(sql, halos, idx, numHalos);
		timingSpecific += timer_getTime();
		io_logging_msg(global_io.log, INT32_C(0),
		               "  Used %12.5gs to write halos", timingSpecific);
		timingSpecific = - timer_getTime();
		ahf_io_sql_writeProfiles(sql, halos, idx, numHalos);
		timingSpecific += timer_getTime();
		io_logging_msg(global_io.log, INT32_C(0),
		               "  Used %12.5gs to write profiles", timingSpecific);
		timingSpecific = - timer_getTime();
#  ifndef AHF_NO_PARTICLES
		ahf_io_sql_writeParticles(sql, halos, idx, numHalos);
		timingSpecific += timer_getTime();
		io_logging_msg(global_io.log, INT32_C(0),
		               "  Used %12.5gs to write particles",
		               timingSpecific);
#  endif /* AHF_NO_PARTICLES */
		timingTotal += timer_getTime();
		io_logging_msg(global_io.log, INT32_C(0),
		               "Used %12.5gs to write data to SQL",
		               timingTotal);

		ahf_io_sql_del(&sql);
   
#ifdef AHF_SQL_ADD_ASCII_FILES
   /*----------------------------------------------------------------------
    * here you can write whatever ASCII file you fancy, but be warned:
    *==================================================================
    * we do not support this and hence do not come crying in case you
    * end up with a segmentation fault or worse because the file you want
    * to write cannot be written given the current DEFINEFLAG combination!
    *----------------------------------------------------------------------*/
   
   
   // nothing so far as the SQL feature has not been 100% implemented yet...
   
   
#endif // AHF_SQL_ADD_ASCII_FILES
	}
#else // AHF_SQL
 {
  double timingTotal, timingSpecific;
  
  
  /* AHF_halos */
  timingTotal = timingSpecific = - timer_getTime();
  ahf_io_WriteHalos(fprefix, halos, idx, numHalos);
  timingSpecific += timer_getTime();
  io_logging_msg(global_io.log, INT32_C(0),
                 "  Used %12.5gs to write halos", timingSpecific);
  
  /* AHF_profiles */
  timingSpecific = -timer_getTime();
  ahf_io_WriteProfiles(fprefix, halos, idx, numHalos);
  timingSpecific += timer_getTime();
  io_logging_msg(global_io.log, INT32_C(0),
                 "  Used %12.5gs to write profiles", timingSpecific);
  
#ifdef AHFdisks
  timingSpecific = -timer_getTime();
  ahf_io_WriteDisks(fprefix, halos, idx, numHalos);
  timingSpecific += timer_getTime();
  io_logging_msg(global_io.log, INT32_C(0),
                 "  Used %12.5gs to write disks", timingSpecific);
#endif
  
#if ((defined AHFsubstructure) && !defined AHFrestart)
  /* AHF_substructure */
  timingSpecific = - timer_getTime();
  ahf_io_WriteSubstructure(fprefix, halos, idx, numHalos);
  timingSpecific += timer_getTime();
  io_logging_msg(global_io.log, INT32_C(0),
                 "  Used %12.5gs to write substructures", timingSpecific);
#endif 
  
#  ifndef AHF_NO_PARTICLES
  /* AHF_particles */
  timingSpecific = - timer_getTime();
  ahf_io_WriteParticles(fprefix, halos, idx, numHalos);
  timingSpecific += timer_getTime();
  io_logging_msg(global_io.log, INT32_C(0),
                 "  Used %12.5gs to write particles", timingSpecific);
#if (defined METALHACK && !defined AHFbinary)
  /* AHF_particlesSTARDUST */
  timingSpecific = - timer_getTime();
  ahf_io_WriteParticlesSTARDUST(fprefix, halos, idx, numHalos);
  timingSpecific += timer_getTime();
  io_logging_msg(global_io.log, INT32_C(0),"  Used %12.5gs to write particlesSTARDUST", timingSpecific);
#endif /* METALHACK */
#  endif /* AHF_NO_PARTICLES */   
  
  timingTotal += timer_getTime();
  io_logging_msg(global_io.log, INT32_C(0),"Used %12.5gs to write data to ASCII files",timingTotal);
 }
#endif // AHF_SQL

#ifdef AHFgeom
	/***************************************************************************/

	/*   FOR TESTING IN STEREO2:
	 *
	 * -DAHFgeom write one addition file in GEOM format
	 *
	 *  this file includes all halos as spheres ... as well as ...
	 *  -DAHFgeom_SIMUPART: a selection of all simulation particles as points
	 *  -DAHFgeom_HALOPART: all particles in every halo as points
	 */

	/**************************************************************************
	 * Printing the _halos.geom file */
	ahf_io_WriteHalosGeom(fprefix, halos, idx, numHalos);
#endif /* AHFgeom*/

  timing.ahf_io += time(NULL);

	/* remove halos[], if we are just interested in a snapshot analysis */
	for (i = 0; i < numHalos; i++)
   {
		if(halos[i].subStruct) free(halos[i].subStruct);
		if(halos[i].ipart)     free(halos[i].ipart);
#ifdef AHFexciseSubhaloStars
    if(halos[i].ipart_uniquestars) free(halos[i].ipart_uniquestars);
#endif
    dest_profile(&(halos[i]));
   }
	free(halos);

	/* free all sorts of things */
	free(gridl1dim);
	free(numIsoRef);
	if (idx) free(idx);
	if (numDensZero > 0)
		free(densZero);


#ifdef VERBOSE
	fprintf(io.logfile,
	        "################## ahf_halos finished ##################\n");
	fflush(io.logfile);
#endif

#ifdef EXTRAE_API_USAGE
  Extrae_user_function(0);
#endif

} /* ahf_halos */

/*
 *
 *===============================================================================
 * END OF MAIN
 *
 *===============================================================================
 */


int
RefCentre(gridls *grid_list, int num_refgrids, SPATIALREF **spatialRef)
{
	gridls  *cur_grid;
	pqptr   cur_pquad;
	cqptr   cur_cquad, icur_cquad;
	nqptr   cur_nquad, icur_nquad;
	nptr    cur_node;
	nptr    tsc_nodes[3][3][3];

	int     i, j, k;
	int     numNodes, tmpNum;
	int     refLevel, isoRefIndex;

	MINMAX  tmpMinMax;
	long    x, y, z;
	double  xx, yy, zz;
  double  xp, yp, zp;
	double  tmpDens, tmpPot;

	partptr current, previous, tmpll;

	int     boundRefIndex;

	double  boxLen, boxVol;
	double  px, py, pz, tmpRad, a, b, c, alpha;
	int     colour;
	double  fl1dim;

	double  xmin, xmax, ymin, ymax, zmin, zmax;

	int     iterate;

	double  cur_shift;


	/***************************************************************************/

	/* Calculating the spatial resolution (boxsize) on each refinement level
	 */

	gridl1dim = NULL;
	if ((gridl1dim = calloc(num_refgrids, sizeof(double))) == NULL) {
		fprintf(stderr, "Error in allocating the memory for halo array\n");
		exit(0);
	}
	cur_grid = global.dom_grid + ahf.min_ref;
	for (i = 0; i < num_refgrids; i++) {
		fl1dim       = ((double)(cur_grid->l1dim));
		gridl1dim[i] = fl1dim;
		cur_grid++;
	}


	/***************************************************************************/

	/* Collecting information about isolated refinements
	 *
	 * number of particles
	 * ll of particles
	 * number of nodes
	 * density centre of the refinement
	 */

	for (iterate = 0, cur_grid = global.dom_grid + ahf.min_ref;
	     iterate < ahf.no_grids;
	     iterate++, cur_grid++) {
		/* shift of cell centre as compared to edge of box [grid units] */
		cur_shift = 0.5 / (double)cur_grid->l1dim;

		for (cur_pquad = cur_grid->pquad; cur_pquad != NULL; cur_pquad = cur_pquad->next) {
			z = cur_pquad->z;
			for (cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++) {
				for (icur_cquad = cur_cquad; icur_cquad != NULL; icur_cquad = icur_cquad->next) {
					y = icur_cquad->y;
					for (cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) {
						for (icur_nquad = cur_nquad; icur_nquad != NULL; icur_nquad = icur_nquad->next) {
							x = icur_nquad->x;
							for (cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++) {
                
								/* Locate the spatial refinement */
								refLevel    = spatialRefIndex[cur_node->force.colour].refLevel;
								isoRefIndex = spatialRefIndex[cur_node->force.colour].isoRefIndex;

								/* Colour of the isolated refinement */
								spatialRef[refLevel][isoRefIndex].colour = cur_node->force.colour;

								/* Count the nodes for the refinement */
								spatialRef[refLevel][isoRefIndex].numNodes++;

								/* Position of the node */
								xx = (double)x / (double)cur_grid->l1dim + cur_shift;
								xx = (double)fmod(xx + 1.0, 1.0);
								yy = (double)y / (double)cur_grid->l1dim + cur_shift;
								yy = (double)fmod(yy + 1.0, 1.0);
								zz = (double)z / (double)cur_grid->l1dim + cur_shift;
								zz = (double)fmod(zz + 1.0, 1.0);

								/* account for periodic boundary conditions */
								px = spatialRefIndex[cur_node->force.colour].periodic.x;
								py = spatialRefIndex[cur_node->force.colour].periodic.y;
								pz = spatialRefIndex[cur_node->force.colour].periodic.z;

								if ((px == 1) && (xx < 0.5)) xx += 1.0;
								if ((py == 1) && (yy < 0.5)) yy += 1.0;
								if ((pz == 1) && (zz < 0.5)) zz += 1.0;

								/*-------------------------------------------
								 * Geometrical centre of the refinement
								 *-------------------------------------------*/
								/* simple spatial centre (no weighing) */
								spatialRef[refLevel][isoRefIndex].centreGEOM.x    += xx;
								spatialRef[refLevel][isoRefIndex].centreGEOM.y    += yy;
								spatialRef[refLevel][isoRefIndex].centreGEOM.z    += zz;
								spatialRef[refLevel][isoRefIndex].centreGEOM.norm += 1.0;


								/*-------------------------------------------
								 * Density centre of the refinement
								 *-------------------------------------------*/
								/* density at node */
								tmpDens = cur_node->dens + simu.mean_dens;

								/* double-check density value */
								if (tmpDens < 0.0) {
#ifdef VERBOSE2
									fprintf(stderr,"RefCentre(): how can we have negative densites?!  l1dim=%ld  x=%g y=%g z=%g  dens=%g  mean_dens=%g\n",
									    cur_grid->l1dim,
									    (x + 0.5) / (double)cur_grid->l1dim,
									    (y + 0.5) / (double)cur_grid->l1dim,
									    (z + 0.5) / (double)cur_grid->l1dim,
									    cur_node->dens,
									    simu.mean_dens);
									tsc_nodes[1][1][1] = cur_node;
									get_TSCnodes(cur_grid,  cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
									if (test_tsc(tsc_nodes) == FALSE)
										fprintf(stderr, "              -> boundary node!\n");
#endif
									tmpDens = 0.0;
								}

#ifdef AHFmaxdenscentre
                
								/* halo centre = position of max. density node **/
								if (tmpDens > spatialRef[refLevel][isoRefIndex].maxDens) {
									spatialRef[refLevel][isoRefIndex].centreDens.x    = xx;
									spatialRef[refLevel][isoRefIndex].centreDens.y    = yy;
									spatialRef[refLevel][isoRefIndex].centreDens.z    = zz;
									spatialRef[refLevel][isoRefIndex].centreDens.norm = 1.0;

									/* store maximum density */
									spatialRef[refLevel][isoRefIndex].maxDens = tmpDens;
								}

#else   /* AHFmaxdenscentre */
      
                /* halo centre = density weighted centre of isolated refinement */
								spatialRef[refLevel][isoRefIndex].centreDens.x    += xx * tmpDens;
								spatialRef[refLevel][isoRefIndex].centreDens.y    += yy * tmpDens;
								spatialRef[refLevel][isoRefIndex].centreDens.z    += zz * tmpDens;
								spatialRef[refLevel][isoRefIndex].centreDens.norm += tmpDens;

								/* store maximum density */
								if (tmpDens > spatialRef[refLevel][isoRefIndex].maxDens)
									spatialRef[refLevel][isoRefIndex].maxDens = tmpDens;
                
#endif   /* AHFmaxdenscentre */


								/*-------------------------------------------
								 * Potential centre of the refinement
								 *-------------------------------------------*/

								/* potential at node (switch on -DAHFpotcentre to get meaningful pot-values!) */
#ifndef AHFlean
								tmpPot = cur_node->pot;
#else
								tmpPot = 0.0;
#endif
								if ((tmpPot > 0.0) || (tmpDens < 0.0))
									tmpPot = 0.0;
								else
									tmpPot = fabs(tmpPot);

								/* halo centre = potential weighted centre of isolated refinement */
								spatialRef[refLevel][isoRefIndex].centrePot.x    += xx * tmpPot;
								spatialRef[refLevel][isoRefIndex].centrePot.y    += yy * tmpPot;
								spatialRef[refLevel][isoRefIndex].centrePot.z    += zz * tmpPot;
								spatialRef[refLevel][isoRefIndex].centrePot.norm += tmpPot;


								/*-------------------------------------------
								 * Particles in the refinement (Number and linked list)
                 *-------------------------------------------*/
								if (cur_node->ll != NULL) {
									/* It is not the first node in this spatial ref to have particels of it */
									if (spatialRef[refLevel][isoRefIndex].ll != NULL) {
										tmpll = spatialRef[refLevel][isoRefIndex].ll;
										spatialRef[refLevel][isoRefIndex].ll = cur_node->ll;

										/* Count the number of particles and find the last particle */
										current  = cur_node->ll;
										previous = current;
										while (current != NULL) {
											/*-------------------------------------------
											 * Centre-Of-Mass of particles on refinement
											 *-------------------------------------------*/
                      xp = current->pos[X];
                      yp = current->pos[Y];
                      zp = current->pos[Z];
                      if ((px == 1) && (xp < 0.5)) xp += 1.0;
                      if ((py == 1) && (yp < 0.5)) yp += 1.0;
                      if ((pz == 1) && (zp < 0.5)) zp += 1.0;
                      
											spatialRef[refLevel][isoRefIndex].centreCMpart.x += xp;
											spatialRef[refLevel][isoRefIndex].centreCMpart.y += yp;
											spatialRef[refLevel][isoRefIndex].centreCMpart.z += zp;
											spatialRef[refLevel][isoRefIndex].centreCMpart.norm += 1.0;

											spatialRef[refLevel][isoRefIndex].numParts++;
											previous = current;
											current  = current->ll;
										}

										/* Now current is the last partical in this nodes linked list */
										previous->ll = tmpll;
										tmpll        = NULL;
									} else {
										current = cur_node->ll;
										while (current != NULL) {
											/*-------------------------------------------
											 * Centre-Of-Mass of particles on refinement
											 *-------------------------------------------*/
                      xp = current->pos[X];
                      yp = current->pos[Y];
                      zp = current->pos[Z];
                      if ((px == 1) && (xp < 0.5)) xp += 1.0;
                      if ((py == 1) && (yp < 0.5)) yp += 1.0;
                      if ((pz == 1) && (zp < 0.5)) zp += 1.0;

											spatialRef[refLevel][isoRefIndex].centreCMpart.x += xp;
											spatialRef[refLevel][isoRefIndex].centreCMpart.y += yp;
											spatialRef[refLevel][isoRefIndex].centreCMpart.z += zp;
											spatialRef[refLevel][isoRefIndex].centreCMpart.norm += 1.0;

											spatialRef[refLevel][isoRefIndex].numParts++;
											current = current->ll;
										}

										spatialRef[refLevel][isoRefIndex].ll = cur_node->ll;
									}
								}
                /* cur_node->ll == NULL */
                else {
#ifdef VERBOSE3 // this leads to a lot of output to the screen that is not really helpful
                  fprintf(stderr,"RefCentre(): cur_node does not contain any particles!\n");
                  fprintf(stderr,"             spatialRef[%d][%d].numNodes=%d\n",refLevel,isoRefIndex,spatialRef[refLevel][isoRefIndex].numNodes);
                  fprintf(stderr,"             spatialRef[%d][%d].numParts=%d\n",refLevel,isoRefIndex,spatialRef[refLevel][isoRefIndex].numParts);
                  fprintf(stderr,"             cur_node->dens = %g\n",cur_node->dens);
#endif
                }
                
							}
						}
					}
				}
			}
		}
	}

	/***************************************************************************/

	/* Finishing the Spatial centre calculation and catching funny refinements
	 */
	densZero = NULL;
#ifdef WITH_OPENMP
	/* TODO: this loop may be parallelized? */
#endif
	for (i = 0; i < num_refgrids; i++) {
		for (j = 0; j < numIsoRef[i]; j++) {
			/*-------------------------------------------------
			 * normalisation of geometrical centre
			 *-------------------------------------------------*/
			if (spatialRef[i][j].centreGEOM.norm > 0) {
				spatialRef[i][j].centreGEOM.x = f1mod(spatialRef[i][j].centreGEOM.x / spatialRef[i][j].centreGEOM.norm + 1.0, 1.0);
				spatialRef[i][j].centreGEOM.y = f1mod(spatialRef[i][j].centreGEOM.y / spatialRef[i][j].centreGEOM.norm + 1.0, 1.0);
				spatialRef[i][j].centreGEOM.z = f1mod(spatialRef[i][j].centreGEOM.z / spatialRef[i][j].centreGEOM.norm + 1.0, 1.0);
			} else {
#ifdef VERBOSE2
				fprintf(
				    stderr,
				    "RefCentre():  centreGEOM  -> spatialRef[%d][%d] numNodes=%d numParts=%d geom.norm=%g\n",
				    i,
				    j,
				    spatialRef[i][j].numNodes,
				    spatialRef[i][j].numParts,
				    spatialRef[i][j].centreGEOM.norm);
#endif
			}

			/*-------------------------------------------------
			 * normalisation of density weighted centre
			 *-------------------------------------------------*/
			if (spatialRef[i][j].centreDens.norm > 0) {
				spatialRef[i][j].centreDens.x = f1mod(spatialRef[i][j].centreDens.x / spatialRef[i][j].centreDens.norm + 1.0, 1.0);
				spatialRef[i][j].centreDens.y = f1mod(spatialRef[i][j].centreDens.y / spatialRef[i][j].centreDens.norm + 1.0, 1.0);
				spatialRef[i][j].centreDens.z = f1mod(spatialRef[i][j].centreDens.z / spatialRef[i][j].centreDens.norm + 1.0, 1.0);
			} else {
#ifdef VERBOSE2
				fprintf(
				    stderr,
				    "RefCentre():  centreDens   -> spatialRef[%d][%d] numNodes=%d numParts=%d dens.norm=%g\n",
				    i,
				    j,
				    spatialRef[i][j].numNodes,
				    spatialRef[i][j].numParts,
				    spatialRef[i][j].centreDens.norm);
        fprintf(
                stderr,
                "              using GEOMcentre = %g %g %g\n",
                spatialRef[i][j].centreGEOM.x*x_fac*1000.,
                spatialRef[i][j].centreGEOM.y*x_fac*1000.,
                spatialRef[i][j].centreGEOM.z*x_fac*1000.);
#endif

				/* assign the only meaningful centre */
				spatialRef[i][j].centreDens.x = spatialRef[i][j].centreGEOM.x;
				spatialRef[i][j].centreDens.y = spatialRef[i][j].centreGEOM.y;
				spatialRef[i][j].centreDens.z = spatialRef[i][j].centreGEOM.z;
			}

			/*-------------------------------------------------
			 * normalisation of particles' centre-of-mass
			 *-------------------------------------------------*/
			if (spatialRef[i][j].centreCMpart.norm > 0) {
				spatialRef[i][j].centreCMpart.x = f1mod(spatialRef[i][j].centreCMpart.x / spatialRef[i][j].centreCMpart.norm + 1.0, 1.0);
				spatialRef[i][j].centreCMpart.y = f1mod(spatialRef[i][j].centreCMpart.y / spatialRef[i][j].centreCMpart.norm + 1.0, 1.0);
				spatialRef[i][j].centreCMpart.z = f1mod(spatialRef[i][j].centreCMpart.z / spatialRef[i][j].centreCMpart.norm + 1.0, 1.0);
			} else {
#ifdef VERBOSE2
				fprintf(
				    stderr,
				    "RefCentre():  centreCMpart -> spatialRef[%d][%d] numNodes=%d numParts=%ld CMpart.norm=%g\n",
				    i,
				    j,
				    spatialRef[i][j].numNodes,
				    spatialRef[i][j].numParts,
				    spatialRef[i][j].centreCMpart.norm);
        fprintf(
            stderr,
            "              using GEOMcentre = %g %g %g\n",
                spatialRef[i][j].centreGEOM.x*x_fac*1000.,
                spatialRef[i][j].centreGEOM.y*x_fac*1000.,
                spatialRef[i][j].centreGEOM.z*x_fac*1000.);
#endif

				/* assign the only meaningful centre */
				spatialRef[i][j].centreCMpart.x = spatialRef[i][j].centreGEOM.x;
				spatialRef[i][j].centreCMpart.y = spatialRef[i][j].centreGEOM.y;
				spatialRef[i][j].centreCMpart.z = spatialRef[i][j].centreGEOM.z;
			}


			/*-------------------------------------------------
			 * normalisation of potential centre
			 *-------------------------------------------------*/
			if (spatialRef[i][j].centrePot.norm > 0) {
				spatialRef[i][j].centrePot.x = f1mod(spatialRef[i][j].centrePot.x / spatialRef[i][j].centrePot.norm + 1.0, 1.0);
				spatialRef[i][j].centrePot.y = f1mod(spatialRef[i][j].centrePot.y / spatialRef[i][j].centrePot.norm + 1.0, 1.0);
				spatialRef[i][j].centrePot.z = f1mod(spatialRef[i][j].centrePot.z / spatialRef[i][j].centrePot.norm + 1.0, 1.0);
			} else {
#if (defined VERBOSE2 && defined AHFpotcentre)
				fprintf(
				    stderr,
				    "RefCentre():  centrePot    -> spatialRef[%d][%d] numNodes=%d numParts=%d pot.norm=%g\n",
				    i,
				    j,
				    spatialRef[i][j].numNodes,
				    spatialRef[i][j].numParts,
				    spatialRef[i][j].centrePot.norm);
        fprintf(
                stderr,
                "              using GEOMcentre = %g %g %g\n",
                spatialRef[i][j].centreGEOM.x*x_fac*1000.,
                spatialRef[i][j].centreGEOM.y*x_fac*1000.,
                spatialRef[i][j].centreGEOM.z*x_fac*1000.);
#endif

				/* assign the only meaningful centre */
				spatialRef[i][j].centrePot.x = spatialRef[i][j].centreGEOM.x;
				spatialRef[i][j].centrePot.y = spatialRef[i][j].centreGEOM.y;
				spatialRef[i][j].centrePot.z = spatialRef[i][j].centreGEOM.z;
			}


			/*-------------------------------------------------
			 * catch spatialRef's with unphysical densities
			 *-------------------------------------------------*/
			if (spatialRef[i][j].maxDens <= MACHINE_ZERO) {
				/* store these spatialRef's in a separate list */
				numDensZero++;

				if (densZero == NULL) {
					if ((densZero = calloc(numDensZero + 1, sizeof(SRINDEX))) == NULL) {
						fprintf(stderr, "calloc failed for densZero\n");
						exit(-1);
					}
				} else {
					if ((densZero = realloc(densZero, (numDensZero + 1) * sizeof(SRINDEX))) == NULL) {
						fprintf(stderr, "realloc failed for densZero\n");
						exit(-1);
					}
				}

				densZero[numDensZero - 1].refLevel    = i;
				densZero[numDensZero - 1].isoRefIndex = j;

				/* assign the only meaningful centre */
				spatialRef[i][j].centreDens.x = spatialRef[i][j].centreGEOM.x;
				spatialRef[i][j].centreDens.y = spatialRef[i][j].centreGEOM.y;
				spatialRef[i][j].centreDens.z = spatialRef[i][j].centreGEOM.z;
			}
      
      /*-------------------------------------------------
			 * catch spatialRef's with unphysical densities
			 *-------------------------------------------------*/
			if (spatialRef[i][j].numParts == 0) {
        numPartZero++;
      }

			/*-------------------------------------------------
			 * which centre to use as the halo centre?
			 *-------------------------------------------------*/
			/* the default = density weighted centre */
			spatialRef[i][j].centre.x = spatialRef[i][j].centreDens.x;
			spatialRef[i][j].centre.y = spatialRef[i][j].centreDens.y;
			spatialRef[i][j].centre.z = spatialRef[i][j].centreDens.z;

#ifdef AHFpotcentre
			spatialRef[i][j].centre.x = spatialRef[i][j].centrePot.x;
			spatialRef[i][j].centre.y = spatialRef[i][j].centrePot.y;
			spatialRef[i][j].centre.z = spatialRef[i][j].centrePot.z;
#else /* AHFpotcentre */
#  ifdef AHFcomcentre

			spatialRef[i][j].centre.x = spatialRef[i][j].centreCMpart.x;
			spatialRef[i][j].centre.y = spatialRef[i][j].centreCMpart.y;
			spatialRef[i][j].centre.z = spatialRef[i][j].centreCMpart.z;
            
#  else /* AHFcomcentre */

#   ifdef AHFgeomcentre
			spatialRef[i][j].centre.x = spatialRef[i][j].centreGEOM.x;
			spatialRef[i][j].centre.y = spatialRef[i][j].centreGEOM.y;
			spatialRef[i][j].centre.z = spatialRef[i][j].centreGEOM.z;
#   endif
         
#  endif /* AHFcomcentre */
         
#endif /* AHFpotcentre */
		}
	}

	/***************************************************************************/

	/* Gathering information needed to calculate the  boundary of the refinement
	 * Calculating the boundRefDiv */
	alpha = 1.1;
#ifdef WITH_OPENMP
	/* TODO: this loop may be parallelized? */
#endif
	for (i = 0; i < num_refgrids; i++) {
		boxLen = ((1.0) / ((double)(gridl1dim[i])));
		boxVol = boxLen * boxLen * boxLen;

		for (j = 0; j < numIsoRef[i]; j++) {
			colour                  = spatialRef[i][j].colour;

			spatialRef[i][j].volume = spatialRef[i][j].numNodes * boxVol;

			px                      = spatialRefIndex[colour].periodic.x;
			py                      = spatialRefIndex[colour].periodic.y;
			pz                      = spatialRefIndex[colour].periodic.z;

			/* Is this a periodic isolated refinement */
			if ((px == 1) || (py == 1) || (pz == 1)) {
				tmpRad = (3.0 * spatialRef[i][j].volume) / (4 * PI);
				tmpRad = pow(tmpRad, 0.333333333);
				tmpRad = tmpRad * alpha;

				if (px == 1) {
					a = tmpRad + spatialRef[i][j].centreDens.x;
					b = 1.0 - tmpRad + spatialRef[i][j].centreDens.x;

					if (b < a) {
#ifdef VERBOSE2
						fprintf(
						    stderr,
						    "'boundRefDiv' b < a The refinement is too big!! leave set min/max as 0:1\n");
#endif
						spatialRef[i][j].boundRefDiv.x = -1.0;
					}

					c                              = (a + b) / 2.0;
					spatialRef[i][j].boundRefDiv.x = fmod(c, 1.0);
				}

				if (py == 1) {
					a = tmpRad + spatialRef[i][j].centreDens.y;
					b = 1.0 - tmpRad + spatialRef[i][j].centreDens.y;

					if (b < a) {
#ifdef VERBOSE2
						fprintf(
						    stderr,
						    "'boundRefDiv' b < a The refinement is too big!! leave set min/max as 0:1\n");
#endif
						spatialRef[i][j].boundRefDiv.y = -1.0;
					}

					c                              = (a + b) / 2.0;
					spatialRef[i][j].boundRefDiv.y = fmod(c, 1.0);
				}

				if (pz == 1) {
					a = tmpRad + spatialRef[i][j].centreDens.z;
					b = 1.0 - tmpRad + spatialRef[i][j].centreDens.z;

					if (b < a) {
#ifdef VERBOSE2
						fprintf(
						    stderr,
						    "'boundRefDiv' b < a The refinement is too big!! leave set min/max as 0:1\n");
#endif
						spatialRef[i][j].boundRefDiv.z = -1.0;
					}

					c                              = (a + b) / 2.0;
					spatialRef[i][j].boundRefDiv.z = fmod(c, 1.0);
				}
			}
		}
	}


	/***************************************************************************/

	/* Collecting information about isolated refinements
	 *
	 * boundary of the refinement
	 */
	for (iterate = 0, cur_grid = global.dom_grid + ahf.min_ref;
	     iterate < ahf.no_grids;
	     iterate++, cur_grid++) {
		/* shift of cell centre as compared to edge of box [grid units] */
		cur_shift = 0.5 / (double)cur_grid->l1dim;

		for (cur_pquad = cur_grid->pquad; cur_pquad != NULL; cur_pquad = cur_pquad->next) {
			z = cur_pquad->z;
			for (cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++) {
				for (icur_cquad = cur_cquad; icur_cquad != NULL; icur_cquad = icur_cquad->next) {
					y = icur_cquad->y;
					for (cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) {
						for (icur_nquad = cur_nquad; icur_nquad != NULL; icur_nquad = icur_nquad->next) {
							x = icur_nquad->x;
							for (cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++) {
								/* Locate the spatial refinement */
								refLevel    = spatialRefIndex[cur_node->force.colour].refLevel;
								isoRefIndex = spatialRefIndex[cur_node->force.colour].isoRefIndex;

								/* Position of the node */
								xx = (double)x / (double)cur_grid->l1dim + cur_shift;
								xx = (double)fmod(xx + 1.0, 1.0);
								yy = (double)y / (double)cur_grid->l1dim + cur_shift;
								yy = (double)fmod(yy + 1.0, 1.0);
								zz = (double)z / (double)cur_grid->l1dim + cur_shift;
								zz = (double)fmod(zz + 1.0, 1.0);


								/*  The extent of the refinement
								 * NOTE :: For periodic boundary isolate refinements max < min !!!!  */

								/**********************************************************/

								/* It is not a boundary refinement */
								if (spatialRef[refLevel][isoRefIndex].boundRefDiv.x < 0.0) {
									tmpMinMax = MinMax(xx,
                                     spatialRef[refLevel][isoRefIndex].x.min,
                                     spatialRef[refLevel][isoRefIndex].x.max
                                     );
									spatialRef[refLevel][isoRefIndex].x.min = tmpMinMax.min;
									spatialRef[refLevel][isoRefIndex].x.max = tmpMinMax.max;
								}
								/* It is a boundary Refinement  */
								else {
									tmpMinMax = MinMaxBound(spatialRef[refLevel][isoRefIndex].boundRefDiv.x,
                                          xx,
                                          spatialRef[refLevel][isoRefIndex].x.min, spatialRef[refLevel][isoRefIndex].x.max
                                          );
									spatialRef[refLevel][isoRefIndex].x.min = tmpMinMax.min;
									spatialRef[refLevel][isoRefIndex].x.max = tmpMinMax.max;
								}


								/*********************************************************
								 * It is not a boundary refinement */
								if (spatialRef[refLevel][isoRefIndex].boundRefDiv.y < 0.0) {
									tmpMinMax = MinMax(
                                     yy,
                                     spatialRef[refLevel][isoRefIndex].y.min,
                                     spatialRef[refLevel][isoRefIndex].y.max
                                     );
									spatialRef[refLevel][isoRefIndex].y.min = tmpMinMax.min;
									spatialRef[refLevel][isoRefIndex].y.max = tmpMinMax.max;
								}
								/* It is a boundary Refinement  */
								else {
									tmpMinMax = MinMaxBound(
                                          spatialRef[refLevel][isoRefIndex].boundRefDiv.y,
                                          yy,
                                          spatialRef[refLevel][isoRefIndex].y.min,
                                          spatialRef[refLevel][isoRefIndex].y.max
                                          );
									spatialRef[refLevel][isoRefIndex].y.min = tmpMinMax.min;
									spatialRef[refLevel][isoRefIndex].y.max = tmpMinMax.max;
								}

								/*********************************************************
								 * It is not a boundary refinement */
								if (spatialRef[refLevel][isoRefIndex].boundRefDiv.z < 0.0) {
									tmpMinMax = MinMax(
                                     zz,
                                     spatialRef[refLevel][isoRefIndex].z.min,
                                     spatialRef[refLevel][isoRefIndex].z.max
                                     );
									spatialRef[refLevel][isoRefIndex].z.min = tmpMinMax.min;
									spatialRef[refLevel][isoRefIndex].z.max = tmpMinMax.max;
								}
								/* It is a boundary Refinement  */
								else {
									tmpMinMax = MinMaxBound(
                                          spatialRef[refLevel][isoRefIndex].boundRefDiv.z,
                                          zz,
                                          spatialRef[refLevel][isoRefIndex].z.min,
                                          spatialRef[refLevel][isoRefIndex].z.max
                                          );
									spatialRef[refLevel][isoRefIndex].z.min = tmpMinMax.min;
									spatialRef[refLevel][isoRefIndex].z.max = tmpMinMax.max;
								}
							}
						}
					}
				}
			}
		}
	}

	/***************************************************************************/

	/* Catching the case where the refinement is periodic but does not extend
	 * beyond the end domain node
	 * I.e. the periodic refinement is too small!!  refer to the 1000 case
	 *
	 */

	for (i = 0; i < num_refgrids; i++) {
		for (j = 0; j < numIsoRef[i]; j++) {
			if (spatialRef[i][j].x.min == 100000.0)
				spatialRef[i][j].x.min = 0.0;

			if (spatialRef[i][j].x.max == -100000.0)
				spatialRef[i][j].x.max = 1.0;

			if (spatialRef[i][j].y.min == 100000.0)
				spatialRef[i][j].y.min = 0.0;

			if (spatialRef[i][j].y.max == -100000.0)
				spatialRef[i][j].y.max = 1.0;

			if (spatialRef[i][j].z.min == 100000.0)
				spatialRef[i][j].z.min = 0.0;

			if (spatialRef[i][j].z.max == -100000.0)
				spatialRef[i][j].z.max = 1.0;
		}
	}

#ifdef AHFspatialReffile
 {
  FILE *fpsref;
  char filename1[MAXSTRING];
  
  sprintf(filename1,"%sz%.3f.spatialRef",global_io.params->outfile_prefix,global.z);
  fpsref = fopen(filename1,"w");
	for (i = 0; i < num_refgrids; i++) {
		for (j = 0; j < numIsoRef[i]; j++) {
      fprintf(fpsref,"%12d %12d   %f %f %f     %d    %d\n",
              i,j,
              spatialRef[i][j].centre.x*simu.boxsize*1000.,
              spatialRef[i][j].centre.y*simu.boxsize*1000.,
              spatialRef[i][j].centre.z*simu.boxsize*1000.,
              spatialRef[i][j].colour,
              spatialRef[i][j].numNodes
              );
    }
  }
  fclose(fpsref);
 }
#endif

	/* Freeing the spatialRefIndex array */
	free(spatialRefIndex);
	spatialRefIndex = NULL;

	return TRUE;
} /* RefCentre */

int
analyseRef(int num_refgrids, SPATIALREF **spatialRef)
{
	int    i, j, k, l, p;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double x, y, z;

	int    DETAILSWEEP = 0;
	int    isoRefIndex;

	int    isoRefIndexNew, isoRefIndexOLD;
	double dx, dy, dz;
	double dist, tmpMin;
	long   maxParts, maxNodes;

	double alpha; /* caution factor */

	gridls *cur_grid;
	double xx, yy, zz;
	double tmpRad;

	int    EMBEDREF;

	int    count1, count2, count3;
	int    *tmpArray;
	int    tmpCount, totalcount;

	int    numcomplex;

	int    *tmpIsoRefIndex;
	int    tmpc;

	int    hirefLevel, hiisoRefIndex;
	int    tmpisoRefIndex, tmprefLevel;

	int    minNodes = 10000000;

  tmpArray = NULL;
  
	/* Isolated refinement under investiagation */
	hirefLevel    = 3;
	hiisoRefIndex = 134;


	/**************************************************************************
	 * What is the minimum number of nodes in an isolated refinement? */
	for (i = 0; i < num_refgrids; i++) {
		for (j = 0; j < numIsoRef[i]; j++) {
			if (spatialRef[i][j].numNodes < minNodes)
				minNodes = spatialRef[i][j].numNodes;
		}
	}
#ifdef VERBOSE2
	fprintf(stderr, "#### Minimum number of nodes = %d\n", minNodes);
#endif

	/**************************************************************************
	 * Searching for embedded refinements */
	for (i = 0; i < num_refgrids - 1; i++) {
		for (j = 0; j < numIsoRef[i]; j++) {
			xmin = spatialRef[i][j].x.min;
			xmax = spatialRef[i][j].x.max;
			ymin = spatialRef[i][j].y.min;
			ymax = spatialRef[i][j].y.max;
			zmin = spatialRef[i][j].z.min;
			zmax = spatialRef[i][j].z.max;

			/* Is the centre of a finer refinement in this refinement? */
			for (k = 0; k < numIsoRef[i + 1]; k++) {
				/* Note we could have used the CofMass centres */
				x = spatialRef[i + 1][k].centreDens.x;
				y = spatialRef[i + 1][k].centreDens.y;
				z = spatialRef[i + 1][k].centreDens.z;


				/* If yes record the occurance */
				EMBEDREF = 0;

				if (xmin < xmax) {
					if ((x > xmin) && (x < xmax))
						EMBEDREF = 1;
				} else {
					if ((x >= 0) && (x < xmax))
						EMBEDREF = 1;

					if ((x > xmin) && (x <= 1.0))
						EMBEDREF = 1;
				}
				if (EMBEDREF == 1) {
					if (ymin < ymax) {
						if ((y > ymin) && (y < ymax))
							EMBEDREF = 2;
					} else {
						if ((y >= 0) && (y < ymax))
							EMBEDREF = 2;

						if ((y > ymin) && (y <= 1.0))
							EMBEDREF = 2;
					}
				}
				if (EMBEDREF == 2) {
					if (zmin < zmax) {
						if ((z > zmin) && (z < zmax))
							EMBEDREF = 3;
					} else {
						if ((z >= 0) && (z < zmax))
							EMBEDREF = 3;

						if ((z > zmin) && (z <= 1.0))
							EMBEDREF = 3;
					}
				}
				if (EMBEDREF == 3) {
					/***************************************************************************************************
					 * Recording the old substruct info */
					if (spatialRef[i][j].numSubStruct > 0) {
            if(tmpArray) free(tmpArray);
						tmpArray = NULL;
						if ((tmpArray = calloc(spatialRef[i][j].numSubStruct, sizeof(int))) == NULL) {
							fprintf(stderr, "Error in allocating the memory for Recording the old substruct info\n");
							exit(0);
						}

						for (l = 0; l < spatialRef[i][j].numSubStruct; l++) {
							tmpArray[l] = spatialRef[i][j].subStruct[l].isoRefIndex;
						}
					}

					/* Adding the new substructure */
					spatialRef[i][j].numSubStruct++;

					if (spatialRef[i][j].subStruct == NULL) {
						if ((spatialRef[i][j].subStruct = calloc(spatialRef[i][j].numSubStruct + 1, sizeof(INDEX))) == NULL) {
							fprintf(stderr, "calloc failed in connecting nodes\n");
							exit(-1);
						}
					} else {
						if ((spatialRef[i][j].subStruct = realloc(spatialRef[i][j].subStruct, (spatialRef[i][j].numSubStruct + 1) * sizeof(INDEX))) == NULL) {
							fprintf(stderr, "realloc failed in connecting nodes\n");
							exit(-1);
						}
					}

					if (spatialRef[i][j].numSubStruct - 1 > 0) {
						for (l = 0; l < spatialRef[i][j].numSubStruct - 1; l++) {
							spatialRef[i][j].subStruct[l].refLevel     = i + 1;
							spatialRef[i][j].subStruct[l].isoRefIndex  = tmpArray[l];
						}
            free(tmpArray);
            tmpArray = NULL;
					}
          
					spatialRef[i][j].subStruct[spatialRef[i][j].numSubStruct-1].refLevel    = i + 1;
					spatialRef[i][j].subStruct[spatialRef[i][j].numSubStruct-1].isoRefIndex = k;


					/***************************************************************************************************
					 * Recording the Parent domain information
					 * i+1 = finer refinment
					 * k   = looping through the isolated refinements on this level */
					if (spatialRef[i + 1][k].numParDom > 0) {
            if(tmpArray) free(tmpArray);
						tmpArray = NULL;
						if ((tmpArray = calloc(spatialRef[i + 1][k].numParDom, sizeof(int))) == NULL) {
							fprintf(stderr, "Error in allocating the memory for Recording the Parent domain information\n");
							exit(0);
						}

						for (l = 0; l < spatialRef[i + 1][k].numParDom; l++)
							tmpArray[l] = spatialRef[i + 1][k].parDom[l].isoRefIndex;
					}

					/*  Recording the Parent domain */
					spatialRef[i + 1][k].numParDom++;

					if (spatialRef[i + 1][k].parDom == NULL) {
						if ((spatialRef[i + 1][k].parDom = calloc(spatialRef[i + 1][k].numParDom + 1, sizeof(INDEX))) == NULL) {
							fprintf(stderr, "calloc failed in connecting nodes\n");
							exit(-1);
						}
					} else {
						if ((spatialRef[i + 1][k].parDom = realloc(spatialRef[i + 1][k].parDom, (spatialRef[i + 1][k].numParDom + 1) * sizeof(INDEX))) == NULL) {
							fprintf(stderr, "realloc failed in connecting nodes\n");
							exit(-1);
						}
					}

					if (spatialRef[i + 1][k].numParDom - 1 > 0) {
						for (l = 0; l < spatialRef[i + 1][k].numParDom - 1; l++) {
							spatialRef[i + 1][k].parDom[l].refLevel    = i;
							spatialRef[i + 1][k].parDom[l].isoRefIndex = tmpArray[l];
						}
            free(tmpArray);
            tmpArray = NULL;
					}
          
					spatialRef[i + 1][k].parDom[spatialRef[i + 1][k].numParDom - 1].refLevel    = i;
					spatialRef[i + 1][k].parDom[spatialRef[i + 1][k].numParDom - 1].isoRefIndex = j;

					/* Flaging if we need to do a complex sweep */
					if (spatialRef[i + 1][k].numParDom > 1) {
						DETAILSWEEP = 1;
					}
				}
			}
		}
	}

	/**************************************************************************
	 * If the refinements are tricky then the simply algorithm above will not work
	 * but this one will :) */
	numcomplex = 0;
	if (DETAILSWEEP == 1) {
#ifdef VERBOSE2
		fprintf(stderr, "analyseRef(): we need to do a complex sweep to find embedded halos\n");
#endif

		for (i = 1; i < num_refgrids - 1; i++) {
			for (j = 0; j < numIsoRef[i]; j++) {
				if (spatialRef[i][j].numParDom > 1) {
					numcomplex++;

					/***********************************************************/

					/* parent halo = halo on next coarser level
					 *               with centre closest to current spatialRef
					 **/
					tmpMin   = 10000000000000.0;
					maxNodes = 0;
					maxParts = 0;
					for (k = 0; k < spatialRef[i][j].numParDom; k++) {
						isoRefIndex    = spatialRef[i][j].parDom[k].isoRefIndex;
						isoRefIndexNew = -1;

						/* calculate the distance to the parent */
						dx = spatialRef[i][j].centreDens.x - spatialRef[i - 1][isoRefIndex].centreDens.x;
						dx = fabs(dx);

						dy = spatialRef[i][j].centreDens.y - spatialRef[i - 1][isoRefIndex].centreDens.y;
						dy = fabs(dy);

						dz = spatialRef[i][j].centreDens.z - spatialRef[i - 1][isoRefIndex].centreDens.z;
						dz = fabs(dz);

						if (dx > 0.5)							dx = 1.0 - dx;
						if (dy > 0.5)							dy = 1.0 - dy;
						if (dz > 0.5)							dz = 1.0 - dz;

						dist = dx * dx + dy * dy + dz * dz;

						/* isoRefIndexNew :: is the RefIndex of the parent
						 *refinement. I.e. the refinement on the coarser level
						 *that is the closest */
						if (dist < tmpMin) {
							isoRefIndexNew = isoRefIndex;
							tmpMin         = dist;
						}
					}

#ifdef VERBOSE2
					if (isoRefIndex < 0)
						fprintf(
						    stderr,
						    "analyseRef(): WARNING ->  no parent   found for spatialRef[%d][%d]  "
						    "x=%g y=%g z=%g   numPart=%lu  numNodes=%d\n",
						    i,
						    j,
						    spatialRef[i][j].centreDens.x,
						    spatialRef[i][j].centreDens.y,
						    spatialRef[i][j].centreDens.z,
						    spatialRef[i][j].numParts,
						    spatialRef[i][j].numNodes);
#endif

					/******************************************************
					 * Correct the substructure listings */

					/* Run through the parents that list this substructure */
					tmpc = 0;
					for (k = 0; k < spatialRef[i][j].numParDom; k++) {
						/* Name of parent to check */
						isoRefIndex = spatialRef[i][j].parDom[k].isoRefIndex;

						/* If this refinment is not the real parent then remove this sub-halo from it's substructure listing
						 * remove sub-halo I.e. (i,j) from the parent halo I.e. (i-1,isoRefIndex) substructure listing */
						if (isoRefIndex != isoRefIndexNew) {
							/* Create tmp array of the names of the sub halos for this parent*/
							if ((tmpArray = calloc(spatialRef[i - 1][isoRefIndex].numSubStruct, sizeof(int))) == NULL) {
								fprintf(stderr, "Error in allocating the memory for remove substrcture array\n");
								exit(0);
							}

							for (p = 0; p < spatialRef[i - 1][isoRefIndex].numSubStruct; p++)
								tmpArray[p] = spatialRef[i - 1][isoRefIndex].subStruct[p].isoRefIndex;

							/* removing the substructure (i,j) from the parent substructure listing */
							spatialRef[i - 1][isoRefIndex].numSubStruct = spatialRef[i - 1][isoRefIndex].numSubStruct - 1;
              if(spatialRef[i - 1][isoRefIndex].subStruct) free(spatialRef[i - 1][isoRefIndex].subStruct);
							spatialRef[i - 1][isoRefIndex].subStruct = NULL;

							/* QUICK AND DIRTY FIX ADDED BY AK ON 30/09/2005 */
							if (spatialRef[i - 1][isoRefIndex].numSubStruct == 0) {
#ifdef VERBOSE2
								fprintf(stderr, "   -> NO MORE SUBSTRUCTURE LEFT IN PARENT HALO (%d, %d)\n", i - 1, isoRefIndex);
#endif

								/* however, the loop over p right below requires at least on INDEX structure */
								if ((spatialRef[i - 1][isoRefIndex].subStruct = calloc(1, sizeof(INDEX))) == NULL) {
									fprintf(stderr, "calloc failed in substructure\n");
									exit(-1);
								}
							} else {
								if ((spatialRef[i - 1][isoRefIndex].subStruct = calloc(spatialRef[i - 1][isoRefIndex].numSubStruct, sizeof(INDEX))) == NULL) {
									fprintf(stderr, "calloc failed in substructure\n");
									exit(-1);
								}
							}

							tmpCount = 0;
							for (p = 0; p < spatialRef[i - 1][isoRefIndex].numSubStruct + 1; p++) {
								/* I.e. include all except our substructure */
								if (tmpArray[p] != j) {
									spatialRef[i - 1][isoRefIndex].subStruct[tmpCount].isoRefIndex = tmpArray[p];
									spatialRef[i - 1][isoRefIndex].subStruct[tmpCount].refLevel    = i;
									tmpCount++;
								}
							}
              
              free(tmpArray);
              tmpArray = NULL;
						}

					}

					/******************************************************
					 * Finalise the parent halo */
					spatialRef[i][j].numParDom = 1;

					if (spatialRef[i][j].parDom == NULL) {
						if ((spatialRef[i][j].parDom = calloc(spatialRef[i][j].numParDom + 1, sizeof(INDEX))) == NULL) {
							fprintf(stderr, "calloc failed in connecting nodes\n");
							exit(-1);
						}
					} else {
						if ((spatialRef[i][j].parDom = realloc(spatialRef[i][j].parDom, (spatialRef[i][j].numParDom + 1) * sizeof(INDEX))) == NULL) {
							fprintf(stderr, "realloc failed in connecting nodes\n");
							exit(-1);
						}
					}

					/* only store parent if we actually found one */
					if (isoRefIndexNew > 0) {
						spatialRef[i][j].parDom[0].isoRefIndex = isoRefIndexNew;
						spatialRef[i][j].parDom[0].refLevel    = i - 1;
					}
					/* otherwise we will try a second time just below... */
					else {
						if (spatialRef[i][j].parDom != NULL)
							free(spatialRef[i][j].parDom);

						spatialRef[i][j].numParDom = 0;
						spatialRef[i][j].parDom    = NULL;
					}
				}
			}
		}
	}


	/********************************************************************************************************************************************************
	 *******************************************************************************************/
        
	/* Double checking
	 * Picking up halos that don't have parents */
	for (i = 1; i < num_refgrids; i++) {
		for (j = 0; j < numIsoRef[i]; j++) {
			if (spatialRef[i][j].parDom == NULL) {
#ifdef AHFDEBUG
				fprintf(
				    stderr,
				    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
				fprintf(
				    stderr,
				    "analyseRef(): no parent for spatialRef[%d][%d].numParDom = %d\n",
				    i,
				    j,
				    spatialRef[i][j].numParDom);
#endif

				/* What is the closest potential parent refinement? */
				x              = spatialRef[i][j].centreDens.x;
				y              = spatialRef[i][j].centreDens.y;
				z              = spatialRef[i][j].centreDens.z;
				tmpMin         = 10000000000000.0;
				tmpisoRefIndex = -1.0;
				tmprefLevel    = -1.0;

				/* loop over all spatial refinements on next coarser level */
				for (p = 0; p < numIsoRef[i - 1]; p++) {
          
					/* irrespective of "#ifdef PARDAU_*" we search for the refinement closest in distance */

					dx = x - spatialRef[i - 1][p].centreDens.x;
					dx = fabs(dx);

					dy = y - spatialRef[i - 1][p].centreDens.y;
					dy = fabs(dy);

					dz = z - spatialRef[i - 1][p].centreDens.z;
					dz = fabs(dz);

					if (dx > 0.5)						dx = 1.0 - dx;
					if (dy > 0.5)						dy = 1.0 - dy;
					if (dz > 0.5)						dz = 1.0 - dz;

					dist = dx * dx + dy * dy + dz * dz;

					if (dist < tmpMin) {
						tmpMin         = dist;
						tmpisoRefIndex = p;
						tmprefLevel    = i - 1;
					}
				}

#ifdef AHFDEBUG
				fprintf(stderr, "PARENT [%d][%d]\n", tmprefLevel, tmpisoRefIndex);
				fprintf(stderr, "numSubStruct=%d daughter=[%d][%d]\n",
                spatialRef[i - 1][tmpisoRefIndex].numSubStruct,
                spatialRef[i - 1][tmpisoRefIndex].daughter.isoRefIndex,
                spatialRef[i - 1][tmpisoRefIndex].daughter.refLevel);
#endif

				/* We have found the lost parent refinement
				 * Name the parent */
				spatialRef[i][j].numParDom   = 1;
        if(spatialRef[i][j].parDom) free(spatialRef[i][j].parDom);
				spatialRef[i][j].parDom      = NULL;
				if ((spatialRef[i][j].parDom = calloc(1, sizeof(INDEX))) == NULL) {
					fprintf(stderr, "Error in allocating the memory for Recording the Parent domain information\n");
					exit(0);
				}
				spatialRef[i][j].parDom[0].isoRefIndex = tmpisoRefIndex;
				spatialRef[i][j].parDom[0].refLevel    = tmprefLevel;


#ifdef AHFDEBUG
				/* Update the parents substructure list */
				fprintf(stderr, "spatialRef[%d][%d].numSubStruct = %d\n", tmprefLevel, tmpisoRefIndex, spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct);
#endif

				/* Create tmp array of the names of the sub halos for this parent */
        if(tmpArray) free(tmpArray);
				tmpArray = NULL;
        
        /* we need to add +1 to make space for the new substructure */
				if ((tmpArray = calloc(spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct+1, sizeof(int))) == NULL) {
					fprintf(stderr, "Error in allocating the memory for remove substrcture array\n");
					exit(0);
				}

        /* loop over the already existing substructures */
				for (p = 0; p < spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct; p++)
					tmpArray[p] = spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].isoRefIndex;

				/* adding the substructure to the end of the parent substructure listing */
				spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct += 1;
        
        /* generate a new .subStruct[] array able to hold one more entry at the end */
				if(spatialRef[tmprefLevel][tmpisoRefIndex].subStruct) free(spatialRef[tmprefLevel][tmpisoRefIndex].subStruct);
				spatialRef[tmprefLevel][tmpisoRefIndex].subStruct = NULL;
				if ((spatialRef[tmprefLevel][tmpisoRefIndex].subStruct = calloc(spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct, sizeof(INDEX))) == NULL) {
					fprintf(stderr, "calloc failed in substrcture\n");
					exit(-1);
				}

        /* copy old substructure over to new .subStruct[] array */
				for (p = 0; p < spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct - 1; p++) {
					spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].isoRefIndex = tmpArray[p];
					spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].refLevel    = i;
				}

				if(tmpArray) free(tmpArray);
				tmpArray = NULL;

        /* fill in new substructure */
				spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct - 1].refLevel    = i;
				spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct - 1].isoRefIndex = j;

#ifdef AHFDEBUG
				fprintf(stderr, "spatialRef[%d][%d].numSubStruct = %d\n",
				    tmprefLevel,
				    tmpisoRefIndex,
				    spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct);
				for (p = 0; p < spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct; p++)
					fprintf(stderr, "[%d][%d], ",
					    spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].refLevel,
					    spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].isoRefIndex);

				fprintf(stderr, "\n");
#endif


				/* Giving this lost refinement the centre from the refinement below
				 * We do this because it might not actually be related.
				 * If it is related then it will sort itself out - if not the
				 * partilces we be removed by rem_unbound() */
				spatialRef[i][j].centreDens.x = spatialRef[tmprefLevel][tmpisoRefIndex].centreDens.x;
				spatialRef[i][j].centreDens.y = spatialRef[tmprefLevel][tmpisoRefIndex].centreDens.y;
				spatialRef[i][j].centreDens.z = spatialRef[tmprefLevel][tmpisoRefIndex].centreDens.z;
				spatialRef[i][j].centreGEOM.x = spatialRef[tmprefLevel][tmpisoRefIndex].centreGEOM.x;
				spatialRef[i][j].centreGEOM.y = spatialRef[tmprefLevel][tmpisoRefIndex].centreGEOM.y;
				spatialRef[i][j].centreGEOM.z = spatialRef[tmprefLevel][tmpisoRefIndex].centreGEOM.z;
			} /* parDom == NULL */
		}
	}

	/********************************************************************************************************************************************************
	 **************************************************************************
	 * Finding the daughter refinement */

	/* just to make the collecting radii a little smaller */
	alpha = 0.9;
	for (i = 0; i < num_refgrids - 1; i++) {
		for (j = 0; j < numIsoRef[i]; j++) {
			/*************************************************
			 * Finding the daughter refinement */
			if (spatialRef[i][j].numSubStruct > 1) {
				tmpMin         = 10000000000000.0;
				maxNodes       = 0;
				maxParts       = 0;
				isoRefIndexNew = -1;

				for (k = 0; k < spatialRef[i][j].numSubStruct; k++) {
					isoRefIndex = spatialRef[i][j].subStruct[k].isoRefIndex;

#ifdef PARDAU_DISTANCE
					dx = spatialRef[i][j].centreDens.x - spatialRef[i + 1][isoRefIndex].centreDens.x;
					dx = fabs(dx);

					dy = spatialRef[i][j].centreDens.y - spatialRef[i + 1][isoRefIndex].centreDens.y;
					dy = fabs(dy);

					dz = spatialRef[i][j].centreDens.z - spatialRef[i + 1][isoRefIndex].centreDens.z;
					dz = fabs(dz);

					if (dx > 0.5)						dx = 1.0 - dx;
					if (dy > 0.5)						dy = 1.0 - dy;
					if (dz > 0.5)						dz = 1.0 - dz;


					dist = dx * dx + dy * dy + dz * dz;

					if (dist < tmpMin) {
						isoRefIndexNew = isoRefIndex;
						tmpMin         = dist;
					}
#endif
#ifdef PARDAU_NODES
					if (spatialRef[i + 1][isoRefIndex].numNodes > maxNodes) {
						isoRefIndexNew = isoRefIndex;
						maxNodes       = spatialRef[i + 1][isoRefIndex].numNodes;
					}
#endif
#ifdef PARDAU_PARTS
					if (spatialRef[i + 1][isoRefIndex].numParts > maxParts) {
						isoRefIndexNew = isoRefIndex;
						maxParts       = spatialRef[i + 1][isoRefIndex].numParts;
					}
#endif
				}


				spatialRef[i][j].daughter.isoRefIndex = isoRefIndexNew;
				spatialRef[i][j].daughter.refLevel    = i + 1;

#ifdef VERBOSE2
				if (isoRefIndexNew < 0)
					fprintf(
					    stderr,
					    "analyseRef(): WARNING ->  no daughter found for spatialRef[%d][%d]  "
					    "x=%g y=%g z=%g   numPart=%lu  numNodes=%d\n",
					    i,
					    j,
					    spatialRef[i][j].centreDens.x,
					    spatialRef[i][j].centreDens.y,
					    spatialRef[i][j].centreDens.z,
					    spatialRef[i][j].numParts,
					    spatialRef[i][j].numNodes);
#endif
			} else if (spatialRef[i][j].numSubStruct == 1) {
				spatialRef[i][j].daughter.isoRefIndex = spatialRef[i][j].subStruct[0].isoRefIndex;
				spatialRef[i][j].daughter.refLevel    = i + 1;
			} else {
				/* no daughter at all */
				spatialRef[i][j].daughter.isoRefIndex = -1;
				spatialRef[i][j].daughter.refLevel    = -1;
			}


			/*************************************************
			 * Calculating the distance to the closest refinement */
			if (spatialRef[i][j].numSubStruct > 1)
#ifdef AHFnewCloseRefDist
      { /* IF loop over all subhaloes within this host [i][j] */
				for (k = 0; k < spatialRef[i][j].numSubStruct; k++) {
					isoRefIndex    = spatialRef[i][j].subStruct[k].isoRefIndex;
					isoRefIndexOLD = isoRefIndex;
          
          // do not try to find closeRefDist for the daughter
          if(isoRefIndex != spatialRef[i][j].daughter.isoRefIndex)
           {
            x = spatialRef[i+1][isoRefIndex].centreDens.x;
            y = spatialRef[i+1][isoRefIndex].centreDens.y;
            z = spatialRef[i+1][isoRefIndex].centreDens.z;
            
            tmpMin = 10000000000000.0;
            /* loop over all subhaloes other than k */
            for (l = 0; l < spatialRef[i][j].numSubStruct; l++) {
              
              // avoid distance to itself
              if (k != l) {
                isoRefIndex = spatialRef[i][j].subStruct[l].isoRefIndex;
                
                dx = x - spatialRef[i+1][isoRefIndex].centreDens.x;
                dx = fabs(dx);
                
                dy = y - spatialRef[i+1][isoRefIndex].centreDens.y;
                dy = fabs(dy);
                
                dz = z - spatialRef[i+1][isoRefIndex].centreDens.z;
                dz = fabs(dz);
                
                if (dx > 0.5)	dx = 1.0 - dx;
                if (dy > 0.5)	dy = 1.0 - dy;
                if (dz > 0.5) dz = 1.0 - dz;
                
                dist = pow2(dx) + pow2(dy) + pow2(dz);
                
                if (dist < tmpMin && spatialRef[i+1][isoRefIndex].numParts > spatialRef[i+1][isoRefIndexOLD].numParts)
                  tmpMin = dist;
              }
            }
            
            // this is highly tunable, but so far the half-distance worked best!
            spatialRef[i+1][isoRefIndexOLD].closeRefDist = 0.75 * sqrt(tmpMin);
          }
        }
			}
#else // AHFnewCloseRefDist
      { /* IF loop over all subhaloes within this host [i][j] */
				for (k = 0; k < spatialRef[i][j].numSubStruct; k++) {
					isoRefIndex    = spatialRef[i][j].subStruct[k].isoRefIndex;
					isoRefIndexOLD = isoRefIndex;
          
          // do not try to find closeRefDist for the daughter
          if(isoRefIndex != spatialRef[i][j].daughter.isoRefIndex)
           {
            x = spatialRef[i+1][isoRefIndex].centreDens.x;
            y = spatialRef[i+1][isoRefIndex].centreDens.y;
            z = spatialRef[i+1][isoRefIndex].centreDens.z;
            
            tmpMin = 10000000000000.0;
            /* loop over all subhaloes other than k */
            for (l = 0; l < spatialRef[i][j].numSubStruct; l++) {
              
              // avoid distance to itself
              if (k != l) {
                isoRefIndex = spatialRef[i][j].subStruct[l].isoRefIndex;
                
                dx  = x - spatialRef[i+1][isoRefIndex].centreDens.x;
                dx = fabs(dx);
                
                dy = y  - spatialRef[i+1][isoRefIndex].centreDens.y;
                dy = fabs(dy);
                
                dz = z - spatialRef[i+1][isoRefIndex].centreDens.z;
                dz = fabs(dz);
                
                if (dx > 0.5)	dx = 1.0 - dx;
                if (dy > 0.5)	dy = 1.0 - dy;
                if (dz > 0.5) dz = 1.0 - dz;
                
                dist = pow2(dx) + pow2(dy) + pow2(dz);
                
                //if (dist < tmpMin && spatialRef[i+1][isoRefIndex].numParts > spatialRef[i+1][isoRefIndexOLD].numParts)
                if (dist < tmpMin)
                  tmpMin = dist;
              }
            }
            
            // this is highly tunable, but so far the half-distance worked best!
            spatialRef[i+1][isoRefIndexOLD].closeRefDist = 0.5 * sqrt(tmpMin);
          }
        }
			} /* IF loop over all subhaloes within this host [i][j] */
#endif // AHFnewCloseRefDist
      
		} // loop over all spatialRef[i][j]
	}


	/********************************************************************************************************************************************************
	 ******************************************************************************************
	 * Checking the refinement reconstructuion */
	totalcount = 0;
	for (i = 0; i < num_refgrids; i++) {
		count1 = 0;
		count2 = 0;
		count3 = 0;

		/* Count the number of daughters */
		if (i != num_refgrids - 1) {
			for (j = 0; j < numIsoRef[i]; j++) {
				if (spatialRef[i][j].daughter.isoRefIndex != -1)
					count1++;
			}
		}

		/* Count the number of substructure */
		if (i != num_refgrids - 1) {
			for (j = 0; j < numIsoRef[i]; j++)
				count2 = count2 + spatialRef[i][j].numSubStruct;
		}

		/* Count the number of parents */
		if (i != 0) {
			for (j = 0; j < numIsoRef[i]; j++) {
				count3 = count3 + spatialRef[i][j].numParDom;
			}
		}
#ifdef VERBOSE
		fprintf(
		    stderr,
		    "%3d || numSubStruct(%12d) numParDom(%12d) numDaughter(%12d) newHalos(%12d)\n",
		    i,
		    count2,
		    count3,
		    count1,
		    count2 - count1);
#endif
		totalcount = totalcount + count2 - count1;
	}
#ifdef VERBOSE
	fprintf(stderr, "Number of Dark Matter Halos = %d\n", totalcount + numIsoRef[0]);
#endif


	/******************************************************************************************
	 * Checking for duplication of subhalos */

#ifdef VERBOSE2
	for (i = 0; i < num_refgrids; i++) {
		tmpIsoRefIndex = NULL;
		tmpCount       = 0;
		for (j = 0; j < numIsoRef[i]; j++) {
      
			for (k = 0; k < spatialRef[i][j].numSubStruct; k++) {
        
				for (l = 0; l < tmpCount; l++) {
					if (tmpIsoRefIndex[l] == spatialRef[i][j].subStruct[k].isoRefIndex)
						fprintf(stderr, "#######  [%d][%d]  (while checking for duplication of halos...)\n",
						    spatialRef[i][j].subStruct[k].refLevel,
						    spatialRef[i][j].subStruct[k].isoRefIndex);
				}

				if (tmpIsoRefIndex == NULL) {
					if ((tmpIsoRefIndex = calloc(tmpCount + 1, sizeof(INDEX))) == NULL) {
						fprintf(stderr, "calloc failed in substructure %d %ld\n", tmpCount, sizeof(INDEX));
						exit(-1);
					}
				} else {
					if ((tmpIsoRefIndex = realloc(tmpIsoRefIndex, (tmpCount + 1) * sizeof(INDEX))) == NULL) {
						fprintf(stderr, "realloc failed in substructure %d %ld\n", tmpCount, sizeof(INDEX));
						exit(-1);
					}
				}
				tmpIsoRefIndex[tmpCount] = spatialRef[i][j].subStruct[k].isoRefIndex;
				tmpCount++;
			}
      
		}
		if (tmpIsoRefIndex != NULL) {
			free(tmpIsoRefIndex);
			tmpIsoRefIndex = NULL;
		}
	}
#endif   /* VERBOSE2 */


	/*******************************************************************************************************************************************************
	 *******************************************************************************************************************************************************
	 ********************************************************************************************************************************************************/


	return TRUE;
} /* analyseRef */


int spatialRef2halos(int num_refgrids, SPATIALREF **spatialRef)
{
	int     i, j, k, f, ii, jj;
	int     count, numNewHalos;
	int     isoRefIndex, refLevel;
	int     SSisoRefIndex, SSrefLevel;
	double  dx, dy, dz, tmpRad;
	int     primHaloIndex, haloIndex;
	partptr current, previous, subCurrent;
	int     countBC, OKcountBC, tmpCount;
	int     countTMP;

	int     kcount;
	int     *tmpSubStruct;

	int     refgrid_start;
	double  oldRad;
	int     jnumpart, inumpart;
	double  maxGathRad, gatherRad2;

	int     tmp;
  
  long unsigned *idxtmp, *idx;
	double        *fsort;

  time_t  tdummy;
  
  idx = NULL;

  /* this version follows the same logic as the loops below filling halos[] with meaning */
	count = 0;
  
  /* treat 1st grid separately */
  i = 0;
  for (j = 0; j < numIsoRef[i]; j++)
   {
    if(spatialRef[i][j].numSubStruct == 0)
      count++;
    else if(spatialRef[i][j].numSubStruct == 1)
      count++;
    else
      count += spatialRef[i][j].numSubStruct;
   }
  
  /* the refinment levels */
	for (i = 1; i < num_refgrids; i++)
   {
		for (j = 0; j < numIsoRef[i]; j++)
     {
      if(spatialRef[i][j].numSubStruct > 1)
        //count += spatialRef[i][j].numSubStruct;
        count += (spatialRef[i][j].numSubStruct-1);   // AK: otherwise we count the main-branch again and again and again
     }
   }
	numHalos = count;

#ifdef VERBOSE
	fprintf(stderr, " spatialRef2halos():\n");
	fprintf(stderr, "  number of isolated refinements  = %d\n", totnumIsoRef);
	fprintf(stderr, "  first guess for number of halos = %d\n", numHalos);
	fprintf(io.logfile, " spatialRef2halos():\n");
	fprintf(io.logfile, "  number of isolated refinements  = %d\n", totnumIsoRef);
	fprintf(io.logfile, "  first guess for number of halos = %d\n", numHalos);
	fflush(io.logfile);
#endif


	/* allocate memory for halos:
	 *----------------------------
	 * rather use MAX(totnumIsoRef,numHalos) than numHalos
	 * because there is still this strange bug leading to "count != numHalos"
	 */
	halos = NULL;
	if ((halos = calloc(MAX(totnumIsoRef+1,numHalos), sizeof(HALO))) == NULL) {
		fprintf(io.logfile, "Error in allocating the memory for halo array\n");
		exit(0);
	}

	/* initialize halo properties */
	for (i = 0; i < MAX(totnumIsoRef+1,numHalos); i++) {
		halos[i].npart      = 0;
		halos[i].nll        = 0;
		halos[i].ipart      = NULL;
		halos[i].ll         = NULL;

		halos[i].pos.x      = 0;
		halos[i].pos.y      = 0;
		halos[i].pos.z      = 0;
		halos[i].vel.x      = 0;
		halos[i].vel.y      = 0;
		halos[i].vel.z      = 0;
		halos[i].M_vir      = -1.0;
		halos[i].R_vir      = -1.0;
		halos[i].sigV       = 0.0;
		halos[i].lambda     = 0.0;
		halos[i].lambdaE    = 0.0;
		halos[i].R_max      = 0.0;
		halos[i].V2_max     = 0.0;
		halos[i].ovdens     = 0.0;
		halos[i].R_edge     = 0.0;
		halos[i].Ekin       = 0.0;
		halos[i].Epot       = 0.0;
		halos[i].Phi0       = 0.0;
		halos[i].axis.x     = 0.0;
		halos[i].axis.y     = 0.0;
		halos[i].axis.z     = 0.0;
		halos[i].E1.x       = 0.0;
		halos[i].E1.y       = 0.0;
		halos[i].E1.z       = 0.0;
		halos[i].E2.x       = 0.0;
		halos[i].E2.y       = 0.0;
		halos[i].E2.z       = 0.0;
		halos[i].E3.x       = 0.0;
		halos[i].E3.y       = 0.0;
		halos[i].E3.z       = 0.0;
		halos[i].AngMom.x   = 0.0;
		halos[i].AngMom.y   = 0.0;
		halos[i].AngMom.z   = 0.0;
		halos[i].com_offset = 0.0;
		halos[i].mbp_offset = 0.0;
		halos[i].r2         = 0.0;
    halos[i].Phi0       = 0.0;
    halos[i].v_esc2     = 0.0;
    halos[i].cNFW       = 0.0;
    halos[i].SurfP      = 0.0;


		halos[i].hostHaloLevel    = -1;
		halos[i].hostHalo         = -1;
    
#ifdef AHFnewHaloIDs
		halos[i].haloID           = 0;
		halos[i].hostHaloID       = 0;
#endif

		halos[i].gatherRad        = 100000000000.0;

		halos[i].numSubStruct     = 0;
		halos[i].subStruct        = NULL;

		halos[i].spaRes           = 0.0;
		halos[i].refLev           = 0;

		halos[i].numNodes         = 0;
	}

	/*********************************************************************************************************
	 *********************************************************************************************************
	 * connecting the halos
	 * NOTE :: IF 'halos[i].hostHalo = i' then it is not the substrcture of any other halo */
	count    = 0;
	tmpCount = 0;
	countBC  = 0;
#ifdef VERBOSE
	fprintf(stderr, "  constructing %d (potential) halos from all %d grid levels...\n", numHalos, num_refgrids);
	fprintf(io.logfile, "  constructing %d (potential) halos from all %d grid levels...\n", numHalos, num_refgrids);
	fflush(io.logfile);
#endif
  
  /*----------------------------------------------------------------------------------------
   * loop over all refinement levels building the halos[] tree from the spatialRef[][] tree
   *----------------------------------------------------------------------------------------*/
  //tdummy  = 0;
  //tdummy -= time(NULL);
  
  /* we might not start on level 0 as this is the level outside the ovlim criterion
   * it makes more sense to start the halo-tree on a level "inside" the virial radius of haloes? */
  refgrid_start = 0;
  
	for (i = refgrid_start; i < num_refgrids; i++) {
#ifdef VERBOSE
		fprintf(stderr, "      grid level %8d (%12d) -> %10d isolated refinements ... ", i, (int)gridl1dim[i], numIsoRef[i]);
		fprintf(io.logfile, "      grid level %8d (%12d) -> %10d isolated refinements ... ", i, (int)gridl1dim[i], numIsoRef[i]);
		fflush(io.logfile);
#endif

		/* treat 1st grid separately
     *
     * the 1st grid encompasses ovlim and hence no halo on this level should be a subalo
     *    -> therefore, all halos[].hostHalo will be initialized to -1
     */
   
		if (i == refgrid_start) {
			for (j = 0; j < numIsoRef[i]; j++) {
				/******************************************************/
				if (spatialRef[i][j].numSubStruct == 0) {   // the halo ends on this level
#ifdef VERBOSE2
          fprintf(stderr,"%d(0),count=%d ",j,count);
#endif
					/* there is no host halo */
					halos[count].hostHalo = -1;
          
					/* halo centre (there is no finer spatialRef[][] and hence assign centre...) */
					halos[count].pos.x = spatialRef[i][j].centre.x;
					halos[count].pos.y = spatialRef[i][j].centre.y;
					halos[count].pos.z = spatialRef[i][j].centre.z;

					/* The spatial resolution of this halo */
					halos[count].spaRes = 1.0 / ((double)gridl1dim[i]);
					halos[count].refLev = i;

					/* particles */
					halos[count].npart  = spatialRef[i][j].numParts;
          
					/*  Substructure */
					halos[count].numSubStruct = 0;
					halos[count].subStruct    = NULL;

					/* nodes */
					halos[count].numNodes = spatialRef[i][j].numNodes;

					count++;
					tmpCount++;

					/******************************************************/
				} else if (spatialRef[i][j].numSubStruct == 1) {  // there is only the halo that continues downwards
#ifdef VERBOSE2
          fprintf(stderr,"%d(1) ",j);
#endif
					/* there is no host halo */
					halos[count].hostHalo = -1;
          
					/* do not assign centre as there is a finer refinement... */

					/* particles */
					halos[count].npart  = spatialRef[i][j].numParts;

					/* nodes */
					halos[count].numNodes = spatialRef[i][j].numNodes;

					/*  Substructure (the one substructure on spatialRef[][] is in fact the host!) */
					halos[count].numSubStruct = 0;
					halos[count].subStruct    = NULL;

					/* Telling daughter about host halo */
					refLevel    = spatialRef[i][j].daughter.refLevel;
					isoRefIndex = spatialRef[i][j].daughter.isoRefIndex;

					if ((refLevel != -1) && (isoRefIndex != -1))
						spatialRef[refLevel][isoRefIndex].haloIndex = count;
#ifdef VERBOSE2
					else
           {
						fprintf(stderr, "      spatialRef2halos():  i==refgrid_start, numSubStruct==1\n", i, j);
						fprintf(stderr, "                           WARNING -> no daughter for spatialRef[%d][%d]\n", i, j);
           }
#endif
					count++;

					/******************************************************/
				} else {     // the spatialRef[][] tree spreads into multiple leaves
#ifdef VERBOSE2
          fprintf(stderr,"%d(%d)\n",j,spatialRef[i][j].numSubStruct);
#endif
					/* there is no host halo */
					halos[count].hostHalo = -1;

					/* count;  If the host halo name is the same as the iterator then it is a base halo */
					primHaloIndex = count;

					/* Telling daughter about host halo */
					refLevel    = spatialRef[i][j].daughter.refLevel;
					isoRefIndex = spatialRef[i][j].daughter.isoRefIndex;

					if ((refLevel != -1) && (isoRefIndex != -1))
						spatialRef[refLevel][isoRefIndex].haloIndex = count;
#ifdef VERBOSE2
					else
           {
						fprintf(stderr, "      spatialRef2halos():  i==refgrid_start, numSubStruct>1\n", i, j);
						fprintf(stderr, "                           WARNING -> no daughter for spatialRef[%d][%d]\n", i, j);
           }
#endif

					/* nodes */
					halos[count].numNodes = spatialRef[i][j].numNodes;

					/*  Substructure */   
					numNewHalos               = spatialRef[i][j].numSubStruct - 1; // -1 because the daughter is also in the SubStruct list
					halos[count].numSubStruct = numNewHalos;
          if(halos[count].subStruct) free(halos[count].subStruct);
					halos[count].subStruct    = NULL;
					if ((halos[count].subStruct = calloc(numNewHalos+1, sizeof(int))) == NULL) // +1 to be on the safe side
           {
						fprintf(stderr,"Error in allocating the memory for halos[count].subStruct array\n");
						exit(0);
           }

					/* increment halo counter as we now... */
					count++;

					/* ...generate the substructure halos of "primHaloIndex" (=the host) */
					kcount = 0; // we loop over all SubStruct, but this also includes the daughter not to be counted in the end
					for (k = 0; k < numNewHalos + 1; k++) {  // loop over "+1" as numNewHalos did not count the daughter
						SSrefLevel    = spatialRef[i][j].subStruct[k].refLevel;
						SSisoRefIndex = spatialRef[i][j].subStruct[k].isoRefIndex;

						/* Make sure not the daughter */
						if (SSisoRefIndex != isoRefIndex) {
							/* What is my host halo ? */
							halos[count].hostHalo      = primHaloIndex;
              halos[count].hostHaloLevel = i;

							/* Telling the subHalos who they are */
							spatialRef[SSrefLevel][SSisoRefIndex].haloIndex = count;

							/* Adding that name to the list of substrcture */
							halos[primHaloIndex].subStruct[kcount] = count;

							/* initial guess for virial radius */
							halos[count].R_vir = spatialRef[SSrefLevel][SSisoRefIndex].closeRefDist;

							/* subhalo has no particles yet */
							halos[count].ll = NULL;

							count++;  // move to next (sub-)halo
							kcount++; // kcount should add up to (spatialRef[i][i].numSubStruct-1)
						}
					}

					/* Assign the remaining particles to the host halo */
					halos[primHaloIndex].npart = spatialRef[i][j].numParts;

        }
			} /* for(j<numIsoRef[i]) */
		}
    
    
    /*------------------------------------------------
		 * the refinement levels (i.e. i > refgrid_start)
     *------------------------------------------------*/
		else {
			for (j = 0; j < numIsoRef[i]; j++) {
        
				/* there is a parent refinement level */
				if (spatialRef[i][j].numParDom != 0) {

					if (spatialRef[i][j].numSubStruct == 0) {
            
						/* haloIndex is the ID of the halo correspondong to the parent refinement level,
             * i.e. we adjust its properties as we have refined information from the present refinement level */
						haloIndex = spatialRef[i][j].haloIndex;

						if (haloIndex == -1) {
							fprintf(stderr,"  HELP haloIndex HELP :: HALO[%d][%d] \n", i, j);
							fprintf(stderr,"  HALO-dau[%d][%d] \n", spatialRef[i][j].daughter.refLevel, spatialRef[i][j].daughter.isoRefIndex);
							fprintf(stderr,"  HALO-par[%d][%d] \n", spatialRef[i][j].parDom[0].isoRefIndex, spatialRef[i][j].parDom[0].refLevel);
              exit(-1);
						}

						/* The spatial resolution of this halo :- to finish ALEX */
						halos[haloIndex].spaRes = 1.0 / ((double)gridl1dim[i]);
						halos[haloIndex].refLev = i;

						/* halo centre -> change to new centre of this finer refinement! */
						halos[haloIndex].pos.x = spatialRef[i][j].centre.x;
						halos[haloIndex].pos.y = spatialRef[i][j].centre.y;
						halos[haloIndex].pos.z = spatialRef[i][j].centre.z;

						/* nodes       -> record number of nodes centre is based upon */
						halos[haloIndex].numNodes = spatialRef[i][j].numNodes;

						/* particles   -> add additional particles */
						halos[haloIndex].npart += spatialRef[i][j].numParts;

						tmpCount++; /* Closing of the halo */


					} else if (spatialRef[i][j].numSubStruct == 1) {

						/* haloIndex is the ID of the halo correspondong to the parent refinement level,
             * i.e. we adjust its properties as we have refined information from the present refinement level */
						haloIndex = spatialRef[i][j].haloIndex;

						if (haloIndex == -1)
             {
							fprintf(stderr, "  HELP2 haloIndex HELP2\n");
              exit(-1);
             }

						/* particles */
						halos[haloIndex].npart += spatialRef[i][j].numParts;

						/* nodes */
						halos[haloIndex].numNodes = spatialRef[i][j].numNodes;

						/* Telling daughter about host halo */
						refLevel    = spatialRef[i][j].daughter.refLevel;
						isoRefIndex = spatialRef[i][j].daughter.isoRefIndex;

						if ((refLevel != -1) && (isoRefIndex != -1))
							spatialRef[refLevel][isoRefIndex].haloIndex = haloIndex;
#ifdef VERBOSE2
						else
							fprintf(stderr,"      spatialRef2halos(1):  WARNING -> wrong daughter for spatialRef[%d][%d]\n",i,j);
#endif

					} else {

						/* haloIndex is the ID of the halo correspondong to the parent refinement level,
             * i.e. we adjust its properties as we have refined information from the present refinement level */
						haloIndex     = spatialRef[i][j].haloIndex;
						primHaloIndex = haloIndex;

						if (haloIndex == -1)
             {
							fprintf(stderr, "  HELP3 haloIndex HELP3\n");
              exit(-1);
             }

						/* nodes */
						halos[haloIndex].numNodes = spatialRef[i][j].numNodes;
            
            /* BUILD 039: quick-and-dirty fix in case this halo did not ever get a position assigned */
            // BUILD 084: we should actually always update the halo position to the best resolution
            //if(halos[haloIndex].pos.x < MACHINE_ZERO && halos[haloIndex].pos.y < MACHINE_ZERO && halos[haloIndex].pos.z < MACHINE_ZERO)
             {
              halos[haloIndex].pos.x = spatialRef[i][j].centre.x;
              halos[haloIndex].pos.y = spatialRef[i][j].centre.y;
              halos[haloIndex].pos.z = spatialRef[i][j].centre.z;
             }

						/* Telling daughter about host halo */
						refLevel    = spatialRef[i][j].daughter.refLevel;
						isoRefIndex = spatialRef[i][j].daughter.isoRefIndex;

						if ((refLevel != -1) && (isoRefIndex != -1))
							spatialRef[refLevel][isoRefIndex].haloIndex = haloIndex;
#ifdef VERBOSE2
						else
							fprintf(stderr,"      spatialRef2halos(2):  WARNING -> wrong daughter for spatialRef[%d][%d]\n",i,j);
#endif


						/*  Substructure */
						numNewHalos  = spatialRef[i][j].numSubStruct - 1;
						kcount       = halos[primHaloIndex].numSubStruct; // this is the number of old substructure halos serving as a loop counter below
          
						halos[primHaloIndex].numSubStruct = kcount + numNewHalos;

						if (halos[primHaloIndex].subStruct == NULL) {
							if ((halos[primHaloIndex].subStruct = calloc(halos[primHaloIndex].numSubStruct + 1, sizeof(int))) == NULL) {
								fprintf(stderr,"Error in allocating the memory for halos[count].subStruct array\n");
								exit(0);
							}
						} else {
							if ((halos[primHaloIndex].subStruct = realloc(halos[primHaloIndex].subStruct,(halos[primHaloIndex].numSubStruct + 1) * sizeof(int))) == NULL) {
								fprintf(stderr,"Error in reallocating the memory for halos[count].subStruct array\n");
								exit(0);
							}
						}
            
						/* Collecting information for the 'new' substructure */
						for (k = 0; k < spatialRef[i][j].numSubStruct; k++) {    // note that the substructure loop counter is kcount!
							SSrefLevel    = spatialRef[i][j].subStruct[k].refLevel;
							SSisoRefIndex = spatialRef[i][j].subStruct[k].isoRefIndex;

							if (SSisoRefIndex != isoRefIndex) {
								/* What is my host halo ? */
								halos[count].hostHalo      = primHaloIndex;
                halos[count].hostHaloLevel = i;

								/* Telling the subHalos who they are */
								spatialRef[SSrefLevel][SSisoRefIndex].haloIndex = count;

								/*  Substructure */
								halos[count].subStruct = NULL;

								/*  Substructure of primary index */
								halos[primHaloIndex].subStruct[kcount] = count;

								/* first guess for subhalo's virial radius */
								halos[count].R_vir = spatialRef[SSrefLevel][SSisoRefIndex].closeRefDist;

								/* subhalo has no particles yet */
								halos[count].ll = NULL;

								count++;
								kcount++;
							}
						}

						/* Assign the remaining particles to the host halo */
						halos[haloIndex].npart += spatialRef[i][j].numParts;            

          }
				} else {
					OKcountBC = 1;
					for (f = 0; f < numDensZero; f++) {
						if ((i == densZero[f].refLevel) && (j == densZero[f].isoRefIndex))
							OKcountBC = 0;
					}
#ifdef VERBOSE
					fprintf(stderr, "      s %g %g %g 0.08 0.0 1.0 1.0\n",
					        spatialRef[i][j].centreGEOM.x * simu.boxsize,
					        spatialRef[i][j].centreGEOM.y * simu.boxsize,
					        spatialRef[i][j].centreGEOM.z * simu.boxsize);
					fprintf(stderr, "      countBC[%d][%d]\n", i, j);
					fprintf(stderr, "      OKcountBC = %d (0 is good)\n", OKcountBC);
#endif

					countBC++;
				}
			} /* for(j<numIsoRef[i]) */
		} /* i !=0 */
    
#ifdef VERBOSE
		fprintf(stderr, "finished\n");
		fprintf(io.logfile, "finished\n");
		fflush(io.logfile);
#endif

	} /* for(i<num_refgrids) */
  //tdummy += time(NULL);
  //fprintf(stderr,"timing: %ld\n",tdummy);
  //exit(0);

#ifdef AHFDEBUG
	fprintf(stderr, "done\n");
#endif

#ifdef VERBOSE
	if (countBC > 0)
   {
		fprintf(stderr,"      countBC = %d => isolated refinements do not have a parent refinement\n",countBC);
		fprintf(io.logfile,"      countBC = %d => isolated refinements do not have a parent refinement\n",countBC);

   }
	fprintf(stderr,"      number of halos: found=%d (max=%d) expected=%d (tmpCount=%d)\n",count,totnumIsoRef,numHalos,tmpCount);
	fprintf(stderr, "      numDensZero = %d, numPartZero = %d\n", numDensZero, numPartZero);

	fprintf(io.logfile,"      number of halos: found=%d (max=%d) expected=%d (tmpCount=%d)\n",count,totnumIsoRef,numHalos,tmpCount);
	fprintf(io.logfile, "      numDensZero = %d, numPartZero = %d\n", numDensZero, numPartZero);
	fflush(io.logfile);
#endif /* VERBOSE */

	/* did we find more halos than initially expected? */
	if (count > numHalos) {
		numHalos = count;
#ifdef VERBOSE
		fprintf(stderr, "      => adjusted number of halos: found=%d (totnumIsoRef=%d) expected=%d (tmpCount=%d)\n", count, totnumIsoRef, numHalos,	tmpCount);
#endif
		fprintf(io.logfile, "      => adjusted number of halos: found=%d (totnumIsoRef=%d) expected=%d (tmpCount=%d)\n", count, totnumIsoRef, numHalos, tmpCount);
		fflush(io.logfile);
	}

  /* make this first guess of haloes globally accessible */
  simu.no_halos = numHalos;
  
	/*********************************************************************************************************
	 *********************************************************************************************************
	 * Ordering the halos wrt mass
	 * The assumption is that sub-halos always have less mass than their host */
  /* Order halos with respects to their mass */
	if (numHalos > 0) {
		idx    = (long unsigned *)calloc(numHalos, sizeof(long unsigned));
		idxtmp = (long unsigned *)calloc(numHalos + 1, sizeof(long unsigned));
		fsort  = (double *)calloc(numHalos + 1, sizeof(double));
		for (i = 0; i < numHalos; i++)
			fsort[i+1] = (double)halos[i].npart;
		indexx(numHalos, fsort, idxtmp);
    
		/* indexx sorts ascending and gives indizes starting at 1 */
		for (i = 0; i < numHalos; i++)
			idx[numHalos - i - 1] = idxtmp[i + 1] - 1;
    
		free(idxtmp);
		free(fsort);
	} else {
		/* If there are no halos, have a proper empty idx array */
		idx = NULL;
	}
 
	/*********************************************************************************************************
	 *********************************************************************************************************
	 * Calculate the Gathering radius:
	 *
	 * the gathering radius is the distance to the closest halo that is more massive than the current halo...
   * we therefore require the halos to be ordered by npart for this part!
   */
   
	maxGathRad = MIN(simu.MaxGatherRad / simu.boxsize, 1. / 4.);

#ifdef WITH_OPENMP
#      pragma omp parallel private(ii,jj,i,j,count,inumpart,jnumpart,dx,dy,dz,tmpRad) shared(idx,halos,numHalos,maxGathRad)
#      pragma omp for schedule(dynamic)
#endif
	for (ii = numHalos - 1; ii >= 0; ii--) {
    
    /* pick halo from ordered list */
    i = idx[ii];
    
		/* root halos do not need to gather anything from larger objects... */
    count    = 0;

    inumpart = halos[i].npart;
    
    /* loop over all halos that are more massive... */
    //for (jj = numHalos - 1; jj >= 0; jj--) {
    for (jj = ii; jj >= 0; jj--) {
      
      /* pick halo from ordered list */
      j = idx[jj];
      
      jnumpart = halos[j].npart;
      
      if ((i != j) && (jnumpart > inumpart)) {
        /* Calculate the distance to this more massive halo */
        dx = fabs(halos[i].pos.x - halos[j].pos.x);
        dy = fabs(halos[i].pos.y - halos[j].pos.y);
        dz = fabs(halos[i].pos.z - halos[j].pos.z);
        if (dx > 0.5)		dx = 1.0 - dx;
        if (dy > 0.5) 	dy = 1.0 - dy;
        if (dz > 0.5) 	dz = 1.0 - dz;
        
        tmpRad = dx * dx + dy * dy + dz * dz;
        
        /* Is it the smallest distance */
        if (tmpRad < halos[i].gatherRad) {
          halos[i].gatherRad = tmpRad;
        }
        
        count++;
      }
    }
    
    /* Finally calculating the gathering radius */
    halos[i].gatherRad = (sqrt(halos[i].gatherRad)) * 0.5;
    
    if (count == 0) {
#ifdef AHFDEBUG
      fprintf(stderr, "There are no other halos bigger than this halo - should see this just once\n");
#endif
      halos[i].gatherRad = maxGathRad;
    }
    
#ifdef AHFmaxGatherRadTest
    /* for host haloes we simply estimate the radius according to our desired virial overdensity criterion */
    if(halos[i].hostHalo < 0)
      halos[i].gatherRad = pow(halos[i].npart*simu.pmass*rho_fac/(4.*PI/3. * global.ovlim*global.rho_vir),0.333333333333333333);
#endif

		/* restrict the gathering radius (Note: R_vir at this stage is closeRefDist) */
		halos[i].gatherRad = MAX(halos[i].gatherRad, halos[i].R_vir);
		halos[i].gatherRad = MIN(halos[i].gatherRad, maxGathRad);
	} /* for ( i=numHalos-1; i>=0; i-- ) */

#ifdef DPhalos
	/*********************************************************************************************************
	 *********************************************************************************************************
	 * Dump preliminary halo properties to file -> these objects are our DPhaloes (Avila Perez et al., in prep.) */
 {
  FILE *fptmp;
  char  prename[MAXSTRING];
  
  
#ifdef WITH_MPI
  fprintf(stderr,"\n MPI rank %d writing DPhalos ... ",global_mpi.rank);
	sprintf(prename,"%s.%04d.z%.3f.AHF_DPhalos", global_io.params->outfile_prefix,global_mpi.rank,fabs(global.z));
  
  if( (fptmp = fopen(prename,"w")) == NULL ) {
    fprintf(stderr,"Could not open DPhalos file %s\nABORTING\n",prename);
    common_terminate(EXIT_FAILURE);
  }

  if(global_mpi.rank == 0)
    fprintf(fptmp,"# x(1) y(2) z(3) npart(4) ncells(5)\n");

  /* Flagging the haloes in the boundary zone */
	flag_boundary_haloes();
#else
  fprintf(stderr,"\n Writing DPhalos ... ");
  sprintf(prename,"%s.z%.3f.AHF_DPhalos",      global_io.params->outfile_prefix,                fabs(global.z));
  if( (fptmp = fopen(prename,"w")) == NULL ) {
    fprintf(stderr,"Could not open DPhalos file %s\nABORTING\n",prename);
    common_terminate(EXIT_FAILURE);
  }
  fprintf(fptmp,"# x(1) y(2) z(3) npart(4) ncells(5)\n");
#endif
  
              
#ifdef DPhalos_WITHCOMMENTS
	fprintf(fptmp, "\n\n   preliminary halo properties:\n");
	fprintf(fptmp, "   ============================\n");
#endif
  for (ii = numHalos - 1; ii >= 0; ii--) {
    /* pick halo from ordered list */
    i = idx[ii];

#ifdef WITH_MPI
    if(halos[i].ignoreme == FALSE) {
#endif
      
#ifdef DPhalos_WITHCOMMENTS
      fprintf(fptmp, "   halos[%d].numNodes         = %16d\n",                          i, halos[i].numNodes);
      fprintf(fptmp, "   halos[%d].npart            = %16ld\n",                         i, halos[i].npart);
      fprintf(fptmp, "   halos[%d].R_vir            = %16.8g kpc/h (=closeRefDist)\n",  i, halos[i].R_vir * x_fac * 1000.);
      fprintf(fptmp, "   halos[%d].gatherRad        = %16.8g kpc/h (=closeHostDist)\n", i, halos[i].gatherRad * x_fac * 1000.);
      fprintf(fptmp, "   halos[%d].spaRes           = %16.8g\n",                        i, halos[i].spaRes);
      fprintf(fptmp, "   halos[%d].pos.x            = %16.8g\n",	                      i, halos[i].pos.x * x_fac);
      fprintf(fptmp, "   halos[%d].pos.y            = %16.8g\n",		                    i, halos[i].pos.y * x_fac);
      fprintf(fptmp, "   halos[%d].pos.z            = %16.8g\n",		                    i, halos[i].pos.z * x_fac);
      fprintf(fptmp, "   halos[%d].hostHalo         = %16d\n",		                      i, halos[i].hostHalo);
      fprintf(fptmp, "   halos[%d].hostHaloLevel    = %16d\n",		                      i, halos[i].hostHaloLevel);
      fprintf(fptmp, "   halos[%d].numSubStruct     = %16d\n",		                      i, halos[i].numSubStruct);
      fprintf(fptmp, "            gatherRad/R_vir   = %16.8g\n",		                    halos[i].gatherRad / halos[i].R_vir);
      fprintf(fptmp, "   #########################################################################\n");
#else
      if(halos[i].npart > 0) {
        fprintf(fptmp,"%16.8f %16.8f %16.8f %16ld %16d\n",
                halos[i].pos.x     * x_fac * 1000,
                halos[i].pos.y     * x_fac * 1000,
                halos[i].pos.z     * x_fac * 1000,
                halos[i].npart,
                halos[i].numNodes);
        fflush(fptmp);
      }
#endif
      
#ifdef WITH_MPI
    } // if(ignoreme == FALSE)
#endif
	} // for(ii=numHalos)
  fclose(fptmp);
 
  // we are not interested in anything else, right?
  fprintf(stderr,"done and exiting now.\n");
  common_terminate(EXIT_SUCCESS);
 }
#endif   /* DPhalos */

  /* remove tempory index array again */
  if(idx != NULL)
    free(idx);
  
	return TRUE;
} /* spatialRef2halos */


/*
 ************************************************************
 ************************************************************
 * Writes information about the centres of the whole tree
 */
int WriteGridtreefile(const char *fprefix, int num_refgrids, SPATIALREF **spatialRef)
{
  char filename[MAXSTRING];
  FILE *fout;
  long i,j,k;
  
  /* Generate the filename */
  strcpy(filename, fprefix);
  strcat(filename, ".AHF_gridtree");
#  ifdef VERBOSE
  fprintf(stderr, "%s\n", filename);
#  endif
  
  /* Open file */
  if ((fout = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "could not open %s\n", filename);
    exit(1);
  }
  
  fprintf(fout,"%d %d\n",ahf.min_ref, ahf.no_grids);
  
  for(i=0; i<num_refgrids; i++) {
    
    fprintf(fout,"%d %d\n",
            i + ahf.min_ref,  // write the actual level (and not the internal count)
            numIsoRef[i]);
    
    for (j=0; j<numIsoRef[i]; j++) {
      fprintf(fout, "%18.14lf %18.14lf %18.14lf %18.14lf %d %ld %d %d %d\n",
              spatialRef[i][j].centre.x,
              spatialRef[i][j].centre.y,
              spatialRef[i][j].centre.z,
              spatialRef[i][j].closeRefDist,
              spatialRef[i][j].numNodes,
              spatialRef[i][j].numParts,
              spatialRef[i][j].daughter.refLevel            + ahf.min_ref,  // write the actual level (and not the internal count)
              spatialRef[i][j].daughter.isoRefIndex,
              spatialRef[i][j].numSubStruct);
      
      if(spatialRef[i][j].daughter.refLevel != i+1 && spatialRef[i][j].daughter.isoRefIndex != -1) {
        fprintf(stderr,"spatialRef[%d][%d].daughter.refLevel=%d != %d+1 (isoRefIndex=%d)\n",i,j,spatialRef[i][j].daughter.refLevel,i,spatialRef[i][j].daughter.isoRefIndex);
        exit(0);
      }

      for(k=0; k<spatialRef[i][j].numSubStruct; k++) {
        fprintf(fout,"   %d %d\n",
                spatialRef[i][j].subStruct[k].refLevel      + ahf.min_ref,  // write the actual level (and not the internal count)
                spatialRef[i][j].subStruct[k].isoRefIndex);
        
        if(spatialRef[i][j].subStruct[k].refLevel != i+1 && spatialRef[i][j].subStruct[k].isoRefIndex != -1) {
          fprintf(stderr,"spatialRef[%d][%d].subStruct[%d].refLevel=%d != %d+1 (isoRefIndex=%d)\n",i,j,k,spatialRef[i][j].subStruct[k].refLevel,i,spatialRef[i][j].daughter.isoRefIndex);
          exit(0);
        }
      }
    }
  }
  
  /* Clean up */
  fclose(fout);
  
  /* Done */
  return(1);  
}  

/*
 ************************************************************
 ************************************************************
 * Orders the particels within the halos
 */

int
compare(struct particle *p, struct particle *q, XYZ *pointer)
{
	double x1, y1, z1, dist1;
	double x2, y2, z2, dist2;
	double xx, yy, zz;
	double dx1, dy1, dz1;
	double dx2, dy2, dz2;

	xx  = pointer->x;
	yy  = pointer->y;
	zz  = pointer->z;

	x1  = p->pos[0];
	y1  = p->pos[1];
	z1  = p->pos[2];

	dx1 = xx - x1;
	dx1 = fabs(dx1);
	if (dx1 > 0.5)
		dx1 = 1.0 - dx1;

	dy1 = yy - y1;
	dy1 = fabs(dy1);
	if (dy1 > 0.5)
		dy1 = 1.0 - dy1;

	dz1 = zz - z1;
	dz1 = fabs(dz1);
	if (dz1 > 0.5)
		dz1 = 1.0 - dz1;

	x2  = q->pos[0];
	y2  = q->pos[1];
	z2  = q->pos[2];

	dx2 = xx - x2;
	dx2 = fabs(dx2);
	if (dx2 > 0.5)
		dx2 = 1.0 - dx2;

	dy2 = yy - y2;
	dy2 = fabs(dy2);
	if (dy2 > 0.5)
		dy2 = 1.0 - dy2;

	dz2 = zz - z2;
	dz2 = fabs(dz2);
	if (dz2 > 0.5)
		dz2 = 1.0 - dz2;


	dist1
	    = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1) * simu.boxsize * simu.boxsize;
	dist2
	    = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2) * simu.boxsize * simu.boxsize;


	if (dist1 < dist2)
		return -1;
	else if (dist1 > dist2)
		return 1;
	else
		return 0;
} /* compare */



/*
 ***********************************************************************
 ***********************************************************************
 * fills those values in that normally rem_unbound() would fill in...
 */
void
rem_nothing(HALO *halo)
{
	long unsigned jpart, npart;
	double        Xc, Yc, Zc, Xp, Yp, Zp, dX, dY, dZ, dist2;
	partptr       cur_part;
	double        R_vir, M_vir, weight;


	if (halo->npart < simu.AHF_MINPART)
		return;

#ifdef VERBOSE
	fprintf(io.logfile,
	        "    rem_nothing:      npart=%12ld -> ",
	        halo->npart);
	fflush(io.logfile);
#endif


	Xc    = halo->pos.x;
	Yc    = halo->pos.y;
	Zc    = halo->pos.z;
	R_vir = -1.0;
	M_vir = 0.0;
	npart = 0;

	/* loop over all particles... */
	for (jpart = 0; jpart < halo->npart; jpart++) {
		/* access particle */
		cur_part = global.fst_part + halo->ipart[jpart];

#ifdef MULTIMASS
		weight = (double)cur_part->weight;
#else
		weight = (double)1.0;
#endif

		/* cumulative mass */
		M_vir += weight;
		npart++;

		/* particle position */
		Xp = (double)cur_part->pos[X];
		Yp = (double)cur_part->pos[Y];
		Zp = (double)cur_part->pos[Z];

		/* put particle into halo rest frame */
		dX = fabs(Xp - Xc);
		dY = fabs(Yp - Yc);
		dZ = fabs(Zp - Zc);

		/* take care of periodic boundary conditions */
		if (dX > 0.5)			dX -= 1.0;
		if (dY > 0.5)			dY -= 1.0;
		if (dZ > 0.5)			dZ -= 1.0;

		/* finally calculate distance squared */
		dist2 = pow2(dX) + pow2(dY) + pow2(dZ);

		/* remember distance of current particle (particles are ordered
		 *distancewise!) */
		if (dist2 > R_vir)
			R_vir = dist2;
	} /* particle for loop */

	halo->npart = npart;
	halo->M_vir = M_vir;
	halo->R_vir = sqrt(R_vir);
	halo->Phi0  = 0.0;

#ifdef VERBOSE
	fprintf(io.logfile, "%12ld\n", halo->npart);
	fflush(io.logfile);
#endif
} /* rem_nothing */

/*
 ************************************************************
 ************************************************************
 * Removes the unbound particles from the halos
 */
void
rem_unbound(HALO *halo)
{
	long unsigned nremove, npart_old;
	double        Xc, Yc, Zc, VXc, VYc, VZc;
	double        Xp, Yp, Zp, VXp, VYp, VZp;
	double        dX, dY, dZ, dVX, dVY, dVZ;
	double        v2_tune, weight;
	double        Phi, Phi0, Phi_infty, M_r, M_vir, M_vel, vel2, v_esc2,
	              d_prev, R_vir;
	double        I_now, I_mid, I_prev, dr;
	partptr       cur_part, pre_part, host_part, tmp_part;
	gasptr        cur_gas;
	double        x, y, z, dx, dy, dz, distA, distB, dist;
	double        xx, yy, zz;
	int           niter;
	int           hostHalo;
	long          jpart;
	long unsigned *bound_ipart, bound_npart;
	double        *mom2;
	long unsigned *idx, no_vbulk;

	if (halo->npart < simu.AHF_MINPART)
		return;

#ifdef VERBOSE2
	fprintf(io.logfile,
	        "    rem_unbound:         npart=%12ld -> ",
	        halo->npart);
	fflush(io.logfile);
#endif

	/* velocity tune parameter */
	v2_tune = pow2(simu.AHF_VTUNE);

	/* how many central particles to use for the initial bulk velocity guess */
	no_vbulk = (long unsigned)(simu.AHF_MINPART / 2);

	/* remember initial number of particles
	 *  npart_old = halos[i].npart; */

	/* halo centre in AMIGA units */
	Xc = halo->pos.x;
	Yc = halo->pos.y;
	Zc = halo->pos.z;

	nremove = 4;
	niter   = 0;

	/*----------------------------------------------------------------------------
	 * remove particles until all particles are bound or halo mass gets too
	 *small
	 *----------------------------------------------------------------------------*/
	while ((nremove > 3)) {
		/* iteration counter */
		niter++;

		/*---------------------------------------------------------------
		 * determine Phi0  (the zero point of the potential is infinity)
		 *---------------------------------------------------------------*/
		/* reset values */
    I_now  = 0.0;
		I_prev = 0.0;
		d_prev = 0.0;
		M_r    = 0.0;
		Phi0   = 0.0;

		/***********************************************************
		 * loop over all sorted particles from inside out to calculate Phi0 */
		for (jpart = 0; jpart < halo->npart; jpart++) {
			/* access particle */
			cur_part = global.fst_part + halo->ipart[jpart];

#ifdef MULTIMASS
			weight = (double)cur_part->weight;
#else
			weight = (double)1.0;
#endif

			/* cumulative mass */
			M_r += weight;

			/* particle position */
			Xp = (double)cur_part->pos[X];
			Yp = (double)cur_part->pos[Y];
			Zp = (double)cur_part->pos[Z];

			/* put particle into halo rest frame */
			dX = fabs(Xp - Xc);
			dY = fabs(Yp - Yc);
			dZ = fabs(Zp - Zc);

			/* take care of periodic boundary conditions */
			if (dX > 0.5)	dX -= 1.0;
			if (dY > 0.5)	dY -= 1.0;
			if (dZ > 0.5)	dZ -= 1.0;

			/* finally calculate distance (conversion to dist_phys=a*dist via phi_fac!)  */
			dist = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));

#ifdef AHFDEBUG
      if(dist<d_prev)
        fprintf(io.logfile,"FATAL in rem_unbound(): particles do not appear to be ordered by distance! dist=%g < d_prev=%g\n",
                dist*x_fac*1000.,d_prev*x_fac*1000.);
      //exit(0);
#endif
      
			/* accumulate Phi0 */
			if (dist > MACHINE_ZERO) {
				/* mid-point integration */
				I_now = M_r / pow2(dist);
#ifdef AHF_LRSI
				/* add the effect of scalar forces */
				I_now*=(1.0+ simu.lrsi_beta*exp(-1.0*x_fac*dist/simu.lrsi_r_s));
#endif
				I_mid = (I_now + I_prev) / 2.;
				dr    = dist - d_prev;

				/* accumulate Phi0 */
				Phi0 += I_mid * dr;
			}

			d_prev = dist;
			I_prev = I_now;
		} /* particle loop for Phi0 determination */

		/* finally calculate Phi0 */
#ifdef AHF_LRSI
		/* add the effect of scalar forces */
		Phi0 += M_r/dist*(1.0+ simu.lrsi_beta*exp(-1.0*x_fac*dist/simu.lrsi_r_s));
#endif
    
    /* do not forget the GM/R term */
		Phi0 += M_r / dist;

		/* remember Phi0 as it will be used when calculating halos[].prof.Epot */
		halo->Phi0 = Phi0;

		/*--------------------------------------------------------------------
		 * determine Phi             ( v_esc^2 = 2 * |Phi| )
		 *--------------------------------------------------------------------*/
		/* reset values */
		nremove = 0;
    I_now   = 0.0;
		I_prev  = 0.0;
		d_prev  = 0.0;
		M_r     = 0.0;
		M_vir   = 0.0;
		Phi     = 0.0;

		/* initial bulk velocity of halo */
		if (niter == 1) {
			mom2 = (double *)calloc(no_vbulk + 1, sizeof(double));
			idx  = (long unsigned *)calloc(no_vbulk + 1, sizeof(long unsigned));

			for (jpart = 0; jpart < no_vbulk; jpart++) {
				cur_part        = global.fst_part + halo->ipart[jpart];
				mom2[jpart + 1] = pow2(cur_part->mom[X]) + pow2(cur_part->mom[Y]) + pow2(cur_part->mom[Z]);
			}
			indexx((long unsigned)(no_vbulk), mom2, idx);

			/* use the median of the innermost particle's velocity */
			jpart = idx[(int)(no_vbulk / 2)] - 1;

			free(mom2);
			free(idx);
		} else {
			/* use the most central particle's velocity */
			jpart = 0;
		}

		cur_part = global.fst_part + halo->ipart[jpart];
#ifdef MULTIMASS
		weight   = (double)cur_part->weight;
#else
		weight   = (double)1.0;
#endif
		M_vel    = weight;
		VXc      = weight * cur_part->mom[X];
		VYc      = weight * cur_part->mom[Y];
		VZc      = weight * cur_part->mom[Z];

		/***********************************************************
		 * loop over all sorted particles from inside out */
		bound_npart = 0;
		bound_ipart = calloc(1, sizeof(long unsigned));  // some realloc()'s do
		                                                 // not like NULL
		                                                 // pointers...
		for (jpart = 0; jpart < halo->npart; jpart++) {
			/* access particle */
			cur_part = global.fst_part + halo->ipart[jpart];

#ifdef MULTIMASS
			weight = (double)cur_part->weight;
#else
			weight = (double)1.0;
#endif
			/* cumulative mass */
			M_r += weight;

			/* particle position */
			Xp = (double)cur_part->pos[X];
			Yp = (double)cur_part->pos[Y];
			Zp = (double)cur_part->pos[Z];

			/* put particle into halo rest frame :: *no* fabs() this time ! */
			dX = (Xp - Xc);
			dY = (Yp - Yc);
			dZ = (Zp - Zc);

			/* take care of periodic boundary conditions */
			if (dX > 0.5)				dX -= 1.0;
			if (dY > 0.5)				dY -= 1.0;
			if (dZ > 0.5)				dZ -= 1.0;
			if (dX < -0.5)			dX += 1.0;
			if (dY < -0.5)			dY += 1.0;
			if (dZ < -0.5)			dZ += 1.0;

			/* finally calculate distance (conversion to dist_phys=a*dist via phi_fac!) */
			dist = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));

			/* get potential escape velocity */
			if (dist > MACHINE_ZERO) {
				/* mid-point integration */
				I_now = M_r / pow2(dist);
#ifdef AHF_LRSI
				I_now*=(1.0+ simu.lrsi_beta*exp(-1.0*x_fac*dist/simu.lrsi_r_s));
#endif
				I_mid = (I_now + I_prev) / 2.;
				dr    = dist - d_prev;

				/* accumulate potential */
				Phi += I_mid * dr;

				/* get escape velocity */
				v_esc2 = (2 * fabs(Phi - Phi0) * phi_fac);
			}
			/* potential="inf" for dist=0 -> set v_esc manually */
			else {
				v_esc2 = 1e30;
			}

			/* get particle velocity in halo rest frame velocity plus Hubble flow */
			VXp = (double)cur_part->mom[X];
			VYp = (double)cur_part->mom[Y];
			VZp = (double)cur_part->mom[Z];

#ifdef AHFnoHubbleDrag
			dVX = (VXp - VXc / M_vel) * v_fac;
			dVY = (VYp - VYc / M_vel) * v_fac;
			dVZ = (VZp - VZc / M_vel) * v_fac;
#else
			dVX = (VXp - VXc / M_vel) * v_fac + Hubble * (dX) * r_fac;
			dVY = (VYp - VYc / M_vel) * v_fac + Hubble * (dY) * r_fac;
			dVZ = (VZp - VZc / M_vel) * v_fac + Hubble * (dZ) * r_fac;
#endif

			/* absolute velocity of current particle */
			vel2 = pow2(dVX) + pow2(dVY) + pow2(dVZ);

#ifndef AHFignore_ugas
#ifdef GAS_PARTICLES
			/* u = 3/2  kT/m  = 1/2  v_therm^2 ?! */
			//vel2 += (cur_part->u < PGAS ? 0.0 : 2 * cur_part->u / pow2(v_fac));
			vel2 += (cur_part->u < PGAS ? 0.0 : 2 * cur_part->u ); // we actually convert all velocities to km/sec when comparing to v_esc2
#endif   /* GAS_PARTICLES */
#endif /* AHFignore_ugas */

			/* unbound particle? */
			if (vel2 > v2_tune * v_esc2) {
				/* count number of particles to be removed */
				nremove++;
			} else {
				/* let particle contribute to bulk/mean velocity of halo */
				M_vel += weight;
				VXc   += weight * cur_part->mom[X];
				VYc   += weight * cur_part->mom[Y];
				VZc   += weight * cur_part->mom[Z];

				/* store bound particles temporarily in bound_ipart[] */
				bound_npart++;
				bound_ipart = (long unsigned *)realloc(bound_ipart, bound_npart * sizeof(long unsigned));
				bound_ipart[bound_npart - 1] = cur_part - global.fst_part;

				/* accumulate virial values (NOTE: for v_esc2 we are accumulating M_r and not M_vir!) */
				M_vir += weight;
				R_vir  = dist;
			}

			I_prev = I_now;
			d_prev = dist;
		} /* particle loop */

		/* double-check new number of bound particles */
		if (bound_npart != (halo->npart - nremove))
			fprintf(stderr, "rem_unbound: better check the unbinding procedure! bound_part=%ld vs. halos[num].npart-nremove=%ld\n", bound_npart, halo->npart - nremove);


		/* update number of particles in halo */
		free(halo->ipart);
		halo->ipart = bound_ipart;
		halo->npart = bound_npart;
		halo->M_vir = M_vir;
		halo->R_vir = R_vir;


		/* there is no point in removing unbound particles once the halo is too small */
		if (((double)(halo->npart)) < ((double)(simu.AHF_MINPART)))
			break;
	} /* while( (nremove > n)  ) */

#ifdef VERBOSE2
	fprintf(io.logfile, "%12ld (Rvir=%g kpc/h)\n", halo->npart, halo->R_vir * x_fac * 1000.);
	fflush(io.logfile);
#endif
} /* rem_unbound */


/* qsort stuff */
int
haloCompare(const void *p1, const void *p2)
{
	if (((HALO *)p1)->npart < ((HALO *)p2)->npart)
		return 1;
	else if (((HALO *)p1)->npart > ((HALO *)p2)->npart)
		return -1;
	else
		return 0;
}

/*
 ************************************************************
 ************************************************************
 *  qsort stuff
 */
int
refCompare(const void *p1, const void *p2)
{
	if (((SPATIALREF *)p1)->volume < ((SPATIALREF *)p2)->volume)
		return 1;
	else if (((SPATIALREF *)p1)->volume > ((SPATIALREF *)p2)->volume)
		return -1;
	else
		return 0;
}

void
sortd(length, a, ind)
int length;


double a[];
int    ind[];
{
	int    i, ii, ij, j, m, m1, n2;
	double t;

	for (i = 0; i < length; i++)
		ind[i] = i;
	m  = 1;
	n2 = length / 2;
	m  = 2 * m;
	while (m <= n2)
		m = 2 * m;
	m  = m - 1;
three:;
	m1 = m + 1;
	for (j = m1 - 1; j < length; j++) {
		t  = a[ind[j]];
		ij = ind[j];
		i  = j - m;
		ii = ind[i];
four:;
		if (t > a[ii]) {
			ind[i + m] = ii;
			i          = i - m;
			if (i >= 0) {
				ii = ind[i];
				goto four;
			}
		}
		ind[i + m] = ij;
	}
	m = m / 2;
	if (m > 0)
		goto three;
	return;
}

/*
 ************************************************************
 ************************************************************
 * Remove particles outside virial radius
 */
void
rem_outsideRvir(HALO *halo, int icall)
{
	partptr       cur_part, host_part, first_part;
	double        M_sphere, ovdens_sphere, lR_sphere, R_sphere, V_sphere;
	double        Xc, Yc, Zc, Xp, Yp, Zp, dX, dY, dZ, weight, dummy;
	long unsigned npart_sphere, jpart;
	int           hostHalo, nbins, ibin;
	double        dist_min, dist_max, ldist_min, ldist_max, ldr, *ovdens,	*r, R_edge, ovdens_min, ovdens_vir, cur_dist;

	if (halo->npart < simu.AHF_MINPART)
		return;

#ifdef VERBOSE2
	fprintf(io.logfile, "    rem_outsideRvir(%d):  npart=%12ld -> ", icall,halo->npart);
	fflush(io.logfile);
#endif

	/* halo position in AMIGA units */
	Xc = halo->pos.x;
	Yc = halo->pos.y;
	Zc = halo->pos.z;
  
  /* before unbinding possibly allow for a bit lower overdensity criterion and hence collecting more particles */
  if (icall == 0)
    ovdens_vir = global.ovlim;
  else
    ovdens_vir = global.ovlim;
  

#ifdef AHFprofilerise
	/*============================================================================================
	 * this part checks for an upturn in the density profile and chops the halo
	 * off at that point
	 *
	 * NOTE: there are two calls to rem_outsideRvir()!
	 *       - one before rem_unbound() and
	 *       - one after  rem_unbound()
	 *
	 *       we should only look for the upturn before rem_unbound() as only
	 *       then the current halo may or may not be embedded within the background
	 *       of the host further, it saves time to not do it...
	 *============================================================================================*/

	if (icall == 0) {
		/* how many bins from where to where for density profile? */
		binning_parameter(*halo, &nbins, &dist_min, &dist_max);

		fprintf(io.logfile, "(nbins=%d, dist_min=%g dist_max=%g kpc/h) ", nbins, dist_min * x_fac * 1000., dist_max * x_fac * 1000.);
		fflush(io.logfile);

		/* logarithmical binning from dist_min to dist_max */
		ldist_min = log10(dist_min);
		ldist_max = log10(dist_max);
		ldr       = (ldist_max - ldist_min) / (double)nbins;

		/* create profile arrays */
		ovdens = (double *)calloc(nbins, sizeof(double));
		r      = (double *)calloc(nbins, sizeof(double));

		/* first particle */
		jpart    = 0;
		cur_dist = -1.0;

		/* calculate binned density profile */
		for (ibin = 0; ibin < nbins; ibin++) {
			/* get current outer radius using logarithmic radial bins */
			lR_sphere = ldist_min + ((double)ibin + 1) * ldr;
			R_sphere  = pow(10., lR_sphere);

			/* this heavily assumes that particles are ordered distancewise... */
			while (cur_dist < R_sphere && jpart < halo->npart) {
				/* access particle */
				cur_part = global.fst_part + halo->ipart[jpart];

#  ifdef MULTIMASS
				weight = (double)cur_part->weight;
#  else
				weight = (double)1.0;
#  endif
				/* accumulate mass */
				M_sphere += weight;

				/* particle position */
				Xp = (double)cur_part->pos[X];
				Yp = (double)cur_part->pos[Y];
				Zp = (double)cur_part->pos[Z];

				/* put particle into halo rest frame */
				dX = (Xp - Xc);
				dY = (Yp - Yc);
				dZ = (Zp - Zc);

				/* take care of periodic boundary conditions */
				if (dX > 0.5)					dX -= 1.0;
				if (dY > 0.5)					dY -= 1.0;
				if (dZ > 0.5)					dZ -= 1.0;
				if (dX < -0.5)  			dX += 1.0;
				if (dY < -0.5)				dY += 1.0;
				if (dZ < -0.5)				dZ += 1.0;

				/* distance of current particle (shouldn't sqrt(), but
				 *whatever...) */
				cur_dist = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));

				/* jump to next particle */
				jpart++;
			} /* while(cur_dist < R_sphere && jpart < halos[i].npart) */

			/* volume of sphere [0, cur_rad] */
			V_sphere     = 4. * PI / 3. * pow3(R_sphere);

			r[ibin]      = R_sphere;
			ovdens[ibin] = M_sphere / V_sphere * rho_fac / global.rho_vir;
		} /* ibin */

		/* we have ovdens(r) available now... */
		R_edge = get_haloedge(r, ovdens, nbins, 3);
	} /* icall==0 */
	else {
		R_edge = 2 * halo->R_vir;
	}
#endif   /* AHFprofilerise */

	/* (re-)set values... */
	jpart         = 0;
	npart_sphere  = 0;
	M_sphere      = 0.0;
	R_sphere      = -1.0;
	ovdens_sphere = 2 * ovdens_vir;


	while (jpart < halo->npart && ovdens_sphere >= ovdens_vir
#ifdef AHFprofilerise
	       && R_sphere <= R_edge
#endif
	       ) {
		/* access particle */
		cur_part = global.fst_part + halo->ipart[jpart];

#ifdef MULTIMASS
		weight = (double)cur_part->weight;
#else
		weight = (double)1.0;
#endif
		/* mass in current sphere */
		M_sphere += weight;

		/* particle position */
		Xp = (double)cur_part->pos[X];
		Yp = (double)cur_part->pos[Y];
		Zp = (double)cur_part->pos[Z];

		/* put particle into halo rest frame */
		dX = fabs(Xp - Xc);
		dY = fabs(Yp - Yc);
		dZ = fabs(Zp - Zc);

		/* take care of periodic boundary conditions */
		if (dX > 0.5)	dX -= 1.0;
		if (dY > 0.5)	dY -= 1.0;
		if (dZ > 0.5)	dZ -= 1.0;

		/* radius of current sphere */
		R_sphere = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));

		/* volume of current sphere */
		V_sphere = 4. * PI / 3. * pow3(R_sphere);

		/* overdensity in current sphere */
		ovdens_sphere = M_sphere / V_sphere * rho_fac / global.rho_vir;

		/* move to next particle */
		npart_sphere++;
		jpart++;
	}

	/* update halo properties (a realloc() nicely chops the ipart[] array...) */
  halo->npart = npart_sphere;
    
	halo->ipart = (long unsigned *)realloc(halo->ipart, halo->npart * sizeof(long unsigned));
	halo->M_vir  = M_sphere;
	halo->R_vir  = R_sphere;
	halo->ovdens = ovdens_sphere;
  if(icall==0)
    halo->Phi0 = 0.0;

	/*****************************************************************************************************************
	 *****************************************************************************************************************
	 ******************************************************************************************************************/


#ifdef VERBOSE2
	fprintf(io.logfile, "%12ld (Rvir=%g kpc/h)\n", halo->npart, halo->R_vir * x_fac * 1000.);
	fflush(io.logfile);
#endif

#ifdef AHFprofilerise
	if (icall == 0) {
		free(ovdens);
		free(r);
	}
#endif
} /* rem_outsideRvir */

/*
 ************************************************************
 ************************************************************
 * a little wrapper for the resetting the properties of DM/gas/star block properties
 */
void
reset_SPECIESPROP(SPECIESPROP *type)
{
	type->npart     = 0;
	type->Mass      = 0.0;
	type->pos_com.x = 0.0;
	type->pos_com.y = 0.0;
	type->pos_com.z = 0.0;
	type->pos_mbp.x = 0.0;
	type->pos_mbp.y = 0.0;
	type->pos_mbp.z = 0.0;
	type->vel.x     = 0.0;
	type->vel.y     = 0.0;
	type->vel.z     = 0.0;
	type->lambda    = 0.0;
	type->lambdaE   = 0.0;
	type->AngMom.x  = 0.0;
	type->AngMom.y  = 0.0;
	type->AngMom.z  = 0.0;
	type->axis.x    = 0.0;
	type->axis.y    = 0.0;
	type->axis.z    = 0.0;
	type->E1.x      = 0.0;
	type->E1.y      = 0.0;
	type->E1.z      = 0.0;
	type->E2.x      = 0.0;
	type->E2.y      = 0.0;
	type->E2.z      = 0.0;
	type->E3.x      = 0.0;
	type->E3.y      = 0.0;
	type->E3.z      = 0.0;
	type->Ekin      = 0.0;
	type->Epot      = 0.0;
} /* reset_SPECIESPROP */

/*
 ************************************************************
 ************************************************************
 * a little wrapper for the cumbersome calculation of the energy based spin parameter
 */
double
calc_lambdaE(double absAngMom, double Mhalo, double Mass, double Ekin, double Epot)
{
	double tmp1, tmp2, tmp3;

	tmp1 = sqrt(m_fac * Mhalo);                 /* tmp1 = sqrt(Mhalo)         */
	tmp1 = pow3(tmp1);                          /* tmp1 = Mhalo^3/2           */
	tmp2 = Ekin * m_fac * pow2(v_fac);          /* tmp2 = Ekin                */
	tmp3 = Epot * m_fac * phi_fac;              /* tmp3 = Epot                */
	tmp2 = fabs(tmp2 + tmp3);                   /* tmp2 = E                   */
	tmp2 = sqrt(tmp2);                          /* tmp2 = sqrt(|E|)           */
	tmp1 = tmp2 / tmp1;                         /* tmp1 = sqrt(|E|)/Mhalo^3/2 */
	tmp2 = m_fac * r_fac * v_fac * absAngMom;   /* tmp2 = L                   */
	tmp2 = tmp2 / (m_fac * Mass);               /* tmp2 = L/Mass              */
	tmp3 = tmp1 * tmp2;

	return tmp3 / Grav;
}

/*
 ************************************************************
 ************************************************************
 * Calculate profiles of halos
 */
int
HaloProfiles(HALO *halo)
{
	long unsigned ihalo, npart, n_prev, jpart, nprofile;
	double        Volume, V_prev, dV, cur_rad, lcur_rad, rad_prev, cur_dist,
	              rad_mid;
	double        M_sphere, M_prev, dM, weight, M_vir, r_vir, Mb_prev;
	partptr       cur_part, mb_part;
	gasptr        cur_gas;
	double        Xc, Yc, Zc, VXc, VYc, VZc;
	double        Xp, Yp, Zp, VXp, VYp, VZp;
	double        dX, dY, dZ, dVX, dVY, dVZ;
	double        R_edge, Lx, Ly, Lz, absAngMom, sig_v;
	double        itensor[3][3];
	double        a11, a22, a33, a12, a13, a23, axis1, axis2, axis3;
	double        pre_dist, I_mid, I_now, I_prev, Phi0, Phi, Epot, Ekin,
	              Tpart, Upart, Epart, Emin, v_esc2;
	double        dist_min, dist_max, ldist_min, ldist_max, dr, ldr;
	int           binpart, nbins, ibin, nignore;
	double        dummy;
	double        *Vcirc2, *dens_r2, *ovdens, *radius, M_sphere_prev, prev_dist;
	double        r2, R_max, V_max, M_max, x_max, y_max;
	double        CoM[NDIM];
  double        M_hires, M_lores;
  double        Ts, frad;
  double        FourPiThird = 4.*PI/3.;

#ifdef GAS_PARTICLES
	double        v_therm;

	long unsigned ndark;
	double        Xcom_dark, Ycom_dark, Zcom_dark, VX_dark, VY_dark,
	              VZ_dark;
	double        Xmbp_dark, Ymbp_dark, Zmbp_dark;
	double        M_dark;
	double        Lx_dark, Ly_dark, Lz_dark;
	double        Ekin_dark, Epot_dark;
	double        a11_dark, a22_dark, a33_dark, a12_dark, a13_dark,
	              a23_dark;
	partptr       mb_dark;
	double        Emin_dark;

	long unsigned ngas, ngas_prev;
	double        Xcom_gas, Ycom_gas, Zcom_gas, VX_gas, VY_gas, VZ_gas;
	double        Xmbp_gas, Ymbp_gas, Zmbp_gas;
	double        M_gas;
	double        Lx_gas, Ly_gas, Lz_gas;
	double        Ekin_gas, Epot_gas;
	double        a11_gas, a22_gas, a33_gas, a12_gas, a13_gas, a23_gas;
	partptr       mb_gas;
	double        Emin_gas;

	long unsigned nstar, nstar_prev;
	double        Xcom_star, Ycom_star, Zcom_star, VX_star, VY_star,
	              VZ_star;
	double        Xmbp_star, Ymbp_star, Zmbp_star;
	double        M_star;
	double        Lx_star, Ly_star, Lz_star;
	double        Ekin_star, Epot_star;
	double        a11_star, a22_star, a33_star, a12_star, a13_star,
	              a23_star;
	partptr       mb_star;
	double        Emin_star;
	double        u_shell_gas;
#endif

#ifdef METALHACK
	double z_star;
	double z_gas;
	double z_star_cum;
	double z_gas_cum;
	double M_shell_gas;
	double M_shell_star;
#endif
  

//  double VVVx=0., VVVy=0., VVVz=0.;
//  jpart = 0;
//  while (jpart < halo->npart) {
//    cur_part = global.fst_part + halo->ipart[jpart];
//    VVVx += cur_part->mom[X];
//    VVVy += cur_part->mom[Y];
//    VVVz += cur_part->mom[Z];
//    jpart++;
//  }
//  VVVx /= (double)halo->npart;
//  VVVy /= (double)halo->npart;
//  VVVz /= (double)halo->npart;
  
  
	/* only consider decent halos
	 * NOTE: npart may have changed since the last if-statement in ConstructHalos() */
	if (halo->npart < simu.AHF_MINPART)
		return TRUE;

#ifdef VERBOSE2
	fprintf(io.logfile, "    HaloProfiles:     npart=%12ld\n", halo->npart);
	fflush(io.logfile);
#endif

	/* how many bins from where to where for density profile? */
	binning_parameter(*halo, &nbins, &dist_min, &dist_max); // halos[i]

	/* logarithmical binning from dist_min to dist_max */
	ldist_min = log10(dist_min);
	ldist_max = log10(dist_max);
	ldr       = (ldist_max - ldist_min) / (double)nbins;

	/* create profile arrays */
	c_profile(halo, nbins);     // &halos[i]
  
#ifdef AHFparticle_Rmax_r2
  nprofile = halo->npart;
  nignore  = AHF_Rmax_r2_NIGNORE;
#else
  nprofile = nbins;
  nignore  = 0;
#endif
	Vcirc2   = (double *)calloc(nprofile, sizeof(double));
	dens_r2  = (double *)calloc(nprofile, sizeof(double));
	ovdens   = (double *)calloc(nprofile, sizeof(double));
  radius   = (double *)calloc(nprofile, sizeof(double));

	/* reset values */
	Phi0     = halo->Phi0;
	Phi      = 0.0;
	pre_dist = 0.0;
	I_prev   = 0.0;
  I_now    = 0.0;

	npart    = 0;
	n_prev   = 0;
	M_prev   = 0;
  Mb_prev  = 0;
	V_prev   = 0.0;
	rad_prev = dist_min;
  
  M_sphere_prev = 0.0;
  prev_dist     = 0.0;
  
	/* reset (cumulative) values based upon particles [0,cur_rad] */
	M_sphere = 0.0;
	VXc      = 0.0;
	VYc      = 0.0;
	VZc      = 0.0;
	a11      = 0.0;
	a22      = 0.0;
	a33      = 0.0;
	a12      = 0.0;
	a13      = 0.0;
	a23      = 0.0;
	sig_v    = 0.0;
	Lx       = 0.0;
	Ly       = 0.0;
	Lz       = 0.0;
	CoM[X]   = 0.0;
	CoM[Y]   = 0.0;
	CoM[Z]   = 0.0;
	Epot     = 0.0;
	Ekin     = 0.0;
	Emin     = 1e30;
	mb_part  = NULL;
  M_hires  = 0.0;
  M_lores  = 0.0;

#ifdef GAS_PARTICLES
	ndark     = 0;
	Xcom_dark = 0.0;
	Ycom_dark = 0.0;
	Zcom_dark = 0.0;
	VX_dark   = 0.0;
	VY_dark   = 0.0;
	VZ_dark   = 0.0;
	Xmbp_dark = 0.0;
	Ymbp_dark = 0.0;
	Zmbp_dark = 0.0;
	M_dark    = 0.0;
	Lx_dark   = 0.0;
	Ly_dark   = 0.0;
	Lz_dark   = 0.0;
	Ekin_dark = 0.0;
	Epot_dark = 0.0;
	a11_dark  = 0.0;
	a22_dark  = 0.0;
	a33_dark  = 0.0;
	a12_dark  = 0.0;
	a13_dark  = 0.0;
	a23_dark  = 0.0;
	mb_dark   = NULL;
	Emin_dark = 1e30;

	ngas      = 0;
  ngas_prev = 0;
	Xcom_gas  = 0.0;
	Ycom_gas  = 0.0;
	Zcom_gas  = 0.0;
	VX_gas    = 0.0;
	VY_gas    = 0.0;
	VZ_gas    = 0.0;
	Xmbp_gas  = 0.0;
	Ymbp_gas  = 0.0;
	Zmbp_gas  = 0.0;
	M_gas     = 0.0;
	Lx_gas    = 0.0;
	Ly_gas    = 0.0;
	Lz_gas    = 0.0;
	Ekin_gas  = 0.0;
	Epot_gas  = 0.0;
	a11_gas   = 0.0;
	a22_gas   = 0.0;
	a33_gas   = 0.0;
	a12_gas   = 0.0;
	a13_gas   = 0.0;
	a23_gas   = 0.0;
	mb_gas    = NULL;
	Emin_gas  = 1e30;

	nstar     = 0;
  nstar_prev = 0;
	Xcom_star = 0.0;
	Ycom_star = 0.0;
	Zcom_star = 0.0;
	VX_star   = 0.0;
	VY_star   = 0.0;
	VZ_star   = 0.0;
	Xmbp_star = 0.0;
	Ymbp_star = 0.0;
	Zmbp_star = 0.0;
	M_star    = 0.0;
	Lx_star   = 0.0;
	Ly_star   = 0.0;
	Lz_star   = 0.0;
	Ekin_star = 0.0;
	Epot_star = 0.0;
	a11_star  = 0.0;
	a22_star  = 0.0;
	a33_star  = 0.0;
	a12_star  = 0.0;
	a13_star  = 0.0;
	a23_star  = 0.0;
	mb_star   = NULL;
	Emin_star = 1e30;
	u_shell_gas = 0.0;
#endif

#ifdef METALHACK
	z_gas      = 0.0;
	z_gas_cum  = 0.0;
	z_star     = 0.0;
	z_star_cum = 0.0;
#endif

	/* halo position in AMIGA units */
	Xc = halo->pos.x;
	Yc = halo->pos.y;
	Zc = halo->pos.z;

	/* first particle */
	jpart    = 0;
	cur_dist = -1.0;
  
	/* loop over all bins */
	for (ibin = 0; ibin < nbins; ibin++) {
		/* get current outer radius using logarithmic radial bins */
		lcur_rad = ldist_min + ((double)ibin + 1) * ldr;
		cur_rad  = pow(10., lcur_rad);
    if (ibin == nbins-1) {
			/*
			 * This will make sure that the last bin is large enough to
			 * count up to the last particle: the calculation above can
			 * miss the last particle which can then in turn mess up the
			 * calculation of the profiled quantities.  The worst effect
			 * that missing a particle could have is that the number
			 * count of particle types is wrong which might then mess up
			 * the writing of the output files.
			 * The ZERO is added because the comparison is 
			 * cur_dist < cur_rad  and we do want the last bin to
			 * include the last particle (whereas for all other bins, a
			 * particle sitting directly at a bin border would be
			 * counted in the more inward bin).
			 */
      cur_rad = dist_max + ZERO;
		}

		/* Reset the shellular values */
#ifdef AHFshellshape
    a11       = 0.0;
    a22       = 0.0;
    a33       = 0.0;
    a12       = 0.0;
    a13       = 0.0;
    a23       = 0.0;
#ifdef GAS_PARTICLES
    a11_dark  = 0.0;
    a22_dark  = 0.0;
    a33_dark  = 0.0;
    a12_dark  = 0.0;
    a13_dark  = 0.0;
    a23_dark  = 0.0;
    a11_gas   = 0.0;
    a22_gas   = 0.0;
    a33_gas   = 0.0;
    a12_gas   = 0.0;
    a13_gas   = 0.0;
    a23_gas   = 0.0;
    a11_star  = 0.0;
    a22_star  = 0.0;
    a33_star  = 0.0;
    a12_star  = 0.0;
    a13_star  = 0.0;
    a23_star  = 0.0;
#endif /* GAS_PARTICLES */
#endif /* AHFshellshape */
#ifdef GAS_PARTICLES
		u_shell_gas = 0.0;
#endif
#ifdef METALHACK
		/* Reset the shellular values */
		z_star      = z_gas = 0.0;
		M_shell_gas = M_shell_star = 0.0;
#endif

		/* this heavily assumes that particles are ordered distancewise... */
		while (cur_dist < cur_rad && jpart < halo->npart) {
			/* access particle */
			cur_part = global.fst_part + halo->ipart[jpart];

			/* increment number of particles counter */
			npart++;

			/*----------------------
			 * calculate properties
			 *----------------------*/
#ifdef MULTIMASS
			weight = (double)cur_part->weight;
      
      /* this expects the high resolution DM particles to have weight=AHF_HIRES_DM_WEIGHT */
      if(fabs(weight-AHF_HIRES_DM_WEIGHT)<ZERO)
        M_hires += weight;
      else if (weight > AHF_HIRES_DM_WEIGHT)
        M_lores += weight;
#else
      weight = (double)1.0;
      M_hires += weight;
#endif

			/* accumulate mass */
			M_sphere += weight;

			/* particle position */
			Xp = (double)cur_part->pos[X];
			Yp = (double)cur_part->pos[Y];
			Zp = (double)cur_part->pos[Z];

			/* put particle into halo rest frame */
			dX = (Xp - Xc);
			dY = (Yp - Yc);
			dZ = (Zp - Zc);

			/* take care of periodic boundary conditions */
			if (dX > 0.5)				dX -= 1.0;
			if (dY > 0.5)				dY -= 1.0;
			if (dZ > 0.5)				dZ -= 1.0;
			if (dX < -0.5)			dX += 1.0;
			if (dY < -0.5)			dY += 1.0;
			if (dZ < -0.5)			dZ += 1.0;

			/* distance of current particle (shouldn't sqrt(), but whatever...) */
			cur_dist = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));

			/* centre of mass with correct boundaries... */
			CoM[X] += weight * (Xc + dX);
			CoM[Y] += weight * (Yc + dY);
			CoM[Z] += weight * (Zc + dZ);

#ifdef AHFreducedinertiatensor
			/* reduced inertia tensor of all particles within sphere(!) */
      if (cur_dist > MACHINE_ZERO) {
        a11 += weight * dX * dX / pow2(cur_dist);
        a22 += weight * dY * dY / pow2(cur_dist);
        a33 += weight * dZ * dZ / pow2(cur_dist);
        a12 += weight * dX * dY / pow2(cur_dist);
        a13 += weight * dX * dZ / pow2(cur_dist);
        a23 += weight * dY * dZ / pow2(cur_dist);
      }
#else
			/* inertia tensor of all particles within sphere(!) */
			a11 += weight * dX * dX;
			a22 += weight * dY * dY;
			a33 += weight * dZ * dZ;
			a12 += weight * dX * dY;
			a13 += weight * dX * dZ;
			a23 += weight * dY * dZ;
#endif
			/* mean halo velocity of particles in sphere [0,cur_rad] */
			VXc += weight * cur_part->mom[X];
			VYc += weight * cur_part->mom[Y];
			VZc += weight * cur_part->mom[Z];

			/* particle velocity in halo rest frame (Hubble correction not needed for L as r x r = 0) */
			VXp = (double)cur_part->mom[X];
			VYp = (double)cur_part->mom[Y];
			VZp = (double)cur_part->mom[Z];
			dVX = (VXp - VXc / M_sphere);
			dVY = (VYp - VYc / M_sphere);
			dVZ = (VZp - VZc / M_sphere);
      

			/* angular momentum of particles within sphere [0, cur_rad] */
			Lx += weight * (dY * dVZ - dZ * dVY);
			Ly += weight * (dZ * dVX - dX * dVZ);
			Lz += weight * (dX * dVY - dY * dVX);

			/* add Hubble flow to velocities */
			dVX += Hubble * dX * r_fac / v_fac;
			dVY += Hubble * dY * r_fac / v_fac;
			dVZ += Hubble * dZ * r_fac / v_fac;

      /* kinetic energy of particle */
			Tpart = weight * (pow2(dVX) + pow2(dVY) + pow2(dVZ));

			/* get potential energy of particle */
			if (cur_dist > MACHINE_ZERO) {
				I_now = M_sphere / pow2(cur_dist);
#ifdef AHF_LRSI
				I_now*=(1.0+ simu.lrsi_beta*exp(-1.0*x_fac*cur_dist/simu.lrsi_r_s));
#endif
				I_mid = (I_now + I_prev) / 2.;
				Phi  += I_mid * (cur_dist - pre_dist);

				Upart = (Phi - Phi0) * weight;
			} else {
				/* simply use old Phi value... */
				Upart = (Phi - Phi0) * weight;
			}
      
#ifdef SUSSING2013
      cur_part->E = 0.5*(Tpart*pow2(v_fac)+Upart*phi_fac)*m_fac;
      //cur_part->d = cur_dist*x_fac*1000.;
#endif
      
#ifdef AHFDEBUG
      if(cur_dist<pre_dist)
        fprintf(stderr,"HaloProfiles: cur_dist<pre_dist   %g<%g (%ld of %ld)\n",cur_dist*x_fac*1000.,pre_dist*x_fac*1000.,npart,halo->npart);
      if(Upart>0.)
        fprintf(stderr,"HaloProfiles: positive potential  %g (%ld of %ld)\n",Upart,npart,halo->npart);
#endif
        
			/* keep track of escape velocity */
			v_esc2 = 2 * fabs(Upart) / weight;      
      
			/* take care of potential energy integration... */
			I_prev   = I_now;
			pre_dist = cur_dist;

			/* potential energy of all particles within sphere [0, cur_rad] */
			Epot += Upart;


#ifdef GAS_PARTICLES
			/*==============
			 * GAS PARTICLE
			 *==============*/
			if (isgreaterequal(cur_part->u, PGAS))
			{
#ifndef AHFignore_ugas
				/* account for thermal energy u = 3/2  kT/m  = 1/2  v_therm^2 ?! */
				v_therm = 2 * cur_part->u / pow2(v_fac);
        //we use v_fac here as we transform 'u' to a velocity which internally requires the additional 'aexp' factor
				Tpart  += weight * v_therm;
#endif

				/* calculate integral properties of gas alone */
				ngas++;
        M_gas    += weight;
				Xcom_gas += weight * (Xc + dX);
				Ycom_gas += weight * (Yc + dY); // this ensures capturing periodic boundary conditions
				Zcom_gas += weight * (Zc + dZ);
				VX_gas   += weight * cur_part->mom[X];
				VY_gas   += weight * cur_part->mom[Y];
				VZ_gas   += weight * cur_part->mom[Z];
				Lx_gas   += weight * (dY * dVZ - dZ * dVY);
				Ly_gas   += weight * (dZ * dVX - dX * dVZ);
				Lz_gas   += weight * (dX * dVY - dY * dVX);
#  ifdef AHFreducedinertiatensor
       if (cur_dist > MACHINE_ZERO){
         a11_gas  += weight * dX * dX / pow2(cur_dist);
         a22_gas  += weight * dY * dY / pow2(cur_dist);
         a33_gas  += weight * dZ * dZ / pow2(cur_dist);
         a12_gas  += weight * dX * dY / pow2(cur_dist);
         a13_gas  += weight * dX * dZ / pow2(cur_dist);
         a23_gas  += weight * dY * dZ / pow2(cur_dist);
       }
#  else
				a11_gas  += weight * dX * dX;
				a22_gas  += weight * dY * dY;
				a33_gas  += weight * dZ * dZ;
				a12_gas  += weight * dX * dY;
				a13_gas  += weight * dX * dZ;
				a23_gas  += weight * dY * dZ;
#  endif
				/* potential energy of gas alone */
				Epot_gas += Upart;

				/* calculate kinetic+internal gas energy */
				Ekin_gas += Tpart;

				/* keep track of most bound gas particle */
				Epart = (0.5 * Tpart + Upart);
				if (Epart < Emin_gas) {
					Emin_gas = Epart;
					mb_gas   = cur_part;
				}
        u_shell_gas += weight * cur_part->u / u_fac; // for whatever reasons we are converting cur_part->u to code internal units at this point!
#  ifdef METALHACK
				z_gas       += weight * (cur_part->z);
				z_gas_cum   += weight * (cur_part->z);
				M_shell_gas += weight;
#  endif
			}

			/*===============
			 * STAR PARTICLE
			 *===============*/
			if (fabs(cur_part->u - PSTAR) < ZERO)
       {
				/* calculate integral properties of stars alone */
				nstar++;
				M_star    += weight;
				Xcom_star += weight * (Xc + dX);
				Ycom_star += weight * (Yc + dY); // this ensures capturing periodic boundary conditions
				Zcom_star += weight * (Zc + dZ);
				VX_star   += weight * cur_part->mom[X];
				VY_star   += weight * cur_part->mom[Y];
				VZ_star   += weight * cur_part->mom[Z];
				Lx_star   += weight * (dY * dVZ - dZ * dVY);
				Ly_star   += weight * (dZ * dVX - dX * dVZ);
				Lz_star   += weight * (dX * dVY - dY * dVX);
#  ifdef AHFreducedinertiatensor
        if (cur_dist > MACHINE_ZERO) {
          a11_star  += weight * dX * dX / pow2(cur_dist);
          a22_star  += weight * dY * dY / pow2(cur_dist);
          a33_star  += weight * dZ * dZ / pow2(cur_dist);
          a12_star  += weight * dX * dY / pow2(cur_dist);
          a13_star  += weight * dX * dZ / pow2(cur_dist);
          a23_star  += weight * dY * dZ / pow2(cur_dist);
        }
#  else
				a11_star  += weight * dX * dX;
				a22_star  += weight * dY * dY;
				a33_star  += weight * dZ * dZ;
				a12_star  += weight * dX * dY;
				a13_star  += weight * dX * dZ;
				a23_star  += weight * dY * dZ;
#  endif
				/* potential energy of stars alone */
				Epot_star += Upart;

				/* calculate kinetic+internal star energy */
				Ekin_star += Tpart;

				/* keep track of most bound star particle */
				Epart = (0.5 * Tpart + Upart);
				if (Epart < Emin_star) {
					Emin_star = Epart;
					mb_star   = cur_part;
				}
#  ifdef METALHACK
				z_star       += weight * (cur_part->z);
				z_star_cum   += weight * (cur_part->z);
				M_shell_star += weight;
#  endif
       }

			/*======================
			 * DARK MATTER PARTICLE
			 *======================*/
			if (fabs(cur_part->u - PDM) < ZERO)
       {
				/* calculate integral properties of dark matter alone */
				ndark++;
				M_dark    += weight;
				Xcom_dark += weight * (Xc + dX);
				Ycom_dark += weight * (Yc + dY); // this ensures capturing periodic boundary conditions
				Zcom_dark += weight * (Zc + dZ);
				VX_dark   += weight * cur_part->mom[X];
				VY_dark   += weight * cur_part->mom[Y];
				VZ_dark   += weight * cur_part->mom[Z];
				Lx_dark   += weight * (dY * dVZ - dZ * dVY);
				Ly_dark   += weight * (dZ * dVX - dX * dVZ);
				Lz_dark   += weight * (dX * dVY - dY * dVX);
#  ifdef AHFreducedinertiatensor
        if (cur_dist > MACHINE_ZERO) {
          a11_dark  += weight * dX * dX / pow2(cur_dist);
          a22_dark  += weight * dY * dY / pow2(cur_dist);
          a33_dark  += weight * dZ * dZ / pow2(cur_dist);
          a12_dark  += weight * dX * dY / pow2(cur_dist);
          a13_dark  += weight * dX * dZ / pow2(cur_dist);
          a23_dark  += weight * dY * dZ / pow2(cur_dist);
        }
#  else
				a11_dark  += weight * dX * dX;
				a22_dark  += weight * dY * dY;
				a33_dark  += weight * dZ * dZ;
				a12_dark  += weight * dX * dY;
				a13_dark  += weight * dX * dZ;
				a23_dark  += weight * dY * dZ;
#  endif
				/* potential energy of DM alone */
				Epot_dark += Upart;

				/* calculate kinetic+internal DM energy */
				Ekin_dark += Tpart;

				/* keep track of most bound DM particle */
				Epart = (0.5 * Tpart + Upart);
				if (Epart < Emin_dark) {
					Emin_dark = Epart;
					mb_dark   = cur_part;
				}
			}
#endif   /* GAS_PARTICLES */


			/* velocity dispersion and kinetic energy of all particles within sphere [0, cur_rad] */
			sig_v += Tpart;
			Ekin  += Tpart;

			/* total energy of current particle */
			Epart = (0.5 * Tpart + Upart);

			/* remember most bound particle */
			if (Epart < Emin) {
				Emin    = Epart;
				mb_part = cur_part;
			}
      
#ifdef AHFparticle_Rmax_r2
      /* accumulate radius[], Vmax2[], ovdens[], and dens_r2[] arrays */
      radius[jpart]  = cur_dist;
      ovdens[jpart]  = M_sphere / (FourPiThird * pow3(cur_dist));
      dV             = FourPiThird*(pow3(cur_dist)-pow3(prev_dist));
#	if (defined AHFdmonly_Rmax_r2 && defined GAS_PARTICLES)
      dM             = M_sphere-M_gas-M_star - M_sphere_prev;
      Vcirc2[jpart]  = (M_sphere-M_gas-M_star) / cur_dist;
      M_sphere_prev  = M_sphere-M_gas-M_star;
# else
      dM             = M_sphere - M_sphere_prev;
      Vcirc2[jpart]  = M_sphere / cur_dist;
      M_sphere_prev  = M_sphere;
#endif
      dens_r2[jpart] = dM/dV * pow2((cur_dist+prev_dist)/2.);


      prev_dist      = cur_dist;
#endif /* AHFparticle_Rmax_r2 */
      

			/* jump to next particle */
			jpart++;
      
		} /* while(cur_dist < cur_rad && jpart < halo->npart) */

    
    /*
     * cur_rad is our set right edge of the bin, but maybe the last particle considered does not have cur_dist==cur_rad
     * -> then it makes more sense to use the actual distance of the last particle used for the calculation! 
     */
    //cur_rad = cur_dist; // but this leads to nan's in _profiles when there are no particles in an anticipated shell :-(
    
    
		/* volume of sphere [0, cur_rad] */
		Volume = FourPiThird * pow3(cur_rad);

		/* differential values */
		dM = M_sphere - M_prev; /* mass   in shell [prev_rad, cur_rad] */
		dV = Volume   - V_prev; /* volume of shell [prev_rad, cur_rad] */

		/* get eigenavalues of inertia tensor */
#ifdef AHFshellshape
    if(npart-n_prev > AHF_MINPART_SHELL)
#else
    if(npart > AHF_MINPART_SHELL)
#endif /* AHFshellshape */
     {
      itensor[0][0] = a11;
      itensor[1][1] = a22;
      itensor[2][2] = a33;
      itensor[0][1] = a12;
      itensor[1][0] = a12;
      itensor[0][2] = a13;
      itensor[2][0] = a13;
      itensor[1][2] = a23;
      itensor[2][1] = a23;
      get_axes(itensor, &axis1, &axis2, &axis3);
     }
    else
     {
      itensor[0][0] = 0;
      itensor[1][1] = 0;
      itensor[2][2] = 0;
      itensor[0][1] = 0;
      itensor[1][0] = 0;
      itensor[0][2] = 0;
      itensor[2][0] = 0;
      itensor[1][2] = 0;
      itensor[2][1] = 0;
      axis1         = 1;
      axis2         = 0;
      axis3         = 0;
     }

		rad_mid                  = (cur_rad + rad_prev) / 2.;
		halo->prof.r[ibin]       = cur_rad;
		halo->prof.npart[ibin]   = npart;
		halo->prof.nvpart[ibin]  = M_sphere;
		halo->prof.ovdens[ibin]  = M_sphere / Volume;
		halo->prof.dens[ibin]    = dM       / dV;
		halo->prof.v2_circ[ibin] = M_sphere / cur_rad;
#ifdef AHF_LRSI
		halo->prof.v2_circ[ibin]*= (1.0+(1.0 + x_fac*cur_rad/simu.lrsi_r_s)*simu.lrsi_beta*exp(-x_fac*cur_rad/simu.lrsi_r_s));
#endif
		halo->prof.v_esc2[ibin]  = v_esc2;
		halo->prof.sig_v[ibin]   = sqrt(sig_v / M_sphere);
		halo->prof.Ekin[ibin]    = 0.5 * Ekin;
		halo->prof.Epot[ibin]    = 0.5 * Epot;
		halo->prof.Lx[ibin]      = Lx;
		halo->prof.Ly[ibin]      = Ly;
		halo->prof.Lz[ibin]      = Lz;
		halo->prof.E1x[ibin]     = itensor[0][0];
		halo->prof.E1y[ibin]     = itensor[1][0];
		halo->prof.E1z[ibin]     = itensor[2][0];
		halo->prof.E2x[ibin]     = itensor[0][1];
		halo->prof.E2y[ibin]     = itensor[1][1];
		halo->prof.E2z[ibin]     = itensor[2][1];
		halo->prof.E3x[ibin]     = itensor[0][2];
		halo->prof.E3y[ibin]     = itensor[1][2];
		halo->prof.E3z[ibin]     = itensor[2][2];
		halo->prof.axis1[ibin]   = 1.0;
    if(axis1 > 0.) {
      halo->prof.axis2[ibin]   = sqrt(axis2 / axis1);  // cf. PhD thesis of Kristin Riebe Eqs. A.7-A.9
      halo->prof.axis3[ibin]   = sqrt(axis3 / axis1);
    }
    else {
      halo->prof.axis2[ibin]   = 0.0;
      halo->prof.axis3[ibin]   = 0.0;
    }
    
#ifdef GAS_PARTICLES
		halo->prof.M_gas[ibin]   = M_gas;
		halo->prof.M_star[ibin]  = M_star;
		halo->prof.u_gas[ibin]   = u_shell_gas;
#ifdef AHFdisks
		halo->prof.Ekin_gas[ibin]    = 0.5 * Ekin_gas;
		halo->prof.Lx_gas[ibin]      = Lx_gas;
		halo->prof.Ly_gas[ibin]      = Ly_gas;
		halo->prof.Lz_gas[ibin]      = Lz_gas;
		halo->prof.Ekin_star[ibin]   = 0.5 * Ekin_star;
		halo->prof.Lx_star[ibin]     = Lx_star;
		halo->prof.Ly_star[ibin]     = Ly_star;
		halo->prof.Lz_star[ibin]     = Lz_star;
    
    /* gas shape */
#ifdef AHFshellshape
    if(ngas-ngas_prev > AHF_MINPART_SHELL)
#else
    if(ngas > AHF_MINPART_GAS)
#endif /* AHFshellshape */
     {
      itensor[0][0] = a11_gas;
      itensor[1][1] = a22_gas;
      itensor[2][2] = a33_gas;
      itensor[0][1] = a12_gas;
      itensor[1][0] = a12_gas;
      itensor[0][2] = a13_gas;
      itensor[2][0] = a13_gas;
      itensor[1][2] = a23_gas;
      itensor[2][1] = a23_gas;
      get_axes(itensor, &axis1, &axis2, &axis3);
     }
    else
     {
      itensor[0][0] = 0;
      itensor[1][1] = 0;
      itensor[2][2] = 0;
      itensor[0][1] = 0;
      itensor[1][0] = 0;
      itensor[0][2] = 0;
      itensor[2][0] = 0;
      itensor[1][2] = 0;
      itensor[2][1] = 0;
      axis1         = 1;
      axis2         = 0;
      axis3         = 0;
     }
		halo->prof.E1x_gas[ibin]     = itensor[0][0];
		halo->prof.E1y_gas[ibin]     = itensor[1][0];
		halo->prof.E1z_gas[ibin]     = itensor[2][0];
		halo->prof.E2x_gas[ibin]     = itensor[0][1];
		halo->prof.E2y_gas[ibin]     = itensor[1][1];
		halo->prof.E2z_gas[ibin]     = itensor[2][1];
		halo->prof.E3x_gas[ibin]     = itensor[0][2];
		halo->prof.E3y_gas[ibin]     = itensor[1][2];
		halo->prof.E3z_gas[ibin]     = itensor[2][2];
		halo->prof.axis1_gas[ibin]   = 1.0;
    if(axis1 > 0.) {
      halo->prof.axis2_gas[ibin]   = sqrt(axis2 / axis1);  // cf. PhD thesis of Kristin Riebe Eqs. A.7-A.9
      halo->prof.axis3_gas[ibin]   = sqrt(axis3 / axis1);
    }
    else {
      halo->prof.axis2_gas[ibin]   = 0.0;
      halo->prof.axis3_gas[ibin]   = 0.0;
    }
    
    /* stellar shape */
#ifdef AHFshellshape
    if(nstar-nstar_prev > AHF_MINPART_SHELL)
#else
    if(nstar > AHF_MINPART_STARS)
#endif /* AHFshellshape */
     {
      itensor[0][0] = a11_star;
      itensor[1][1] = a22_star;
      itensor[2][2] = a33_star;
      itensor[0][1] = a12_star;
      itensor[1][0] = a12_star;
      itensor[0][2] = a13_star;
      itensor[2][0] = a13_star;
      itensor[1][2] = a23_star;
      itensor[2][1] = a23_star;
      get_axes(itensor, &axis1, &axis2, &axis3);
     }
    else
     {
      itensor[0][0] = 0;
      itensor[1][1] = 0;
      itensor[2][2] = 0;
      itensor[0][1] = 0;
      itensor[1][0] = 0;
      itensor[0][2] = 0;
      itensor[2][0] = 0;
      itensor[1][2] = 0;
      itensor[2][1] = 0;
      axis1         = 1;
      axis2         = 0;
      axis3         = 0;
     }
		halo->prof.E1x_star[ibin]     = itensor[0][0];
		halo->prof.E1y_star[ibin]     = itensor[1][0];
		halo->prof.E1z_star[ibin]     = itensor[2][0];
		halo->prof.E2x_star[ibin]     = itensor[0][1];
		halo->prof.E2y_star[ibin]     = itensor[1][1];
		halo->prof.E2z_star[ibin]     = itensor[2][1];
		halo->prof.E3x_star[ibin]     = itensor[0][2];
		halo->prof.E3y_star[ibin]     = itensor[1][2];
		halo->prof.E3z_star[ibin]     = itensor[2][2];
		halo->prof.axis1_star[ibin]   = 1.0;
    if(axis1 > 0.) {
      halo->prof.axis2_star[ibin]   = sqrt(axis2 / axis1);  // cf. PhD thesis of Kristin Riebe Eqs. A.7-A.9
      halo->prof.axis3_star[ibin]   = sqrt(axis3 / axis1);
    }
    else {
      halo->prof.axis2_star[ibin]   = 0.0;
      halo->prof.axis3_star[ibin]   = 0.0;
    }
#endif /* AHFdisks */
#endif /* GAS_PARTICLES */

#ifdef METALHACK
		halo->prof.z_gas[ibin]   = (isgreater(M_shell_gas, 0.0)  ? z_gas  / M_shell_gas  : 0.0);
		halo->prof.z_star[ibin]  = (isgreater(M_shell_star, 0.0) ? z_star / M_shell_star : 0.0);
#endif


    
    
#ifndef AHFparticle_Rmax_r2
    /* accumulate radius[], ovdens[], Vcirc2[], and dens_r2[] arrays using binned profile */
    radius[ibin]  = halo->prof.r[ibin];    
		ovdens[ibin]  = halo->prof.ovdens[ibin];
#	if (defined AHFdmonly_Rmax_r2 && defined GAS_PARTICLES)
    Vcirc2[ibin]  = (M_sphere-M_gas-M_star)/cur_rad;
    dens_r2[ibin] = ((M_sphere-M_gas-M_star)-(M_prev-Mb_prev))/dV * pow2(rad_mid);
#else
		Vcirc2[ibin]  = halo->prof.v2_circ[ibin];
		dens_r2[ibin] = halo->prof.dens[ibin] * pow2(rad_mid);
#endif
#endif /* AHFparticle_Rmax_r2 */

    
    
    
		/* store values from present bin */
		M_prev   = M_sphere;
#ifdef GAS_PARTICLES
    Mb_prev  = M_gas+M_star;
#endif
		V_prev     = Volume;
		n_prev     = npart;
#ifdef GAS_PARTICLES
    ngas_prev  = ngas;
    nstar_prev = nstar;
#endif
		rad_prev   = cur_rad;
	} /* ibin */

	/* centre of mass */
	CoM[X] = fmod(CoM[X] / M_sphere + 1., 1.);
	CoM[Y] = fmod(CoM[Y] / M_sphere + 1., 1.);
	CoM[Z] = fmod(CoM[Z] / M_sphere + 1., 1.);

	/* get scale radius r2 */
#ifdef AHFsplinefit
	find_max_spline(radius+nignore, dens_r2+nignore, nprofile-nignore, 3, &x_max, &y_max, 10000);
#else
	find_max(radius+nignore, dens_r2+nignore, nprofile-nignore, 3, &x_max, &y_max);
#endif
	r2 = x_max;

	/* get peak position R_max of rotation curve */
#ifdef AHFsplinefit
	find_max_spline(radius+nignore, Vcirc2+nignore, nprofile-nignore, 1, &x_max, &y_max, 10000);
#else
	find_max(radius+nignore, Vcirc2+nignore, nprofile-nignore, 1, &x_max, &y_max);
#endif
	R_max = x_max;
  
	/* get peak height V_max of rotation curve using all matter */
  ibin=0;
  while(radius[ibin] < x_max && ibin < nprofile-1)
    ibin++;
  
#ifdef AHFparticle_Rmax_r2
  M_max = ovdens[ibin] * FourPiThird * pow3(radius[ibin]);
#else
  M_max = halo->prof.nvpart[ibin]/radius[ibin] * R_max; //TODO: use more sophisticated guess for mass in-between two bins
#endif
	V_max = M_max/R_max;

	/* get physical outer radius: either upturn or virial radius  (NOTE: ovdens[] will be modified/smoothed) */
#ifdef AHFfindRedge
	R_edge = get_haloedge(radius, ovdens, nprofile, 1);
#endif
  
  /* free all temporary profile arrays */
	free(ovdens);
	free(Vcirc2);
  free(radius);
	free(dens_r2);
  
  
	/* get angular momentum */
	Lx        = halo->prof.Lx[nbins - 1];
	Ly        = halo->prof.Ly[nbins - 1];
	Lz        = halo->prof.Lz[nbins - 1];
	absAngMom = sqrt(pow2(Lx) + pow2(Ly) + pow2(Lz));

	/* store final (integral) halo properties (we assume that profile stores cumulative values...) */
	halo->M_vir    = M_sphere;
	halo->vel.x    = VXc / M_sphere;
	halo->vel.y    = VYc / M_sphere;
	halo->vel.z    = VZc / M_sphere;
	halo->R_edge   = R_edge; // NOTE, we do neither calculate nor write this information anymore
	halo->R_max    = R_max;
	halo->V2_max   = V_max;
	halo->r2       = r2;
	halo->sigV     = halo->prof.sig_v[nbins - 1];
	halo->v_esc2   = halo->prof.v_esc2[nbins - 1];
	halo->Ekin     = halo->prof.Ekin[nbins - 1];
	halo->Epot     = halo->prof.Epot[nbins - 1];
	halo->axis.x   = halo->prof.axis1[nbins - 1];
	halo->axis.y   = halo->prof.axis2[nbins - 1];
	halo->axis.z   = halo->prof.axis3[nbins - 1];
	halo->E1.x     = halo->prof.E1x[nbins - 1];
	halo->E1.y     = halo->prof.E1y[nbins - 1];
	halo->E1.z     = halo->prof.E1z[nbins - 1];
	halo->E2.x     = halo->prof.E2x[nbins - 1];
	halo->E2.y     = halo->prof.E2y[nbins - 1];
	halo->E2.z     = halo->prof.E2z[nbins - 1];
	halo->E3.x     = halo->prof.E3x[nbins - 1];
	halo->E3.y     = halo->prof.E3y[nbins - 1];
	halo->E3.z     = halo->prof.E3z[nbins - 1];
	halo->AngMom.x = Lx / absAngMom;
	halo->AngMom.y = Ly / absAngMom;
	halo->AngMom.z = Lz / absAngMom;
  
  /* Bullock et al. (2001) spin parameter */
	halo->lambda   = absAngMom / halo->M_vir / sqrt(2. * halo->M_vir * halo->R_vir);
	halo->lambda  *= v_fac * sqrt(r_fac / (Grav * m_fac));
  
  /* energy based spin parameter (ala Peebles) */
	halo->lambdaE  = calc_lambdaE(absAngMom, halo->M_vir, halo->M_vir, halo->Ekin, halo->Epot);
  
  /* NFW concentration ala Prada et al. (2012) */
  halo->cNFW     = calc_cNFW(halo->V2_max, halo->M_vir/halo->R_vir);
  
  /* surface pressure term ala Shaw et al. (2006) */
  Ts    = 2.0*(halo->prof.Ekin[nbins-1]-halo->prof.Ekin[nbins-2]);
  frad  = fabs(halo->prof.r[nbins-2]/halo->prof.r[nbins-1]);
  halo->SurfP    = -0.125*pow3(1.+frad)/(1.-pow3(frad))*Ts;;
  
  /* mass fraction of high-resolution particles */
  if(M_hires>0)
    halo->fMhires  = M_hires/(M_hires+M_lores); 
  else
    halo->fMhires  = 0.0; 
  
  /* offset of most-bound particle from halo centre */
	if (mb_part != NULL) {
		dX = fabs(mb_part->pos[X] - Xc);
		dY = fabs(mb_part->pos[Y] - Yc);
		dZ = fabs(mb_part->pos[Z] - Zc);
		if (dX > 0.5)	dX -= 1.0;
		if (dY > 0.5)	dY -= 1.0;
		if (dZ > 0.5)	dZ -= 1.0;
		halo->mbp_offset = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
		halo->pos_mbp.x  = mb_part->pos[X];
		halo->pos_mbp.y  = mb_part->pos[Y];
		halo->pos_mbp.z  = mb_part->pos[Z];
		halo->vel_mbp.x  = mb_part->mom[X];
		halo->vel_mbp.y  = mb_part->mom[Y];
		halo->vel_mbp.z  = mb_part->mom[Z];
	} else {
		halo->mbp_offset = -1.0;
		halo->pos_mbp.x  = 0.0;
		halo->pos_mbp.y  = 0.0;
		halo->pos_mbp.z  = 0.0;
		halo->vel_mbp.x  = 0.0;
		halo->vel_mbp.y  = 0.0;
		halo->vel_mbp.z  = 0.0;
	}

  /* offset of particles' centre-of-mass from halo centre */
	dX = fabs(CoM[X] - Xc);
	dY = fabs(CoM[Y] - Yc);
	dZ = fabs(CoM[Z] - Zc);
	if (dX > 0.5)	dX -= 1.0;
	if (dY > 0.5)	dY -= 1.0;
	if (dZ > 0.5)	dZ -= 1.0;
	halo->com_offset = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
	halo->pos_com.x  = CoM[X];
	halo->pos_com.y  = CoM[Y];
	halo->pos_com.z  = CoM[Z];

#ifdef GAS_PARTICLES
	/*-----------------
	 * GAS properties
	 *-----------------*/
	if (ngas > 0) {
		/* mass of gas alone */
		halo->gas_only.npart = ngas;
		halo->gas_only.Mass  = M_gas;

		/* centre-of-mass properties */
		halo->gas_only.pos_com.x = fmod(Xcom_gas / M_gas + 1.0, 1.0);
		halo->gas_only.pos_com.y = fmod(Ycom_gas / M_gas + 1.0, 1.0);
		halo->gas_only.pos_com.z = fmod(Zcom_gas / M_gas + 1.0, 1.0);
		halo->gas_only.vel.x     = VX_gas / M_gas;
		halo->gas_only.vel.y     = VY_gas / M_gas;
		halo->gas_only.vel.z     = VZ_gas / M_gas;

		/* energy of gas alone */
		halo->gas_only.Ekin = 0.5 * Ekin_gas;
		halo->gas_only.Epot = 0.5 * Epot_gas;

    if(ngas > AHF_MINPART_GAS)
     {
      /* get spin parameter of gas particles (all gas part's within sphere!) */
      absAngMom = sqrt(pow2(Lx_gas) + pow2(Ly_gas) + pow2(Lz_gas));
      halo->gas_only.AngMom.x = Lx_gas / absAngMom;
      halo->gas_only.AngMom.y = Ly_gas / absAngMom;
      halo->gas_only.AngMom.z = Lz_gas / absAngMom;
      halo->gas_only.lambda   = absAngMom / M_gas / sqrt(2. * halo->M_vir * halo->R_vir);
      halo->gas_only.lambda  *= v_fac * sqrt(r_fac / (Grav * m_fac));
      halo->gas_only.lambdaE  = calc_lambdaE(absAngMom, halo->M_vir, halo->gas_only.Mass, halo->gas_only.Ekin, halo->gas_only.Epot);
      
      /* get shape of gas particles (all gas part's within sphere!) */
      itensor[0][0] = a11_gas;
      itensor[1][1] = a22_gas;
      itensor[2][2] = a33_gas;
      itensor[0][1] = a12_gas;
      itensor[1][0] = a12_gas;
      itensor[0][2] = a13_gas;
      itensor[2][0] = a13_gas;
      itensor[1][2] = a23_gas;
      itensor[2][1] = a23_gas;
      get_axes(itensor, &axis1, &axis2, &axis3);
      
      halo->gas_only.axis.x = 1.0;
      if(axis1 > 0.) {
        halo->gas_only.axis.y = sqrt(axis2 / axis1);
        halo->gas_only.axis.z = sqrt(axis3 / axis1);
      }
      else {
        halo->gas_only.axis.y = 0.0;
        halo->gas_only.axis.z = 0.0;
      }
      halo->gas_only.E1.x   = itensor[0][0];
      halo->gas_only.E1.y   = itensor[1][0];
      halo->gas_only.E1.z   = itensor[2][0];
      halo->gas_only.E2.x   = itensor[0][1];
      halo->gas_only.E2.y   = itensor[1][1];
      halo->gas_only.E2.z   = itensor[2][1];
      halo->gas_only.E3.x   = itensor[0][2];
      halo->gas_only.E3.y   = itensor[1][2];
      halo->gas_only.E3.z   = itensor[2][2];
     }
    
		/* position of most bound gas particle */
		if (mb_gas != NULL) {
			halo->gas_only.pos_mbp.x = mb_gas->pos[X];
			halo->gas_only.pos_mbp.y = mb_gas->pos[Y];
			halo->gas_only.pos_mbp.z = mb_gas->pos[Z];
		} else {
			halo->gas_only.pos_mbp.x = 0.0;
			halo->gas_only.pos_mbp.y = 0.0;
			halo->gas_only.pos_mbp.z = 0.0;
		}
#  ifdef METALHACK
		/* Get the mean metallicity of the stars */
		halo->mean_z_gas = z_gas_cum / M_gas;
#  endif
	} else {
		reset_SPECIESPROP(&(halo->gas_only));
	}

	/*-----------------
	 * STAR properties
	 *-----------------*/
	if (nstar > 0) {
		/* mass of stars alone */
		halo->stars_only.npart = nstar;
		halo->stars_only.Mass  = M_star;

		/* centre-of-mass properties */
		halo->stars_only.pos_com.x = fmod(Xcom_star / M_star + 1.0, 1.0);
		halo->stars_only.pos_com.y = fmod(Ycom_star / M_star + 1.0, 1.0);
		halo->stars_only.pos_com.z = fmod(Zcom_star / M_star + 1.0, 1.0);
		halo->stars_only.vel.x     = VX_star / M_star;
		halo->stars_only.vel.y     = VY_star / M_star;
		halo->stars_only.vel.z     = VZ_star / M_star;

		/* energy of stars alone */
		halo->stars_only.Ekin = 0.5 * Ekin_star;
		halo->stars_only.Epot = 0.5 * Epot_star;

    if(nstar > AHF_MINPART_STARS)
     {
      /* get spin parameter of star particles (all star part's within sphere!) */
      absAngMom                 = sqrt(pow2(Lx_star) + pow2(Ly_star) + pow2(Lz_star));
      halo->stars_only.AngMom.x = Lx_star / absAngMom;
      halo->stars_only.AngMom.y = Ly_star / absAngMom;
      halo->stars_only.AngMom.z = Lz_star / absAngMom;
      halo->stars_only.lambda   = absAngMom / M_star / sqrt(2. * halo->M_vir * halo->R_vir);
      halo->stars_only.lambda  *= v_fac * sqrt(r_fac / (Grav * m_fac));
      halo->stars_only.lambdaE  = calc_lambdaE(absAngMom, halo->M_vir, halo->stars_only.Mass, halo->stars_only.Ekin, halo->stars_only.Epot);
      
      /* get shape of gas particles (all gas part's within sphere!) */
      itensor[0][0] = a11_star;
      itensor[1][1] = a22_star;
      itensor[2][2] = a33_star;
      itensor[0][1] = a12_star;
      itensor[1][0] = a12_star;
      itensor[0][2] = a13_star;
      itensor[2][0] = a13_star;
      itensor[1][2] = a23_star;
      itensor[2][1] = a23_star;
      get_axes(itensor, &axis1, &axis2, &axis3);
      
      halo->stars_only.axis.x = 1.0;
      if(axis1 > 0.) {
        halo->stars_only.axis.y = sqrt(axis2 / axis1);
        halo->stars_only.axis.z = sqrt(axis3 / axis1);
      }
      else {
        halo->stars_only.axis.y = 0.0;
        halo->stars_only.axis.z = 0.0;
      }
      halo->stars_only.E1.x   = itensor[0][0];
      halo->stars_only.E1.y   = itensor[1][0];
      halo->stars_only.E1.z   = itensor[2][0];
      halo->stars_only.E2.x   = itensor[0][1];
      halo->stars_only.E2.y   = itensor[1][1];
      halo->stars_only.E2.z   = itensor[2][1];
      halo->stars_only.E3.x   = itensor[0][2];
      halo->stars_only.E3.y   = itensor[1][2];
      halo->stars_only.E3.z   = itensor[2][2];
     }

		/* position of most bound star particle */
		if (mb_star != NULL) {
			halo->stars_only.pos_mbp.x = mb_star->pos[X];
			halo->stars_only.pos_mbp.y = mb_star->pos[Y];
			halo->stars_only.pos_mbp.z = mb_star->pos[Z];
		} else {
			halo->stars_only.pos_mbp.x = 0.0;
			halo->stars_only.pos_mbp.y = 0.0;
			halo->stars_only.pos_mbp.z = 0.0;
		}
#  ifdef METALHACK
		/* Get the mean metallicity of the stars */
		halo->mean_z_star = z_star_cum / M_star;
#  endif
	} else {
		reset_SPECIESPROP(&(halo->stars_only));
	}

#ifdef AHFaddDMonlyproperties
	/*-----------------
	 * DM properties
	 *-----------------*/
	if (ndark > 0) {
		/* mass of DM alone */
		halo->DM_only.npart = ndark;
		halo->DM_only.Mass  = M_dark;

		/* centre-of-mass properties */
		halo->DM_only.pos_com.x = fmod(Xcom_dark / M_dark + 1.0, 1.0);
		halo->DM_only.pos_com.y = fmod(Ycom_dark / M_dark + 1.0, 1.0);
		halo->DM_only.pos_com.z = fmod(Zcom_dark / M_dark + 1.0, 1.0);
		halo->DM_only.vel.x     = VX_dark / M_dark;
		halo->DM_only.vel.y     = VY_dark / M_dark;
		halo->DM_only.vel.z     = VZ_dark / M_dark;

		/* energy of DM alone */
		halo->DM_only.Ekin = 0.5 * Ekin_dark;
		halo->DM_only.Epot = 0.5 * Epot_dark;

		/* get spin parameter of dark matter particles (all dark matter part's within sphere!) */
		absAngMom = sqrt(pow2(Lx_dark) + pow2(Ly_dark) + pow2(Lz_dark));
		halo->DM_only.AngMom.x = Lx_dark / absAngMom;
		halo->DM_only.AngMom.y = Ly_dark / absAngMom;
		halo->DM_only.AngMom.z = Lz_dark / absAngMom;
		halo->DM_only.lambda   = absAngMom / M_dark / sqrt(2. * halo->M_vir * halo->R_vir);
		halo->DM_only.lambda  *= v_fac * sqrt(r_fac / (Grav * m_fac));
		halo->DM_only.lambdaE  = calc_lambdaE(absAngMom,
		                                      halo->M_vir,
		                                      halo->DM_only.Mass,
		                                      halo->DM_only.Ekin,
		                                      halo->DM_only.Epot);

		/* get shape of DM particles (all DM part's within sphere!) */
		itensor[0][0] = a11_dark;
		itensor[1][1] = a22_dark;
		itensor[2][2] = a33_dark;
		itensor[0][1] = a12_dark;
		itensor[1][0] = a12_dark;
		itensor[0][2] = a13_dark;
		itensor[2][0] = a13_dark;
		itensor[1][2] = a23_dark;
		itensor[2][1] = a23_dark;
		get_axes(itensor, &axis1, &axis2, &axis3);

		halo->DM_only.axis.x = 1.0;
    if(axis1 > 0.) {
      halo->DM_only.axis.y = sqrt(axis2 / axis1);
      halo->DM_only.axis.z = sqrt(axis3 / axis1);
    }
    else {
      halo->DM_only.axis.y = 0.0;
      halo->DM_only.axis.z = 0.0;
    }
		halo->DM_only.E1.x   = itensor[0][0];
		halo->DM_only.E1.y   = itensor[1][0];
		halo->DM_only.E1.z   = itensor[2][0];
		halo->DM_only.E2.x   = itensor[0][1];
		halo->DM_only.E2.y   = itensor[1][1];
		halo->DM_only.E2.z   = itensor[2][1];
		halo->DM_only.E3.x   = itensor[0][2];
		halo->DM_only.E3.y   = itensor[1][2];
		halo->DM_only.E3.z   = itensor[2][2];

		/* position of most bound DM particle */
		if (mb_dark != NULL) {
			halo->DM_only.pos_mbp.x = mb_dark->pos[X];
			halo->DM_only.pos_mbp.y = mb_dark->pos[Y];
			halo->DM_only.pos_mbp.z = mb_dark->pos[Z];
		} else {
			halo->DM_only.pos_mbp.x = 0.0;
			halo->DM_only.pos_mbp.y = 0.0;
			halo->DM_only.pos_mbp.z = 0.0;
		}
	} else {
		reset_SPECIESPROP(&(halo->DM_only));
	}
#endif /* AHFaddDMonlyproperties */
#endif /* GAS_PARTICLES */

	return TRUE;
} /* HaloProfiles */

#ifdef AHFphspdens
int
HaloProfilesPhaseSpace(HALO *halo)
{
	double  halo_pos[3];
	double  halo_vel[3];
	int     nbins, ibin;
	int     npart_sp, npart_sh;
	double  binpart;
	double  lcur_rad, cur_rad, cur_dist;
	double  ldist_min, ldist_max, ldr;
	double  dist_min, dist_max;
	partptr cur_part, fst_part, lst_part;
	double  x, y, z, vx, vy, vz;
	double  r, theta, phi, vr, vtheta, vphi;
	double  tmp;
	double  mean_vx_sp, mean_vy_sp, mean_vz_sp;
	double  mean_vr_sp, mean_vtheta_sp, mean_vphi_sp;
	double  sum_vx_sp, sum_vy_sp, sum_vz_sp;
	double  sum_vr_sp, sum_vtheta_sp, sum_vphi_sp;
	double  mean_vx_sh, mean_vy_sh, mean_vz_sh;
	double  mean_vr_sh, mean_vtheta_sh, mean_vphi_sh;
	double  sigma2_vx_sh, sigma2_vy_sh, sigma2_vz_sh;
	double  sigma2_vr_sh, sigma2_vtheta_sh, sigma2_vphi_sh;
	long    jpart, fst_jpart, lst_jpart;

	/* Don't do it for haloes too small */
	if (halo->npart < simu.AHF_MINPART)
		return TRUE;

	/* Get the halo restframe */
	halo_pos[0] = halo->pos.x;
	halo_pos[1] = halo->pos.y;
	halo_pos[2] = halo->pos.z;
	halo_vel[0] = halo->vel.x;
	halo_vel[1] = halo->vel.y;
	halo_vel[2] = halo->vel.z;


	/* how many bins from where to where for density profile? */
	binning_parameter(*halo, &nbins, &dist_min, &dist_max);

	/* logarithmical binning from dist_min to dist_max */
	ldist_min = log10(dist_min);
	ldist_max = log10(dist_max);
	ldr       = (ldist_max - ldist_min) / (double)nbins;

	/* Reset all quantities in the sphere */
	cur_dist   = -1.0;
	jpart      = 0;
	npart_sp   = 0;
	mean_vx_sp = mean_vy_sp = mean_vz_sp = 0.0;
	mean_vr_sp = mean_vtheta_sp = mean_vphi_sp = 0.0;
	sum_vx_sp  = sum_vy_sp = sum_vz_sp = 0.0;
	sum_vr_sp  = sum_vtheta_sp = sum_vphi_sp = 0.0;

	/* Loop over all bins */
	for (ibin = 0; ibin < nbins; ibin++) {
		/* Current outer radius */
		lcur_rad = ldist_min + ((double)ibin + 1) * ldr;
		cur_rad  = pow(10., lcur_rad);

		/* Figure out first and last particle */
		cur_part  = global.fst_part + halo->ipart[jpart];
		fst_part  = lst_part = cur_part;
		fst_jpart = lst_jpart = jpart;
		npart_sh  = 0;
		while (jpart < halo->npart && cur_dist < cur_rad) {
			/* Access particle*/
			cur_part = global.fst_part + halo->ipart[jpart];

			/* Get position */
			x = ((double)cur_part->pos[0]) - halo_pos[0];
			y = ((double)cur_part->pos[1]) - halo_pos[1];
			z = ((double)cur_part->pos[2]) - halo_pos[2];

			/* Correct for periodic boundary */
			if (x < -0.5)				x += 1.0;
			if (y < -0.5)				y += 1.0;
			if (z < -0.5)				z += 1.0;
			if (x > 0.5)				x -= 1.0;
			if (y > 0.5)				y -= 1.0;
			if (z > 0.5)				z -= 1.0;

			/* Get the distance from the center */
			cur_dist = sqrt(x * x + y * y + z * z);

			/* Check whether that particle counts */
			if (cur_dist < cur_rad) {
				lst_part  = cur_part;
				lst_jpart = jpart;
				jpart++;
				npart_sp++;
				npart_sh++;
			}
		}

		/* Loop over all particles in the bin to get mean velocities */
		mean_vx_sh = mean_vy_sh = mean_vz_sh = 0.0;
		mean_vr_sh = mean_vtheta_sh = mean_vphi_sh = 0.0;
		if (npart_sh > 0) {
			for (jpart = fst_jpart;
			     jpart <= lst_jpart;
			     jpart++) {
				/* Access particle */
				cur_part = global.fst_part + halo->ipart[jpart];

				/* Particle position in halo rest frame */
				x = ((double)cur_part->pos[0]) - halo_pos[0];
				y = ((double)cur_part->pos[1]) - halo_pos[1];
				z = ((double)cur_part->pos[2]) - halo_pos[2];

				/* Correct for periodic boundary */
				if (x < -0.5)					x += 1.0;
				if (y < -0.5)					y += 1.0;
				if (z < -0.5)					z += 1.0;
				if (x > 0.5)					x -= 1.0;
				if (y > 0.5)					y -= 1.0;
				if (z > 0.5)					z -= 1.0;

				/* Used for the Hubble flow correction */
				tmp = Hubble * r_fac / v_fac;

				/* Particle velocity in halo rest frame w/ Hubble */
				vx = ((double)cur_part->mom[0]) - halo_vel[0] + x * tmp;
				vy = ((double)cur_part->mom[1]) - halo_vel[1] + y * tmp;
				vz = ((double)cur_part->mom[2]) - halo_vel[2] + z * tmp;

				/* Put positions to spherical coordinates */
				r     = sqrt(x * x + y * y + z * z);
				theta = acos(z / r);
				phi   = atan2(y, x);

				/* Put velocities to spherical coordinates */
				vr     = 1. / r * (x * vx + y * vy + z * vz);
				vtheta = -1.0
				         / sqrt(1 - z * z
				                / (r * r)) * (vz / r - z * vr / (r * r));
				vphi = 1.
				       / (1 + y * y
				          / (x * x)) * (vy / x - y * vx / (x * x));

				/* Keep track of the mean velocities in the shell
				 * Here it is only the sum, division by npart takes place
				 * after the loop */
				mean_vx_sh     += vx;
				mean_vy_sh     += vy;
				mean_vz_sh     += vz;
				mean_vr_sh     += vr;
				mean_vtheta_sh += vtheta;
				mean_vphi_sh   += vphi;

				/* Keep track of the sum of the velocities in the sphere */
				sum_vx_sp     += vx;
				sum_vy_sp     += vy;
				sum_vz_sp     += vz;
				sum_vr_sp     += vr;
				sum_vtheta_sp += vtheta;
				sum_vphi_sp   += vphi;
			}
			/* Calculate real mean velocity in the shell */
			tmp             = 1. / npart_sh;
			mean_vx_sh     *= tmp;
			mean_vy_sh     *= tmp;
			mean_vz_sh     *= tmp;
			mean_vr_sh     *= tmp;
			mean_vtheta_sh *= tmp;
			mean_vphi_sh   *= tmp;
		} /* if (npart_sh > 0)*/

		/* Calculate the real mean velocity within the sphere */
		tmp            = 1. / npart_sp;
		mean_vx_sp     = sum_vx_sp * tmp;
		mean_vy_sp     = sum_vy_sp * tmp;
		mean_vz_sp     = sum_vz_sp * tmp;
		mean_vr_sp     = sum_vr_sp * tmp;
		mean_vtheta_sp = sum_vtheta_sp * tmp;
		mean_vphi_sp   = sum_vphi_sp * tmp;

		/* Loop again over all particles in the bin to get dispersion */
		sigma2_vx_sh = sigma2_vy_sh = sigma2_vz_sh = 0.0;
		sigma2_vr_sh = sigma2_vtheta_sh = sigma2_vphi_sh = 0.0;
		if (npart_sh > 0) {
			for (jpart = fst_jpart;
			     jpart <= lst_jpart;
			     jpart++) {
				/* Access particle */
				cur_part = global.fst_part + halo->ipart[jpart];

				/* Particle position in halo rest frame */
				x = ((double)cur_part->pos[0]) - halo_pos[0];
				y = ((double)cur_part->pos[1]) - halo_pos[1];
				z = ((double)cur_part->pos[2]) - halo_pos[2];

				/* Correct for periodic boundary */
				if (x < -0.5)					x += 1.0;
				if (y < -0.5)					y += 1.0;
				if (z < -0.5)					z += 1.0;
				if (x > 0.5)					x -= 1.0;
				if (y > 0.5)					y -= 1.0;
				if (z > 0.5)					z -= 1.0;

				/* Used for the Hubble flow correction */
				tmp = Hubble * r_fac / v_fac;

				/* Particle velocity in halo rest frame w/ Hubble */
				vx = ((double)cur_part->mom[0]) - halo_vel[0] + x * tmp;
				vy = ((double)cur_part->mom[1]) - halo_vel[1] + y * tmp;
				vz = ((double)cur_part->mom[2]) - halo_vel[2] + z * tmp;

				/* Put positions to spherical coordinates */
				r     = sqrt(x * x + y * y + z * z);
				theta = acos(z / r);
				phi   = atan2(y, x);

				/* Put velocities to spherical coordinates */
				vr     = 1. / r * (x * vx + y * vy + z * vz);
				vtheta = -1.0
				         / sqrt(1 - z * z
				                / (r * r)) * (vz / r - z * vr / (r * r));
				vphi = 1.
				       / (1 + y * y
				          / (x * x)) * (vy / x - y * vx / (x * x));

				/* Keep track of shell dispersions*(npart_sh-1) */
				sigma2_vx_sh     += pow2(vx - mean_vx_sh);
				sigma2_vy_sh     += pow2(vy - mean_vy_sh);
				sigma2_vz_sh     += pow2(vz - mean_vz_sh);
				sigma2_vr_sh     += pow2(vr - mean_vr_sh);
				sigma2_vtheta_sh += pow2(vtheta - mean_vtheta_sh);
				sigma2_vphi_sh   += pow2(vphi - mean_vphi_sh);
			}
			/* Get the real dispersions in the shell */
			tmp               = 1. / (npart_sh - 1);
			sigma2_vx_sh     *= tmp;
			sigma2_vy_sh     *= tmp;
			sigma2_vz_sh     *= tmp;
			sigma2_vr_sh     *= tmp;
			sigma2_vtheta_sh *= tmp;
			sigma2_vphi_sh   *= tmp;
		} /* if (npart > 0) */

		/* Plug the quantities into the profile array */
		halo->prof.sigma2_vx_sh[ibin]     = sigma2_vx_sh;
		halo->prof.sigma2_vy_sh[ibin]     = sigma2_vy_sh;
		halo->prof.sigma2_vz_sh[ibin]     = sigma2_vz_sh;
		halo->prof.sigma2_vr_sh[ibin]     = sigma2_vr_sh;
		halo->prof.sigma2_vtheta_sh[ibin] = sigma2_vtheta_sh;
		halo->prof.sigma2_vphi_sh[ibin]   = sigma2_vphi_sh;
#  ifdef AHFmeanvelocities
		halo->prof.mean_vx_sh[ibin]       = mean_vx_sh;
		halo->prof.mean_vy_sh[ibin]       = mean_vy_sh;
		halo->prof.mean_vz_sh[ibin]       = mean_vz_sh;
		halo->prof.mean_vr_sh[ibin]       = mean_vr_sh;
		halo->prof.mean_vtheta_sh[ibin]   = mean_vtheta_sh;
		halo->prof.mean_vphi_sh[ibin]     = mean_vphi_sh;
		halo->prof.mean_vx_sp[ibin]       = mean_vx_sp;
		halo->prof.mean_vy_sp[ibin]       = mean_vy_sp;
		halo->prof.mean_vz_sp[ibin]       = mean_vz_sp;
		halo->prof.mean_vr_sp[ibin]       = mean_vr_sp;
		halo->prof.mean_vtheta_sp[ibin]   = mean_vtheta_sp;
		halo->prof.mean_vphi_sp[ibin]     = mean_vphi_sp;
#  endif
	} /* End of for-loop over the bins */

	return TRUE;
} /* HaloProfilesPhaseSpace */

#endif   /* AHFphspdens */

#ifdef AHFdisks
int
HaloProfilesDisk(HALO *halo)
{
	double  Xc, Yc, Zc;
  double  VXc, VYc, VZc;
	int     nbins, ibin;
	double  lcur_rad, cur_rad, cur_dist;
	double  ldist_min, ldist_max, ldr;
	double  dist_min, dist_max;
  double  weight;
	partptr cur_part;
	double  dx, dy, dz;
	double  dvx, dvy, dvz;
	long    jpart;
  
  double  k_gas, k_star;
  double  Lx, Ly, Lz, Lzz;
  
  double  Lx_gas_bin, Ly_gas_bin, Lz_gas_bin, absL_gas_bin;
  double  Lx_star_bin, Ly_star_bin, Lz_star_bin, absL_star_bin;
  
	/* Don't do it for haloes too small */
	if (halo->npart < simu.AHF_MINPART)
		return TRUE;
  
	/* Get the halo restframe */
	Xc  = halo->pos.x;
	Yc  = halo->pos.y;
	Zc  = halo->pos.z;
	VXc = halo->vel.x;
	VYc = halo->vel.y;
	VZc = halo->vel.z;
  
  
	/* how many bins from where to where for density profile? */
	binning_parameter(*halo, &nbins, &dist_min, &dist_max);
  if(nbins != halo->prof.nbins) {
    fprintf(io.logfile, "    HaloProfilesDisk():     nbins=%d vs. halo->nbins=%d THIS SHOULD NOT HAPPEN! INVESTIGATE!\n",nbins,halo->prof.nbins);
    fflush(io.logfile);    
  }
  
	/* logarithmical binning from dist_min to dist_max */
	ldist_min = log10(dist_min);
	ldist_max = log10(dist_max);
	ldr       = (ldist_max - ldist_min) / (double)nbins;
  
	/* Reset all quantities in the sphere */
  k_gas    = 0.0;
  k_star   = 0.0;
	cur_dist = -1.0;
	jpart    = 0;
  
  /*====================
   * loop over all bins
   *====================*/
	for (ibin = 0; ibin < nbins; ibin++) {

		lcur_rad  = ldist_min + ((double)ibin + 1) * ldr;
		cur_rad   = pow(10., lcur_rad);    
    
    Lx_gas_bin    = halo->prof.Lx_gas[ibin];
    Ly_gas_bin    = halo->prof.Ly_gas[ibin];
    Lz_gas_bin    = halo->prof.Lz_gas[ibin];
    absL_gas_bin  = sqrt(pow2(Lx_gas_bin)+pow2(Ly_gas_bin)+pow2(Lz_gas_bin));
    
    Lx_star_bin   = halo->prof.Lx_star[ibin];
    Ly_star_bin   = halo->prof.Ly_star[ibin];
    Lz_star_bin   = halo->prof.Lz_star[ibin];
    absL_star_bin = sqrt(pow2(Lx_star_bin)+pow2(Ly_star_bin)+pow2(Lz_star_bin));
    
    /*========================================
     * loop over all particles in present bin
     *========================================*/
		while (cur_dist < cur_rad && jpart < halo->npart) {
      
			/* Access particle */
			cur_part = global.fst_part + halo->ipart[jpart];
      
      /* Get particle mass */
#ifdef MULTIMASS
      weight = cur_part->weight;
#else
      weight = 1.0;
#endif
      
			/* Get particle position in halo rest-frame */
			dx = ((double)cur_part->pos[X]) - Xc;
			dy = ((double)cur_part->pos[Y]) - Yc;
			dz = ((double)cur_part->pos[Z]) - Zc;
      
			/* Correct for periodic boundary */
			if (dx < -0.5)			dx += 1.0;
			if (dy < -0.5)			dy += 1.0;
			if (dz < -0.5)			dz += 1.0;
			if (dx > 0.5)				dx -= 1.0;
			if (dy > 0.5)				dy -= 1.0;
			if (dz > 0.5)				dz -= 1.0;
      
      /* Get particle velocity in halo rest-frame */
      dvx = ((double)cur_part->mom[X]) - VXc;
      dvy = ((double)cur_part->mom[Y]) - VYc;
      dvz = ((double)cur_part->mom[Z]) - VZc;
      
			/* Get the distance from the center */
			cur_dist = sqrt(pow2(dx) + pow2(dy) + pow2(dz));
      
			/*==============
			 * GAS PARTICLE
			 *==============*/
			if (isgreaterequal(cur_part->u, PGAS)) {

        Lx = weight * (dy * dvz - dz * dvy);
        Ly = weight * (dz * dvx - dx * dvz);
        Lz = weight * (dx * dvy - dy * dvx);
        
        if(absL_gas_bin > 0.) {
          Lzz = (Lx*Lx_gas_bin + Ly*Ly_gas_bin + Lz*Lz_gas_bin)/absL_gas_bin;
        }
        else {
          Lzz = 0.;
        }

        if (cur_dist > MACHINE_ZERO)
          k_gas += pow2(Lzz/cur_dist) / weight;
       }
      
      /*===============
			 * STAR PARTICLE
			 *===============*/
			else if (fabs(cur_part->u - PSTAR) < ZERO) {

        Lx = weight * (dy * dvz - dz * dvy);
        Ly = weight * (dz * dvx - dx * dvz);
        Lz = weight * (dx * dvy - dy * dvx);

        if(absL_star_bin > 0.) {
          Lzz = (Lx*Lx_star_bin + Ly*Ly_star_bin + Lz*Lz_star_bin)/absL_star_bin;
        }
        else {
          Lzz = 0.;
        }
        
        if (cur_dist > MACHINE_ZERO)
          k_star += pow2(Lzz/cur_dist) / weight;
      }
      /*===============
			 * DM PARTICLE
			 *===============*/
      /* we are not performing any operations for DM particles */

      /* move to next particle */
      jpart++;
		}
    
    /* store k(<R_ibin) values in halo profile */
    halo->prof.k_gas[ibin]  = 0.5 * k_gas;
    halo->prof.k_star[ibin] = 0.5 * k_star;
    
	} /* End of for-loop over the bins */
  
	return TRUE;
}
#endif

/*==============================================================================
 *  check Rho[] for rise and return position of rise
 *  in addition, check if Rho[] declines like r[]^AHF_SLOPE out to rising point
 *==============================================================================*/
double
get_haloedge(double *r, double *Rho, int iRadOut, int ismooth)
{
	int    ir, iRadIn, iRadEnd, nrise;
	int    irm1, irp1, irp2;
	double gradlf, gradrt, R_edge;
	double Rholeft, Rhoright, rleft, rright, Rhomin, rmin;

	/* is there a well pronounced minimum in the density profile? */
#ifdef AHFsplinefit
	find_min_spline(r, Rho, iRadOut, ismooth, &rmin, &Rhomin, 10000);
#else
	find_min(r, Rho, iRadOut, ismooth, &rmin, &Rhomin);
#endif

	if (rmin < r[iRadOut - 1]) {
		return rmin;
	} else {
		/* loop inwards until above virial overdensity */
		for (ir = iRadOut - 2; ir > 0; ir--) {
			if (Rho[ir] > global.ovlim) {
				rleft    = r[ir];
				rright   = r[ir + 1];
				Rholeft  = Rho[ir];
				Rhoright = Rho[ir + 1];

				/* get interpolated virial radius */
				R_edge
				    = (rleft
				       * (global.ovlim
				          - Rhoright) + rright * (Rholeft - global.ovlim))
				      /                  (Rholeft - Rhoright);
				goto found_edge;
			}
		}

		/* Rho never below virial overdensity -> search for upturn radius */
		iRadIn = 0;

		/* does the profile ever rise indicative of a sub-halo or mis-placed centre */
		nrise = 0;
		for (ir = iRadIn + 1; ir < iRadOut; ir++) {
			/* profile is rising */
			if (Rho[ir] >= AHF_RISE * Rho[ir - 1])
				nrise++;
			else
				nrise = 0;

			if (nrise > AHF_MAXNRISE) {
				goto found_rise;
			}
		}

found_rise:
		iRadEnd = MIN(ir, iRadOut - 1);

		/* reduce iRadEnd to the point where profiles declines faster than r^-slope */
		irm1   = MAX(iRadEnd - 1, iRadIn);
		irp1   = MIN(iRadEnd + 1, iRadOut - 1);
		gradlf = Rho[irm1] * pow(r[irm1], AHF_SLOPE);
		gradrt = Rho[irp1] * pow(r[irp1], AHF_SLOPE);

		while (gradlf < gradrt) {
			iRadEnd--;

			if (iRadEnd <= iRadIn) {
				return r[iRadIn];
			}

			irm1   = MAX(iRadEnd - 1, iRadIn);
			irp1   = MIN(iRadEnd + 1, iRadOut - 1);
			gradlf = Rho[irm1] * pow(r[irm1], AHF_SLOPE);
			gradrt = Rho[irp1] * pow(r[irp1], AHF_SLOPE);
		}

		R_edge = r[iRadEnd];

found_edge:
		return R_edge;
	}
} /* get_haloedge */

/*==============================================================================
 * sort the particles belonging to a halo with respects to distance
 *==============================================================================*/
void
sort_halo_particles(HALO *halo)
{
	double        Xc, Yc, Zc, Xp, Yp, Zp;
	long unsigned npart, i;
	double        *dist, dX, dY, dZ, dR;
	long unsigned *idx;
	long unsigned *ipart;
	partptr       cur_part;
    
	if (halo->npart < simu.AHF_MINPART)
      return;
  
#ifdef VERBOSE2
	fprintf(io.logfile,
	        "    sort_particles:      npart=%12ld => ", halo->npart);
	fflush(io.logfile);
#endif

	Xc    = halo->pos.x;
	Yc    = halo->pos.y;
	Zc    = halo->pos.z;
	npart = halo->npart;
  

	ipart = (long unsigned *)calloc(npart, sizeof(long unsigned));
	idx   = (long unsigned *)calloc(npart + 1, sizeof(long unsigned));
	dist  = (double *)calloc(npart + 1, sizeof(double));
  

	/* fill dist[] array with dR^2 */
	for (i = 0; i < npart; i++) {
		
    /* access to particle */
		cur_part = global.fst_part + halo->ipart[i];
    
		/* particle position */
		Xp = (double)cur_part->pos[X];
		Yp = (double)cur_part->pos[Y];
		Zp = (double)cur_part->pos[Z];
    
		/* put particle into halo rest frame */
		dX = (Xp - Xc);
		dY = (Yp - Yc);
		dZ = (Zp - Zc);
    
		/* take care of periodic boundary conditions */
		if (dX > 0.5)			dX -= 1.0;
		if (dY > 0.5)			dY -= 1.0;
		if (dZ > 0.5)			dZ -= 1.0;
		if (dX < -0.5)		dX += 1.0;
		if (dY < -0.5)		dY += 1.0;
		if (dZ < -0.5)		dZ += 1.0;
    
		/* distance^2 of current particle */
		dR          = (pow2(dX) + pow2(dY) + pow2(dZ));
    
		dist[i + 1] = dR;
	}
  
	/* sort particles according to dist[] */
	indexx(npart, dist, idx);
  
  
	/* generate an ordered array ... */
	for (i = 0; i < npart; i++)
		ipart[i] = halo->ipart[idx[i + 1] - 1];
  
  
	/* ... and replace halo.ipart[] */
	for (i = 0; i < npart; i++)
		halo->ipart[i] = ipart[i];
  
	free(idx);
	free(dist);
	free(ipart);
  
#ifdef VERBOSE2
	Xc = halo->pos.x;
	Yc = halo->pos.y;
	Zc = halo->pos.z;
  
	for (i = 0; i < 1; i++) {
		cur_part = global.fst_part + halo->ipart[i];
		Xp       = (double)cur_part->pos[X];
		Yp       = (double)cur_part->pos[Y];
		Zp       = (double)cur_part->pos[Z];
		dX       = (Xp - Xc);
		dY       = (Yp - Yc);
		dZ       = (Zp - Zc);
		if (dX > 0.5)	 dX -= 1.0;
		if (dY > 0.5)	 dY -= 1.0;
		if (dZ > 0.5)	 dZ -= 1.0;
		if (dX < -0.5) dX += 1.0;
		if (dY < -0.5) dY += 1.0;
		if (dZ < -0.5) dZ += 1.0;
		dR = (pow2(dX) + pow2(dY) + pow2(dZ));
    
#ifdef VERBOSE2
		fprintf(io.logfile, "dR_min = %16.8g kpc/h\n", sqrt(dR) * x_fac * 1000.);
		fflush(io.logfile);
#endif
		/*  fprintf(stderr,"dR[%ld] = %16.8g kpc/h\n",i,sqrt(dR)*x_fac*1000.); */
	}
#endif
  
	Xc = halo->pos.x;
	Yc = halo->pos.y;
	Zc = halo->pos.z;
} /* sort_halo_particles */


/*================================================================
 * check whether the subhalo lies within the virial radius or not
 *================================================================*/
boolean check_subhalo(HALO *host, HALO *sub)
{
  double dx,dy,dz;
  
  dx = fabs(host->pos.x - sub->pos.x);
  dy = fabs(host->pos.y - sub->pos.y);
  dz = fabs(host->pos.z - sub->pos.z);
  if(dx > 0.5) dx=1.-dx;
  if(dy > 0.5) dy=1.-dy;
  if(dz > 0.5) dz=1.-dz;
  
  /* we add it to the subhalo list as soon as the two radii overlap */
  if(pow2(dx)+pow2(dy)+pow2(dz) < pow2(host->R_vir + AHF_HOSTSUBOVERLAP*sub->R_vir) && sub->npart >= simu.AHF_MINPART && host->npart >= simu.AHF_MINPART)
    return(TRUE);
  else
    return(FALSE);
}

#if (defined AHFrestart || defined WITH_MPI)
void
rem_boundary_haloes(void)
{
	int         i;
	sfc_key_t   key;
	sfc_curve_t ctype;
	sfc_key_t   minkey;
	sfc_key_t   maxkey;
	int         level;
	HALO        *tmpHalos;
	int         k;
	int         new_numHalos = 0;

	io_logging_msg(global_io.log, INT32_C(0),
	               "\nFiguring out which halo can be ignored.");

#  ifdef AHFrestart
	minkey = global_info.minkey;
	maxkey = global_info.maxkey;
	ctype  = global_info.ctype;
	level  = global_info.level;
#  else
	minkey = global_info.loadbal->fstkey[global_mpi.rank];
	maxkey = global_info.loadbal->lstkey[global_mpi.rank];
	ctype  = global_info.loadbal->ctype;
	level  = global_info.loadbal->level;
#  endif

	/* Loop over all halos and figure out which should be used */
#ifdef WITH_OPENMP
#  pragma omp parallel for                \
	schedule (dynamic)                    \
	shared(numHalos, halos, level, ctype, \
	global_io, minkey, maxkey)            \
	private(key)                          \
	reduction (+:new_numHalos)
#endif
	for (i = 0; i < numHalos; i++) {
		/* First get the key*/
		key = sfc_curve_calcKey(ctype,
		                        (double)(halos[i].pos.x),
		                        (double)(halos[i].pos.y),
		                        (double)(halos[i].pos.z),
		                        BITS_PER_DIMENSION);
		/* Reduce it */
		key = sfc_curve_contract(level,
		                         BITS_PER_DIMENSION,
		                         ctype,
		                         key);

		/* Check if it is with our key range */
		if ((key >= minkey) && (key <= maxkey)) {
			halos[i].ignoreme = FALSE;
			new_numHalos++;
		} else {
			halos[i].ignoreme = TRUE;
		}
	}
	io_logging_msg(global_io.log, INT32_C(0),
	               "Done with identifying halos to be ignored.");

	/* Now remove them */
	io_logging_msg(global_io.log, INT32_C(0),
	               "Now really getting rid of outsider halos.");

	/* That's going to be the new Halo list */
	tmpHalos = (HALO *)malloc(sizeof(HALO) * new_numHalos);
	if (tmpHalos == NULL) {
    fprintf(stderr,
            "Not enough memory to allocate tmpHalos array "
            "in %s.  Aborting.\n\n", __func__);
		common_terminate(EXIT_FAILURE);
  }

	/* Loop over all halos and only pick the good ones */
	new_numHalos = 0;
	for (k = 0; k < numHalos; k++) {
		if (halos[k].ignoreme == FALSE) {
			memcpy((void *)(tmpHalos + new_numHalos),
			       (const void *)(halos + k),
			       sizeof(HALO));
			new_numHalos++;
		}
    else {
      /* free() those ipart[] and SubStruct[] arrays as well as the profile for haloes not kept */
      if(halos[k].ipart)     free(halos[k].ipart);
      if(halos[k].subStruct) free(halos[k].subStruct);
      dest_profile(&(halos[k]));
    }
	}

	/* Status information */
	io_logging_msg(global_io.log, INT32_C(0),
	               "Old number of halos: %i   "
	               "New number of halos: %i",
	               (int)numHalos, (int)new_numHalos);

	/* Replace the old halo list with the new one */
	free(halos);
	halos    = tmpHalos;
	numHalos = new_numHalos;
	io_logging_msg(global_io.log, INT32_C(0),
	               "Done with getting rid of halos.");

	/* Done */
	return;
} /* rem_boundary_haloes */

void
flag_boundary_haloes(void)
{
	int         i;
	sfc_key_t   key;
	sfc_curve_t ctype;
	sfc_key_t   minkey;
	sfc_key_t   maxkey;
	int         level;
	int         k;
	int         new_numHalos = 0;
  
	io_logging_msg(global_io.log, INT32_C(0), "\nFiguring out which halo can be ignored.");
  
#  ifdef AHFrestart
	minkey = global_info.minkey;
	maxkey = global_info.maxkey;
	ctype  = global_info.ctype;
	level  = global_info.level;
#  else
	minkey = global_info.loadbal->fstkey[global_mpi.rank];
	maxkey = global_info.loadbal->lstkey[global_mpi.rank];
	ctype  = global_info.loadbal->ctype;
	level  = global_info.loadbal->level;
#  endif
  
	/* Loop over all halos and figure out which should be used */
#ifdef WITH_OPENMP
#  pragma omp parallel for                \
schedule (dynamic)                    \
shared(numHalos, halos, level, ctype, \
global_io, minkey, maxkey)            \
private(key)                          \
reduction (+:new_numHalos)
#endif
	for (i = 0; i < numHalos; i++) {
		/* First get the key*/
		key = sfc_curve_calcKey(ctype,
		                        (double)(halos[i].pos.x),
		                        (double)(halos[i].pos.y),
		                        (double)(halos[i].pos.z),
		                        BITS_PER_DIMENSION);
		/* Reduce it */
		key = sfc_curve_contract(level,
		                         BITS_PER_DIMENSION,
		                         ctype,
		                         key);
    
		/* Check if it is with our key range */
		if ((key >= minkey) && (key <= maxkey)) {
			halos[i].ignoreme = FALSE;
			new_numHalos++;
		} else {
			halos[i].ignoreme = TRUE;
		}
	}
	io_logging_msg(global_io.log, INT32_C(0), "Done with flagging halos to be ignored: numHalos=%ld new_numHalos=%ld",numHalos,new_numHalos);

	/* Done */
	return;
} /* flag_boundary_haloes */
#endif   /* (defined AHFrestart || defined WITH_MPI) */

#ifdef AHFexciseSubhaloStars
/*=====================================================================
 * extract all those stars from the host halo that belong to subhaloes
 *=====================================================================*/
int id_comp(const void *id1, const void *id2)
{
	if      ( *((long *)id1) > *((long*)id2) )
		return 1;
	else if ( *((long *)id1) < *((long*)id2) )
		return -1;
	else
		return 0;
}

void exciseSubhaloStars(HALO *halos, long ihost)
{
  long *idsubstar;
  long *idhoststar;
  long *ipart_uniquestars;
  long unsigned nhoststar, nsubstar, k, khoststar, ksubstar, isub, i, npart_uniquestars;
	partptr       cur_part;
  
  /* only excise stars if there is substructure present */
  if(halos[ihost].numSubStruct > 0) {
    
    /* allocate temporary arrays (one malloc() for an array big enough to hold everything) */
    nhoststar  = halos[ihost].stars_only.npart;
    
    idhoststar = (long *) malloc(nhoststar*sizeof(long));
    idsubstar  = (long *) malloc(nhoststar*sizeof(long));
    if(idsubstar == NULL || idhoststar == NULL) {
      fprintf(stderr, "Not enough memory to allocate arrays in %s.  Aborting.\n\n", __func__);
      common_terminate(EXIT_FAILURE);
    }
    
    /* copy host star particle IDs over to temporary array */
    khoststar = 0;
    for(k=0; k<halos[ihost].npart; k++) {
      cur_part = global.fst_part + halos[ihost].ipart[k];
      if(fabs(cur_part->u - PSTAR) < ZERO) {
        idhoststar[khoststar] = cur_part - global.fst_part;
        khoststar++;
      }
    }
    /* double-check number of star particles */
    if(khoststar != nhoststar) {
      fprintf(stderr,"exciseSubhaloStars(1): number of stars in host does not add up: khoststar=%ld nhoststar=%ld\n",khoststar,nhoststar);
      common_terminate(EXIT_FAILURE);
    }
    
    
    /* copy all subhalo star particle IDs over to temporary array */
    ksubstar = 0;
    for(i=0; i<halos[ihost].numSubStruct; i++) {
      isub = halos[ihost].subStruct[i];
      
      for(k=0; k<halos[isub].npart; k++) {
        cur_part = global.fst_part + halos[isub].ipart[k];
        if(fabs(cur_part->u - PSTAR) < ZERO) {
          idsubstar[ksubstar] = cur_part - global.fst_part;
          ksubstar++;
          
          // just in case all subhaloes combined contain more stars than the host (this can happen if there are sub-subhaloes present!)
          if(ksubstar >= nhoststar)
            idsubstar = (long *) realloc(idsubstar, (ksubstar+1)*sizeof(long));
        }
      }
    }
    nsubstar = ksubstar;
    
    /* the subhaloes contain stars */
    if(nsubstar > 0)
     {
      /* remove excess array elements */
      //idhoststar = (long *) realloc(idhoststar, nhoststar*sizeof(long));
      idsubstar  = (long *) realloc(idsubstar,  nsubstar *sizeof(long));
      if(idsubstar == NULL || idhoststar == NULL) {
        fprintf(stderr, "Not enough memory to re-allocate arrays in %s (nhoststar=%ld, nsubstar=%ld).  Aborting.\n\n", __func__,nhoststar,nsubstar);
        common_terminate(EXIT_FAILURE);
      }
      
      /* qsort both idhoststar[] and idsubstar[] arrays */
      qsort((void *)idhoststar, nhoststar, sizeof(long), id_comp);
      qsort((void *)idsubstar,  nsubstar,  sizeof(long), id_comp);
      
      /* check every single host star for being in a subhalo */
      npart_uniquestars  = 0;
      ipart_uniquestars  = (long *) malloc(nhoststar*sizeof(long));
      for(k=0; k<nhoststar; k++) {
        
        /* host star id cannot be found in subhalo list -> keep it! */
        if(bsearch(&(idhoststar[k]), idsubstar, nsubstar, sizeof(long), id_comp) == NULL) {
          ipart_uniquestars[npart_uniquestars] = idhoststar[k];
          npart_uniquestars++;
        }
        
      }
      /* free excess memory in ipart_uniquestars[] array */
      ipart_uniquestars = (long *) realloc(ipart_uniquestars, npart_uniquestars*sizeof(long));
      
      /* store informaion in host halo structure */
      halos[ihost].npart_uniquestars = npart_uniquestars;
      halos[ihost].ipart_uniquestars = ipart_uniquestars;
      
      /* free both temporary arrays again */
      free(idsubstar);
      free(idhoststar);
     }
    
    /* the subhaloes do not contain any stars */
    else
     {
      /* simply use all host stars */
      halos[ihost].npart_uniquestars = nhoststar;
      halos[ihost].ipart_uniquestars = idhoststar;
      
      /* free temporary array again (do not free idhoststar as it has to be kept alive!) */
      free(idsubstar);
     }
    
  } // if(numSubStruct>0)
  
  /* if there is no substructure, simply copy all stars over to 'uniquestars' */
  else {
    /* there is no substructure but we need to copy the star particles over to the "unique" array anyways! */
    halos[ihost].npart_uniquestars = halos[ihost].stars_only.npart;
    halos[ihost].ipart_uniquestars = (long *) malloc(halos[ihost].stars_only.npart*sizeof(long));
    
    if(halos[ihost].ipart_uniquestars == NULL) {
      fprintf(stderr, "Not enough memory to allocate array in %s (halos[ihost].stars_only.npart=%ld).  Aborting.\n\n", __func__,halos[ihost].stars_only.npart);
      common_terminate(EXIT_FAILURE);
    }
    
    khoststar = 0;
    for(k=0; k<halos[ihost].npart; k++) {
      cur_part = global.fst_part + halos[ihost].ipart[k];
      if(fabs(cur_part->u - PSTAR) < ZERO) {
        halos[ihost].ipart_uniquestars[khoststar] = halos[ihost].ipart[k];
        khoststar++;
      }
    }
    
    /* double-check number of star particles */
    if(khoststar != halos[ihost].stars_only.npart) {
      fprintf(stderr,"exciseSubhaloStars(2): number of stars in host does not add up: khoststar=%ld halos[ihost].stars_only.npart=%ld\n",khoststar,halos[ihost].stars_only.npart);
      common_terminate(EXIT_FAILURE);
    }

  }
  
  // update Mstar_excised
  halos[ihost].Mstar_excised = 0.0;
  for(i=0; i<halos[ihost].npart_uniquestars; i++) {
    cur_part = global.fst_part + halos[ihost].ipart_uniquestars[i];
    halos[ihost].mean_z_star_excised   += cur_part->weight * (cur_part->z);
    halos[ihost].Mstar_excised         += cur_part->weight;
  }
  halos[ihost].mean_z_star_excised /= halos[ihost].Mstar_excised;
  
}
#endif /* AHFexciseSubhaloStars */

#endif // AHF
