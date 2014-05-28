/*==========================================================================
 *
 * This is an adaption of AHF to the new libtree.a
 *
 *==========================================================================*/

#ifdef AHF2

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
#include "ahf_halos_sfc.h"
#include "ahf_io.h"
#include "patchtree2halos.h"
#include "../libutility/utility.h"
#include "../libutility/timer.h"
#include "../libamr_serial/amr_serial.h"


#ifdef AHF_SQL
#  include "ahf_io_sql.h"
#endif


/***************************************************************************
 *   Macros and assorted definitions                                       *
 ***************************************************************************/


/***************************************************************************
 *   Structure and type definitions                                        *
 ***************************************************************************/



/***************************************************************************
 *   Prototypes of local functions                                         *
 ***************************************************************************/
double  get_haloedge                (double *, double *, int, int);
boolean check_subhalo               (HALO *host, HALO *sub);
void    exciseSubhaloStars          (HALO *halos, long ihost);
int     compare                     (struct particle *, struct particle *, XYZ *);

/***************************************************************************
 *   Global variables                                                      *
 ***************************************************************************/
static int     numHalos;
static HALO    *halos;

/* Arrays */
static double *gridl1dim;

/* conversion factors etc. */
double r_fac, x_fac, v_fac, m_fac, rho_fac, phi_fac, u_fac, Hubble;


/***************************************************************************
 *   Implementation of exported functions                                  *
 ***************************************************************************/
void
ahf_halos(ahf2_patches_t *patches)
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
  int           numSubStruct, *SubStruct, ihost, isub;
  
  halo_s       *ptree2halos;

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
	fprintf(stderr, "\n#################### ahf_halos ####################\n");
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
	fprintf(io.logfile,"\n#################### ahf_halos ####################\n");
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


	/**************************************************************************/
  /* Convert the patch_tree to a halo tree to be stored in halos[]          */
	/**************************************************************************/
#ifdef VERBOSE
	fprintf(stderr, "\nConverting the patch_tree[][] to a halos[] array:\n");
#endif
	fprintf(io.logfile, "\nConverting the patch_tree[][] to a halos[] array:\n");
	fflush(io.logfile);
  timing.patchtree2halos -= time(NULL);
  ptree2halos = patchtree2halos(patches->tree, patches->n_patches);
  halos    = ptree2halos->halos;
  numHalos = ptree2halos->nhalos;
  timing.patchtree2halos += time(NULL);
  
	if (halos == NULL) {
#ifdef VERBOSE
		fprintf(stderr, "Stuffed up converting patch_tree[][] to a halos[]\n");
#endif
		fprintf(io.logfile, "Stuffed up converting patch_tree[][] to a halos[]\n");
    fflush(io.logfile);
		exit(-1);
	}
  
  //======================================================================================
  // remove patchtree[][] from memory
  //======================================================================================
#ifdef VERBOSE
  fprintf(stderr,"Freeing ahf2_halos structure ... ");
#endif
  ahf2_patches_free(&patches);
  patches=NULL;
#ifdef VERBOSE
  fprintf(stderr,"done\n");
#endif
  //fprintf(stderr,"[%s:%d] MEGA-WARNING!! We are NOT freeing ahf2_halos structure!!!\n", __FILE__, __LINE__);
  
  
#ifdef AHF2_read_preliminary_halos
 {
  FILE *fp;
  char infile[MAXSTRING];
  
  // remove all halos[]
  for(i=0; i<numHalos; i++) {
    if(halos[i].subStruct) free(halos[i].subStruct);
  }
  if(halos) free(halos);
  
  
  sprintf(infile,"%s.AHF_preliminaryhalos",fprefix);
  fprintf(stderr,"Reading all halos[] array from %s ... ",infile);
  fp = fopen(infile,"r");
  if(fp == NULL) {
    fprintf(stderr,"could not open %s\n",infile);
    exit(0);
  }
  
  // number of halos (allocating memory for all of them...)
  fscanf(fp,"%d",&numHalos);
  fprintf(stderr,"(numHalos=%d) ",numHalos);
  halos = (HALO *) calloc(numHalos, sizeof(HALO));
  
  // loop over all halos
  for(i=0; i<numHalos; i++) {
    
    fscanf(fp,"%lu %lf %lf %lf %lf %lf %lf %d %d %d %d %d",
           &(halos[i].npart),
           &(halos[i].pos.x),
           &(halos[i].pos.y),
           &(halos[i].pos.z),
           &(halos[i].gatherRad),
           &(halos[i].R_vir),
           &(halos[i].spaRes),
           &(halos[i].refLev),
           &(halos[i].numNodes),
           &(halos[i].hostHalo),
           &(halos[i].hostHaloLevel),
           &(halos[i].numSubStruct));
    
    // none of this is really necessary to set in patchtree2halos() (but we might as well as we have the knowledge...)
    halos[i].R_vir    = -1.0;
    halos[i].spaRes   = -1.0;
    halos[i].refLev   = -1;
    halos[i].numNodes = -1;

    // substructure
    if(halos[i].numSubStruct > 0) {
      halos[i].subStruct = (int *) calloc(halos[i].numSubStruct, sizeof(int));
      if(halos[i].subStruct == NULL) {
        fprintf(stderr,"could not allocate halos[%d].subStruct for halos[%d].numSubStruct=%d\n",i,i,halos[i].numSubStruct);
      }
      
      for(k=0; k<halos[i].numSubStruct; k++) {
        fscanf(fp,"%d",&(halos[i].subStruct[k]));
      }
    }
    else {
      halos[i].subStruct = NULL;
    }
  }
  
  fclose(fp);
  fprintf(stderr,"done\n");
 }
#endif

#ifdef AHF2_write_preliminary_halos
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
    for(k=0; k<halos[i].numSubStruct; k++) {
      /* quick-and-easy access to host and subhalo in halos[] */
      isub  = halos[i].subStruct[k];
      ihost = halos[isub].hostHalo; // this should be identical to i !?
      
#ifdef DEBUG_AHF2_substructuretree
      if(ihost != i) {
        int kk;
        fprintf(stderr,"ihost=%d i=%d isub=%d\n",ihost,i,halos[i].numSubStruct,isub);
        fprintf(stderr,"   numSubStruct=%d of %d = :\n",halos[i].numSubStruct,i);
        for(kk=0; kk<halos[i].numSubStruct; kk++) {
          fprintf(stderr,"       isub=%d\n",halos[i].subStruct[kk]);
          fprintf(stderr,"         ihost=%d\n",halos[halos[i].subStruct[kk]].hostHalo);
        }
      }
      if(isub < 0) {
        fprintf(stderr,"ihost=%d i=%d isub=%d\n",ihost,i,isub);
      }
#endif
      
      /* if subhalo spawned from AHF_HOSTHALOLEVEL (or below) we check for distance (and mass!) */
      if(halos[isub].hostHaloLevel >= AHF_HOSTHALOLEVEL) {
        if(check_subhalo(halos+ihost, halos+isub) == TRUE) {
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
        else {
          /* if it lies outside mark it as field halo */
          halos[isub].hostHalo = -1;
        }
      }
      
      /* if it spawned from above AHF_HOSTHALOLEVEL mark it as field halo */
      else {
        halos[isub].hostHalo = -1;
      }
    }
    
    /* copy the new substructure list over to host halo structure */
    halos[i].numSubStruct = numSubStruct;
    if(numSubStruct>0) {
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
	if (idx) free(idx);


#ifdef VERBOSE
	fprintf(io.logfile,
	        "################## ahf_halos finished ##################\n");
	fflush(io.logfile);
#endif
} /* ahf_halos */

/*
 *
 *===============================================================================
 * END OF MAIN
 *
 *===============================================================================
 */


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

			/* particle velocity in halo rest frame (Hubble correction not
			 *needed for L as r x r = 0) */
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
  Ts    = (halo->prof.Ekin[nbins-1]-halo->prof.Ekin[nbins-2]);
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
	fprintf(io.logfile, "    sort_particles:      npart=%12ld => ", halo->npart);
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
  
}
#endif /* AHFexciseSubhaloStars */

#endif // AHF2
