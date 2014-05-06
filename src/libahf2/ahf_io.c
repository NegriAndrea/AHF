#ifdef AHF2

/***********************************************************************
 *    Includes                                                         *
 ***********************************************************************/
/* Standard includes */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "ahf_io.h"
#include "ahf_halos.h"

/** Okay, we need those to catch all the defines etc. */
#include "../define.h"
#include "../param.h"
#include "../common.h"

/** This is required for the HALO structure */
#include "../tdef.h"

/** Get the needed utility functions */
#include "../libutility/cosmology.h"
#include "../libutility/utility.h"


/***********************************************************************
 *    Gloval variables                                                  *
 ***********************************************************************/
extern double r_fac, x_fac, v_fac, m_fac, rho_fac, phi_fac, u_fac, Hubble;


/***********************************************************************
 *    Definitions of local functions                                   *
 ***********************************************************************/
static void
WriteProfilesLegacy(FILE          *fout,
                    HALO          *halos,
                    unsigned long *idx,
                    int           numHalos);
#ifdef AHFdisks
static void
WriteDisksLegacy(FILE          *fout,
                 HALO          *halos,
                 unsigned long *idx,
                 int           numHalos);
#endif

static void
WriteHalosLegacy(FILE          *fout,
                 HALO          *halos,
                 unsigned long *idx,
                 int           numHalos);

static void
WriteParticlesLegacy(FILE          *fout,
                     HALO          *halos,
                     unsigned long *idx,
                     int           numHalos);

static void
WriteSTARDUSTParticlesLegacy(FILE          *fout,
                             HALO          *halos,
                             unsigned long *idx,
                             int           numHalos);
static void
WriteSTARDUSTexcisedParticlesLegacy(FILE          *fout,
                                    HALO          *halos,
                                    unsigned long *idx,
                                    int           numHalos);

void write_halos_line(FILE *fout, HALO *halos, long unsigned j, long unsigned i, unsigned long *idx, int numHalos);

/***********************************************************************
 *    Implemenation of exported functions                              *
 ***********************************************************************/
extern void
ahf_io_WriteProfiles(const char    *fprefix,
                     HALO          *halos,
                     unsigned long *idx,
                     int           numHalos)
{
#ifdef AHFbinary
  ahf_binwrite_profiles(fprefix, halos, idx, numHalos);
#else /* AHFbinary */
	char filename[MAXSTRING];
	FILE *fout;

	/* Generate the filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_profiles");
#ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#endif

	/* Open file*/
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}
	WriteProfilesLegacy(fout, halos, idx, numHalos);

	/* Clean up */
	fclose(fout);
#endif /* AHFbinary */
  
	return;

} /* ahf_io_WriteProfiles */

extern void
ahf_io_WriteHalos(const char    *fprefix,
                  HALO          *halos,
                  unsigned long *idx,
                  int           numHalos)
{
#ifdef AHFbinary
  ahf_binwrite_halos(fprefix,halos,idx,numHalos);
#else /* AHFbinary */
	char filename[MAXSTRING];
	FILE *fout;

	/* Generate the filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_halos");
#ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#endif

	/* Open file*/
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}
	WriteHalosLegacy(fout, halos, idx, numHalos);

	/* Clean up */
	fclose(fout);
#endif /* AHFbinary */
  
	return;

} /* ahf_io_WriteHalos */

#ifdef AHFdisks
extern void
ahf_io_WriteDisks(const char    *fprefix,
                  HALO          *halos,
                  unsigned long *idx,
                  int           numHalos)
{
	char filename[MAXSTRING];
	FILE *fout;
  
	/* Generate the filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_disks");
#ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#endif
  
	/* Open file*/
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}
	WriteDisksLegacy(fout, halos, idx, numHalos);
  
	/* Clean up */
	fclose(fout);
  
	return;  
} /* ahf_io_WriteDisks */
#endif

extern void
ahf_io_WriteParticles(const char    *fprefix,
                      HALO          *halos,
                      unsigned long *idx,
                      int           numHalos)
{
#ifdef AHFbinary
  ahf_binwrite_particles(fprefix, halos, idx, numHalos);
#else /* AHFbinary */

	char filename[MAXSTRING];
	FILE *fout;

	/* Generate the filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_particles");
#ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#endif

	/* Open file*/
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}
	WriteParticlesLegacy(fout, halos, idx, numHalos);
  
	/* Clean up */
	fclose(fout);
#endif /* AHFbinary */
  
	/* Done */
	return;
} /* ahf_io_WriteParticles */

#ifdef METALHACK
extern void
ahf_io_WriteParticlesSTARDUST(const char    *fprefix,
                              HALO          *halos,
                              unsigned long *idx,
                              int           numHalos)
{
  
  
	char filename[MAXSTRING];
	FILE *fout;
  
	/* Generate the filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_particlesSTARDUST");
#ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#endif
  
	/* Open file*/
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}
  
	WriteSTARDUSTParticlesLegacy(fout, halos, idx, numHalos);
  
	/* Clean up */
	fclose(fout);
    
  
#ifdef AHFexciseSubhaloStars
	/* Generate the filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_particlesSTARDUSTexcised");
#ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#endif
  
	/* Open file*/
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}
  
	WriteSTARDUSTexcisedParticlesLegacy(fout, halos, idx, numHalos);
  
	/* Clean up */
	fclose(fout);  
#endif /* AHFexciseSubhaloStars */
  
	/* Done */
	return;
} /* ahf_io_WriteParticlesSTARDUST */
#endif /* METALHACK */

#ifdef AHFcentrefile
extern void
ahf_io_WriteCenterfile(const char *fprefix,
                       HALO       *halos,
                       int        numHalos)
{
	char filename[MAXSTRING];
	FILE *fout;
	long  i,j;

	/* Generate the filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_centres");
#  ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#  endif

	/* Open file */
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}

#ifdef STEREO2
	/* Write centres */
	for (i = 0; i < numHalos; i++) {
		fprintf(fout, "P %g %g %g 0 1 0 6\n",
		        halos[i].pos.x * simu.boxsize,
		        halos[i].pos.y * simu.boxsize,
		        halos[i].pos.z * simu.boxsize);
	}
#else
	/* Write centres */
	for (i = 0; i < numHalos; i++)
#ifdef AHFcentrefileBASIC
   {
		fprintf(fout, "%16.8f %16.8f %16.8f\n",
		        halos[i].pos.x * simu.boxsize,
		        halos[i].pos.y * simu.boxsize,
		        halos[i].pos.z * simu.boxsize);
   }
#else
   {
    fprintf(fout, "%10ld",			          i);
    fprintf(fout, "\t%10d",               halos[i].hostHalo);
    fprintf(fout, "\t%10d",               halos[i].numSubStruct);
    fprintf(fout, "\t%12.6g",			        halos[i].M_vir * m_fac);
    fprintf(fout, "\t%10ld",		          halos[i].npart);
#ifdef AHFundoPositionShiftAndScale
    fprintf(fout, "\t%16.8f",			        ((halos[i].pos.x/simu.pos_scale)-simu.pos_shift[X]));
    fprintf(fout, "\t%16.8f",			        ((halos[i].pos.y/simu.pos_scale)-simu.pos_shift[Y]));
    fprintf(fout, "\t%16.8f",			        ((halos[i].pos.z/simu.pos_scale)-simu.pos_shift[Z]));
#else
    fprintf(fout, "\t%16.8f",			        halos[i].pos.x * x_fac*1000.);
    fprintf(fout, "\t%16.8f",			        halos[i].pos.y * x_fac*1000.);
    fprintf(fout, "\t%16.8f",			        halos[i].pos.z * x_fac*1000.);
#endif
    fprintf(fout, "\t%8.2f",			        halos[i].vel.x * v_fac);
    fprintf(fout, "\t%8.2f",			        halos[i].vel.y * v_fac);
    fprintf(fout, "\t%8.2f",			        halos[i].vel.z * v_fac);
    fprintf(fout, "\t%10.2f",			        halos[i].R_vir * x_fac * 1000.);
    fprintf(fout, "\t%10.2f",			        halos[i].R_max * x_fac * 1000.);
    fprintf(fout, "\t%10.5f",			        halos[i].r2 * x_fac * 1000.);
    fprintf(fout, "\t%10.5f",			        halos[i].mbp_offset * x_fac * 1000.);
    fprintf(fout, "\t%10.5f",			        halos[i].com_offset * x_fac * 1000.);
    fprintf(fout, "\t%8.2f",			        sqrt(halos[i].V2_max * phi_fac));
    fprintf(fout, "\t%12.6f",			        sqrt(halos[i].v_esc2 * phi_fac));
    fprintf(fout, "\t%8.2f",			        halos[i].sigV * v_fac);
    fprintf(fout, "\t%10.6f",			        halos[i].lambda);
    fprintf(fout, "\t%12.6f",			        halos[i].lambdaE);
    fprintf(fout, "\t%12.4g",			        halos[i].AngMom.x);
    fprintf(fout, "\t%12.4g",			        halos[i].AngMom.y);
    fprintf(fout, "\t%12.4g",			        halos[i].AngMom.z);
    fprintf(fout, "\t%10.6f",			        halos[i].axis.y);
    fprintf(fout, "\t%10.6f",			        halos[i].axis.z);
    fprintf(fout, "\t%10.6f",			        halos[i].E1.x);
    fprintf(fout, "\t%10.6f",			        halos[i].E1.y);
    fprintf(fout, "\t%10.6f",			        halos[i].E1.z);
    fprintf(fout, "\t%10.6f",			        halos[i].E2.x);
    fprintf(fout, "\t%10.6f",			        halos[i].E2.y);
    fprintf(fout, "\t%10.6f",			        halos[i].E2.z);
    fprintf(fout, "\t%10.6f",			        halos[i].E3.x);
    fprintf(fout, "\t%10.6f",			        halos[i].E3.y);
    fprintf(fout, "\t%10.6f",			        halos[i].E3.z);
    fprintf(fout, "\t%8.2f",			        halos[i].ovdens);
    fprintf(fout, "\t%6d",			          halos[i].prof.nbins);
    fprintf(fout, "\t%8.6f",              halos[i].fMhires);
    fprintf(fout, "\t%12.6g",			        halos[i].Ekin * m_fac * pow2(v_fac));
    fprintf(fout, "\t%12.6g",			        halos[i].Epot * m_fac * phi_fac);
    fprintf(fout, "\t%12.6g",			        halos[i].SurfP * m_fac * pow2(v_fac));
    fprintf(fout, "\t%12.6g",			        halos[i].Phi0 * phi_fac);
    fprintf(fout, "\t%12.6g",			        halos[i].cNFW);
#ifdef AHFvmbp
    fprintf(fout, "\t%8.2f",			        halos[i].vel_mbp.x * v_fac);
    fprintf(fout, "\t%8.2f",			        halos[i].vel_mbp.y * v_fac);
    fprintf(fout, "\t%8.2f",			        halos[i].vel_mbp.z * v_fac);
#endif
#  ifdef GAS_PARTICLES
    fprintf(fout, "\t%10ld",			        halos[i].gas_only.npart);
    fprintf(fout, "\t%12.6g",			        halos[i].gas_only.Mass * m_fac);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.lambda);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.lambdaE);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.AngMom.x);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.AngMom.y);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.AngMom.z);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.axis.y);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.axis.z);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E1.x);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E1.y);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E1.z);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E2.x);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E2.y);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E2.z);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E3.x);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E3.y);
    fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E3.z);
    fprintf(fout, "\t%12.6g",			        halos[i].gas_only.Ekin * m_fac * pow2(v_fac));
    fprintf(fout, "\t%12.6g",			        halos[i].gas_only.Epot * m_fac * phi_fac);
    fprintf(fout, "\t%10ld",			        halos[i].stars_only.npart);
    fprintf(fout, "\t%12.6g",			        halos[i].stars_only.Mass * m_fac);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.lambda);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.lambdaE);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.AngMom.x);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.AngMom.y);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.AngMom.z);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.axis.y);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.axis.z);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E1.x);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E1.y);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E1.z);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E2.x);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E2.y);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E2.z);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E3.x);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E3.y);
    fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E3.z);
    fprintf(fout, "\t%12.6g",			        halos[i].stars_only.Ekin * m_fac * pow2(v_fac));
    fprintf(fout, "\t%12.6g",			        halos[i].stars_only.Epot * m_fac * phi_fac);
#  endif
    
#  ifdef METALHACK
    fprintf(fout, "\t%e",			        halos[i].mean_z_gas);
    fprintf(fout, "\t%e",			        halos[i].mean_z_star);
#  endif
    fprintf(fout, "\n");
   }
#endif // AHFcentrefileBASIC
    
	
#endif
  
	/* Clean up */
	fclose(fout);

	/* Done */
	return;
} /* ahf_io_WriteCenterfile */

#endif /* AHFcentrefile */

#if (defined AHFsubstructure && !defined AHFrestart)
extern void
ahf_io_WriteSubstructure(const char    *fprefix,
                         HALO          *halos,
                         unsigned long *idx,
                         int           numHalos)
{
#ifdef AHFbinary
#else /* AHFbinary */
	FILE          *fout;
	char          filename[MAXSTRING];
	unsigned long i;
	int           j, k, numGoodHalos=0;

	/* Generate filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_substructure");
#  ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#  endif

	/* Open file */
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}

	/* Write the substructure file */
	//fprintf(fout, "%d\n", 0); // just a placeholder for the time being
	for (j = 0; j < numHalos; j++)
   {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      /* only write information if there is something to write */
      if(halos[i].npart >= simu.AHF_MINPART && halos[i].numSubStruct>0)
       {
        numGoodHalos++;
        
#ifdef AHFnewHaloIDs
        fprintf(fout, "%"PRIu64" %12d\n", halos[i].haloID, halos[i].numSubStruct);
        for (k = 0; k < halos[i].numSubStruct; k++)
          fprintf(fout, "\t%"PRIu64, halos[halos[i].subStruct[k]].haloID);
#else
        fprintf(fout, "%10d %12d\n", j, halos[i].numSubStruct);
        for (k = 0; k < halos[i].numSubStruct; k++)
          fprintf(fout, "%10d ", idx_inv(idx, numHalos, halos[i].subStruct[k]));
#endif
        
        fprintf(fout, "\n");
       }
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
   }
  
  /* readjust number of halos */
//  rewind(fout);
//  fprintf(fout, "%d\n", numGoodHalos);
//	fseek(fout, 0L, SEEK_END);

	/* Clean up*/
	fclose(fout);
#endif /* AHFbinary */
  
	/* Done */
	return;
} /* ahf_io_WriteSubstructure */

#endif /* Substructure selection */

#ifdef AHFgeom
extern void
ahf_io_WriteHalosGeom(const char    *fprefix,
                      HALO          *halos,
                      unsigned long *idx,
                      int           numHalos)
{
	char          filename[MAXSTRING];
	FILE          *fout;
	unsigned long ifrac;
	unsigned long ipart;
	unsigned long j;
	partptr       cur_part, current;
	int           i;
	double        r, g, b;


	/* Generate filename */
	strcpy(filename, fprefix);
	strcat(filename, ".AHF_halos.geom");
#  ifdef VERBOSE
	fprintf(stderr, "%s\n", filename);
#  endif

	/* Open file */
	if ((fout = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "could not open %s\n", filename);
		exit(1);
	}

#  ifdef AHFgeom_SIMUPART
	/* 1. write about 500000 simulation particles to AHF_geom file */
#    ifdef WITH_MPI
	if (global_info.no_part > 500000)
		ifrac = (long unsigned)((double)global_info.no_part
		                        / ((double)(500000.)));
#    else
	if (global.no_part > 500000)
		ifrac = (long unsigned)((double)global.no_part / (double)(500000.));
#    endif
	else
		ifrac = (long unsigned)1;

#    if (defined WITH_MPI || defined AHFrestart)
	for (ipart = 0; ipart < global_info.no_part; ipart += ifrac)
#    else
	for (ipart = 0; ipart < global.no_part; ipart += ifrac)
#    endif
	{
		cur_part = global.fst_part + ipart;
#    ifdef MULTIMASS
		if (cur_part->weight < simu.max_weight)
			fprintf(fout,
			        "p %12.6g %12.6g %12.6g 0 0 1\n",
			        cur_part->pos[X] * x_fac,
			        cur_part->pos[Y] * x_fac,
			        cur_part->pos[Z] * x_fac);
		else
			fprintf(fout,
			        "p %12.6g %12.6g %12.6g 0.3 0 0\n",
			        cur_part->pos[X] * x_fac,
			        cur_part->pos[Y] * x_fac,
			        cur_part->pos[Z] * x_fac);
#    else
		fprintf(fout,
		        "p %12.6g %12.6g %12.6g 0 0 1\n",
		        cur_part->pos[X] * x_fac,
		        cur_part->pos[Y] * x_fac,
		        cur_part->pos[Z] * x_fac);
#    endif
	}
#  endif     /* AHFgeom_SIMUPART */

	/* 2. write halos */
	for (j = 0; j < numHalos; j++) {
		i = idx[j];

		r = 1.0 - (double)i / (double)numHalos;
		g = (double)i / (double)numHalos;
		b = 0.0;

		r = 1.0;
		g = 0.0;
		b = 0.0;

		if ((halos[i].npart >= simu.AHF_MINPART)) {
#  ifdef AHFgeom_HALOPART
			current = halos[i].ll;
			while (current != NULL) {
				fprintf(fout,
				        "p %12.6g %12.6g %12.6g  %g %g %g\n",
				        current->pos[X] * x_fac,
				        current->pos[Y] * x_fac,
				        current->pos[Z] * x_fac,
				        r,
				        g,
				        b);
				current = current->ll;
			}
#  else
			fprintf(fout,
			        "s %12.6g %12.6g %12.6g %12.6g  %g %g %g\n",
			        halos[i].pos.x * x_fac,
			        halos[i].pos.y * x_fac,
			        halos[i].pos.z * x_fac,
			        halos[i].R_vir * x_fac,
			        r,
			        g,
			        b);
#  endif
		}
	}

	/* Clean */
	fclose(fout);

	/* Done */
	return;
} /* ahf_io_WriteHalosGeom */

#endif /* AHFgeom */


/***********************************************************************
 *    Implementation of local functions                                *
 ***********************************************************************/

/* This is used to print the helpful notes about which columns holds what */
#  define COLUMN_INFO(f, info, counter) fprintf(f, info "(%i)\t", counter++);

static void
WriteProfilesLegacy(FILE          *fout,
                    HALO          *halos,
                    unsigned long *idx,
                    int           numHalos)
{
	int           column = 1;
	unsigned long i, j;
	int           r_conv_i, ibin;
	double        rad, t_relax, t0;
  
#  if (defined WITH_MPI || defined AHFrestart)
#    ifdef WITH_MPI
	if (global_mpi.rank == 0)
#    else
  if (global_info.rank == 0)
#    endif
#  endif
     {
      fprintf(fout, "#");
      COLUMN_INFO(fout, "r", column);
      COLUMN_INFO(fout, "npart", column);
      COLUMN_INFO(fout, "M_in_r", column);
      COLUMN_INFO(fout, "ovdens", column);
      COLUMN_INFO(fout, "dens", column);
      COLUMN_INFO(fout, "vcirc", column);
      COLUMN_INFO(fout, "vesc", column);
      COLUMN_INFO(fout, "sigv", column);
      COLUMN_INFO(fout, "Lx", column);
      COLUMN_INFO(fout, "Ly", column);
      COLUMN_INFO(fout, "Lz", column);
      COLUMN_INFO(fout, "b", column);
      COLUMN_INFO(fout, "c", column);
      COLUMN_INFO(fout, "Eax", column);
      COLUMN_INFO(fout, "Eay", column);
      COLUMN_INFO(fout, "Eaz", column);
      COLUMN_INFO(fout, "Ebx", column);
      COLUMN_INFO(fout, "Eby", column);
      COLUMN_INFO(fout, "Ebz", column);
      COLUMN_INFO(fout, "Ecx", column);
      COLUMN_INFO(fout, "Ecy", column);
      COLUMN_INFO(fout, "Ecz", column);
      COLUMN_INFO(fout, "Ekin", column);
      COLUMN_INFO(fout, "Epot", column);
#  ifdef GAS_PARTICLES
      COLUMN_INFO(fout, "M_gas", column);
      COLUMN_INFO(fout, "M_star", column);
      COLUMN_INFO(fout, "u_gas", column);
#  endif
#  ifdef AHFphspdens
      COLUMN_INFO(fout, "sigma2_vx_sh", column);
      COLUMN_INFO(fout, "sigma2_vy_sh", column);
      COLUMN_INFO(fout, "sigma2_vz_sh", column);
      COLUMN_INFO(fout, "sigma2_vr_sh", column);
      COLUMN_INFO(fout, "sigma2_vtheta_sh", column);
      COLUMN_INFO(fout, "sigma2_vphi_sh", column);
#    ifdef AHFmeanvelocities
      COLUMN_INFO(fout, "mean_vx_sh", column);
      COLUMN_INFO(fout, "mean_vy_sh", column);
      COLUMN_INFO(fout, "mean_vx_sh", column);
      COLUMN_INFO(fout, "mean_vr_sh", column);
      COLUMN_INFO(fout, "mean_vtheta_sh", column);
      COLUMN_INFO(fout, "mean_vphi_sh", column);
      COLUMN_INFO(fout, "mean_vx_sp", column);
      COLUMN_INFO(fout, "mean_vy_sp", column);
      COLUMN_INFO(fout, "mean_vx_sp", column);
      COLUMN_INFO(fout, "mean_vr_sp", column);
      COLUMN_INFO(fout, "mean_vtheta_sp", column);
      COLUMN_INFO(fout, "mean_vphi_sp", column);
#    endif
#  endif
#  ifdef METALHACK
      COLUMN_INFO(fout, "Z_gas_sh", column);
      COLUMN_INFO(fout, "Z_star_sh", column);
#  endif
      fprintf(fout, "\n");
     }
  
	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART)) {
        for (r_conv_i = 0, ibin = 0; ibin < halos[i].prof.nbins; ibin++) {
          rad     = halos[i].prof.r[ibin];
          
          /* check for converged radius (Power et al. 2003) */
          //t_relax = halos[i].prof.npart[ibin] / log(rad / halos[i].spaRes) * rad / sqrt(halos[i].prof.v2_circ[ibin]);           // Eq.(4)
          //t0      = 0.9*calc_t(global.a) * simu.t_unit * Mpc / Gyr; /* age of the Universe in Gyr/h */
          
          t_relax = halos[i].prof.npart[ibin] / log(halos[i].prof.npart[ibin]) / 8. * rad / sqrt(halos[i].prof.v2_circ[ibin]);  // Eq.(20)
          t0      = 0.6*halos[i].prof.r[halos[i].prof.nbins-1] / sqrt(halos[i].prof.v2_circ[halos[i].prof.nbins-1]) * r_fac/sqrt(phi_fac) * 1E3;
          
          
          /* convert to (Mpc/h) / (km/sec) ~ (Mpc/h) / (kpc/Gyr) */
          t_relax *= r_fac / sqrt(phi_fac);
          
          /* convert to Gyr/h */
          t_relax *= 1E3;
          
          
          /* if not converged, write negative radius into *.AHF_profiles */
          if (t_relax < t0)
           /* The +1 is merely to keep the name consistent:
            * r_conv_i should give the smallest converged radius,
            * not the bin before it. Hence to check for converged
            * bin or not, i<r_conv_i is the thing to do. */
            r_conv_i = ibin + 1;
        }
        
        for (ibin = 0; ibin < halos[i].prof.nbins; ibin++) {
          /* If the radius is not converged, print the radius negative to the .AHF_profiles */
          rad = (ibin < r_conv_i) ? -halos[i].prof.r[ibin] : halos[i].prof.r[ibin];
          fprintf(fout, "%8.4f",				      rad * x_fac * 1000.);
          fprintf(fout, "\t%10ld",			      halos[i].prof.npart[ibin]);
          fprintf(fout, "\t%e",				        halos[i].prof.nvpart[ibin] * m_fac);
          fprintf(fout, "\t%12.2f",			      halos[i].prof.ovdens[ibin] * rho_fac / global.rho_vir);
          fprintf(fout, "\t%12.2f",			      halos[i].prof.dens[ibin] * rho_fac / global.rho_vir);
          fprintf(fout, "\t%8.2f",			      sqrt(halos[i].prof.v2_circ[ibin] * phi_fac));
          fprintf(fout, "\t%10.6f",		        sqrt(halos[i].prof.v_esc2[ibin] * phi_fac));
          fprintf(fout, "\t%8.2f",			      halos[i].prof.sig_v[ibin] * v_fac);
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Lx[ibin] * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Ly[ibin] * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Lz[ibin] * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.axis2[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.axis3[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1x[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1y[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1z[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2x[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2y[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2z[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3x[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3y[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3z[ibin]);
          fprintf(fout, "\t%e",				        halos[i].prof.Ekin[ibin] * m_fac * pow2(v_fac));
          fprintf(fout, "\t%e",				        halos[i].prof.Epot[ibin] * m_fac * phi_fac);
#  ifdef GAS_PARTICLES
          fprintf(fout, "\t%e",				        halos[i].prof.M_gas[ibin] * m_fac);
          fprintf(fout, "\t%e",				        halos[i].prof.M_star[ibin] * m_fac);
          fprintf(fout, "\t%e",				        halos[i].prof.u_gas[ibin] * u_fac);  // remember, in HaloProfiles() we converted cur_part->u to code internal units, here we need to undo that!
#  endif
#  ifdef AHFphspdens
          fprintf(fout, "\t%g",				        halos[i].prof.sigma2_vx_sh[ibin] * v_fac * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.sigma2_vy_sh[ibin] * v_fac * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.sigma2_vz_sh[ibin] * v_fac * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.sigma2_vr_sh[ibin] * v_fac * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.sigma2_vtheta_sh[ibin] * v_fac * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.sigma2_vphi_sh[ibin] * v_fac * v_fac);
#    ifdef AHFmeanvelocities
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vx_sh[ibin] * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vy_sh[ibin] * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vz_sh[ibin] * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vr_sh[ibin] * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vtheta_sh[ibin]);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vphi_sh[ibin]);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vx_sp[ibin] * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vy_sp[ibin] * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vz_sp[ibin] * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vr_sp[ibin] * v_fac);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vtheta_sp[ibin]);
          fprintf(fout, "\t%g",				        halos[i].prof.mean_vphi_sp[ibin]);
#    endif
#  endif
#  ifdef METALHACK
          fprintf(fout, "\t%e",				        halos[i].prof.z_gas[ibin]);
          fprintf(fout, "\t%e",				        halos[i].prof.z_star[ibin]);
#  endif
          fprintf(fout, "\n");
        }
      }
      
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
	}
  
	/* Done */
	return;
} /* WriteProfilesLegacy */

static void
WriteHalosLegacy(FILE          *fout,
                 HALO          *halos,
                 unsigned long *idx,
                 int           numHalos)
{
	int           column = 1;
	unsigned long i, j;


#  if (defined WITH_MPI || defined AHFrestart)
#    ifdef WITH_MPI
	if (global_mpi.rank == 0)
#    else
	if (global_info.rank == 0)
#    endif
#  endif
   {
		fprintf(fout, "#");
    COLUMN_INFO(fout, "ID", column);
    COLUMN_INFO(fout, "hostHalo", column);
    COLUMN_INFO(fout, "numSubStruct", column);
		COLUMN_INFO(fout, "Mvir", column);
 		COLUMN_INFO(fout, "npart", column);
		COLUMN_INFO(fout, "Xc", column);
		COLUMN_INFO(fout, "Yc", column);
		COLUMN_INFO(fout, "Zc", column);
		COLUMN_INFO(fout, "VXc", column);
		COLUMN_INFO(fout, "VYc", column);
		COLUMN_INFO(fout, "VZc", column);
		COLUMN_INFO(fout, "Rvir", column);
		COLUMN_INFO(fout, "Rmax", column);
		COLUMN_INFO(fout, "r2", column);
		COLUMN_INFO(fout, "mbp_offset", column);
		COLUMN_INFO(fout, "com_offset", column);
		COLUMN_INFO(fout, "Vmax", column);
		COLUMN_INFO(fout, "v_esc", column);
		COLUMN_INFO(fout, "sigV", column);
		COLUMN_INFO(fout, "lambda", column);
		COLUMN_INFO(fout, "lambdaE", column);
		COLUMN_INFO(fout, "Lx", column);
		COLUMN_INFO(fout, "Ly", column);
		COLUMN_INFO(fout, "Lz", column);
    COLUMN_INFO(fout, "b", column);
    COLUMN_INFO(fout, "c", column);
		COLUMN_INFO(fout, "Eax", column);
		COLUMN_INFO(fout, "Eay", column);
		COLUMN_INFO(fout, "Eaz", column);
		COLUMN_INFO(fout, "Ebx", column);
		COLUMN_INFO(fout, "Eby", column);
		COLUMN_INFO(fout, "Ebz", column);
		COLUMN_INFO(fout, "Ecx", column);
		COLUMN_INFO(fout, "Ecy", column);
		COLUMN_INFO(fout, "Ecz", column);
		COLUMN_INFO(fout, "ovdens", column);
		COLUMN_INFO(fout, "nbins", column);
		COLUMN_INFO(fout, "fMhires", column);
		COLUMN_INFO(fout, "Ekin", column);
		COLUMN_INFO(fout, "Epot", column);
		COLUMN_INFO(fout, "SurfP", column);
		COLUMN_INFO(fout, "Phi0", column);
    COLUMN_INFO(fout, "cNFW", column);
#ifdef AHFvmbp
    COLUMN_INFO(fout, "mbp_Vx", column);
    COLUMN_INFO(fout, "mbp_Vy", column);
    COLUMN_INFO(fout, "mbp_Vz", column);
#endif
#  ifdef GAS_PARTICLES
		COLUMN_INFO(fout, "n_gas", column);
		COLUMN_INFO(fout, "M_gas", column);
		COLUMN_INFO(fout, "lambda_gas", column);
		COLUMN_INFO(fout, "lambdaE_gas", column);
		COLUMN_INFO(fout, "Lx_gas", column);
		COLUMN_INFO(fout, "Ly_gas", column);
		COLUMN_INFO(fout, "Lz_gas", column);
    COLUMN_INFO(fout, "b_gas", column);
    COLUMN_INFO(fout, "c_gas", column);
		COLUMN_INFO(fout, "Eax_gas", column);
		COLUMN_INFO(fout, "Eay_gas", column);
		COLUMN_INFO(fout, "Eaz_gas", column);
		COLUMN_INFO(fout, "Ebx_gas", column);
		COLUMN_INFO(fout, "Eby_gas", column);
		COLUMN_INFO(fout, "Ebz_gas", column);
		COLUMN_INFO(fout, "Ecx_gas", column);
		COLUMN_INFO(fout, "Ecy_gas", column);
		COLUMN_INFO(fout, "Ecz_gas", column);
		COLUMN_INFO(fout, "Ekin_gas", column);
		COLUMN_INFO(fout, "Epot_gas", column);
		COLUMN_INFO(fout, "n_star", column);
		COLUMN_INFO(fout, "M_star", column);
		COLUMN_INFO(fout, "lambda_star", column);
		COLUMN_INFO(fout, "lambdaE_star", column);
		COLUMN_INFO(fout, "Lx_star", column);
		COLUMN_INFO(fout, "Ly_star", column);
		COLUMN_INFO(fout, "Lz_star", column);
    COLUMN_INFO(fout, "b_star", column);
    COLUMN_INFO(fout, "c_star", column);
		COLUMN_INFO(fout, "Eax_star", column);
		COLUMN_INFO(fout, "Eay_star", column);
		COLUMN_INFO(fout, "Eaz_star", column);
		COLUMN_INFO(fout, "Ebx_star", column);
		COLUMN_INFO(fout, "Eby_star", column);
		COLUMN_INFO(fout, "Ebz_star", column);
		COLUMN_INFO(fout, "Ecx_star", column);
		COLUMN_INFO(fout, "Ecy_star", column);
		COLUMN_INFO(fout, "Ecz_star", column);
		COLUMN_INFO(fout, "Ekin_star", column);
		COLUMN_INFO(fout, "Epot_star", column);
#    ifdef METALHACK
		COLUMN_INFO(fout, "mean_z_gas", column);
		COLUMN_INFO(fout, "mean_z_star", column);
#    endif
#  endif /* GAS_PARTICLES */
		fprintf(fout, "\n");
   }

	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART)) {
        write_halos_line(fout, halos, j, i, idx, numHalos);
      }
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
	}

	/* Done */
	return;
} /* WriteHalosLegacy */

#ifdef AHFdisks
static void
WriteDisksLegacy(FILE          *fout,
                 HALO          *halos,
                 unsigned long *idx,
                 int           numHalos)
{
	int           column = 1;
	unsigned long i, j;
	int           r_conv_i, ibin;
	double        rad, t_relax, t0;
  
#  if (defined WITH_MPI || defined AHFrestart)
#    ifdef WITH_MPI
	if (global_mpi.rank == 0)
#    else
    if (global_info.rank == 0)
#    endif
#  endif
     {
      fprintf(fout, "#");
      COLUMN_INFO(fout, "r", column);
      COLUMN_INFO(fout, "M_in_r", column);
      COLUMN_INFO(fout, "M_gas", column);
      COLUMN_INFO(fout, "Ekin_gas", column);
      COLUMN_INFO(fout, "k_gas", column);
      COLUMN_INFO(fout, "Lx_gas", column);
      COLUMN_INFO(fout, "Ly_gas", column);
      COLUMN_INFO(fout, "Lz_gas", column);
      COLUMN_INFO(fout, "b_gas", column);
      COLUMN_INFO(fout, "c_gas", column);
      COLUMN_INFO(fout, "Eax_gas", column);
      COLUMN_INFO(fout, "Eay_gas", column);
      COLUMN_INFO(fout, "Eaz_gas", column);
      COLUMN_INFO(fout, "Ebx_gas", column);
      COLUMN_INFO(fout, "Eby_gas", column);
      COLUMN_INFO(fout, "Ebz_gas", column);
      COLUMN_INFO(fout, "Ecx_gas", column);
      COLUMN_INFO(fout, "Ecy_gas", column);
      COLUMN_INFO(fout, "Ecz_gas", column);
      COLUMN_INFO(fout, "M_star", column);
      COLUMN_INFO(fout, "Ekin_star", column);
      COLUMN_INFO(fout, "k_star", column);
      COLUMN_INFO(fout, "Lx_star", column);
      COLUMN_INFO(fout, "Ly_star", column);
      COLUMN_INFO(fout, "Lz_star", column);
      COLUMN_INFO(fout, "b_star", column);
      COLUMN_INFO(fout, "c_star", column);
      COLUMN_INFO(fout, "Eax_star", column);
      COLUMN_INFO(fout, "Eay_star", column);
      COLUMN_INFO(fout, "Eaz_star", column);
      COLUMN_INFO(fout, "Ebx_star", column);
      COLUMN_INFO(fout, "Eby_star", column);
      COLUMN_INFO(fout, "Ebz_star", column);
      COLUMN_INFO(fout, "Ecx_star", column);
      COLUMN_INFO(fout, "Ecy_star", column);
      COLUMN_INFO(fout, "Ecz_star", column);
      fprintf(fout, "\n");
     }
  
	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART) && (halos[i].stars_only.npart >= AHF_MINPART_STARS)) {
        
        /* write the _halos line */
        write_halos_line(fout, halos, j, i, idx, numHalos);
        
        /* write the _disks profile */
        for (ibin = 0; ibin < halos[i].prof.nbins; ibin++) {
          fprintf(fout, "%8.4f",				 fabs(halos[i].prof.r[ibin])         * x_fac * 1000.);
          fprintf(fout, "\t%e",				        halos[i].prof.nvpart[ibin]     * m_fac);
          fprintf(fout, "\t%e",				        halos[i].prof.M_gas[ibin]      * m_fac);
          fprintf(fout, "\t%e",				        halos[i].prof.Ekin_gas[ibin]   * m_fac * pow2(v_fac));
          fprintf(fout, "\t%e",				        halos[i].prof.k_gas[ibin]      * m_fac * pow2(v_fac));
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Lx_gas[ibin]     * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Ly_gas[ibin]     * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Lz_gas[ibin]     * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.axis2_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.axis3_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1x_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1y_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1z_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2x_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2y_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2z_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3x_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3y_gas[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3z_gas[ibin]);
          
          fprintf(fout, "\t%e",				        halos[i].prof.M_star[ibin]      * m_fac);
          fprintf(fout, "\t%e",				        halos[i].prof.Ekin_star[ibin]   * m_fac * pow2(v_fac));
          fprintf(fout, "\t%e",				        halos[i].prof.k_star[ibin]      * m_fac * pow2(v_fac));
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Lx_star[ibin]     * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Ly_star[ibin]     * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%12.4g",			      halos[i].prof.Lz_star[ibin]     * m_fac * r_fac * v_fac);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.axis2_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.axis3_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1x_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1y_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E1z_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2x_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2y_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E2z_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3x_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3y_star[ibin]);
          fprintf(fout, "\t%10.6f",			      halos[i].prof.E3z_star[ibin]);
          fprintf(fout, "\n");
        }
      }
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
  }
  
  /* Done */
  return;
} /* WriteDisksLegacy */
#endif /* AHFdisks */

static void
WriteParticlesLegacy(FILE          *fout,
                     HALO          *halos,
                     unsigned long *idx,
                     int           numHalos)
{
	unsigned long i, j;
	unsigned long numGoodHalos = 0;
	int           ipart;
	partptr       cur_part, current;

	/* This is the number of halos, currently 0, will be filled correctly at the end */
	fprintf(fout, "%15lu\n", numGoodHalos);

	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART)) {
        /* Count this halo as a good one */
        numGoodHalos++;
        
        /* start with dumping the actual number of particle IDs to follow alongside the actual haloID... */
        //fprintf(fout, "%lu\n", (unsigned long)(halos[i].npart));
#ifdef AHFnewHaloIDs
        fprintf(fout, "%lu %22"PRIu64"\n", (unsigned long)(halos[i].npart),halos[i].haloID);
#else
#ifdef SUSSING2013
        fprintf(fout, "%lu %10"PRIu64"\n",		 (unsigned long)(halos[i].npart), getSussing2013ID(simu.isnap,j));
#else
        fprintf(fout, "%lu %10ld\n",		 (unsigned long)(halos[i].npart), (long)j);
#endif
#endif
        
        for (ipart = 0; ipart < halos[i].npart; ipart++) {
          current = global.fst_part + halos[i].ipart[ipart];
        
#ifdef SUSSING2013
          fprintf(fout, "%" PRIpartid " %g %lf %lf %lf %lf %lf %lf\n",
                  current->id,
                  current->E,
                  current->pos[X]*x_fac*1000.,
                  current->pos[Y]*x_fac*1000.,
                  current->pos[Z]*x_fac*1000.,
                  current->mom[X]*v_fac,
                  current->mom[Y]*v_fac,
                  current->mom[Z]*v_fac
                  );
#else /* SUSSING2013 */

#if (defined GAS_PARTICLES)
          fprintf(fout, "%" PRIpartid "\t%d\n", current->id, (current->u >= PGAS) ? (int)PGAS : (int)(-current->u));
#else // GAS_PARTICLES
#		if (!(defined AHF_NO_PARTICLES && defined AHFlean))
          fprintf(fout, "%" PRIpartid "\t%d\n", current->id,(int)(-PDM));
#		endif
#endif     /* GAS_PARTICLES */

#endif /* SUSSING2013 */
          
        }
        
      }
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
	}

	/* Rewriting the number of halos */
	rewind(fout);
	fprintf(fout, "%15lu", numGoodHalos);
	fseek(fout, 0L, SEEK_END);

	/* Done */
	return;
} /* WriteParticlesLegacy */

#ifdef METALHACK
static void
WriteSTARDUSTParticlesLegacy(FILE          *fout,
                             HALO          *halos,
                             unsigned long *idx,
                             int           numHalos)
{
	unsigned long i, j, nstars;
	unsigned long numGoodHalos = 0;
	int           ipart;
	partptr       cur_part, current;
  double        Tnow, Tformation;
  
	/* This is the number of halos, currently 0, will be filled
	 * correctly at the end */
	fprintf(fout, "%15lu\n", numGoodHalos);
  
	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART)) {
        /* Count this halo as a good one */
        numGoodHalos++;
        
        /* start with dumping the actual number of particle IDs ... */
        fprintf(fout, "%lu  ", (unsigned long)(halos[i].stars_only.npart));
        
        /* ...and the haloid */
#ifdef AHFnewHaloIDs
        fprintf(fout, "%22"PRIu64, halos[i].haloID);
#else
        fprintf(fout, "%10ld\n",	 (long)j);
#endif

        /* double checking the number of stars */
        nstars = 0;
        
        /* THIS IS THE STANDARD FORMAT */
        for (ipart = 0; ipart < halos[i].npart; ipart++) {
          current = global.fst_part + halos[i].ipart[ipart];
          
          /* filter star particles */
          if(fabs(current->u - PSTAR) < ZERO)
           {
            nstars++;
            Tnow       = calc_t(global.a);
            Tformation = calc_t(current->age);
            fprintf(fout, "%ld\t%g\t%g\t%g\n",
                    (long int)current->id,
                    current->weight   * m_fac,
                    (Tnow-Tformation) * simu.t_unit * Mpc / Gyr,
                    current->z);
           }
        }
        
        if (nstars != halos[i].stars_only.npart)
         {
          fprintf(stderr,"WriteSTARDUSTParticlesLegacy: mismatch of number of star particles in halo[%ld]   stars_only.npart=%ld, nstars=%ld\n",
                  i,halos[i].stars_only.npart,nstars);
          exit(EXIT_FAILURE);
         }
      }
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
	}
  
	/* Rewriting the number of halos */
	rewind(fout);
	fprintf(fout, "%15lu", numGoodHalos);
	fseek(fout, 0L, SEEK_END);
  
	/* Done */
	return;
} /* WriteSTARDUSTParticlesLegacy */

#ifdef AHFexciseSubhaloStars
static void
WriteSTARDUSTexcisedParticlesLegacy(FILE          *fout,
                                    HALO          *halos,
                                    unsigned long *idx,
                                    int           numHalos)
{
	unsigned long i, j, nstars;
	unsigned long numGoodHalos = 0;
	int           ipart;
	partptr       cur_part, current;
  double        Tnow, Tformation;
  
	/* This is the number of halos, currently 0, will be filled correctly at the end */
	fprintf(fout, "%15lu\n", numGoodHalos);
  
	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART)) {
        /* Count this halo as a good one */
        numGoodHalos++;
        
        /* start with dumping the actual number of particle IDs ... */
        fprintf(fout, "%lu  ", (unsigned long)(halos[i].npart_uniquestars));
        
        /* ...and the haloid */
#ifdef AHFnewHaloIDs
        fprintf(fout, "%22"PRIu64, halos[i].haloID);
#else
        fprintf(fout, "%10ld\n",	 (long)j);
#endif

        /* double checking the number of stars */
        nstars = 0;
        
        /* THIS IS THE STANDARD FORMAT */
        for (ipart = 0; ipart < halos[i].npart_uniquestars; ipart++) {
          current = global.fst_part + halos[i].ipart_uniquestars[ipart];
          
          /* filter star particles */
          if(fabs(current->u - PSTAR) < ZERO)
           {
            nstars++;
            Tnow       = calc_t(global.a);
            Tformation = calc_t(current->age);
            fprintf(fout, "%ld\t%g\t%g\t%g\n",
                    (long int)current->id,
                    current->weight   * m_fac,
                    (Tnow-Tformation) * simu.t_unit * Mpc / Gyr,
                    current->z);
           }
        }
        
        if (nstars != halos[i].npart_uniquestars)
         {
          fprintf(stderr,"WriteSTARDUSTexcisedParticlesLegacy: mismatch of number of star particles in halo[%ld]   npart_uniquestars=%ld, nstars=%ld\n",
                  i,halos[i].npart_uniquestars,nstars);
          exit(EXIT_FAILURE);
         }
      }
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
	}
  
	/* Rewriting the number of halos */
	rewind(fout);
	fprintf(fout, "%15lu", numGoodHalos);
	fseek(fout, 0L, SEEK_END);
  
	/* Done */
	return;
} /* WriteSTARDUSTexcisedParticlesLegacy */
#endif /* AHFexciseSubhaloStars */ 

#endif /* METALHACK */


#ifdef AHFbinary

#  define LEN_TYPEFIELD 128
#  define WRITE_I       fwrite(&tmp_int, sizeof(bin_int_t), 1, f)
#  define WRITE_F       fwrite(&tmp_float, sizeof(bin_float_t), 1, f)
#  define WRITE_LU      fwrite(&tmp_ulong, sizeof(bin_ulong_t), 1, f) 

typedef uint32_t bin_int_t;
typedef uint64_t bin_ulong_t;
typedef float    bin_float_t;


void
ahf_binwrite_open_files(FILE **f,
                     FILE **f_info,
                     char *prefix,
                     char *suffix)
{
	char fname[2048];
	char fname_info[2048];
  
	/* Generate the filenames */
	sprintf(fname, "%s.%s_bin", prefix, suffix);
	sprintf(fname_info, "%s.%s_info", prefix, suffix);
  
	/* Open binary output file */
	*f = fopen(fname, "wb");
	if (*f == NULL) {
		fprintf(stderr, "could not open %s with mode `wb'\n", fname);
		exit(EXIT_FAILURE);
	}
  
	/* Open the info file (in MPI mode, only one does that) */
	*f_info = NULL;
#  if (defined WITH_MPI || defined AHFrestart)
#    ifdef WITH_MPI
	if (global_mpi.rank == 0)
#    else
    if (global_info.rank == 0)
#    endif
#  endif
     {
      *f_info = fopen(fname_info, "w");
      if (*f_info == NULL) {
        fprintf(stderr, "could not open %s with mode `w'\n",
                fname_info);
        exit(EXIT_FAILURE);
      }
     }
  
#  if (defined VERBOSE)
	fprintf(stderr, "%s\t", fname);
#  endif
  
	/* Done */
	return;
} /* ahf_binwrite_open_files */

void
ahf_binwrite_profiles(char          *prefix,
                   HALO          *halos,
                   unsigned long *idx,
                   int           numHalos)
{
	FILE        *f;
	FILE        *f_info;
	int         i = 0, j = 0;
	int         num_columns;
	uint32_t    sizes[2];
	bin_int_t   tmp_int;
	bin_float_t tmp_float;
	int         r_conv_i, ibin;
	double      age, rad, t_relax, t0;
	uint64_t    real_num_halos;
	uint64_t    total_num_lines;
	uint64_t    tmp;
	int8_t      typefield[LEN_TYPEFIELD];
  int32_t     one=1;
  
	/* Store the number of bytes used for integer and float values */
	sizes[0] = sizeof(bin_int_t);
	sizes[1] = sizeof(bin_float_t);
  
	/* Initialize typefield to float */
	for (i = 0; i < LEN_TYPEFIELD; i++)
		typefield[i] = INT8_C(1);
  
	/* Open the files */
	ahf_binwrite_open_files(&f, &f_info, prefix, "AHF_profiles");
  
	/* Write info, if required */
	i = 0;
	if (f_info != NULL)
		fprintf(f_info,
            "one (used to check for byteswap)\n"
            "numHalos\n"
            "numColumns\n"
		        "r(%i) npart(%i) M_in_r(%i) ovdens(%i) dens(%i) "
		        "vcirc(%i) vesc(%i) sigv(%i) Lx(%i) Ly(%i) Lz(%i) "
		        "b(%i) c(%i) Eax(%i) Eay(%i) Eaz(%i) "
		        "Ebx(%i) Eby(%i) Ebz(%i) "
		        "Ecx(%i) Ecy(%i) Ecz(%i) "
		        "Ekin(%i) Epot(%i) ",
		        i + 1, i + 2, i + 3, i + 4, i + 5,
		        i + 6, i + 7, i + 8, i + 9, i + 10,
		        i + 11, i + 12, i + 13, i + 14,
		        i + 15, i + 16, i + 17, i + 18,
		        i + 19, i + 20, i + 21, i + 22,
            i + 23, i + 24);
	i += 24;
	typefield[1] = INT8_C(0);
#  if (defined GAS_PARTICLES)
	if (f_info != NULL)
		fprintf(f_info,
		        "M_gas(%i) M_star(%i) u_gas(%i) ",
		        i + 1, i + 2, i + 3);
	i += 3;
#ifdef METALHACK
	if (f_info != NULL)
		fprintf(f_info,
		        "Z_gas_sh(%i) Z_star_sh(%i)",
		        i + 1, i + 2);
	i += 2;
#endif
  if (f_info != NULL)
    fprintf(f_info,"\n");
#  endif
#  if (defined AHFphspdens)
	if (f_info != NULL)
		fprintf(
            f_info,
            "sigma2_vx_sh(%i) sigma2_vy_sh(%i) sigma2_vz_sh(%i) "
            "sigma2_vr_sh(%i) sigma2_vtheta_sh(%i) sigma2_vphi_sh(%i) ",
            i + 1,
            i + 2,
            i + 3,
            i + 4,
            i + 5,
            i + 6);
	i += 6;
#    if (defined AHFmeanvelocities)
	if (f_info != NULL)
		fprintf(f_info,
		        "mean_vx_sh(%i) mean_vy_sh(%i) mean_vz_sh(%i) "
		        "mean_vr_sh(%i) mean_vtheta_sh(%i) mean_vphi_sh(%i) "
		        "mean_vx_sp(%i) mean_vy_sp(%i) mean_vz_sp(%i)\n"
		        "mean_vr_sp(%i) mean_vtheta_sp(%i) mean_vphi_sp(%i) ",
		        i + 1, i + 2, i + 3,
		        i + 4, i + 5, i + 6,
		        i + 7, i + 8, i + 9,
		        i + 10, i + 11, i + 12);
	i += 12;
#    endif
#  endif
	if (f_info != NULL)
		fprintf(f_info,"\n");
	num_columns = i;
  
  /* Just write the number "1" into the file */
  fwrite(&one, sizeof(int32_t), 1, f);

	/* Store the sizes */
	//fwrite(sizes, sizeof(uint32_t), 2, f);
  
	/* Write a dummy value for the number of halos and total number of lines */
	real_num_halos  = UINT64_C(0);
	total_num_lines = UINT64_C(0);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
	//fwrite(&total_num_lines, sizeof(uint64_t), 1, f);
  
	/* Store the number of columns */
	fwrite(&num_columns, sizeof(uint32_t), 1, f);
  
	/* Store the typefield */
	//fwrite(typefield, sizeof(int8_t), num_columns, f);
  
	/* Now loop over all haloes for writing */
	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART)) {
        /* First identify the converged radial bins */
        for (r_conv_i = 0, ibin = 0; ibin < halos[i].prof.nbins; ibin++) {
          rad     = halos[i].prof.r[ibin];
          
          /* check for converged radius (Power et al. 2003) */
          //t_relax = halos[i].prof.npart[ibin] / log(rad / halos[i].spaRes) * rad / sqrt(halos[i].prof.v2_circ[ibin]);           // Eq.(4)
          //t0      = 0.9*calc_t(global.a) * simu.t_unit * Mpc / Gyr; /* age of the Universe in Gyr/h */
          
          t_relax = halos[i].prof.npart[ibin] / log(halos[i].prof.npart[ibin]) / 8. * rad / sqrt(halos[i].prof.v2_circ[ibin]);  // Eq.(20)
          t0      = 0.6*halos[i].prof.r[halos[i].prof.nbins-1] / sqrt(halos[i].prof.v2_circ[halos[i].prof.nbins-1]) * r_fac/sqrt(phi_fac) * 1E3;
          
          
          /* convert to (Mpc/h) / (km/sec) ~ (Mpc/h) / (kpc/Gyr) */
          t_relax *= r_fac / sqrt(phi_fac);
          
          /* convert to Gyr/h */
          t_relax *= 1E3;
          
          
          /* if not converged, write negative radius into *.AHF_profiles */
          if (t_relax < t0)
           /* The +1 is merely to keep the name consistent:
            * r_conv_i should give the smallest converged radius,
            * not the bin before it. Hence to check for converged
            * bin or not, i<r_conv_i is the thing to do. */
            r_conv_i = ibin + 1;
        }

        /* Write the number of profile lines and update total number */
        tmp_int = (bin_int_t)(halos[i].prof.nbins);
        fwrite(&tmp_int, sizeof(bin_int_t), 1, f);
        total_num_lines += tmp;
        
        /* Now we actually write the stuff to the file */
        for (ibin = 0; ibin < halos[i].prof.nbins; ibin++) {
          tmp_float  = (bin_float_t)(halos[i].prof.r[ibin]);
          tmp_float  = (ibin < r_conv_i) ? -tmp_float : tmp_float;
          tmp_float *= x_fac * 1000.;
          WRITE_F;
          tmp_int    = (bin_int_t)(halos[i].prof.npart[ibin]);
          WRITE_I;
          tmp_float  = (bin_float_t)(halos[i].prof.nvpart[ibin] * m_fac);
          WRITE_F;
          tmp_float  = (bin_float_t)(halos[i].prof.ovdens[ibin] * rho_fac / global.rho_vir);
          WRITE_F;
          tmp_float  = (bin_float_t)(halos[i].prof.dens[ibin] * rho_fac   / global.rho_vir);
          WRITE_F;
          tmp_float  = (bin_float_t)(sqrt(halos[i].prof.v2_circ[ibin] * phi_fac));
          WRITE_F;
          tmp_float  = (bin_float_t)(sqrt(halos[i].prof.v_esc2[ibin] * phi_fac));
          WRITE_F;
          tmp_float  = (bin_float_t)(halos[i].prof.sig_v[ibin] * v_fac);
          WRITE_F;
          tmp_float  = (bin_float_t)(halos[i].prof.Lx[ibin] * m_fac * r_fac * v_fac);
          WRITE_F;
          tmp_float  = (bin_float_t)(halos[i].prof.Ly[ibin] * m_fac * r_fac * v_fac);
          WRITE_F;
          tmp_float  = (bin_float_t)(halos[i].prof.Lz[ibin] * m_fac * r_fac * v_fac);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.axis2[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.axis3[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E1x[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E1y[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E1z[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E2x[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E2y[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E2z[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E3x[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E3y[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.E3z[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.Ekin[ibin] * m_fac * pow2(v_fac));
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.Epot[ibin] * m_fac * phi_fac);
          WRITE_F;
#  if (defined GAS_PARTICLES)
          tmp_float = (bin_float_t)(halos[i].prof.M_gas[ibin] * m_fac);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.M_star[ibin] * m_fac);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.u_gas[ibin] * u_fac); // remember, in HaloProfiles() we converted cur_part->u to code internal units, here we need to undo that!
          WRITE_F;
#ifdef METALHACK
          tmp_float = (bin_float_t)(halos[i].prof.z_gas[ibin]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.z_star[ibin]);
          WRITE_F;
#endif
#  endif
#  if (defined AHFphspdens)
          tmp_float
          = (bin_float_t)(halos[i].prof.sigma2_vx_sh[ibin]
                          * v_fac
                          * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.sigma2_vy_sh[ibin]
                          * v_fac
                          * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.sigma2_vz_sh[ibin]
                          * v_fac
                          * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.sigma2_vr_sh[ibin]
                          * v_fac
                          * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.sigma2_vtheta_sh[ibin]);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.sigma2_vphi_sh[ibin]);
          WRITE_F;
#    if (defined AHFmeanvelocities)
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vx_sh[ibin] * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vy_sh[ibin] * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vz_sh[ibin] * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vr_sh[ibin] * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vtheta_sh[ibin
                                                       ]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.mean_vphi_sh[ibin]);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vx_sp[ibin] * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vy_sp[ibin] * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vz_sp[ibin] * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vr_sp[ibin] * v_fac);
          WRITE_F;
          tmp_float
          = (bin_float_t)(halos[i].prof.mean_vtheta_sp[ibin
                                                       ]);
          WRITE_F;
          tmp_float = (bin_float_t)(halos[i].prof.mean_vphi_sp[ibin]);
          WRITE_F;
#    endif
#  endif
        } /* End of for-loop over profile lines */
        
        /* And finally increment the halo counter */
        real_num_halos++;
      } /* End of if for suitable halo */
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
	}
  
	/* Rewind and put the right numbers in the front (after the sizes) */
	fseek(f, (long)sizeof(int32_t), SEEK_SET); // skip the "one"
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
	//fwrite(&total_num_lines, sizeof(uint64_t), 1, f);
  
	/* Close the files */
	fclose(f);
	if (f_info != NULL)
		fclose(f_info);
  
#  if (defined VERBOSE)
	/* End the 'Writing file..' statement started when opening the file */
	fprintf(stderr, "done\n");
#  endif
  
	/* Done */
	return;
} /* ahf_binwrite_profiles */

void
ahf_binwrite_particles(char          *prefix,
                    HALO          *halos,
                    unsigned long *idx,
                    int           numHalos)
{
	FILE        *f;
	FILE        *f_info;
	uint32_t    sizes[2];
	bin_int_t   tmp_int;
	bin_float_t tmp_float;
  bin_ulong_t tmp_ulong;
	int         i = 0, j = 0, k = 0;
	uint32_t    num_columns;
	uint64_t    real_num_halos;
	uint64_t    total_num_particles;
	uint64_t    tmp;
	int32_t     particle_type;
  int32_t     one=1;
	double      particle_mass;
	partptr     cur_part;
	int8_t      typefield[LEN_TYPEFIELD];
  
	/* Store the number of bytes used for integer and float values */
	sizes[0] = sizeof(bin_int_t);
	sizes[1] = sizeof(bin_float_t);
  
	/* Initialize typefield to int */
	for (i = 0; i < LEN_TYPEFIELD; i++)
		typefield[i] = INT8_C(0);
  
	/* Open the files */
	ahf_binwrite_open_files(&f, &f_info, prefix, "AHF_particles");
  
	/* Write info, if required */
	if (f_info != NULL)
    fprintf(f_info, "one (used to check for byteswap)\nnumHalos\nnumColumns\n numPart\n");
	i = 0;
	if (f_info != NULL)
		fprintf(f_info, "  ID(%i) ", i + 1);
	i += 1;
	if (f_info != NULL)
		fprintf(f_info, "pType(%i)\n", i + 1);
	i          += 1;
	num_columns = i;
  
  /* Just write the number "1" into the file */
  fwrite(&one, sizeof(int32_t), 1, f);

	/* Store the sizes */
	//fwrite(sizes, sizeof(uint32_t), 2, f);
  
	/* Write dummy values */
	real_num_halos      = UINT64_C(0);
	total_num_particles = UINT64_C(0);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
	//fwrite(&total_num_particles, sizeof(uint64_t), 1, f);
  
	/* Write number of columns */
	fwrite(&num_columns, sizeof(uint32_t), 1, f);
  
	/* Store the typefield */
	//fwrite(typefield, sizeof(int8_t), num_columns, f);
  
	/* Now loop over all haloes for writing */
	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART)) {
        /* Write the number of particles in this halo and update total number */
        tmp                  = (uint64_t)(halos[i].npart);
        fwrite(&tmp, sizeof(uint64_t), 1, f);
        total_num_particles += tmp;
        
        /* Write the haloID */
#ifdef AHFnewHaloIDs
        tmp = (uint64_t)halos[i].haloID;
#else
        tmp = (uint64_t)j;
#endif
        fwrite(&tmp, sizeof(uint64_t), 1, f);
        
        /* Loop over all particles in this halo */
        for (k = 0; k < halos[i].npart; k++) {
          cur_part = global.fst_part + halos[i].ipart[k];
          tmp_ulong  = (bin_ulong_t)(cur_part->id);
          WRITE_LU;
#if (defined GAS_PARTICLES)
          particle_type = (cur_part->u >= PGAS) ? (int32_t)PGAS : (int32_t)(-cur_part->u);
#else
          particle_type = 1;
#endif
          fwrite(&particle_type, sizeof(int32_t), 1, f);
        } /* End of loop over particles */
        
        /* Count this halo */
        real_num_halos++;
      } /* end of if selecting 'proper' haloes */
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
	} /* end of for looping over all haloes*/
  
	/* Rewind and put the right numbers in the front */
	fseek(f, (long)sizeof(int32_t), SEEK_SET); // skip the "one"
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
	//fwrite(&total_num_particles, sizeof(uint64_t), 1, f);
  
	/* Close the files */
	fclose(f);
	if (f_info != NULL)
		fclose(f_info);
  
#  if (defined VERBOSE)
	/* End the 'Writing file..' statement started when opening the file */
	fprintf(stderr, "done\n");
#  endif
  
	/* Done */
	return;
} /* ahf_binwrite_particles */

void
ahf_binwrite_halos(char          *prefix,
                   HALO          *halos,
                   unsigned long *idx,
                   int           numHalos)
{
	FILE        *f;
	FILE        *f_info;
	uint32_t    num_columns;
	uint32_t    sizes[2];
	int         i, j;
	bin_int_t   tmp_int;
	bin_float_t tmp_float;
  bin_ulong_t tmp_ulong;
	uint64_t    real_num_halos;
	int8_t      typefield[LEN_TYPEFIELD];
  int32_t     one=1;
  
	/* Store the number of bytes used for integer and float values */
	sizes[0] = sizeof(bin_int_t);
	sizes[1] = sizeof(bin_float_t);
  
	/* Initialize typefield to float */
	for (i = 0; i < LEN_TYPEFIELD; i++)
		typefield[i] = INT8_C(1);
  
	/* Open the files */
	ahf_binwrite_open_files(&f, &f_info, prefix, "AHF_halos");
  
	/* Write info, if required */
	i = 0;
	if (f_info != NULL)
		fprintf(f_info,
            "one (used to check for byteswap)\n"
            "numHalos\n"
            "numColumns\n"
		        "ID(%i) hostHalo(%i) numSubStruct(%i) "
		        "Mvir(%i) npart(%i) "
            "Xc(%i) Yc(%i) Zc(%i) "
            "VXc(%i) VYc(%i) VZc(%i) "
            "Rvir(%i) Rmax(%i) r2(%i) mbp_offset(%i) com_offset(%i) "
            "Vmax(%i) v_esc(%i) sigV(%i) "
            "lambda(%i) lambdaE(%i) "
            "Lx(%i) Ly(%i) Lz(%i) "
            "b(%i) c(%i) "
		        "Eax(%i) Eay(%i) Eaz(%i) "
		        "Ebx(%i) Eby(%i) Ebz(%i) "
		        "Ecx(%i) Ecy(%i) Ecz(%i) "
		        "ovdens(%i) nbins(%i) fMhires(%i) Ekin(%i) Epot(%i) SurfP(%i) "
		        "Phi0(%i) cNFW(%i) ",
		        i + 1, i + 2, i + 3,
		        i + 4, i + 5, i + 6, i + 7, i + 8, i + 9,
		        i + 10, i + 11, i + 12, i + 13, i + 14,
		        i + 15, i + 16, i + 17, i + 18,
		        i + 19, i + 20, i + 21, i + 22,
		        i + 23, i + 24, i + 25, i + 26,
		        i + 27, i + 28, i + 29, i + 30,
		        i + 31, i + 32, i + 33, i + 34,
		        i + 35, i + 36, i + 37, i + 38,
            i + 39, i + 40, i + 41, i + 42, i + 43);
	i += 43;
	typefield[0] = INT8_C(0);
	typefield[1] = INT8_C(0);
	typefield[2] = INT8_C(0);
	typefield[4] = INT8_C(0);
	typefield[36] = INT8_C(0);
#  if (defined GAS_PARTICLES)
	if (f_info != NULL)
		fprintf(f_info,
		        "n_gas(%i) M_gas(%i) lambda_gas(%i) lambdaE_gas(%i) "
		        "Lx_gas(%i) Ly_gas(%i) Lz_gas(%i) "
		        "b_gas(%i) c_gas(%i) "
            "Eax_gas(%i) Eay_gas(%i) Eaz_gas(%i) "
		        "Ebx_gas(%i) Eby_gas(%i) Ebz_gas(%i) "
		        "Ecx_gas(%i) Ecy_gas(%i) Ecz_gas(%i) "
		        "Ekin_gas(%i) Epot_gas(%i) "
		        "n_star(%i) M_star(%i) lambda_star(%i) lambdaE_star(%i) "
		        "Lx_star(%i) Ly_star(%i) Lz_star(%i) "
		        "b_star(%i) c_star(%i) "
            "Eax_star(%i) Eay_star(%i) Eaz_star(%i) "
		        "Ebx_star(%i) Eby_star(%i) Ebz_star(%i) "
		        "Ecx_star(%i) Ecy_star(%i) Ecz_star(%i) "
		        "Ekin_star(%i) Epot_star(%i) ",
		        i + 1, i + 2, i + 3,
		        i + 4, i + 5, i + 6,
		        i + 7, i + 8, i + 9, i + 10,
		        i + 11, i + 12, i + 13, i + 14,
		        i + 15, i + 16, i + 17, i + 18,
		        i + 19, i + 20, i + 21, i + 22,
		        i + 23, i + 24, i + 25,
		        i + 26, i + 27, i + 28,
		        i + 29, i + 30, i + 31, i + 32,
		        i + 33, i + 34, i + 35, i + 36,
		        i + 37, i + 38, i + 39, i + 40);
	i += 40;
	typefield[43] = INT8_C(0);
	typefield[63] = INT8_C(0);
#ifdef METALHACK
	if (f_info != NULL)
		fprintf(f_info,
		        "mean_z_gas(%i) mean_z_star(%i) ",
		        i + 1, i + 2);
	i += 2;
#endif
#  endif
	if (f_info != NULL)
		fprintf(f_info,"\n");
	num_columns   = (uint32_t)i;
  
  /* Just write the number "1" into the file */
  fwrite(&one, sizeof(int32_t), 1, f);
  
	/* Store the sizes */
	//fwrite(sizes, sizeof(uint32_t), 2, f);
  
	/* Write a dummy value for the number of halos */
	real_num_halos = UINT64_C(0);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
  
	/* Store the number of columns used */
	fwrite(&num_columns, sizeof(uint32_t), 1, f);
  
	/* Store the typefield */
	//fwrite(typefield, sizeof(int8_t), num_columns, f);
  
	/* Now loop over all haloes for writing */
	for (j = 0; j < numHalos; j++) {
		i = idx[j];
    
#if (defined WITH_MPI || defined AHFrestart)
    if(halos[i].ignoreme == FALSE) {
#endif
      if ((halos[i].npart >= simu.AHF_MINPART)) {
#ifdef AHFnewHaloIDs
        tmp_ulong = (bin_ulong_t)halos[i].haloID;
        WRITE_LU;
        tmp_ulong = (bin_ulong_t)halos[i].hostHaloID;
        WRITE_LU;
#else
        tmp_ulong = (bin_ulong_t)real_num_halos;
        WRITE_LU;
        tmp_ulong = (bin_ulong_t)(idx_inv(idx, numHalos, halos[i].hostHalo));
        WRITE_LU;
#endif
        tmp_int   = (bin_int_t)(halos[i].numSubStruct);
        WRITE_I;
        tmp_float = (bin_float_t)(halos[i].M_vir * m_fac);
        WRITE_F;
        tmp_int   = (bin_int_t)(halos[i].npart);
        WRITE_I;
        tmp_float = (bin_float_t)(halos[i].pos.x * x_fac * 1000.);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].pos.y * x_fac * 1000.);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].pos.z * x_fac * 1000.);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].vel.x * v_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].vel.y * v_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].vel.z * v_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].R_vir * x_fac * 1000.);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].R_max * x_fac * 1000.);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].r2 * x_fac * 1000.);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].mbp_offset * x_fac * 1000.);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].com_offset * x_fac * 1000.);
        WRITE_F;
        tmp_float = (bin_float_t)sqrt(halos[i].V2_max * phi_fac);
        WRITE_F;
        tmp_float = (bin_float_t)sqrt(halos[i].v_esc2 * phi_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].sigV * v_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].lambda);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].lambdaE);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].AngMom.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].AngMom.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].AngMom.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].axis.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].axis.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E1.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E1.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E1.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E2.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E2.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E2.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E3.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E3.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].E3.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].ovdens);
        WRITE_F;
        tmp_int   = (bin_int_t)  (halos[i].prof.nbins);
        WRITE_I;
        tmp_float = (bin_float_t)(halos[i].fMhires);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].Ekin  * m_fac * pow2(v_fac));
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].Epot  * m_fac * phi_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].SurfP * m_fac * pow2(v_fac));
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].Phi0  * m_fac * phi_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].cNFW);
        WRITE_F;
#  if (defined GAS_PARTICLES)
        tmp_int   = (bin_int_t)(halos[i].gas_only.npart);
        WRITE_I;
        tmp_float = (bin_float_t)(halos[i].gas_only.Mass * m_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.lambda);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.lambdaE);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.AngMom.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.AngMom.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.AngMom.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.axis.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.axis.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E1.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E1.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E1.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E2.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E2.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E2.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E3.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E3.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.E3.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.Ekin * m_fac * pow2(v_fac));
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].gas_only.Epot * m_fac * phi_fac);
        WRITE_F;
        tmp_int   = (bin_int_t)(halos[i].stars_only.npart);
        WRITE_I;
        tmp_float = (bin_float_t)(halos[i].stars_only.Mass * m_fac);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.lambda);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.lambdaE);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.AngMom.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.AngMom.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.AngMom.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.axis.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.axis.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E1.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E1.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E1.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E2.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E2.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E2.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E3.x);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E3.y);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.E3.z);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.Ekin * m_fac * pow2(v_fac));
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].stars_only.Epot * m_fac * phi_fac);
        WRITE_F;
#ifdef METALHACK
        tmp_float = (bin_float_t)(halos[i].mean_z_gas);
        WRITE_F;
        tmp_float = (bin_float_t)(halos[i].mean_z_star);
        WRITE_F;
#endif
#  endif
        
        real_num_halos++;
      } /* End of if selecting proper haloes */
#if (defined WITH_MPI || defined AHFrestart)
    }
#endif
	} /* End of halo loop */
  
	/* Rewind and put the right numbers in the front */
  fseek(f, (long)sizeof(int32_t),SEEK_SET); // skip the "one"
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
  
	/* Close the files */
	fclose(f);
	if (f_info != NULL)
		fclose(f_info);
  
#  if (defined VERBOSE)
	/* End the 'Writing file..' statement started when opening the file */
	fprintf(stderr, "done\n");
#  endif
  
	/* Done */
	return;
} /* ahf_binwrite_halos */

#  undef BUFFER_TOO_SMALL_ERROR_PROFILES

#endif /* AHFbinary */


void write_halos_line(FILE *fout, HALO *halos, long unsigned j, long unsigned i, unsigned long *idx, int numHalos)
{
#ifdef AHFnewHaloIDs
  fprintf(fout, "%22"PRIu64,		        halos[i].haloID);
  fprintf(fout, "%22"PRIu64,            halos[i].hostHaloID);
#else
#ifdef SUSSING2013
  fprintf(fout, "%10"PRIu64,	          getSussing2013ID(simu.isnap, j));
  if(halos[i].hostHalo < 0)
    fprintf(fout, "\t%10ld",            (long)0);
  else
    fprintf(fout, "\t%10"PRIu64,        getSussing2013ID(simu.isnap, idx_inv(idx, numHalos, halos[i].hostHalo)));
#else
  fprintf(fout, "%10lu",			          (long unsigned)j);
  fprintf(fout, "\t%10d",               idx_inv(idx, numHalos, halos[i].hostHalo));
#endif
#endif
  fprintf(fout, "\t%10d",               halos[i].numSubStruct);
  fprintf(fout, "\t%12.6g",			        halos[i].M_vir * m_fac);
  fprintf(fout, "\t%10ld",		          halos[i].npart);
#ifdef AHFundoPositionShiftAndScale
  fprintf(fout, "\t%16.8f",			        ((halos[i].pos.x/simu.pos_scale)-simu.pos_shift[X]));
  fprintf(fout, "\t%16.8f",			        ((halos[i].pos.y/simu.pos_scale)-simu.pos_shift[Y]));
  fprintf(fout, "\t%16.8f",			        ((halos[i].pos.z/simu.pos_scale)-simu.pos_shift[Z]));
#else
  fprintf(fout, "\t%16.8f",			        halos[i].pos.x * x_fac*1000.);
  fprintf(fout, "\t%16.8f",			        halos[i].pos.y * x_fac*1000.);
  fprintf(fout, "\t%16.8f",			        halos[i].pos.z * x_fac*1000.);
#endif
  fprintf(fout, "\t%8.2f",			        halos[i].vel.x * v_fac);
  fprintf(fout, "\t%8.2f",			        halos[i].vel.y * v_fac);
  fprintf(fout, "\t%8.2f",			        halos[i].vel.z * v_fac);
  fprintf(fout, "\t%10.2f",			        halos[i].R_vir * x_fac * 1000.);
  fprintf(fout, "\t%10.2f",			        halos[i].R_max * x_fac * 1000.);
  fprintf(fout, "\t%10.5f",			        halos[i].r2 * x_fac * 1000.);
  fprintf(fout, "\t%10.5f",			        halos[i].mbp_offset * x_fac * 1000.);
  fprintf(fout, "\t%10.5f",			        halos[i].com_offset * x_fac * 1000.);
  fprintf(fout, "\t%8.2f",			        sqrt(halos[i].V2_max * phi_fac));
  fprintf(fout, "\t%12.6f",			        sqrt(halos[i].v_esc2 * phi_fac));
  fprintf(fout, "\t%8.2f",			        halos[i].sigV * v_fac);
  fprintf(fout, "\t%10.6f",			        halos[i].lambda);
  fprintf(fout, "\t%12.6f",			        halos[i].lambdaE);
  fprintf(fout, "\t%12.4g",			        halos[i].AngMom.x);
  fprintf(fout, "\t%12.4g",			        halos[i].AngMom.y);
  fprintf(fout, "\t%12.4g",			        halos[i].AngMom.z);
  fprintf(fout, "\t%10.6f",			        halos[i].axis.y);
  fprintf(fout, "\t%10.6f",			        halos[i].axis.z);
  fprintf(fout, "\t%10.6f",			        halos[i].E1.x);
  fprintf(fout, "\t%10.6f",			        halos[i].E1.y);
  fprintf(fout, "\t%10.6f",			        halos[i].E1.z);
  fprintf(fout, "\t%10.6f",			        halos[i].E2.x);
  fprintf(fout, "\t%10.6f",			        halos[i].E2.y);
  fprintf(fout, "\t%10.6f",			        halos[i].E2.z);
  fprintf(fout, "\t%10.6f",			        halos[i].E3.x);
  fprintf(fout, "\t%10.6f",			        halos[i].E3.y);
  fprintf(fout, "\t%10.6f",			        halos[i].E3.z);
  fprintf(fout, "\t%8.2f",			        halos[i].ovdens);
  fprintf(fout, "\t%6d",			          halos[i].prof.nbins);
  fprintf(fout, "\t%8.6f",              halos[i].fMhires);
  fprintf(fout, "\t%12.6g",			        halos[i].Ekin * m_fac * pow2(v_fac));
  fprintf(fout, "\t%12.6g",			        halos[i].Epot * m_fac * phi_fac);
  fprintf(fout, "\t%12.6g",			        halos[i].SurfP * m_fac * pow2(v_fac));
  fprintf(fout, "\t%12.6g",			        halos[i].Phi0 * phi_fac);
  fprintf(fout, "\t%12.6g",			        halos[i].cNFW);
#ifdef AHFvmbp
  fprintf(fout, "\t%8.2f",			        halos[i].vel_mbp.x * v_fac);
  fprintf(fout, "\t%8.2f",			        halos[i].vel_mbp.y * v_fac);
  fprintf(fout, "\t%8.2f",			        halos[i].vel_mbp.z * v_fac);
#endif
#  ifdef GAS_PARTICLES
  fprintf(fout, "\t%10ld",			        halos[i].gas_only.npart);
  fprintf(fout, "\t%12.6g",			        halos[i].gas_only.Mass * m_fac);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.lambda);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.lambdaE);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.AngMom.x);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.AngMom.y);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.AngMom.z);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.axis.y);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.axis.z);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E1.x);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E1.y);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E1.z);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E2.x);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E2.y);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E2.z);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E3.x);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E3.y);
  fprintf(fout, "\t%10.6f",			        halos[i].gas_only.E3.z);
  fprintf(fout, "\t%12.6g",			        halos[i].gas_only.Ekin * m_fac * pow2(v_fac));
  fprintf(fout, "\t%12.6g",			        halos[i].gas_only.Epot * m_fac * phi_fac);
  fprintf(fout, "\t%10ld",			        halos[i].stars_only.npart);
  fprintf(fout, "\t%12.6g",			        halos[i].stars_only.Mass * m_fac);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.lambda);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.lambdaE);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.AngMom.x);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.AngMom.y);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.AngMom.z);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.axis.y);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.axis.z);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E1.x);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E1.y);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E1.z);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E2.x);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E2.y);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E2.z);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E3.x);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E3.y);
  fprintf(fout, "\t%10.6f",			        halos[i].stars_only.E3.z);
  fprintf(fout, "\t%12.6g",			        halos[i].stars_only.Ekin * m_fac * pow2(v_fac));
  fprintf(fout, "\t%12.6g",			        halos[i].stars_only.Epot * m_fac * phi_fac);
#  endif
  
#  ifdef METALHACK
  fprintf(fout, "\t%e",			        halos[i].mean_z_gas);
  fprintf(fout, "\t%e",			        halos[i].mean_z_star);
#  endif
  fprintf(fout, "\n");
}

#endif /* AHF2 */
