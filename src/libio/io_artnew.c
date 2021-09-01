// Copyright (C) 2011, Steffen Knollmann
// Released under the terms of the GNU General Public License version 3.
// This file is part of `AHF'.


/*--- Doxygen file description ------------------------------------------*/

/**
 * @file  io_artnew.c
 *
 * Implements the functionality to read from ART files.
 */


/*--- Includes ----------------------------------------------------------*/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stddef.h>
#ifdef WITH_MPI
#  include <mpi.h>
#endif
#include "io_artnew.h"
#include "io_file.h"
#include "io_util.h"
#include "art.h"
#include "xmem.h"
#include "xfile.h"
#include "xstring.h"
#include "stai.h"




// ARTNEW funktioniert eh etwas
// anders: Die Datei, die man in der ahf.input angibt, ist nicht die
// ART-Datei, sondern eine kleine Text-Datei, die so aussieht:
// Path = dollerPfad/
// Suffix = _super.DAT
// numFiles = 3
// truncateNrowc = true/false

// Das wuerde dann die Dateien
// dollerPfad/PMcrd_super.DAT
// und
// dollerPfad/PMcrs0_super.DAT
// dollerPfad/PMcrs1_super.DAT
// dollerPfad/PMcrs2_super.DAT
// lesen (in etwa, bin mir gerade nicht 100% sicher, wie genau er den
//        Namen bastelt). truncateNrowc macht irgendwas komischen mit den
// Pagesizes. War irgendwie relevant fuer Multi-File ART die groesser als
// 1024^3 sind (da ist die Pagesize immer 1024).  Format-Voodoo....




/*--- Implementations of the main structure -----------------------------*/
#include "io_artnew_def.h"


/*--- Prototypes of local functions -------------------------------------*/
inline static void
local_getAllFromInfo(const char *fname,
                     char       **path,
                     char       **suffix,
                     int        *numFiles,
                     bool       *truncateNrowc);

inline static void
local_convertStrgToStais(io_file_strg_struct_t *strg, stai_t data[6]);

inline static void
local_convertUnits(io_file_strg_struct_t *s,
                   double                xfac,
                   double                vfac,
                   uint64_t              pr);

inline static void
local_handleStrgWeightIdU(io_file_strg_struct_t *s,
                          uint64_t              ps,
                          uint64_t              pr);


/*--- Implementations of exported functions -----------------------------*/
extern io_artnew_t
io_artnew_open(io_logging_t   log,
               char           *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t       reader)
{
	io_artnew_t f;
	char        *path;
	char        *suffix;
	int         numFiles;
	bool        truncateNrowc;

	f = xmalloc(sizeof(struct io_artnew_struct));

	local_getAllFromInfo(fname, &path, &suffix, &numFiles, &truncateNrowc);
	f->ftype  = IO_FILE_ARTNEW;
	f->art    = art_new(path, suffix, numFiles);
	art_setTruncateNrowc(f->art, truncateNrowc);
	f->header = NULL;
	f->pmass  = 0.0;
#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	MPI_Comm_split(MPI_COMM_WORLD, 1,
	               f->rank, &(f->mycomm));
	MPI_Comm_size(f->mycomm, &(f->size_mycomm));
	MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
#endif
	io_logging_msg(log, INT32_C(1),
	               "New ARTnew reader created.");
	xfree(path);
	xfree(suffix);

	return f;
}

extern void
io_artnew_close(io_logging_t log,
                io_artnew_t  *f)
{
	assert(log != NULL);

	if ((f != NULL) && (*f != NULL))
		return;

	if ((*f)->art != NULL)
		art_del(&((*f)->art));

	xfree(*f);
	*f = NULL;
}

extern void
io_artnew_init(io_logging_t log,
               io_artnew_t  f)
{
	assert(f != NULL);

	art_attachHeaderFromFile(f->art);
	f->header = art_getHeaderHandle(f->art);
	f->pmass  = artHeader_getFactorFileWeightToMsunh(f->header)
	            * artHeader_getWspecies(f->header, 0);

	io_logging_msg(log, INT32_C(1),
	               "ARTnew reader initialized.");
}

extern uint64_t
io_artnew_readpart(io_logging_t          log,
                   io_artnew_t           f,
                   uint64_t              pskip,
                   uint64_t              pread,
                   io_file_strg_struct_t strg)
{
	stai_t   data[6];
	uint64_t numParticlesRead = UINT64_C(0);
	double   xfac, vfac;

	assert(f != NULL);
	assert(f->art != NULL);

	local_convertStrgToStais(&strg, data);

	numParticlesRead = art_read(f->art, pskip, pread, data);

	// Converting from file units to Mpc/h and km/s.
	xfac  = artHeader_getFactorFilePositionToMpch(f->header);
	vfac  = artHeader_getFactorFileVelocityToKms(f->header);
	// Now converting from Mpc/h and km/s to AHF units.
	xfac /= artHeader_getBoxsizeInMpch(f->header);
	vfac *= artHeader_getAexpn(f->header)
	        / (artHeader_getBoxsizeInMpch(f->header) * 100.);

	local_convertUnits(&strg, xfac, vfac, numParticlesRead);
	local_handleStrgWeightIdU(&strg, pskip, numParticlesRead);

	for (int i = 0; i < 6; i++)
		stai_del(data + i);

	return numParticlesRead;
}

extern bool
io_artnew_get(io_logging_t  log,
              io_artnew_t   f,
              io_file_get_t what,
              void          *res)
{
	switch (what) {
	case IO_FILE_GET_NOPART_IN_FILE:
		*((long *)res) = (long)artHeader_getNumParticlesTotal(f->header);
		break;
	case IO_FILE_GET_NOPART:
		*((long *)res) = (long)artHeader_getNumParticlesTotal(f->header);
		break;
	case IO_FILE_GET_NOVPART:
		// This only works if there is only 1 particle species!
		*((double *)res) = (double)artHeader_getNumParticlesTotal(f->header);
		break;
	case IO_FILE_GET_NOSPECIES:
		// This only works if there is only 1 particle species!
		*((int *)res) = artHeader_getNspecies(f->header);
		break;
	case IO_FILE_GET_BOXSIZE:
		*((double *)res) = artHeader_getBoxsizeInMpch(f->header);
		break;
	case IO_FILE_GET_PMASS:
		*((double *)res) = f->pmass;
		break;
	case IO_FILE_GET_ZINITIAL:
		*((double *)res) = 1. / (artHeader_getAexp0(f->header)) - 1.;
		break;
	case IO_FILE_GET_Z:
		*((double *)res) = 1. / (artHeader_getAexpn(f->header)) - 1.;
		break;
	case IO_FILE_GET_AINITIAL:
		*((double *)res) = artHeader_getAexp0(f->header);
		break;
	case IO_FILE_GET_A:
		*((double *)res) = artHeader_getAexpn(f->header);
		break;
	case IO_FILE_GET_OMEGA0:
		*((double *)res) = artHeader_getOm0(f->header);
		break;
	case IO_FILE_GET_OMEGAL:
		*((double *)res) = artHeader_getOml0(f->header);
		break;
	case IO_FILE_GET_H:
		*((double *)res) = artHeader_getHubble(f->header);
		break;
	case IO_FILE_GET_DOUBLE:
		io_logging_warn(log, INT32_C(1),
		                "ART files are always single precisions.");
		*((int *)res) = 0;
		break;
	case IO_FILE_GET_MMASS:
		io_logging_warn(log, INT32_C(1),
		                "Assuming that there is only 1 particle species.");
		*((int *)res) = 0;
		break;
	case IO_FILE_GET_NOTSTEP:
		*((int32_t *)res) = artHeader_getIstep(f->header);
		break;
	case IO_FILE_GET_TSTEP:
		*((double *)res) = artHeader_getAstep(f->header);
		break;
	case IO_FILE_GET_HEADERSTR:
		*((const char **)res) = artHeader_getHeaderString(f->header);
		break;
	case IO_FILE_GET_MINWEIGHT:
	case IO_FILE_GET_MAXWEIGHT:
		// We assume that there is only 1 particle species!
		*((double *)res) = 1.0;
		break;
	default:
		io_logging_fatal(log, "Requesting something unkown in %s.",
		                 __func__);
		return false;
	} /* switch */
	return true;
} /* io_artnew_get */

extern void
io_artnew_log(io_logging_t log, io_artnew_t f)
{
	assert(f != NULL);
	art_prettyPrint(f->art, NULL, log->logfile);
	if (f->header != NULL)
		artHeader_prettyPrint(f->header, NULL, log->logfile);
}

/*--- Prototypes of local functions -------------------------------------*/
static void
local_getAllFromInfo(const char *fname,
                     char       **path,
                     char       **suffix,
                     int        *numFiles,
                     bool       *truncateNrowc)
{
	FILE *f;
	char buffer[2048];

	f = xfopen(fname, "r");

	fscanf(f, "Path = %s \n", buffer);
	*path = xstrdup(buffer);

	fscanf(f, "Suffix = %s \n", buffer);
	*suffix = xstrdup(buffer);

	fscanf(f, "numFiles = %i \n", numFiles);

	fscanf(f, "truncateNrowc = %s \n", buffer);
	if (strncmp(buffer, "true", 4) == 0)
		*truncateNrowc = true;
	else
		*truncateNrowc = false;


	xfclose(&f);
}

inline static void
local_convertStrgToStais(io_file_strg_struct_t *strg, stai_t data[6])
{
	data[0] = stai_new(strg->posx.val, strg->bytes_float,
	                   strg->posx.stride);
	data[1] = stai_new(strg->posy.val, strg->bytes_float,
	                   strg->posy.stride);
	data[2] = stai_new(strg->posz.val, strg->bytes_float,
	                   strg->posz.stride);
	data[3] = stai_new(strg->momx.val, strg->bytes_float,
	                   strg->momx.stride);
	data[4] = stai_new(strg->momy.val, strg->bytes_float,
	                   strg->momy.stride);
	data[5] = stai_new(strg->momz.val, strg->bytes_float,
	                   strg->momz.stride);
}

inline static void
local_convertUnits(io_file_strg_struct_t *s,
                   double                xfac,
                   double                vfac,
                   uint64_t              pr)
{
	for (uint64_t i = UINT64_C(0); i < pr; i++) {
		if (s->bytes_float == sizeof(float)) {
			*((float *)(s->posx.val)) -= 1.0f;
			*((float *)(s->posx.val)) *= xfac;
			*((float *)(s->posy.val)) -= 1.0f;
			*((float *)(s->posy.val)) *= xfac;
			*((float *)(s->posz.val)) -= 1.0f;
			*((float *)(s->posz.val)) *= xfac;
			*((float *)(s->momx.val)) *= vfac;
			*((float *)(s->momy.val)) *= vfac;
			*((float *)(s->momz.val)) *= vfac;
		} else {
			*((double *)(s->posx.val)) -= 1.0f;
			*((double *)(s->posx.val)) *= xfac;
			*((double *)(s->posy.val)) -= 1.0f;
			*((double *)(s->posy.val)) *= xfac;
			*((double *)(s->posz.val)) -= 1.0f;
			*((double *)(s->posz.val)) *= xfac;
			*((double *)(s->momx.val)) *= vfac;
			*((double *)(s->momy.val)) *= vfac;
			*((double *)(s->momz.val)) *= vfac;
		}
		s->posx.val = (char *)(s->posx.val) + s->posx.stride;
		s->posy.val = (char *)(s->posy.val) + s->posy.stride;
		s->posz.val = (char *)(s->posz.val) + s->posz.stride;
		s->momx.val = (char *)(s->momx.val) + s->momx.stride;
		s->momy.val = (char *)(s->momy.val) + s->momy.stride;
		s->momz.val = (char *)(s->momz.val) + s->momz.stride;
	}
}

inline static void
local_handleStrgWeightIdU(io_file_strg_struct_t *s,
                          uint64_t              ps,
                          uint64_t              pr)
{
	if (s->weight.val != NULL) {
		for (uint64_t i = UINT64_C(0); i < pr; i++) {
			if (s->bytes_float == sizeof(float))
				*((float *)(s->weight.val)) = 1.0f;
			else
				*((double *)(s->weight.val)) = 1.0;
			s->weight.val = (void *)(((char *)(s->weight.val))
			                         + s->weight.stride);
		}
	}
	if (s->id.val != NULL) {
		
		for (uint64_t i = UINT64_C(0); i < pr; i++) {
			if (s->bytes_int == sizeof(uint32_t))
				*((uint32_t *)(s->id.val)) = (uint32_t)(ps + i);
			else
				*((uint64_t *)(s->id.val)) = ps + i;
			s->id.val = (void *)(((char *)(s->id.val)) + s->id.stride);
		}
	}
	if (s->u.val != NULL) {
		for (uint64_t i = UINT64_C(0); i < pr; i++) {
			if (s->bytes_float == sizeof(float))
				*((float *)(s->u.val)) = -1.0f;
			else
				*((double *)(s->u.val)) = -1.0;
			s->u.val = (void *)(((char *)(s->u.val)) + s->u.stride);
		}
	}
}
