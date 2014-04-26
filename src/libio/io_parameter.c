/**
 * \file io_parameter.c
 *
 * Provides functions for reading in AMIGA parameter files.
 */


/*--- Includes ----------------------------------------------------------*/
#include <stdlib.h>
#include <assert.h>

#include "io_parameter.h"
#include "io_util.h"
#include "parse_ini.h"
#ifdef WITH_MPI
#	include "../comm.h"
#endif


/*--- Local variables ---------------------------------------------------*/
static const char *local_secName = "AHF";


/*--- Prototypes of local functions -------------------------------------*/
inline static io_parameter_t
local_allocStruct(void);

inline static parse_ini_t
local_getIni(const char *fname);

inline static void
local_initStruct(io_parameter_t params);

static void
local_readRequired(io_parameter_t params, parse_ini_t ini);

static void
local_printStruct(FILE *f, const io_parameter_t params);


/*--- Implementations of exported functions -----------------------------*/
extern io_parameter_t
io_parameter_get(char *fname)
{
	parse_ini_t ini;
	io_parameter_t params = NULL;

	if ((ini = local_getIni(fname)) != NULL) {
		params = local_allocStruct();
		local_initStruct(params);
		local_readRequired(params, ini);
		parse_ini_close(&ini);
		local_printStruct(stderr, params);
	}

	return params;
}

extern void
io_parameter_del(io_parameter_t *params)
{
	if (params == NULL)
		return;

	if (*params == NULL)
		return;

	free((*params)->icfile_name);
	free((*params)->outfile_prefix);

#ifdef DARK_ENERGY
	free((*params)->defile_name);
#endif

	free(*params);
	*params = NULL;

	return;
}

/*--- Implementations of local functions --------------------------------*/
inline static io_parameter_t
local_allocStruct(void)
{
	io_parameter_t dummy;

	dummy = (io_parameter_t)malloc(sizeof(io_parameter_struct_t));
	if (dummy == NULL) {
		fprintf(stderr, "Could not allocate memory for parameters.\n");
		return NULL;
	}

	return dummy;
}

inline static parse_ini_t
local_getIni(const char *fname)
{
	parse_ini_t ini;

	if ((fname == NULL) || (fname[0] == '\0'))
		ini = parse_ini_open(IO_PARAMETER_FNAME);
	else
		ini = parse_ini_open(fname);

	if (ini == NULL) {
		fprintf(stderr,
		        "FATAL: Could not open parameter file (%s)\n",
		        ((fname == NULL) ? IO_PARAMETER_FNAME : fname));
		return NULL;
	}

	return ini;
}

inline static void
local_initStruct(io_parameter_t params)
{
	assert(params != NULL);

	params->icfile_name    = NULL;
	params->ic_filetype    = IO_FILE_AMIGA;
	params->reader         = 1;
	params->outfile_prefix = NULL;
	params->NGRID_DOM      = -1;
	params->Nth_dom        = -1.0;
	params->Nth_ref        = -1.0;
	params->UseRhoBack     = 0;
	params->UserDvir       = -1.0;
	params->MaxGatherRad   = -1.0;
	params->lb_level       = 0;
  params->AHF_VTUNE      = 0;
  params->NGRID_MAX      = 0;
  params->AHF_MINPART    = 0;
#ifdef AHF_LRSI
	params->lrsi_beta      = 0.0;
	params->lrsi_r_s       = 0.0;
#endif
#ifdef DARK_ENERGY
	params->defile_name	   = NULL;
#endif
}

static void
local_readRequired(io_parameter_t params, parse_ini_t ini)
{
	assert(params != NULL);
	assert(ini != NULL);

	int32_t tmpInt;

	/* Get the input file name */
	getFromIni(&(params->icfile_name), parse_ini_get_string,
	           ini, "ic_filename", local_secName);

	/* Get the type of the input file */
	getFromIni(&tmpInt, parse_ini_get_int32,
	           ini, "ic_filetype", local_secName);
	params->ic_filetype = (io_file_type_t)tmpInt;

	/* Get the outfile prefix */
	getFromIni(&(params->outfile_prefix), parse_ini_get_string,
	           ini, "outfile_prefix", local_secName);

	/* Get the domain grid size */
	getFromIni(&(params->NGRID_DOM), parse_ini_get_int32,
	           ini, "LgridDomain", local_secName);

	/* Get the maximum grid size */
	getFromIni(&(params->NGRID_MAX), parse_ini_get_int32,
	           ini, "LgridMax", local_secName);
  
	/* Get the refinment criterion for the domain grid */
	getFromIni(&(params->Nth_dom), parse_ini_get_double,
	           ini, "NperDomCell", local_secName);

	/* Get the refinement criterion for the refinment grids */
	getFromIni(&(params->Nth_ref), parse_ini_get_double,
	           ini, "NperRefCell", local_secName);

	/* Shall we use rho_back or rho_crit for halo edge */
	getFromIni(&(params->UseRhoBack), parse_ini_get_int32,
	           ini, "RhoVir", local_secName);

	/* Get virial overdensity criterion */
	getFromIni(&(params->UserDvir), parse_ini_get_double,
	           ini, "Dvir", local_secName);

	/* Get the maximum gather radius */
	getFromIni(&(params->MaxGatherRad), parse_ini_get_double,
	           ini, "MaxGatherRad", local_secName);

	/* Get the tune parameter for the escape velocity */
	getFromIni(&(params->AHF_VTUNE), parse_ini_get_double,
	           ini, "VescTune", local_secName);
  
	/* Get the maximum gather radius */
	getFromIni(&(params->AHF_MINPART), parse_ini_get_int32,
	           ini, "NminPerHalo", local_secName);
  
#ifdef WITH_MPI
	/* Get the number of readers */
	getFromIni(&(params->reader), parse_ini_get_uint32,
	           ini, "NcpuReading", local_secName);
	/* Get the LOADBALANCE_DOMAIN_GRID */
	getFromIni(&(params->lb_level), parse_ini_get_int32,
	           ini, "LevelDomainDecomp", local_secName);
#else
  params->reader   = 1;
	params->lb_level = 0;
#endif

#ifdef AHF_LRSI
	/* Get LRSI specific things */
	getFromIni(&(params->lrsi_beta), parse_ini_get_double,
	           ini, "lrsi_beta", local_secName);
	getFromIni(&(params->lrsi_r_s), parse_ini_get_double,
	           ini, "lrsi_r_s", local_secName);
#endif
  
#if (defined AHFmixHaloIDandSnapID || defined SUSSING2013)
	getFromIni(&(params->isnap), parse_ini_get_uint64,
	           ini, "snapID", local_secName);
#endif
  
#ifdef DARK_ENERGY
	/* Get name of file containing the dark energy relevant tables */
	getFromIni(&(params->defile_name), parse_ini_get_string,
	           ini, "de_filename", local_secName);
#endif
  
  

  /* read additional information for some input file formats and write the respective IO_FILE.info */
  
  /*-----------------
   *     GADGET
   *-----------------*/
  if(params->ic_filetype == IO_FILE_GADGET || params->ic_filetype == IO_FILE_MGADGET)
   {
    /* Get GADGET_MUNIT */
    getFromIni(&(params->GADGET_m2Msunh), parse_ini_get_double,
               ini, "GADGET_MUNIT", "GADGET");
    /* Get GADGET_LUNIT */
    getFromIni(&(params->GADGET_l2Mpch), parse_ini_get_double,
               ini, "GADGET_LUNIT", "GADGET");
    
#ifdef VERBOSE2
    fprintf(stderr,"GADGET_LUNIT       = %g\n",params->GADGET_l2Mpch);
    fprintf(stderr,"GADGET_MUNIT       = %g\n",params->GADGET_m2Msunh);
    fprintf(stderr,"\n");
#endif
    
    // nothing to write here as GADGET_stuff will be passed on ;-)
    
   }
  
  /*-----------------
   *       ART
   *-----------------*/
  if(params->ic_filetype == IO_FILE_ART)
   {
    FILE   *fpout;
    double ART_BOXSIZE, ART_MUNIT;
    
    getFromIni(&ART_BOXSIZE, parse_ini_get_double,
               ini, "ART_BOXSIZE", "ART");
    getFromIni(&ART_MUNIT, parse_ini_get_double,
               ini, "ART_MUNIT", "ART");
    
#ifdef VERBOSE
    fprintf(stderr,"ART_BOXSIZE       = %20.10lf\n",ART_BOXSIZE);
    fprintf(stderr,"ART_MUNIT         = %20.10lf\n",ART_MUNIT);
    fprintf(stderr,"\n");
#endif

#ifdef WITH_MPI
    if(global_mpi.rank == 0) // only the master should write the file
     {
#endif
      if( (fpout = fopen("art.info","w")) == NULL)
       {
        fprintf(stderr,"Could not open art.info\nAborting\n");
        exit(0);
       }    
      fprintf(fpout,"%20.10lf \t\t ART_BOXSIZE\n",ART_BOXSIZE);
      fprintf(fpout,"%20.10lf \t\t ART_MUNIT\n",ART_MUNIT);
      fclose(fpout);
#ifdef WITH_MPI
     }
#endif
   }
  
  /*-----------------
   *      TIPSY
   *-----------------*/
  if(params->ic_filetype == IO_FILE_TIPSY)
   {
    FILE   *fpout;
    double TIPSY_BOXSIZE, TIPSY_MUNIT, TIPSY_VUNIT, TIPSY_EUNIT, TIPSY_OMEGA0, TIPSY_LAMBDA0;
    
    getFromIni(&TIPSY_BOXSIZE, parse_ini_get_double,
               ini, "TIPSY_BOXSIZE", "TIPSY");
    getFromIni(&TIPSY_MUNIT, parse_ini_get_double,
               ini, "TIPSY_MUNIT", "TIPSY");
    getFromIni(&TIPSY_VUNIT, parse_ini_get_double,
               ini, "TIPSY_VUNIT", "TIPSY");
    getFromIni(&TIPSY_EUNIT, parse_ini_get_double,
               ini, "TIPSY_EUNIT", "TIPSY");
    getFromIni(&TIPSY_OMEGA0, parse_ini_get_double,
               ini, "TIPSY_OMEGA0", "TIPSY");
    getFromIni(&TIPSY_LAMBDA0, parse_ini_get_double,
               ini, "TIPSY_LAMBDA0", "TIPSY");
    
#ifdef VERBOSE
    fprintf(stderr,"TIPSY_OMEGA0      = %20.10lf\n",TIPSY_OMEGA0);
    fprintf(stderr,"TIPSY_LAMBDA0     = %20.10lf\n",TIPSY_LAMBDA0);
    fprintf(stderr,"TIPSY_BOXSIZE     = %20.10lf\n",TIPSY_BOXSIZE);
    fprintf(stderr,"TIPSY_VUNIT       = %20.10lf\n",TIPSY_VUNIT);
    fprintf(stderr,"TIPSY_MUNIT       = %20.10lf\n",TIPSY_MUNIT);
    fprintf(stderr,"TIPSY_EUNIT       = %20.10lf\n",TIPSY_EUNIT);
    fprintf(stderr,"\n");
#endif
    
#ifdef WITH_MPI
    if(global_mpi.rank == 0) // only the master should write the file
     {
#endif
      if( (fpout = fopen("tipsy.info","w")) == NULL)
       {
        fprintf(stderr,"Could not open tipsy.info\nAborting\n");
        exit(0);
       }    
      fprintf(fpout,"%20.10lf \t\t TIPSY_OMEGA0\n",TIPSY_OMEGA0);
      fprintf(fpout,"%20.10lf \t\t TIPSY_LAMBDA0\n",TIPSY_LAMBDA0);
      fprintf(fpout,"%20.10lf \t\t TIPSY_BOXSIZE\n",TIPSY_BOXSIZE);
      fprintf(fpout,"%20.10lf \t\t TIPSY_VUNIT\n",TIPSY_VUNIT);
      fprintf(fpout,"%20.10lf \t\t TIPSY_MUNIT\n",TIPSY_MUNIT);
      fprintf(fpout,"%20.10lf \t\t TIPSY_EUNIT\n",TIPSY_EUNIT);
      fclose(fpout);
#ifdef WITH_MPI
     }
#endif

   }

  /*-----------------
   *    CUBEP3M
   *-----------------*/
  if(params->ic_filetype == IO_FILE_CUBEP3M || params->ic_filetype == IO_FILE_MCUBEP3M)
   {
    FILE   *fpout;
    double CUBEP3M_BOXSIZE, CUBEP3M_NGRID, CUBEP3M_NODES_DIM, CUBEP3M_OMEGA0, CUBEP3M_LAMBDA0;
    
    getFromIni(&CUBEP3M_BOXSIZE, parse_ini_get_double,
               ini, "CUBEP3M_BOXSIZE", "CUBEP3M");
    getFromIni(&CUBEP3M_NGRID, parse_ini_get_double,
               ini, "CUBEP3M_NGRID", "CUBEP3M");
    getFromIni(&CUBEP3M_NODES_DIM, parse_ini_get_double,
               ini, "CUBEP3M_NODES_DIM", "CUBEP3M");
    getFromIni(&CUBEP3M_OMEGA0, parse_ini_get_double,
               ini, "CUBEP3M_OMEGA0", "CUBEP3M");
    getFromIni(&CUBEP3M_LAMBDA0, parse_ini_get_double,
               ini, "CUBEP3M_LAMBDA0", "CUBEP3M");
    
#ifdef VERBOSE
    fprintf(stderr,"CUBEP3M_OMEGA0      = %20.10lf\n",CUBEP3M_OMEGA0);
    fprintf(stderr,"CUBEP3M_LAMBDA0     = %20.10lf\n",CUBEP3M_LAMBDA0);
    fprintf(stderr,"CUBEP3M_BOXSIZE     = %20.10lf\n",CUBEP3M_BOXSIZE);
    fprintf(stderr,"CUBEP3M_NGRID       = %20.10lf\n",CUBEP3M_NGRID);
    fprintf(stderr,"CUBEP3M_NODES_DIM   = %20.10lf\n",CUBEP3M_NODES_DIM);
    fprintf(stderr,"\n");
#endif
    
#ifdef WITH_MPI
    if(global_mpi.rank == 0) // only the master should write the file
     {
#endif
      if( (fpout = fopen("cubep3m.info","w")) == NULL)
       {
        fprintf(stderr,"Could not open cubep3m.info\nAborting\n");
        exit(0);
       }    
      fprintf(fpout,"%20.10lf \t\t CUBEP3M_OMEGA0\n",CUBEP3M_OMEGA0);
      fprintf(fpout,"%20.10lf \t\t CUBEP3M_LAMBDA0\n",CUBEP3M_LAMBDA0);
      fprintf(fpout,"%20.10lf \t\t CUBEP3M_BOXSIZE\n",CUBEP3M_BOXSIZE);
      fprintf(fpout,"%20.10lf \t\t CUBEP3M_NGRID\n",CUBEP3M_NGRID);
      fprintf(fpout,"%20.10lf \t\t CUBEP3M_NODES_DIM\n",CUBEP3M_NODES_DIM);
      fclose(fpout);
#ifdef WITH_MPI
     }
#endif
   }
  
} /* local_readRequired */

static void
local_printStruct(FILE *f, const io_parameter_t params)
{
#ifdef WITH_MPI
  if(global_mpi.rank == 0)
   {
#endif
    assert(f != NULL);
    assert(params != NULL);
    
    fprintf(f, "\n");
    fprintf(f, "These are the input parameters you provided, please check them carefully again:\n");
    fprintf(f, "===============================================================================\n");
    fprintf(f, "ic_filename       = %s\n", params->icfile_name);
    fprintf(f, "ic_filetype       = %i\n", (int)(params->ic_filetype));
    fprintf(f, "outfile_prefix    = %s\n", params->outfile_prefix);
//    fprintf(f, "LgridDomain       = %i\n", params->NGRID_DOM);
//    fprintf(f, "LgridMax          = %i\n", params->NGRID_MAX);
//    fprintf(f, "NminPerHalo       = %i\n", params->AHF_MINPART);
//    fprintf(f, "VescTune          = %g\n", params->AHF_VTUNE);
//    fprintf(f, "NperDomCell       = %g\n", params->Nth_dom);
//    fprintf(f, "NperRefCell       = %g\n", params->Nth_ref);
//    fprintf(f, "RhoVir            = %i\n", params->UseRhoBack);
//    fprintf(f, "Dvir              = %g\n", params->UserDvir);
//    fprintf(f, "MaxGatherRad      = %g Mpc/h\n", params->MaxGatherRad);
#ifdef WITH_MPI
//    fprintf(f, "LevelDomainDecomp = %i\n", params->lb_level);
//    fprintf(f, "NcpuReading       = %" PRIu32 "\n", params->reader);
#endif
#ifdef AHF_LRSI
//    fprintf(f, "lrsi_beta         = %g\n", params->lrsi_beta);
//    fprintf(f, "lrsi_r_s          = %g\n", params->lrsi_r_s);
#endif
#ifdef DARK_ENERGY
//    fprintf(f, "de_filename       = %s\n", params->defile_name);
#endif
//    fprintf(f, "\n");
#ifdef WITH_MPI
   }
#endif
}
