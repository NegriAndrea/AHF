/* Belaid & Doug Nov. 13, 2018 */

/**
 * \file io_gizmo.c
 *
 * Provides functions for reading and writing Gizmo files.
 */

#ifdef WITH_HDF5

/**********************************************************************\
 *    Includes                                                        *
 \**********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif
#include <limits.h>

#include "io_gizmo.h"
#include "io_gizmo_header.h"
#include "io_util.h"

#include "../define.h"
#include "../common.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               *
 \**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   *
 \**********************************************************************/

/**
 * \brief Helper function to open the file
 *
 * This function is supposed to get inlined anyway. Makes
 * io_gizmo_open more readable.
 *
 * \param log   The logging object.
 * \param f     The Gizmo file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_gizmo_t
local_openopen(io_logging_t log, io_gizmo_t f, io_file_mode_t mode);

/**
 * \brief Helper funtion to set the swap-state and write some
 *        messages
 *
 * \param log      The logging object.
 * \param f        The Gizmo file object sofar.
 * \param swapped  The swap state.
 *
 * \return Nothing.
 */
inline static void
local_openswapped(io_logging_t log,
                  io_gizmo_t f,
                  io_file_swap_t swapped);

/**
 * \brief Try to find out swapping status.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Returns the file object or NULL in case of an error.
 */
inline static io_gizmo_t
local_opengetswap(io_logging_t log, io_gizmo_t f);

/**
 * \brief Tries to figure out which Gizmo file version is to be used.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Nothing.
 */
inline static void
local_openversion(io_logging_t log, io_gizmo_t f);

/**
 * \brief Check before writing particles and place file pointer at the
 *        right spot.
 *
 * \param log
 * \param f
 * \param pskip
 * \param bytes
 *
 * \return
 */
inline static int32_t
local_write_common(io_logging_t log,
                   io_gizmo_t f,
                   uint64_t pskip,
                   int32_t bytes);

static uint64_t
local_get_block_pos(io_logging_t log,
                    io_gizmo_t f,
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg,
                    hid_t hdf5_grp[]);

static uint64_t
local_get_block_vel(io_logging_t log,
                    io_gizmo_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg,
                    hid_t hdf5_grp[]);

static uint64_t
local_get_block_id(io_logging_t log,
                   io_gizmo_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg,
                   hid_t hdf5_grp[]);

static uint64_t
local_get_block_mass(io_logging_t log,
                     io_gizmo_t f,
                     uint64_t pskip,
                     uint64_t pread,
                     io_file_strg_struct_t strg,
                     hid_t hdf5_grp[]);

static uint64_t
local_get_block_u(io_logging_t log,
                  io_gizmo_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg,
                  hid_t hdf5_grp[]);

void local_find_block(io_gizmo_t, char *);

#ifdef METALHACK
static uint64_t
local_get_block_z(io_logging_t log,
                  io_gizmo_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg);

static uint64_t
local_get_block_age(io_logging_t log,
                    io_gizmo_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg);

static void
local_skip_GIZMO1_blocks(FILE *, int, int);
#endif



/**********************************************************************\
 *    Implementation of global functions                              *
 \**********************************************************************/
extern io_gizmo_t
io_gizmo_open(io_logging_t log,
              char *fname,
              io_file_swap_t swapped,
              io_file_mode_t mode,
              uint32_t reader)
{
  io_gizmo_t f;
  
  /* Get memory for the structure */
  f = (io_gizmo_t)malloc(sizeof(io_gizmo_struct_t));
  if (f == NULL) {
    io_logging_memfatal(log,  "io_gizmo structure");
    return NULL;
  }
  
  /* Start filling the structure */
  
  /* Store the filename */
  f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
  if (f->fname == NULL) {
    io_logging_memfatal(log, "filename of GizmoFile");
    free(f);
    return NULL;
  }
  strncpy(f->fname, fname, strlen(fname)+1);
  
  /* Okay, we are a Gizmo file */
  f->ftype = IO_FILE_GIZMO;
  
  /* And we can just copy in the parallel information */
#	ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
  MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
  if (f->size >= reader) {
    /* TODO
     * THIS IS JUST A QUICK HACK TO PREVENT MGIZMO FILES TO
     * AGAIN TRY TO SPLIT THE COMMUNICATOR, THAT IS ALREADY DONE
     * IN mgizmo.c
     * TODO
     */
    MPI_Comm_split(MPI_COMM_WORLD, 1, f->rank, &(f->mycomm));
    MPI_Comm_size(f->mycomm, &(f->size_mycomm));
    MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
  } else {
    f->mycomm = MPI_COMM_NULL;
    f->size_mycomm = -1;
    f->rank_mycomm = -1;
  }
  
#	endif
  
  /* Try to open the file and set mode */
  if (local_openopen(log, f, mode) == NULL) {
    free(f->fname);
    free(f);
    return NULL;
  }
  
  /* Set swapping */
  /*local_openswapped(log, f, swapped);
   if (    (f->mode == IO_FILE_READ)
   && (f->swapped == IO_FILE_UNKOWN_SWAPPING)) {
   if (local_opengetswap(log, f) != f) {
			io_logging_fatal(log, "Cannot open this file.");
			free(f->fname);
			free(f);
			return NULL;
   }
   }*/
  
  /* Identify Gizmo format */
  local_openversion(log, f);
  
  /* Nothing for the header for now */
  f->header = NULL;
  
  /* Set some dummy values */
  f->no_part           = UINT64_C(0);
  f->no_part_with_mass = UINT64_C(0);
  f->multimass         = INT8_C(0);
  f->mmass             = 1e40;
  f->minweight         = 1e40;
  f->maxweight         = 0.0;
  f->sumweight         = 0.0;
  f->no_species        = INT32_C(0);
  f->posscale          = 1.0;
  f->weightscale       = 1.0;
  
  return f;
}

extern void
io_gizmo_close(io_logging_t log,
               io_gizmo_t *f)
{
  /* Catch NULLs */
  if (f == NULL || *f == NULL)
    return;
  
  /* Put header to the file if necessary */
  if (    ((*f)->mode == IO_FILE_WRITE)
	     && ((*f)->header != NULL)) {
    io_gizmo_header_write(log, (*f)->header, *f);
  }
  
  /* Close */
  if ((*f)->header != NULL)
    io_gizmo_header_del(log, &((*f)->header));
  if ((*f)->fname != NULL)
    free((*f)->fname);
#	ifdef WITH_MPI
  if ((*f)->mycomm != MPI_COMM_NULL)
    MPI_Comm_free(&((*f)->mycomm));
#	endif
  
  /* Actually close the file */
  if ((*f)->file > 0)
  {
    H5Fclose((*f)->file);
  }
  
  /* Cleaning */
  free(*f);
  *f = NULL;
  
  return;
}

extern void
io_gizmo_init(io_logging_t log,
              io_gizmo_t f)
{
  if (f == NULL)
    return;
  
  if (f->header != NULL) {
    io_logging_warn(log, INT32_C(1),
                    "Already have the header information! Rereading.");
    io_gizmo_header_del(log, &(f->header));
  }
  
  if (f->mode != IO_FILE_READ) {
    io_logging_warn(log, INT32_C(1),
                    "%s is not opened for reading. "
                    "Will do nothing.",
                    f->fname);
    return;
  }
  
  io_logging_msg(log, INT32_C(5),
                 "Starting to initialize file object from %s",
                 f->fname);
  f->header = io_gizmo_header_get(log, f);
  io_logging_msg(log, INT32_C(5),
                 "Done with initializing file object from %s",
                 f->fname);
  
  /* Check for multimass file and also sum up the particle count */
  f->multimass = 0;
  {
    int i;
    for (i=0; i<6; i++) {
      f->no_part += (uint64_t)(f->header->np[i]);
      if (    islessequal(f->header->massarr[i], 0.0)
          && (f->header->np[i]>0)) {
        f->multimass |= (1<<i);
        io_logging_msg(log, INT32_C(5),
                       "Particle type %d requires reading "
                       "of a mass array.", i);
        f->no_part_with_mass += (uint64_t)(f->header->np[i]);
      }
    }
  }
  
  return;
}

extern uint64_t
io_gizmo_readpart(io_logging_t log,
                  io_gizmo_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg)
{
  uint64_t particles_read, tmp;
  double box[3], shift[3];
  double scale_pos, scale_mom, scale_weight;
  
  /*
   * First read the particles unscaled. This will set important
   * scaling information
   */
  particles_read = io_gizmo_readpart_raw(log, f, pskip, pread, strg);
  
  if (particles_read != pread) {
    return UINT64_C(0);
  }
  
  /* And do the scaling */
#ifdef WITH_MPI
  io_gizmo_scale_global(log, f->mycomm,  f->maxpos, f->minpos, &(f->mmass));
#endif
  tmp = io_gizmo_scale_particles(log, f->maxpos, f->minpos,
	                                &(f->header->boxsize),
	                                f->header->expansion,
	                                f->posscale, f->mmass,
	                                particles_read, strg);
  if (tmp != particles_read) {
    return tmp;
  }
  
  /* Wow, we are done! */
  return particles_read;
}

extern uint64_t
io_gizmo_readpart_raw(io_logging_t log,
                      io_gizmo_t f,
                      uint64_t pskip,
                      uint64_t pread,
                      io_file_strg_struct_t strg)
{
  long skipsize;
  char str[5], buf[100];
  uint32_t blocksize, blocksize2;
  uint64_t funcrtn;
  uint32_t nextblocksize;
  int type;
  int rank, pcsum;
  hid_t hdf5_grp[6], hdf5_dataspace_in_file;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_memory, hdf5_dataset;
  hsize_t dims[2], count[2], start[2];
  
#ifdef FOPENCLOSE
  //fprintf(stderr,"FOPENCLOSE: io_gizmo_readpart_raw() opening %s ... ",f->fname);
  //f->file = fopen(f->fname,IO_FILE_MODE_READ);
  f->file = H5Fopen(f->fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (f->file < 0) {
    io_logging_fatal(log,"io_gizmo_readpart_raw(): could not open file %s for reading",f->fname);
    return UINT64_C(0);
  }
#endif
  
  /* Check if we actually have to do something */
  if ( (f == NULL) || (f->header == NULL) )
    return UINT64_C(0);
  
  for (type = 0; type < 6; type++)
  {
    if (f->header->np[type] > 0)
    {
      sprintf(buf, "/PartType%d", type);
      hdf5_grp[type] = H5Gopen(f->file, buf);
    }
  }
  
  /* Coordinates */
  funcrtn = local_get_block_pos(log, f, &pskip, &pread, strg, hdf5_grp);
  if (funcrtn != pread) {
    io_logging_fatal(log,
                     "Expected to read %"PRIu64
                     " particle positions, but only got %"PRIu64
                     ".  Aborting.");
    return UINT64_C(0);
  }
  
  /* Velocities */
  funcrtn = local_get_block_vel(log, f, pskip, pread, strg, hdf5_grp);
  if (funcrtn != pread) {
    io_logging_fatal(log,
                     "Expected to read %"PRIu64
                     " particle velocities, but only got %"PRIu64
                     ".  Aborting.");
    return UINT64_C(0);
  }
  
  /* Identities */
  funcrtn = local_get_block_id(log, f, pskip, pread, strg, hdf5_grp);
  if (funcrtn != pread) {
    io_logging_fatal(log,
                     "Expected to read %"PRIu64
                     " particle identities, but only got %"PRIu64
                     ".  Aborting.");
    return UINT64_C(0);
  }
  
  /* Mass (might not exist!) */
  funcrtn = local_get_block_mass(log, f, pskip, pread, strg, hdf5_grp);
  if (funcrtn != pread) {
    io_logging_fatal(log,
                     "Expected to obtain %"PRIu64
                     " particle masses, but only got %"PRIu64
                     ".  Aborting.");
    return UINT64_C(0);
  }
  
  funcrtn = local_get_block_u(log, f, pskip, pread, strg, hdf5_grp);
  if (funcrtn != pread) {
    io_logging_fatal(log,
                     "Expected to read %"PRIu64
                     " particle identities, but only got %"PRIu64
                     ".  Aborting.");
    return UINT64_C(0);
  }
  
  /* @TODO: D. Rennehan: not implemented yet! */
#	ifdef METALHACK
  local_get_block_age(log, f, pskip, pread, strg);
  local_get_block_z(log, f, pskip, pread, strg);
#	endif
  
#ifdef FOPENCLOSE
  for (type = 5; type >= 0; type--)
  {
    if (f->header->np[type] > 0)
    {
      H5Gclose(hdf5_grp[type]);
    }
  }
  
  H5Fclose(f->file);
  
  f->file = -1; // we are not going to read any more from the file and hence can indicate to close it by setting f->file=NULL
#endif
  
  /* Return the number of particles read */
  return pread;
}


extern uint64_t
io_gizmo_writepart(io_logging_t log,
                   io_gizmo_t f,
                   uint64_t pskip,
                   uint64_t pwrite,
                   io_file_strg_struct_t strg)
{
  return 0;
}

extern uint64_t
io_gizmo_writepart_ord(io_logging_t log,
                       io_gizmo_t f,
                       uint64_t pskip,
                       uint64_t pwrite,
                       void *nxt_part,
                       io_file_strg_struct_t strg)
{
  return UINT64_C(0);
}

extern bool
io_gizmo_get(io_logging_t log,
             io_gizmo_t f,
             io_file_get_t what,
             void *res)
{
  if ( (f == NULL) || (f->header == NULL) )
    return false;
  
  switch (what) {
    case IO_FILE_GET_NOPART_IN_FILE:
    case IO_FILE_GET_NOPART:
      *((long *)res) = (long)f->no_part;
      break;
    case IO_FILE_GET_NOVPART:
      if (f->no_part > UINT64_C(0)) {
        if (    isgreater(f->sumweight, 0.0)
            && isgreater(f->mmass, 0.0))
          *((double *)res) = f->sumweight / f->mmass;
        else
          *((double *)res) = (double)f->no_part;
      } else {
        io_logging_warn(log, INT32_C(0),
                        "Cannot calculate novpart yet. "
                        "You first need to read the "
                        "particles.");
        return false;
      }
      break;
    case IO_FILE_GET_NOSPECIES:
      if (f->multimass) {
        if (f->no_species >= 1)
          *((int *)res) = f->no_species;
        else {
          io_logging_warn(log, INT32_C(0),
                          "Don't know the number of particle "
                          "species yet. You first need to "
                          "read the particles.");
          return false;
        }
      } else {
        *((int *)res) = 1;
      }
      break;
    case IO_FILE_GET_BOXSIZE:
      *((double *)res) =   f->header->boxsize * f->posscale;
      break;
    case IO_FILE_GET_PMASS:
      if (f->multimass) {
        *((double *)res) = f->mmass * f->weightscale;
      } else {
        *((double *)res) =   f->header->massarr[1] * f->weightscale;
      }
      break;
    case IO_FILE_GET_ZINITIAL:
      io_logging_warn(log, INT32_C(1),
                      "zinitial is not set in a Gizmo file, "
                      "using current redshift");
    case IO_FILE_GET_Z:
      *((double *)res) = f->header->redshift;
      break;
    case IO_FILE_GET_AINITIAL:
      io_logging_warn(log, INT32_C(1),
                      "ainitial is not set in a Gizmo file, "
                      "using current expansion.");
    case IO_FILE_GET_A:
      *((double *)res) = f->header->expansion;
      break;
    case IO_FILE_GET_OMEGA0:
      *((double *)res) = f->header->omega0;
      break;
    case IO_FILE_GET_OMEGAL:
      *((double *)res) = f->header->omegalambda;
      break;
    case IO_FILE_GET_H:
      *((double *)res) = f->header->hubbleparameter;
      break;
    case IO_FILE_GET_DOUBLE:
      *((int *)res) = f->header->flagdoubleprecision;
      break;
    case IO_FILE_GET_MMASS:
      if (isgreater(f->header->massarr[1], 0.0))
        *((int *)res) = 0;
      else
        *((int *)res) = 1;
      break;
    case IO_FILE_GET_NOTSTEP:
      io_logging_warn(log, INT32_C(1),
                      "Gizmo files don't store the step number. "
                      "Setting to 0.");
      *((int32_t *)res) = 0;
      break;
    case IO_FILE_GET_TSTEP:
      io_logging_warn(log, INT32_C(1),
                      "Gizmo files don't store the timestep. "
                      "Setting to 0.0");
      *((double *)res) = 0.0;
      break;
    case IO_FILE_GET_HEADERSTR:
      io_logging_warn(log, INT32_C(1),
                      "Gizmo files don't have a header string. "
                      "Using a dummy one.");
      *((char **)res) = "No header string.";
      break;
    case IO_FILE_GET_MINWEIGHT:
      if (isgreater(f->mmass, 0.0))
        *((double *)res) = f->minweight / f->mmass;
      else {
        io_logging_warn(log, INT32_C(1),
                        "Don't know minweight yet, setting to "
                        "0.0.");
        *((double *)res) = 0.0;
      }
      break;
    case IO_FILE_GET_MAXWEIGHT:
      if (isgreater(f->mmass, 0.0))
        *((double *)res) = f->maxweight / f->mmass;
      else {
        io_logging_warn(log, INT32_C(1),
                        "Don't know maxweight yet, setting to "
                        "0.0.");
        *((double *)res) = 0.0;
      }
      break;
    default:
      io_logging_fatal(log, "Requesting something unkown in %s.",
                       __func__);
      return false;
  }
  
  return true;
}

extern bool
io_gizmo_set(io_logging_t log,
             io_gizmo_t f,
             io_file_get_t what,
             void *res)
{
  if ( (f == NULL) || (f->header == NULL) )
    return false;
  
  switch (what) {
    case IO_FILE_GET_BOXSIZE:
      f->header->boxsize = *((double *)res);
      break;
    case IO_FILE_GET_PMASS:
      f->header->massarr[1] = *((double *)res);
      break;
    case IO_FILE_GET_Z:
      f->header->redshift = *((double *)res);
      break;
    case IO_FILE_GET_A:
      f->header->expansion = *((double *)res);
      break;
    case IO_FILE_GET_OMEGA0:
      f->header->omega0 = *((double *)res);
      break;
    case IO_FILE_GET_OMEGAL:
      f->header->omegalambda = *((double *)res);
      break;
    case IO_FILE_GET_H:
      f->header->hubbleparameter = *((double *)res);
      break;
    default:
      io_logging_fatal(log, "Requesting something unkown in %s.",
                       __func__);
      return false;
  }
  
  return true;
}

extern void
io_gizmo_log(io_logging_t log, io_gizmo_t f)
{
  io_logging_msg(log, INT32_C(5),
                 "Fileobject information:");
  io_logging_msg(log, INT32_C(5),
                 "  Filetype:             %s",
                 io_file_typestr(f->ftype));
  io_logging_msg(log, INT32_C(5),
                 "  Filename:             %s",
                 f->fname);
  io_logging_msg(log, INT32_C(5),
                 "  Mode:                 %" PRIi8,
                 f->mode);
  io_logging_msg(log, INT32_C(5),
                 "  Swapping:             %" PRIi8,
                 f->swapped);
  io_logging_msg(log, INT32_C(5),
                 "  File version:         %" PRIi8,
                 f->ver);
  io_logging_msg(log, INT32_C(5),
                 "  Header size:          %" PRIi8,
                 f->ver);
  io_logging_msg(log, INT32_C(5),
                 "  No. particles:        %" PRIu64,
                 f->no_part);
  io_logging_msg(log, INT32_C(5),
                 "  No. particles w/mass: %" PRIu64,
                 f->no_part_with_mass);
  io_logging_msg(log, INT32_C(5),
                 "  Multimass:            %" PRIi8,
                 f->multimass);
  {
    int i;
    for (i=0; i<6; i++)
      io_logging_msg(log, INT32_C(5),
                     "      Part type %d:      %" PRIi8, i,
                     ((int8_t)(f->multimass) & (1<<i)) >> i);
  }
  io_logging_msg(log, INT32_C(5),
                 "  MMass (Halo parts):   %g",
                 f->mmass);
  io_logging_msg(log, INT32_C(5),
                 "  Minimal Weight:       %g",
                 f->minweight);
  io_logging_msg(log, INT32_C(5),
                 "  Maximal Weight:       %g",
                 f->maxweight);
  io_logging_msg(log, INT32_C(5),
                 "  Sum of all weights:   %g",
                 f->sumweight);
  io_logging_msg(log, INT32_C(5),
                 "  No. of species:       %" PRIi32,
                 f->no_species);
  io_logging_msg(log, INT32_C(5),
                 "  Position scale:       %g",
                 f->posscale);
  io_logging_msg(log, INT32_C(5),
                 "  Weight scale:         %g",
                 f->weightscale);
  io_gizmo_header_log(log, f->header);
  
  return;
}

extern void
io_gizmo_resetscale(io_logging_t log,
                    io_gizmo_t f,
                    double posscale,
                    double weightscale) {
  if (f == NULL)
    return;
  
  io_logging_msg(log, INT32_C(8),
                 "Old posscale: %g   New posscale: %g",
                 f->posscale, posscale);
  io_logging_msg(log, INT32_C(8),
                 "Old weightscale: %g   New weightscale: %g",
                 f->weightscale, weightscale);
  f->posscale = posscale;
  f->weightscale = weightscale;
  
  return;
}

extern uint64_t
io_gizmo_scale_particles(io_logging_t log,
                         double maxpos[],
                         double minpos[],
                         double *boxsize,
                         double expansion,
                         double posscale,
                         double mmass,
                         uint64_t particles_read,
                         io_file_strg_struct_t strg)
{
  double box[3], shift[3];
  double scale_pos, scale_mom, scale_weight, scale_u;
  uint64_t i;
  
  /* Now we can do the scaling */
  box[0] = fabs(maxpos[0] - minpos[0]);
  box[1] = fabs(maxpos[1] - minpos[1]);
  box[2] = fabs(maxpos[2] - minpos[2]);
  if (isgreater(box[0], *boxsize)) {
    io_logging_warn(log, INT32_C(1),
                    "x-Separation of particles exceeds boxsize "
                    "(%g > %g), resetting boxsize.",
                    box[0], *boxsize);
    *boxsize = box[0];
  }
  if (isgreater(box[1], *boxsize)) {
    io_logging_warn(log, INT32_C(1),
                    "y-Separation of particles exceeds boxsize "
                    "(%g > %g), resetting boxsize.",
                    box[1], *boxsize);
    *boxsize = box[1];
  }
  if (isgreater(box[2], *boxsize)) {
    io_logging_warn(log, INT32_C(1),
                    "z-Separation of particles exceeds boxsize "
                    "(%g > %g), resetting boxsize.",
                    box[2], *boxsize);
    *boxsize = box[2];
  }
  io_logging_msg(log, INT32_C(4),
                 "Extreme positions: xmin = %g  xmax = %g",
                 minpos[0], maxpos[0]);
  io_logging_msg(log, INT32_C(4),
                 "                   ymin = %g  ymax = %g",
                 minpos[1], maxpos[1]);
  io_logging_msg(log, INT32_C(4),
                 "                   zmin = %g  zmax = %g",
                 minpos[2], maxpos[2]);
  shift[0] = (isless(minpos[0], 0.0) ? -minpos[0] : 0.0);
  shift[1] = (isless(minpos[1], 0.0) ? -minpos[1] : 0.0);
  shift[2] = (isless(minpos[2], 0.0) ? -minpos[2] : 0.0);
  io_logging_msg(log, INT32_C(4),
                 "Applying shift: (%g, %g, %g)",
                 shift[0], shift[1], shift[2]);
  
  
  
  /*================================================================
   *             conversion to AHF internal units
   *
   *    NOTE: the thermal energy is scaled to (km/sec)^2 !!!!!
   *================================================================*/
  scale_pos    = 1.0/(*boxsize);
#ifdef NO_EXPANSION
  scale_mom    = 1.0 / (*boxsize * posscale * 100.);
#else
  scale_mom    = sqrt(expansion) * expansion / (*boxsize * posscale * 100.);
#endif
  scale_weight = 1.0/mmass;
  scale_u      = 1.0; // already in (km/sec)^2 !?
  
  
  
  io_logging_msg(log, INT32_C(3),  "Scaling by:  positions:  %g", scale_pos);
  io_logging_msg(log, INT32_C(3),  "             velocities: %g", scale_mom);
  io_logging_msg(log, INT32_C(3),  "             weights:    %g", scale_weight);
  io_logging_msg(log, INT32_C(3),  "             u:          %g", scale_u);
  
  
  /* keep track of the applied shifts and scales */
  simu.pos_shift[0] = shift[0];
  simu.pos_shift[1] = shift[1];
  simu.pos_shift[2] = shift[2];
  simu.pos_scale    = scale_pos;
  
  
  
  
  /* Define the actual scaling calls type independent */
#	define SCALE_CALL(type) {\
*((type *)strg.posx.val) += (type)(shift[0]); \
*((type *)strg.posx.val) *= (type)(scale_pos); \
*((type *)strg.posy.val) += (type)(shift[1]); \
*((type *)strg.posy.val) *= (type)(scale_pos); \
*((type *)strg.posz.val) += (type)(shift[2]); \
*((type *)strg.posz.val) *= (type)(scale_pos); \
*((type *)strg.momx.val) *= (type)(scale_mom); \
*((type *)strg.momy.val) *= (type)(scale_mom); \
*((type *)strg.momz.val) *= (type)(scale_mom); \
if (strg.weight.val != NULL) *((type *)strg.weight.val) *= (type)(scale_weight); \
if (strg.u.val != NULL) {if(*((type *)strg.u.val) >= 0.) *((type *)strg.u.val) *= (type)(scale_u);} \
strg.posx.val = (void *)(((char *)strg.posx.val) + strg.posx.stride); \
strg.posy.val = (void *)(((char *)strg.posy.val) + strg.posy.stride); \
strg.posz.val = (void *)(((char *)strg.posz.val) + strg.posz.stride); \
strg.momx.val = (void *)(((char *)strg.momx.val) + strg.momx.stride); \
strg.momy.val = (void *)(((char *)strg.momy.val) + strg.momy.stride); \
strg.momz.val = (void *)(((char *)strg.momz.val) + strg.momz.stride); \
if (strg.weight.val != NULL) strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride); \
if (strg.u.val != NULL) strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride); \
}
  
  /* Do the scaling depending on the storage type */
  if (strg.bytes_float == sizeof(float)) {
    for (i=0; i<particles_read; i++) {
      SCALE_CALL(float);
    }
  } else if (strg.bytes_float == sizeof(double)) {
    for (i=0; i<particles_read; i++) {
      SCALE_CALL(double);
    }
  } else if (strg.bytes_float == sizeof(long double)) {
    for (i=0; i<particles_read; i++) {
      SCALE_CALL(long double);
    }
  } else {
    io_logging_fatal(log,
                     "Don't know which floating point types "
                     "has %" PRIi32 " bytes. Aborting read.",
                     strg.bytes_float);
    return UINT64_C(0);
  }
  
  /* Clean the macro away */
#	undef SCALE_CALL
  
  /* And we are done */
  return particles_read;
}

#ifdef WITH_MPI
extern void
io_gizmo_scale_global(io_logging_t log,
                      MPI_Comm comm,
                      double *maxpos,
                      double *minpos,
                      double *mmass)
{
  int size, rank;
  double buffer[3];
  
  io_logging_msg(log, INT32_C(5),
                 "Updating local scale values to global values.");
  MPI_Allreduce((void *)maxpos, (void *)buffer, 3,
                MPI_DOUBLE, MPI_MAX, comm);
  io_logging_msg(log, INT32_C(5),
                 "local : maxpos[0] = %g \t"
                 "maxpos[1] = %g \t"
                 "maxpos[2] = %g",
                 maxpos[0], maxpos[1], maxpos[2]);
  maxpos[0] = buffer[0];
  maxpos[1] = buffer[1];
  maxpos[2] = buffer[2];
  io_logging_msg(log, INT32_C(5),
                 "global: maxpos[0] = %g \t"
                 "maxpos[1] = %g \t"
                 "maxpos[2] = %g",
                 maxpos[0], maxpos[1], maxpos[2]);
  
  MPI_Allreduce((void *)minpos, (void *)buffer, 3,
                MPI_DOUBLE, MPI_MIN, comm);
  io_logging_msg(log, INT32_C(5),
                 "local : minpos[0] = %g \t"
                 "minpos[1] = %g \t"
                 "minpos[2] = %g",
                 minpos[0], minpos[1], minpos[2]);
  minpos[0] = buffer[0];
  minpos[1] = buffer[1];
  minpos[2] = buffer[2];
  io_logging_msg(log, INT32_C(5),
                 "global: minpos[0] = %g \t"
                 "minpos[1] = %g \t"
                 "minpos[2] = %g",
                 minpos[0], minpos[1], minpos[2]);
  
  MPI_Allreduce((void *)mmass, (void *)buffer, 1,
                MPI_DOUBLE, MPI_MIN, comm);
  io_logging_msg(log, INT32_C(5), "local : mmass = %g", *mmass);
  *mmass = buffer[0];
  io_logging_msg(log, INT32_C(5), "global: mmass = %g", *mmass);
  
  return;
}
#endif /* WITH_MPI */


/**********************************************************************\
 *    Implementation of local functions                               *
 \**********************************************************************/

inline static io_gizmo_t
local_openopen(io_logging_t log, io_gizmo_t f, io_file_mode_t mode)
{
  if (mode == IO_FILE_READ) {
    //f->file = fopen(f->fname, IO_FILE_MODE_READ);
    f->file = H5Fopen(f->fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f->file < 0) {
      io_logging_fatal(log, "Could not open '%s' for reading.", f->fname);
      return NULL;
    }
  } else {
    io_logging_fatal(log, "No write implementation for GIZMO HDF5.");
    return NULL;
    //f->file = fopen(f->fname, IO_FILE_MODE_WRITE);
    //if (f->file == NULL) {
    //	io_logging_fatal(log, "Could not open '%s' for writing.", f->fname);
    //	return NULL;
    //}
  }
  
  f->mode = mode;
  
  return f;
}

inline static void
local_openswapped(io_logging_t log,
                  io_gizmo_t f,
                  io_file_swap_t swapped)
{
  if (f->mode == IO_FILE_READ) {
    switch(swapped) {
      case IO_FILE_ISNOT_SWAPPED:
        io_logging_msg(log, INT32_C(3),
                       "Assuming unswapped file");
        break;
      case IO_FILE_IS_SWAPPED:
        io_logging_msg(log, INT32_C(3),
                       "Assuming swapped file");
        break;
      case IO_FILE_UNKOWN_SWAPPING:
      default:
        io_logging_msg(log, INT32_C(3),
                       "Will try to find out swap status");
    }
  }
  
  /* Now set the swapping */
  f->swapped = swapped;
  
  return;
}

inline static void
local_openversion(io_logging_t log, io_gizmo_t f)
{
  f->ver = 1;
  io_logging_msg(log, INT32_C(2), "Assuming a public GIZMO v1 file.");
  
  return;
}

inline static int32_t
local_write_common(io_logging_t log,
                   io_gizmo_t f,
                   uint64_t pskip,
                   int32_t bytes)
{
  return 0;
}

#define CHECK_FLOATBYTES(bfile, bstore) {\
if (bfile == sizeof(float)) { \
io_logging_msg(log, INT32_C(1), \
"Obviously the file uses float for " \
"floating point values (%" PRIi32 " bytes).", \
bfile); \
} else if (bfile == sizeof(double)) { \
io_logging_msg(log, INT32_C(1), \
"Obviously the file uses double for " \
"floating point values (%" PRIi32 " bytes).", \
bfile); \
} else { \
io_logging_fatal(log, \
"No clue what kind of floating point uses " \
"%" PRIi32 " bytes. Aborting reading.", \
bfile); \
return UINT64_C(0); \
}\
if (bfile < bstore) { \
io_logging_msg(log, INT32_C(1), \
"The floating point values in the file have " \
"less precision than the particle storage " \
"(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
"No problem, will upcast.", \
bfile, bstore); \
} else if (bfile > bstore) { \
io_logging_warn(log, INT32_C(1), \
"The floating point values in the file have " \
"a higher precision than the particle storage " \
"(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
"Will downcast, but precision might be lost.", \
bfile, bstore); \
} \
}
static uint64_t
local_get_block_pos(io_logging_t log,
                    io_gizmo_t f,
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg,
                    hid_t hdf5_grp[])
{
  uint32_t blocksize, blocksize2;
  int32_t partsize;
  uint32_t bytes_file;
  double fposx, fposy, fposz;
  float dummy;
  uint64_t i;
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  int num_bytes = 0;
  int num_elements = 3; // 3 Coordinates: x,y,z
  
  /* Set extreme position detectors */
  f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e40;
  f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e40;
  
  io_logging_msg(log, INT32_C(2), "local_get_block_pos(): loop over particle types");
  
  for (type = 0; type < 6; type ++)
  {
    io_logging_msg(log, INT32_C(2), "local_get_block_pos(): read %d of part type %d ", f->header->np[type], type);
    
    if (f->header->np[type] == 0)
    {
      io_logging_msg(log, INT32_C(2), "local_get_block_pos(): no particles of type %d, skip! ", type);
      continue;
    }
    
    if (f->header->flagdoubleprecision)
    {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
      num_bytes = f->header->np[type] * num_elements * sizeof(double);
      CommBuffer = (double *)malloc(num_bytes);
    } else {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
      num_bytes = f->header->np[type] * num_elements * sizeof(float);
      CommBuffer = (float *)malloc(num_bytes);
    }
    
    if (CommBuffer == NULL)
    {
      io_logging_fatal(log, "local_get_block_pos(): could not allocate %d bytes", (int)num_bytes);
      return UINT64_C(0);
    }
    
    io_logging_msg(log, INT32_C(2), "local_get_block_pos(): part type %d: allocated %d bytes for CommBuffer", type, (int)num_bytes);
    
    io_util_readhdf5(log, f, "Coordinates", type, num_elements, hdf5_datatype, hdf5_grp, CommBuffer);
    
    io_logging_msg(log, INT32_C(2), "local_get_block_pos(): part type %d: read particles into CommBuffer", type);
    
    /* Loop over the particle positions */
    for (i = 0; i < num_elements * f->header->np[type]; i += num_elements) {
      //fposx = *((double *)CommBuffer)++;
      //fposy = *((double *)CommBuffer)++;
      //fposz = *((double *)CommBuffer)++;
      
      if (f->header->flagdoubleprecision)
      {
        fposx = ((double *)CommBuffer)[i];
        fposy = ((double *)CommBuffer)[i + 1];
        fposz = ((double *)CommBuffer)[i + 2];
      } else {
        fposx = (double)((float *)CommBuffer)[i];
        fposy = (double)((float *)CommBuffer)[i + 1];
        fposz = (double)((float *)CommBuffer)[i + 2];
      }
      
      //printf("index: %d\n", i);
      
      /* STEP 2:  Detect extreme positions */
      if (isless(fposx, f->minpos[0]))
        f->minpos[0] = fposx;
      if (isless(fposy, f->minpos[1]))
        f->minpos[1] = fposy;
      if (isless(fposz, f->minpos[2]))
        f->minpos[2] = fposz;
      if (isgreater(fposx, f->maxpos[0]))
        f->maxpos[0] = fposx;
      if (isgreater(fposy, f->maxpos[1]))
        f->maxpos[1] = fposy;
      if (isgreater(fposz, f->maxpos[2]))
        f->maxpos[2] = fposz;
      
      /* STEP 3:  Store the particle in the array */
      if (strg.bytes_float != sizeof(float))
      {
        *((double *)strg.posx.val) = fposx;
        *((double *)strg.posy.val) = fposy;
        *((double *)strg.posz.val) = fposz;
      } else {
        *((float *)strg.posx.val) = (float)fposx;
        *((float *)strg.posy.val) = (float)fposy;
        *((float *)strg.posz.val) = (float)fposz;
      }
      
      /* STEP 4:  Increment the pointers to the next particle */
      strg.posx.val = (void *)(((char *)strg.posx.val)
                               + strg.posx.stride);
      strg.posy.val = (void *)(((char *)strg.posy.val)
                               + strg.posy.stride);
      strg.posz.val = (void *)(((char *)strg.posz.val)
                               + strg.posz.stride);
    } /* End of particle position loop */
    
    io_logging_msg(log, INT32_C(2), "local_get_block_pos(): part type %d: now free CommBuffer", type);
    
    free(CommBuffer);
  }
  
  *pread = f->no_part;
  return *pread;
}

static uint64_t
local_get_block_vel(io_logging_t log,
                    io_gizmo_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg,
                    hid_t hdf5_grp[])
{
  uint32_t blocksize, blocksize2;
  int32_t partsize;
  uint32_t bytes_file;
  double fmomx, fmomy, fmomz;
  float dummy;
  uint64_t i;
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  int num_bytes = 0;
  int num_elements = 3; // 3 velocity components: x,y,z
  
  for (type = 0; type < 6; type ++)
  {
    io_logging_msg(log, INT32_C(2), "local_get_block_vel(): read %d of part type %d ", f->header->np[type], type);
    
    if (f->header->np[type] == 0)
    {
      io_logging_msg(log, INT32_C(2), "local_get_block_vel(): no particles of type %d, skip! ", type);
      continue;
    }
    
    if (f->header->flagdoubleprecision)
    {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
      num_bytes = f->header->np[type] * num_elements * sizeof(double);
      CommBuffer = (double *)malloc(num_bytes);
    } else {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
      num_bytes = f->header->np[type] * num_elements * sizeof(float);
      CommBuffer = (float *)malloc(num_bytes);
    }
    
    if (CommBuffer == NULL)
    {
      io_logging_fatal(log, "local_get_block_vel(): could not allocate %d bytes", (int)num_bytes);
      return UINT64_C(0);
    }
    
    io_logging_msg(log, INT32_C(2), "local_get_block_vel(): part type %d: allocated %d bytes for CommBuffer", type, (int)num_bytes);
    
    io_util_readhdf5(log, f, "Velocities", type, num_elements, hdf5_datatype, hdf5_grp, CommBuffer);
    
    /* Loop over the particle positions */
    for (i = 0; i < num_elements * f->header->np[type]; i += num_elements) {
      if (f->header->flagdoubleprecision)
      {
        fmomx = ((double *)CommBuffer)[i];
        fmomy = ((double *)CommBuffer)[i + 1];
        fmomz = ((double *)CommBuffer)[i + 2];
      } else {
        fmomx = (double)((float *)CommBuffer)[i];
        fmomy = (double)((float *)CommBuffer)[i + 1];
        fmomz = (double)((float *)CommBuffer)[i + 2];
      }
      
      if (strg.bytes_float != sizeof(float))
      {
        *((double *)strg.momx.val) = fmomx;
        *((double *)strg.momy.val) = fmomy;
        *((double *)strg.momz.val) = fmomz;
      } else {
        *((float *)strg.momx.val) = (float)fmomx;
        *((float *)strg.momy.val) = (float)fmomy;
        *((float *)strg.momz.val) = (float)fmomz;
      }
      
      /* Increment the pointers to the next particle */
      strg.momx.val = (void *)(((char *)strg.momx.val)
                               + strg.momx.stride);
      strg.momy.val = (void *)(((char *)strg.momy.val)
                               + strg.momy.stride);
      strg.momz.val = (void *)(((char *)strg.momz.val)
                               + strg.momz.stride);
    } /* End of particle velocity loop */
    
    free(CommBuffer);
  }
  
  pread = f->no_part;
  return pread;
}

static uint64_t
local_get_block_id(io_logging_t log,
                   io_gizmo_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg,
                   hid_t hdf5_grp[])
{
  uint32_t blocksize, blocksize2;
  int32_t partsize;
  uint32_t bytes_int_file;
  uint64_t fid;
  uint64_t i;
  uint32_t dummy_int;
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  int num_bytes = 0;
  
  for (type = 0; type < 6; type ++)
  {
    io_logging_msg(log, INT32_C(2), "local_get_block_id(): read %d of part type %d ", f->header->np[type], type);
    
    if (f->header->np[type] == 0)
    {
      io_logging_msg(log, INT32_C(2), "local_get_block_id(): no particles of type %d, skip! ", type);
      continue;
    }
    
    // @TODO: Worry about this
    hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
    num_bytes = f->header->np[type] * sizeof(uint32_t);
    CommBuffer = (uint32_t *)malloc(num_bytes);
    
    if (CommBuffer == NULL)
    {
      io_logging_fatal(log, "local_get_block_id(): could not allocate %d bytes", (int)num_bytes);
      return UINT64_C(0);
    }
    
    io_logging_msg(log, INT32_C(2), "local_get_block_id(): part type %d: allocated %d bytes for CommBuffer", type, (int)num_bytes);
    
    io_util_readhdf5(log, f, "ParticleIDs", type, 1, hdf5_datatype, hdf5_grp, CommBuffer);
    
    /* Loop over the particle IDs */
    for (i = 0; i < f->header->np[type]; i++) {
      fid = ((uint32_t *)CommBuffer)[i];
      *((uint32_t *)strg.id.val) = (uint32_t)fid;
      
      /* Increment the pointers to the next particle */
      strg.id.val = (void *)(((char *)strg.id.val) + strg.id.stride);
      
    } /* End of particle ID loop */
    
    free(CommBuffer);
  }
  
  pread = f->no_part;
  return pread;
}

static uint64_t
local_get_block_mass(io_logging_t log,
                     io_gizmo_t f,
                     uint64_t pskip,
                     uint64_t pread,
                     io_file_strg_struct_t strg,
                     hid_t hdf5_grp[])
{
  uint32_t blocksize=0, blocksize2=0;
  int32_t partsize;
  uint32_t bytes_file;
  double fweight, oldfweight;
  uint32_t curprtt;
  float dummy;
  uint64_t i, j;
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  int num_bytes = 0;
  
  /* Initialize some things */
  f->sumweight = 0.0;
  f->no_species = 0;
  oldfweight = 0.0;
  
  /* We are going to need that a few times for book-keeping */
#       define BOOK_KEEPING {\
if (   isgreater(fweight, oldfweight) \
|| isless(fweight, oldfweight)) { \
f->no_species++; \
oldfweight = fweight; \
if (isless(fweight,f->minweight)) \
f->minweight = fweight; \
if (isgreater(fweight, f->maxweight)) \
f->maxweight = fweight; \
if (type == 1 && isless(fweight, f->mmass)) \
f->mmass = fweight; \
} \
}
  
  if (f->multimass)
  {
    for (type = 0; type < 6; type ++)
    {
      io_logging_msg(log, INT32_C(2), "local_get_block_mass(): read %d of part type %d ", f->header->np[type], type);
      
      if (f->header->np[type] == 0)
      {
        io_logging_msg(log, INT32_C(2), "local_get_block_mass(): no particles of type %d, skip! ", type);
        continue;
      }
      
      if (f->header->flagdoubleprecision)
      {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
        num_bytes = f->header->np[type] * sizeof(double);
        CommBuffer = (double *)malloc(num_bytes);
      } else {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
        num_bytes = f->header->np[type] * sizeof(float);
        CommBuffer = (float *)malloc(num_bytes);
      }
      
      if (CommBuffer == NULL)
      {
        io_logging_fatal(log, "local_get_block_mass(): could not allocate %d bytes", (int)num_bytes);
        return UINT64_C(0);
      }
      
      io_util_readhdf5(log, f, "Masses", type, 1, hdf5_datatype, hdf5_grp, CommBuffer);
      
      /* Loop over the particle massess */
      for (i = 0; i < f->header->np[type]; i++) {
        if (f->header->flagdoubleprecision)
        {
          fweight = ((double *)CommBuffer)[i];
        } else {
          fweight = (double)((float *)CommBuffer)[i];
        }
        
        if (strg.bytes_float != sizeof(float))
        {
          *((double *)strg.weight.val) = fweight;
        } else {
          *((float *)strg.weight.val) = (float)fweight;
        }
        
        /* Increment the pointers to the next particle */
        strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride);
        
        BOOK_KEEPING;
        f->sumweight += fweight;
      } /* End of particle mass loop */
      
      free(CommBuffer);
    }
  }
  else {
    for (type = 0; type < 6; type++)
    {
      for (i = 0; i < f->no_part; i++)
      {
        fweight = f->header->massarr[type];
        
        BOOK_KEEPING;
        f->sumweight += fweight;
        
        if (f->header->flagdoubleprecision)
        {
          *((double *)strg.weight.val) = fweight;
        } else {
          *((float *)strg.weight.val) = (float)fweight;
        }
        
        /* Increment the pointers to the next particle */
        strg.weight.val = (void *)(((char *)strg.weight.val) + strg.weight.stride);
      }
    }
  }
  
  pread = f->no_part;
  return pread;
}

static uint64_t
local_get_block_u(io_logging_t log,
                  io_gizmo_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg,
                  hid_t hdf5_grp[])
{
  uint32_t blocksize, blocksize2;
  int32_t partsize;
  uint32_t bytes_file;
  double fweight, oldfweight;
  uint32_t curprtt;
  uint64_t i, psum;
  float dummy;
  double fu;
  int ptype;
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  int num_bytes = 0;
  
  for (type = 0; type < 6; type ++)
  {
    if (type == 0 && f->header->np[type] != 0) {
      if (f->header->flagdoubleprecision)
      {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
        num_bytes = f->header->np[type] * sizeof(double);
        CommBuffer = (double *)malloc(num_bytes);
      } else {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
        num_bytes = f->header->np[type] * sizeof(float);
        CommBuffer = (float *)malloc(num_bytes);
      }
      
      if (CommBuffer == NULL)
      {
        io_logging_fatal(log, "local_get_block_u(): could not allocate %d bytes", (int)num_bytes);
        return UINT64_C(0);
      }
      
      io_logging_msg(log, INT32_C(2), "local_get_block_u(): part type %d: allocated %d bytes for CommBuffer", type, (int)num_bytes);
      
      io_util_readhdf5(log, f, "InternalEnergy", type, 1, hdf5_datatype, hdf5_grp, CommBuffer);
    }
    
    /* Loop over the particles */
    for (i = 0; i < f->header->np[type]; i++)
    {
      if (type == 0)
      {
        if (f->header->flagdoubleprecision)
        {
          fu = ((double *)CommBuffer)[i];
        } else {
          fu = (double)((float *)CommBuffer)[i];
        }
      } else {
        fu = -1.0 * type;
      }
      
      /* Store the particle in the array */
      if (strg.bytes_float != sizeof(float))
      {
        *((double *)strg.u.val) = fu;
      } else {
        *((float *)strg.u.val) = (float)fu;
      }
      
      /* Increment the pointers to the next particle */
      strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride);
    } /* End of particle velocity loop */
    
    if (type == 0 && f->header->np[type] != 0)
    {
      free(CommBuffer);
    }
  }
  
  pread = f->no_part;
  return pread;
}

#ifdef METALHACK
#include "../define.h"
/* No we define a bunch of macros to make life easier */
#define SKIP  {io_util_readuint32(f->file, &blocksize,  f->swapped);}
#define SKIP2 {io_util_readuint32(f->file, &blocksize2, f->swapped);}
#define CHECK_BLOCK {\
if (blocksize != blocksize2) {\
io_logging_fatal(log,\
"The block boundaries (beginning: %" PRIu32\
" end: %" PRIu32 ") are not identical. "\
"Corrupt file?", blocksize, blocksize2);\
fflush(NULL);\
exit(999);\
return UINT64_C(0);\
} else {\
io_logging_msg(log, INT32_C(5),\
"Block claimed correctly to be %f MB long.", \
(float)(blocksize/1024./1024.));\
}\
}
#define GET_BLOCK {\
SKIP;\
io_util_readstring(f->file, str, (size_t)4);\
io_util_readuint32(f->file, &nextblocksize, f->swapped);\
io_logging_msg(log, INT32_C(1),\
"Arrived at block %s, size of it will be %" \
PRIi32, str, nextblocksize);\
SKIP2;\
CHECK_BLOCK;\
}

static uint64_t
local_get_block_z(io_logging_t log,
                  io_gizmo_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg)
{
  uint32_t blocksize, blocksize2, nextblocksize;
  int32_t partsize, bytes_file;
  char str[5];
  int tries;
  uint64_t i, metalzread, starsoffset, metalzskip;
  double fz;
  float dummy;
  
  /* if this is not a Gizmo 2 file, skip until METAL block */
  /* (assuming the orderting POS,VEL,ID,MASS,UGAS,RHO,NE,NH,HSML,SFR,AGE,Z) */
  if (f->ver != 2) {
    local_skip_GIZMO1_blocks(f->file, f->swapped, 11); // we skip the first 11 blocks
  }
  
  /* we have a Gizmo 2 file and hence use the HEAD to find Z */
  else {
    /* Go to the metal block */
    str[0] = 'A';
    tries = 0;
    nextblocksize = 0;
    while ( (strncmp(str, "Z   ", 4) != 0) && (tries < 100)) {
#ifdef GIZMO_MAGNETICUM
      fprintf(stderr,"skipping block %s of size %"PRIu32"\n",str,nextblocksize);
#endif
      fseek(f->file, nextblocksize, SEEK_CUR);
      GET_BLOCK;
      tries++;
    }
#ifdef GIZMO_MAGNETICUM
    fprintf(stderr,"found block %s of size %"PRIu32"\n",str,nextblocksize);
#endif
  }
  
  /* Start with the block */
  SKIP;
  bytes_file = blocksize / (f->header->np[0]+f->header->np[4]);
  CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
  io_logging_msg(log, INT32_C(1),
                 "A total of %" PRIu64 " metallicities "
                 "with %" PRIi32 " bytes per float (%f MB total) "
                 "are stored.",
                 f->header->np[0]+f->header->np[4], bytes_file,
                 (float)(blocksize/1024./1024.));
  
  /* Set the particle size */
  partsize = bytes_file;
  
  /* Find the right spot */
  metalzskip = UINT64_C(0);
  metalzread = UINT64_C(0);
  starsoffset =   f->header->np[0] + f->header->np[1]
  + f->header->np[2] + f->header->np[3];
  if (pskip < f->header->np[0])
    metalzskip = pskip; // Skip the first gas particles
  else if (pskip < starsoffset)
    metalzskip = f->header->np[0]; // Skip all gas
  else if (pskip < starsoffset + f->header->np[4])
    metalzskip = f->header->np[0] + pskip - starsoffset; // Skip all gas and a few stars
  else
    metalzskip = f->header->np[0] + f->header->np[4]; // Skip all gas and all stars
  
  /* Now we skip to right position in the file */
  fseek(f->file, metalzskip * partsize, SEEK_CUR);
  
  /* Read the gas metalz */
  for (i=0; i<pread && i+pskip<f->header->np[0]; i++) {
    /* STEP 1:  Read the metal from the file */
    if (bytes_file == sizeof(float)) {
      io_util_readfloat(f->file, &dummy, f->swapped);
      fz = (double)dummy;
    } else  {
      io_util_readdouble(f->file, &fz, f->swapped);
    }
    /* STEP 2:  Store the metal in the array */
    if (strg.bytes_float == sizeof(float)) {
      *((float *)strg.z.val) = (float)fz;
    } else {
      *((double *)strg.z.val) = fz;
    }
    /* STEP 3:  Increment the pointers to the next particle */
    strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
    metalzread++;
  }
  
  /* Conjure up metalz for the nonmetal particles */
  for (; i<pread && i+pskip<starsoffset; i++) {
    fz = 0.0;
    if (strg.bytes_float == sizeof(float)) {
      *((float *)strg.z.val) = (float)fz;
    } else {
      *((double *)strg.z.val) = fz;
    }
    strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
  }
  
  /* Read the star metalz */
  for (; i<pread && i+pskip<starsoffset+f->header->np[4]; i++) {
    /* STEP 1:  Read the metal from the file */
    if (bytes_file == sizeof(float)) {
      io_util_readfloat(f->file, &dummy, f->swapped);
      fz = (double)dummy;
    } else  {
      io_util_readdouble(f->file, &fz, f->swapped);
    }
    /* STEP 2:  Store the metal in the array */
    if (strg.bytes_float == sizeof(float)) {
      *((float *)strg.z.val) = (float)fz;
    } else {
      *((double *)strg.z.val) = fz;
    }
    /* STEP 3:  Increment the pointers to the next particle */
    strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
    metalzread++;
  }
  
  /* Conjure up metalz for the nonmetal particles */
  for (; i<pread; i++) {
    fz = 0.0;
    if (strg.bytes_float == sizeof(float)) {
      *((float *)strg.z.val) = (float)fz;
    } else {
      *((double *)strg.z.val) = fz;
    }
    strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
  }
  
  /* Go to the end of the block */
  fseek(f->file,
        partsize*(f->header->np[0]+f->header->np[4]-metalzread-metalzskip),
        SEEK_CUR);
  
  io_logging_msg(log, INT32_C(1),
                 "A total of %" PRIu64 " metals were read.",
                 metalzread);
  
  /* Finish the block */
  SKIP2;
  CHECK_BLOCK;
  
  /* And return the number of read particles for error checking */
  return pread;
}

static uint64_t
local_get_block_age(io_logging_t log,
                    io_gizmo_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg)
{
  uint32_t blocksize, blocksize2, nextblocksize;
  int32_t partsize, bytes_file;
  char str[5];
  int tries;
  uint64_t i, agesread, agesskip, starsoffset;
  double fage;
  float dummy;
  
  /* if this is not a Gizmo 2 file, skip until AGE block */
  /* (assuming the orderting POS,VEL,ID,MASS,UGAS,RHO,NE,NH,HSML,SFR,AGE,METAL) */
  if (f->ver != 2) {
    local_skip_GIZMO1_blocks(f->file, f->swapped, 10); // we skip the first 10 blocks
  }
  
  /* we have a Gizmo 2 file and hence use the HEAD to find AGE */
  else {
    /* Go to the AGE block */
    str[0] = 'X';
    tries = 0;
    nextblocksize = 0;
    while ( (strncmp(str, "AGE ", 4) != 0) && (tries < 100)) {
#ifdef GIZMO_MAGNETICUM
      fprintf(stderr,"skipping block %s of size %"PRIu32"\n",str,nextblocksize);
#endif
      fseek(f->file, nextblocksize, SEEK_CUR);
      GET_BLOCK;
      tries++;
    }
#ifdef GIZMO_MAGNETICUM
    fprintf(stderr,"found block %s of size %"PRIu32"\n",str,nextblocksize);
#endif
  }
  
  /* Start with the block */
  SKIP;
  bytes_file = blocksize / (f->header->np[4]);
  CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
  io_logging_msg(log, INT32_C(1),
                 "A total of %" PRIu64 " star "
                 "with %" PRIi32 " bytes per float (%f MB total) "
                 "are stored.",
                 f->header->np[4], bytes_file,
                 (float)(blocksize/1024./1024.));
  
  /* Set the particle size */
  partsize = bytes_file;
  
  /* Find the right spot in the file */
  agesskip = UINT64_C(0);
  agesread = UINT64_C(0);
  starsoffset =   f->header->np[0] + f->header->np[1]
  + f->header->np[2] + f->header->np[3];
  if (pskip < starsoffset)
    agesskip = 0;
  else if (pskip < starsoffset + f->header->np[4])
    agesskip = pskip - starsoffset;  // Skip the first few stars
  else
    agesskip = f->header->np[4]; // Skip all stars
  
  fseek(f->file, partsize*(agesskip), SEEK_CUR);
  
  /* Set ages of particles type 0-3*/
  for (i=0; i<pread && i+pskip<starsoffset; i++) {
    fage = 0.0;
    if (strg.bytes_float == sizeof(float)) {
      *((float *)strg.age.val) = (float)fage;
    } else {
      *((double *)strg.age.val) = fage;
    }
    strg.age.val = (void *)(((char *)strg.age.val) + strg.age.stride);
  }
  
  /* Read the star ages */
  for (; i<pread && i+pskip<starsoffset+f->header->np[4]; i++) {
    /* STEP 1:  Read the age from the file */
    if (bytes_file == sizeof(float)) {
      io_util_readfloat(f->file, &dummy, f->swapped);
      fage = (double)dummy;
    } else  {
      io_util_readdouble(f->file, &fage, f->swapped);
    }
    /* STEP 2:  Store the age in the array */
    if (strg.bytes_float == sizeof(float)) {
      *((float *)strg.age.val) = (float)fage;
    } else {
      *((double *)strg.age.val) = fage;
    }
    /* STEP 3:  Increment the pointers to the next particle */
    strg.age.val = (void *)(((char *)strg.age.val) + strg.age.stride);
    agesread++;
  }
  
  /* Set ages of particle type 5 */
  for (; i<pread; i++) {
    fage = 0.0;
    if (strg.bytes_float == sizeof(float)) {
      *((float *)strg.age.val) = (float)fage;
    } else {
      *((double *)strg.age.val) = fage;
    }
    strg.age.val = (void *)(((char *)strg.age.val) + strg.age.stride);
  }
  
  io_logging_msg(log, INT32_C(1),
                 "A total of %" PRIu64 " ages were read.",
                 agesread);
  
  /* Go to the end of the block */
  fseek(f->file,
        partsize*(f->header->np[4]-agesread-agesskip),
        SEEK_CUR);
  
  /* Finish the block */
  SKIP2;
  CHECK_BLOCK;
  
  /* And return the number of read particles for error checking */
  return pread;
}

void local_skip_GIZMO1_blocks(FILE *fpgizmo, int swapped, int nblocks)
{
  double       ddummy;
  int          idummy;
  int          nskipblocks;
  unsigned int uidummy;
  unsigned int blocksize1, blocksize2;
  int          iblock;
  char         unused[GIZMO_HEADER_FILLHEADER+1];
  
  /* set the number of blocks to be skipped (excluding the header!) */
  nskipblocks = nblocks;
  
  /* rewind file back to the beginning */
  rewind(fpgizmo);
  
#ifdef VERBOSE_GIZMO1_HEADER
  /* skip header */
  io_util_readuint32(fpgizmo, &blocksize1, swapped); // GIZMO1-SKIP
  
  io_util_readint32(fpgizmo, &idummy, swapped);   //fprintf(stderr,"np0=%d\n",idummy);
  io_util_readint32(fpgizmo, &idummy, swapped);   //fprintf(stderr,"np1=%d\n",idummy);
  io_util_readint32(fpgizmo, &idummy, swapped);   //fprintf(stderr,"np2=%d\n",idummy);
  io_util_readint32(fpgizmo, &idummy, swapped);   //fprintf(stderr,"np3=%d\n",idummy);
  io_util_readint32(fpgizmo, &idummy, swapped);   //fprintf(stderr,"np4=%d\n",idummy);
  io_util_readint32(fpgizmo, &idummy, swapped);   //fprintf(stderr,"np5=%d\n",idummy);
  io_util_readdouble(fpgizmo, &ddummy, swapped);
  io_util_readdouble(fpgizmo, &ddummy, swapped);
  io_util_readdouble(fpgizmo, &ddummy, swapped);
  io_util_readdouble(fpgizmo, &ddummy, swapped);
  io_util_readdouble(fpgizmo, &ddummy, swapped);
  io_util_readdouble(fpgizmo, &ddummy, swapped);
  io_util_readdouble(fpgizmo, &ddummy, swapped);  //fprintf(stderr,"expansion=%g\n",ddummy);
  io_util_readdouble(fpgizmo, &ddummy, swapped);  //fprintf(stderr,"redshift=%g\n",ddummy);
  io_util_readint32(fpgizmo, &idummy, swapped);
  io_util_readint32(fpgizmo, &idummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readint32(fpgizmo, &idummy, swapped);
  io_util_readint32(fpgizmo, &idummy, swapped);
  io_util_readdouble(fpgizmo, &ddummy, swapped);  //fprintf(stderr,"boxsize=%g\n",ddummy);
  io_util_readdouble(fpgizmo, &ddummy, swapped);  //fprintf(stderr,"omega=%g\n",ddummy);
  io_util_readdouble(fpgizmo, &ddummy, swapped);  //fprintf(stderr,"omegalambda=%g\n",ddummy);
  io_util_readdouble(fpgizmo, &ddummy, swapped);  //fprintf(stderr,"hubble=%g\n",ddummy);
  io_util_readint32(fpgizmo, &idummy, swapped);
  io_util_readint32(fpgizmo, &idummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readuint32(fpgizmo, &uidummy, swapped);
  io_util_readint32(fpgizmo, &idummy, swapped);
  io_util_readstring(fpgizmo, unused, GIZMO_HEADER_FILLHEADER);
  
  io_util_readuint32(fpgizmo, &blocksize2, swapped); // GIZMO1-SKIP
  
  //fprintf(stderr,"skipped HEADER of size %u vs. %u bytes\n",blocksize1,blocksize2);
  
#else /* VERBOSE_GIZMO1_HEADER */
  nskipblocks++;
#endif
  
  /* skip GIZMO1 blocks */
  for(iblock=0; iblock<nskipblocks; iblock++)
  {
    io_util_readuint32(fpgizmo, &blocksize1, swapped); // GIZMO1-SKIP
    fseek(fpgizmo, blocksize1, SEEK_CUR);              // block-SKIP
    io_util_readuint32(fpgizmo, &blocksize2, swapped); // GIZMO1-SKIP
    
    //fprintf(stderr,"skipped block #%d (out of %d in total) of size %u vs. %d MB\n",iblock,nblocks,blocksize1/1024/1024,blocksize2/1024/1024);
    
    /* the file appears to be corrupted */
    if(blocksize1 != blocksize2)
    {
      fprintf(stderr,"We are already trying to help you with your ancient GIZMO1 file, but enough is enough!\n");
      exit(-1);
    }
  }
}
#undef GET_BLOCK
#undef CHECK_FLOATBYTES
#undef SKIP
#undef SKIP2
#undef CHECK_BLOCK
#endif // METALHACK

#endif // WITH_HDF5
