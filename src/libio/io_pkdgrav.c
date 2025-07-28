/* Belaid, Doug & Isaac Jul. 5, 2022 */

/**
 * \file io_pkdgrav.c
 *
 * Provides functions for reading and writing PKDGRAV files.
 * It is basically a copy of io_gizmo.c but changing the unit handling and
 * some other small details
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

#include "io_pkdgrav.h"
#include "io_pkdgrav_header.h"
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
 * io_pkdgrav_open more readable.
 *
 * \param log   The logging object.
 * \param f     The Gizmo file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_pkdgrav_t
local_openopen(io_logging_t log, io_pkdgrav_t f, io_file_mode_t mode);

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
                  io_pkdgrav_t f,
                  io_file_swap_t swapped);

/**
 * \brief Try to find out swapping status.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Returns the file object or NULL in case of an error.
 */
inline static io_pkdgrav_t
local_opengetswap(io_logging_t log, io_pkdgrav_t f);

/**
 * \brief Tries to figure out which Gizmo file version is to be used.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Nothing.
 */
inline static void
local_openversion(io_logging_t log, io_pkdgrav_t f);

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
                   io_pkdgrav_t f,
                   uint64_t pskip,
                   int32_t bytes);

static uint64_t
local_get_block_pos(io_logging_t log,
                    io_pkdgrav_t f,
                    uint64_t *pskip,
                    uint64_t *pread,
                    io_file_strg_struct_t strg,
                    hid_t hdf5_grp[]);

static uint64_t
local_get_block_vel(io_logging_t log,
                    io_pkdgrav_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg,
                    hid_t hdf5_grp[]);

static uint64_t
local_get_block_id(io_logging_t log,
                   io_pkdgrav_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg,
                   hid_t hdf5_grp[]);

static uint64_t
local_get_block_mass(io_logging_t log,
                     io_pkdgrav_t f,
                     uint64_t pskip,
                     uint64_t pread,
                     io_file_strg_struct_t strg,
                     hid_t hdf5_grp[]);

static uint64_t
local_get_block_u(io_logging_t log,
                  io_pkdgrav_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg,
                  hid_t hdf5_grp[]);

void local_find_block(io_pkdgrav_t, char *);

#ifdef METALHACK
static uint64_t
local_get_block_z(io_logging_t log,
                  io_pkdgrav_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg,
                  hid_t hdf5_grp[]);

static uint64_t
local_get_block_age(io_logging_t log,
                    io_pkdgrav_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg,
                    hid_t hdf5_grp[]);
#endif



/**********************************************************************\
 *    Implementation of global functions                              *
 \**********************************************************************/
extern io_pkdgrav_t
io_pkdgrav_open(io_logging_t log,
              char *fname,
              io_file_swap_t swapped,
              io_file_mode_t mode,
              uint32_t reader)
{
  io_pkdgrav_t f;
  
  /* Get memory for the structure */
  f = (io_pkdgrav_t)malloc(sizeof(io_pkdgrav_struct_t));
  if (f == NULL) {
    io_logging_memfatal(log,  "io_pkdgrav structure");
    return NULL;
  }
  
  /* Start filling the structure */
  
  /* Store the filename */
  f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
  if (f->fname == NULL) {
    io_logging_memfatal(log, "filename of PkdgravFile");
    free(f);
    return NULL;
  }
  strncpy(f->fname, fname, strlen(fname)+1);
  
  /* Okay, we are a PKDGRAV file */
  f->ftype = IO_FILE_PKDGRAV;
  
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
  
  /* Identify PKDGRAV3 format */
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
io_pkdgrav_close(io_logging_t log,
               io_pkdgrav_t *f)
{
  /* Catch NULLs */
  if (f == NULL || *f == NULL)
    return;
  
  /* Put header to the file if necessary */
  if (    ((*f)->mode == IO_FILE_WRITE)
      && ((*f)->header != NULL)) {
    io_pkdgrav_header_write(log, (*f)->header, *f);
  }
  
  /* Close */
  if ((*f)->header != NULL)
    io_pkdgrav_header_del(log, &((*f)->header));
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
io_pkdgrav_init(io_logging_t log,
              io_pkdgrav_t f)
{
  if (f == NULL)
    return;
  
  if (f->header != NULL) {
    io_logging_warn(log, INT32_C(1),
                    "Already have the header information! Rereading.");
    io_pkdgrav_header_del(log, &(f->header));
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
  f->header = io_pkdgrav_header_get(log, f);
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
io_pkdgrav_readpart(io_logging_t log,
                  io_pkdgrav_t f,
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
  particles_read = io_pkdgrav_readpart_raw(log, f, pskip, pread, strg);
  
  if (particles_read != pread) {
    return UINT64_C(0);
  }
  
  /* And do the scaling */
#ifdef WITH_MPI
  io_pkdgrav_scale_global(log, f->mycomm,  f->maxpos, f->minpos, &(f->mmass));
#endif
  /* I do not understand why we pass the mmass (minimal mass) as the mass
   * scale... so I have changed it.
   * It seems that the normalization by mmass is somewhat arbitrary.
   */
  f->mmass = 1.;

  tmp = io_pkdgrav_scale_particles(log, f->maxpos, f->minpos,
                                 f->header->boxsize,
                                 f->header->expansion,
                                 f->header->KpcUnit,
                                 f->header->MsolUnit,
                                 f->header->ErgPerGmUnit,
                                 f->header->KmPerSecUnit,
                                 particles_read, strg);
  if (tmp != particles_read) {
    return tmp;
  }
  
  /* Wow, we are done! */
  return particles_read;
}

extern uint64_t
io_pkdgrav_readpart_raw(io_logging_t log,
                      io_pkdgrav_t f,
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
  f->file = H5Fopen(f->fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (f->file < 0) {
    io_logging_fatal(log,"io_pkdgrav_readpart_raw(): could not open file %s for reading",f->fname);
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
  
#ifdef METALHACK
  funcrtn = local_get_block_age(log, f, pskip, pread, strg, hdf5_grp);
  if (funcrtn != pread) {
    io_logging_fatal(log,
                     "Expected to read %"PRIu64
                     " particle identities, but only got %"PRIu64
                     ".  Aborting.");
    return UINT64_C(0);
  }
  
  funcrtn = local_get_block_z(log, f, pskip, pread, strg, hdf5_grp);
  if (funcrtn != pread) {
    io_logging_fatal(log,
                     "Expected to read %"PRIu64
                     " particle identities, but only got %"PRIu64
                     ".  Aborting.");
    return UINT64_C(0);
  }
#endif
  
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
io_pkdgrav_writepart(io_logging_t log,
                   io_pkdgrav_t f,
                   uint64_t pskip,
                   uint64_t pwrite,
                   io_file_strg_struct_t strg)
{
  return 0;
}

extern uint64_t
io_pkdgrav_writepart_ord(io_logging_t log,
                       io_pkdgrav_t f,
                       uint64_t pskip,
                       uint64_t pwrite,
                       void *nxt_part,
                       io_file_strg_struct_t strg)
{
  return UINT64_C(0);
}

extern bool
io_pkdgrav_get(io_logging_t log,
             io_pkdgrav_t f,
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
      *((double *)res) =   f->header->boxsize;
      break;
    case IO_FILE_GET_PMASS:
      *((double *)res) = f->header->MsolUnit * f->header->hubbleparameter;
      break;
    case IO_FILE_GET_ZINITIAL:
      io_logging_warn(log, INT32_C(1),
                      "zinitial is not set in a PKDGRAV file, "
                      "using current redshift");
    case IO_FILE_GET_Z:
      *((double *)res) = f->header->redshift;
      break;
    case IO_FILE_GET_AINITIAL:
      io_logging_warn(log, INT32_C(1),
                      "ainitial is not set in a PKDGRAV file, "
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
                      "PKDGRAV files don't store the step number. "
                      "Setting to 0.");
      *((int32_t *)res) = 0;
      break;
    case IO_FILE_GET_TSTEP:
      io_logging_warn(log, INT32_C(1),
                      "PKDGRAV files don't store the timestep. "
                      "Setting to 0.0");
      *((double *)res) = 0.0;
      break;
    case IO_FILE_GET_HEADERSTR:
      io_logging_warn(log, INT32_C(1),
                      "PKDGRAV files don't have a header string. "
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
io_pkdgrav_set(io_logging_t log,
             io_pkdgrav_t f,
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
io_pkdgrav_log(io_logging_t log, io_pkdgrav_t f)
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
  io_pkdgrav_header_log(log, f->header);
  
  return;
}


extern uint64_t
io_pkdgrav_scale_particles(io_logging_t log,
                         double maxpos[],
                         double minpos[],
                         double boxsize,
                         double expansion,
                         double KpcUnit,
                         double MsolUnit,
                         double ErgPerGmUnit,
                         double KmPerSecUnit,
                         uint64_t particles_read,
                         io_file_strg_struct_t strg)
{
  double box[3], shift[3];
  double scale_pos, scale_mom, scale_weight, scale_u;
  uint64_t i;

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
  /* We use directly the unit information in the snapshot to make the convertions
   */
  scale_pos    = 1.0;
  /* Velocity output in PKDGRAV3 is \dot{x} directly, for safety we convert to
   * km/s and then apply the conversion to internal units as in the docs
   */
#ifdef NO_EXPANSION
  scale_mom    = KmPerSecUnit / (boxsize * 100.);
#else
  scale_mom    = KmPerSecUnit * expansion * expansion /
                        ( boxsize * 100.);
#endif

  /* The mass unit seems to be arbitrary, so we left it there */
  scale_weight = 1.0;

  /* Internal energy per unit mass conversion to (km/sec)^2 */
  scale_u      = ErgPerGmUnit / 1e10;
  
  
  
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
io_pkdgrav_scale_global(io_logging_t log,
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

inline static io_pkdgrav_t
local_openopen(io_logging_t log, io_pkdgrav_t f, io_file_mode_t mode)
{
  if (mode == IO_FILE_READ) {
    //f->file = fopen(f->fname, IO_FILE_MODE_READ);
    f->file = H5Fopen(f->fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f->file < 0) {
      io_logging_fatal(log, "Could not open '%s' for reading.", f->fname);
      return NULL;
    }
  } else {
    io_logging_fatal(log, "No write implementation for PKDGRAV HDF5.");
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
                  io_pkdgrav_t f,
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
local_openversion(io_logging_t log, io_pkdgrav_t f)
{
  f->ver = 1;
  io_logging_msg(log, INT32_C(2), "Assuming a public PKDGRAV3 v1 file.");
  
  return;
}

inline static int32_t
local_write_common(io_logging_t log,
                   io_pkdgrav_t f,
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
                    io_pkdgrav_t f,
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
  uint64_t num_bytes = 0;
  int num_elements = 3; // 3 Coordinates: x,y,z
  
  /* Set extreme position detectors */
  f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e40;
  f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e40;
  
  io_logging_msg(log, INT32_C(2), "local_get_block_pos(): loop over particle types");
  
  for (type = 0; type < 6; type ++)
  {
    io_logging_msg(log, INT32_C(2), "local_get_block_pos(): read %" PRIu64 " of part type %d ", f->header->np[type], type);
    
    if (f->header->np[type] == 0)
    {
      io_logging_msg(log, INT32_C(2), "local_get_block_pos(): no particles of type %d, skip! ", type);
      continue;
    }
    
    if (f->header->flagdoubleprecision)
    {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
      num_bytes = (uint64_t) (f->header->np[type]) * num_elements * sizeof(double);
      CommBuffer = (double *)malloc(num_bytes);
    } else {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
      num_bytes = (uint64_t) (f->header->np[type]) * num_elements * sizeof(float);
      CommBuffer = (float *)malloc(num_bytes);
    }
    
    if (CommBuffer == NULL)
    {
      io_logging_fatal(log, "local_get_block_pos(): could not allocate %" PRIu64 " bytes", num_bytes);
      return UINT64_C(0);
    }
    
    io_logging_msg(log, INT32_C(2), "local_get_block_pos(): part type %d: allocated %" PRIu64 " bytes for CommBuffer", type, num_bytes);
    
    io_util_readhdf5_pkdgrav(log, f, "Coordinates", type, num_elements, hdf5_datatype, hdf5_grp, CommBuffer);
    
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
                    io_pkdgrav_t f,
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
  uint64_t num_bytes = 0;
  int num_elements = 3; // 3 velocity components: x,y,z
  
  for (type = 0; type < 6; type ++)
  {
    io_logging_msg(log, INT32_C(2), "local_get_block_vel(): read %" PRIu64 " of part type %d ", f->header->np[type], type);
    
    if (f->header->np[type] == 0)
    {
      io_logging_msg(log, INT32_C(2), "local_get_block_vel(): no particles of type %d, skip! ", type);
      continue;
    }
    
    if (f->header->flagdoubleprecision)
    {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
      num_bytes = (uint64_t) (f->header->np[type]) * num_elements * sizeof(double);
      CommBuffer = (double *)malloc(num_bytes);
    } else {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
      num_bytes = (uint64_t) (f->header->np[type]) * num_elements * sizeof(float);
      CommBuffer = (float *)malloc(num_bytes);
    }
    
    if (CommBuffer == NULL)
    {
      io_logging_fatal(log, "local_get_block_vel(): could not allocate %" PRIu64 " bytes", num_bytes);
      return UINT64_C(0);
    }
    
    io_logging_msg(log, INT32_C(2), "local_get_block_vel(): part type %d: allocated %" PRIu64 " bytes for CommBuffer", type, num_bytes);
    
    io_util_readhdf5_pkdgrav(log, f, "Velocities", type, num_elements, hdf5_datatype, hdf5_grp, CommBuffer);
    
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
                   io_pkdgrav_t f,
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
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  uint64_t num_bytes = 0;
  /*bool long_ids = false;
  
  // D. Rennehan: Going to assume that if there is a high word for any particle we are using long-ids.
  for (type = 0; type < 6; type ++)
  {
    if (f->header->nallhighw[type] != 0)
    {
      long_ids = true;
    }
  }*/

  /* PKDGRAV3 outputs long ids by default */
  bool long_ids = true;
  
  for (type = 0; type < 6; type ++)
  {
    io_logging_msg(log, INT32_C(2), "local_get_block_id(): read %" PRIu64 " of part type %d ", f->header->np[type], type);
    
    if (f->header->np[type] == 0)
    {
      io_logging_msg(log, INT32_C(2), "local_get_block_id(): no particles of type %d, skip! ", type);
      continue;
    }
    
    if (long_ids)
    { 
      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
      num_bytes = (uint64_t) (f->header->np[type]) * sizeof(uint64_t);
      CommBuffer = (uint64_t *)malloc(num_bytes);
    } else {
      hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
      num_bytes = (uint64_t) (f->header->np[type]) * sizeof(uint32_t);
      CommBuffer = (uint32_t *)malloc(num_bytes);
    }
    
    if (CommBuffer == NULL)
    {
      io_logging_fatal(log, "local_get_block_id(): could not allocate %" PRIu64 " bytes", num_bytes);
      return UINT64_C(0);
    }
    
    io_logging_msg(log, INT32_C(2), "local_get_block_id(): part type %d: allocated %" PRIu64 " bytes for CommBuffer", type, num_bytes);
    
    io_util_readhdf5_pkdgrav(log, f, "ParticleIDs", type, 1, hdf5_datatype, hdf5_grp, CommBuffer);
    
    /* Loop over the particle IDs */
    for (i = 0; i < f->header->np[type]; i++) {
      if (long_ids)
      {
        fid = ((uint64_t *)CommBuffer)[i];
        *((uint64_t *)strg.id.val) = (uint64_t)fid;
      } else {
        fid = ((uint32_t *)CommBuffer)[i];
        *((uint32_t *)strg.id.val) = (uint32_t)fid;
      }
      
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
                     io_pkdgrav_t f,
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
  uint64_t num_bytes = 0;
  
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
      io_logging_msg(log, INT32_C(2), "local_get_block_mass(): read %" PRIu64 " of part type %d ", f->header->np[type], type);
      
      if (f->header->np[type] == 0)
      {
        io_logging_msg(log, INT32_C(2), "local_get_block_mass(): no particles of type %d, skip! ", type);
        continue;
      }
      
      if (f->header->flagdoubleprecision)
      {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
        num_bytes = (uint64_t) (f->header->np[type]) * sizeof(double);
        CommBuffer = (double *)malloc(num_bytes);
      } else {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
        num_bytes = (uint64_t) (f->header->np[type]) * sizeof(float);
        CommBuffer = (float *)malloc(num_bytes);
      }
      
      if (CommBuffer == NULL)
      {
        io_logging_fatal(log, "local_get_block_mass(): could not allocate %" PRIu64 " bytes", num_bytes);
        return UINT64_C(0);
      }
      
      io_util_readhdf5_pkdgrav(log, f, "Masses", type, 1, hdf5_datatype, hdf5_grp, CommBuffer);
      
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
                  io_pkdgrav_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg,
                  hid_t hdf5_grp[])
{
  uint64_t i;
  float dummy;
  double fu;
  int ptype;
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  uint64_t num_bytes = 0;
  
  for (type = 0; type < 6; type ++)
  {
    if (type == 0 && f->header->np[type] != 0) {
      if (f->header->flagdoubleprecision)
      {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
        num_bytes = (uint64_t) (f->header->np[type]) * sizeof(double);
        CommBuffer = (double *)malloc(num_bytes);
      } else {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
        num_bytes = (uint64_t) (f->header->np[type]) * sizeof(float);
        CommBuffer = (float *)malloc(num_bytes);
      }
      
      if (CommBuffer == NULL)
      {
        io_logging_fatal(log, "local_get_block_u(): could not allocate %" PRIu64 " bytes", num_bytes);
        return UINT64_C(0);
      }
      
      io_logging_msg(log, INT32_C(2), "local_get_block_u(): part type %d: allocated %" PRIu64 " bytes for CommBuffer", type, num_bytes);
      
      io_util_readhdf5_pkdgrav(log, f, "InternalEnergy", type, 1, hdf5_datatype, hdf5_grp, CommBuffer);
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
static uint64_t
local_get_block_z(io_logging_t log,
                  io_pkdgrav_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg,
                  hid_t hdf5_grp[])
{
  double fz;
  float dummy;
  uint64_t i;
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  uint64_t num_bytes = 0;
  int num_elements = f->header->flagmetals;
  
  for (type = 0; type < 6; type ++)
  {
    io_logging_msg(log, INT32_C(2), "local_get_block_z(): read %" PRIu64 " of part type %d ", f->header->np[type], type);
    
    if (f->header->np[type] == 0 )
    {
      io_logging_msg(log, INT32_C(2), "local_get_block_z(): no particles of type %d, skip! ", type);
      continue;
    }
    
    /* Only gas and stars (pt=0 & pt=4) have metals in cosmo sims */
    if (type == 0 || type == 4)
    {    
      if (f->header->flagdoubleprecision)
      {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
        num_bytes = (uint64_t) (f->header->np[type]) * num_elements * sizeof(double);
        CommBuffer = (double *)malloc(num_bytes);
      } else {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
        num_bytes = (uint64_t) (f->header->np[type]) * num_elements * sizeof(float);
        CommBuffer = (float *)malloc(num_bytes);
      }
      
      if (CommBuffer == NULL)
      {
        io_logging_fatal(log, "local_get_block_z(): could not allocate %" PRIu64 " bytes", num_bytes);
        return UINT64_C(0);
      }
      
      io_logging_msg(log, INT32_C(2), "local_get_block_z(): part type %d: allocated %" PRIu64 " bytes for CommBuffer", type, num_bytes);
      
      io_util_readhdf5_pkdgrav(log, f, "Metallicity", type, num_elements, hdf5_datatype, hdf5_grp, CommBuffer);
      
      /* Loop over the particle positions */
      for (i = 0; i < num_elements * f->header->np[type]; i += num_elements) {
        if (f->header->flagdoubleprecision)
        {
          fz = ((double *)CommBuffer)[i]; // Only take the first element, the mass fraction
        } else {
          fz = (double)((float *)CommBuffer)[i];
        }
        
        if (strg.bytes_float != sizeof(float))
        {
          *((double *)strg.z.val) = fz;
        } else {
          *((float *)strg.z.val) = (float)fz;
        }
        
        /* Increment the pointers to the next particle */
        strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
      } /* End of particle metal loop */
      
      free(CommBuffer); // We allocated this for both pt=0 and pt=4, so free every time
    } else {
      for (i = 0; i < f->header->np[type]; i++) {
        fz = 0.0; // This is already a double
        
        if (strg.bytes_float != sizeof(float))
        {
          *((double *)strg.z.val) = fz;
        } else {
          *((float *)strg.z.val) = (float)fz;
        }
        
        strg.z.val = (void *)(((char *)strg.z.val) + strg.z.stride);
      }
    }
  }
  
  pread = f->no_part;
  return pread;
}

static uint64_t
local_get_block_age(io_logging_t log,
                    io_pkdgrav_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg,
                    hid_t hdf5_grp[])
{
  uint64_t i;
  float dummy;
  double fage;
  int ptype;
  int type;
  hid_t hdf5_datatype = 0;
  void * CommBuffer;
  uint64_t num_bytes = 0;
  
  for (type = 0; type < 6; type ++)
  {
    /* Type==4 is usually stars formed in the simulations, so only treat them */
    if (type == 4 && f->header->np[type] != 0) {
      if (f->header->flagdoubleprecision)
      {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
        num_bytes = (uint64_t) (f->header->np[type]) * sizeof(double);
        CommBuffer = (double *)malloc(num_bytes);
      } else {
        hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
        num_bytes = (uint64_t) (f->header->np[type]) * sizeof(float);
        CommBuffer = (float *)malloc(num_bytes);
      }
      
      if (CommBuffer == NULL)
      {
        io_logging_fatal(log, "local_get_block_age(): could not allocate %" PRIu64 " bytes", num_bytes);
        return UINT64_C(0);
      }
      
      io_logging_msg(log, INT32_C(2), "local_get_block_age(): part type %d: allocated %" PRIu64 " bytes for CommBuffer", type, num_bytes);
      
      /* StellarAge should be the proper tag in the HDF5 file */
      io_util_readhdf5_pkdgrav(log, f, "StellarFormationTime", type, 1, hdf5_datatype, hdf5_grp, CommBuffer);
    }
    
    /* Loop over the particles */
    for (i = 0; i < f->header->np[type]; i++)
    {
      if (type == 4)
      {
        if (f->header->flagdoubleprecision)
        {
          fage = ((double *)CommBuffer)[i];
        } else {
          fage = (double)((float *)CommBuffer)[i];
        }
      } else {
        fage = 0.0;
      }
      
      /* Store the particle in the array */
      if (strg.bytes_float != sizeof(float))
      {
        *((double *)strg.age.val) = fage;
      } else {
        *((float *)strg.age.val) = (float)fage;
      }
      
      /* Increment the pointers to the next particle */
      strg.age.val = (void *)(((char *)strg.age.val) + strg.age.stride);
    } /* End of particle velocity loop */
    
    if (type == 4 && f->header->np[type] != 0)
    {
      free(CommBuffer);
    }
  }
  
  pread = f->no_part;
  return pread;
}

#endif // METALHACK

#endif // WITH_HDF5
