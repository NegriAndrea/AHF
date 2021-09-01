#ifndef IO_ART_HEADER_DEF_H 
#define IO_ART_HEADER_DEF_H

/**
 * \file io_art_header_def.h
 *
 * Provides the structure definition for the ART header
 * structure. Including useful typedefs.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/*
 * The size (in bytes) reserved at the beginning of a file for
 * the header 
 */
#define ART_HEADER_SIZE 529

/*
 * We actually store additional information in the header structure 
 */
#define ART_HEADER_EXTRA 48

/*
 * The header structure itself
 */
struct io_art_header_struct {
   char    header_string[45];
   float   aexpn;
   float   aexp0;
   float   amplt;
   float   astep;
   int     istep;
   float   partw;
   float   tintg;
   float   ekin;
   float   ekin1;
   float   ekin2;
   float   au0;
   float   aeu0;
   int     nrowc;
   int     ngridc;
   int     nspecies;
   int     nseed;
   float   Om0;
   float   Oml0;
   float   hubble;
   float   wp5;
   float   Ocurv;
   float   extras[100];
  
  /* additional information */
  double   boxsize;
  double   munit;
  double   nvpart;
  uint64_t N_particles;
  uint64_t N_pages;
  uint64_t N_in_last;
};

/** Convenient typedef */
typedef struct io_art_header_struct io_art_header_struct_t;

/** Convenient typedef */
typedef io_art_header_struct_t *io_art_header_t;


#endif /* IO_ART_HEADER_DEF_H */
