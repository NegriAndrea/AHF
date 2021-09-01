#ifndef IO_CUBEP3M_HEADER_DEF_H
#define IO_CUBEP3M_HEADER_DEF_H

/**
 * \file io_cubep3m_header_def.h
 *
 * Provides the structure definition for the CUBEP3M header
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
#define CUBEP3M_HEADER_SIZE 48

/* The size (in bytes) of the header in chunked files. */
#define CUBEP3M_HEADER_CHUNK_SIZE 208

/* The size (in bytes) of the file identification magic in chunked files. */
#define CUBEP3M_HEADER_MAGIC_SIZE 20

/*
 * The header structure itself
 */
struct io_cubep3m_header_struct {
	// This is the basic CubePM header information.
	//   48 bytes
	int32_t np_local;
	float   a;
	float   t;
	float   tau;
	int32_t nts;
	float   dt_f_acc;
	float   dt_pp_acc;
	float   dt_c_acc;
	int32_t cur_checkpoint;
	int32_t cur_projection;
	int32_t cur_halofind;
	float   mass_p;

	// This is the chunk header information (not always present)
	//   140 bytes
	int64_t tot_np_in_chunk;
	int64_t np_in_chunk_file;
	int     num_chunk_files;
	double  chunk_offset[3];
	double  chunk_start_full_data[3];
	double  chunk_end_full_data[3];
	double  chunk_start_real_data[3];
	double  chunk_end_real_data[3];

	// this is additional information >>not<< stored in the file!
	double   omega0;   // in rho_crit_0
	double   lambda0;  // in rho_crit_0
	double   boxsize;  // in Mpc/h
	double   lunit;    // conversion factor from internal to Mpc/h
	double   vunit;    // conversion factor from internal to physical km/s
	double   munit;    // mass of a particle with weigth 1. in Msun/h
	uint64_t nptotal;
	uint64_t ngrid;
	int      nodes_dim;
	int      file_number;
	double   offset[3];
};

/** Convenient typedef */
typedef struct io_cubep3m_header_struct io_cubep3m_header_struct_t;

/** Convenient typedef */
typedef io_cubep3m_header_struct_t      *io_cubep3m_header_t;


#endif /* IO_CUBEP3M_HEADER_DEF_H */
