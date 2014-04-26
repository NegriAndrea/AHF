#ifndef IO_DEVA_HEADER_DEF_H 
#define IO_DEVA_HEADER_DEF_H

/**
 * \file io_deva_header_def.h
 *
 * Provides the structure definition for the DEVA header
 * structure. Including useful typedefs.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/** 
 * The size (in bytes) reserved at the beginning of a file for
 * the header 
 */
#define DEVA_HEADER_SIZE 800

/** Defines the total used size (in bytes) of the header */
#define DEVA_HEADER_HEADERSIZE ( 800 )

/** Defines the number of unused bytes in the header */
#define DEVA_HEADER_FILLHEADER (  DEVA_HEADER_SIZE - DEVA_HEADER_HEADERSIZE)


/**
 * The header structure itself
 */
struct io_deva_header_struct {
  int   itime;
  int   itstop;
  int   itdump;
  int   iout;
  int   nsformed;
  int   nsdead;
  int   irun;
  int   nobj;
  int   ngas;
  int   ndark;
  int   L;
  float CHEMEVOL;
  float ANOTHER;
  float COOL;
  float REFINEMENT;
  float DEVA_HYDRO; // HYDRO is a DEFINEFLAGS and hence reserved!
  float GRAVITY;
  int   ISOLATED;
  float EXPAND;
  float COMOVING;
  float STARFORM;
  float GRADH;
  int   INITIALCOND;
  int   nstar;
  int   iseed1;
  int   ispec;
  int   indxsp;
  int   n_neigh;
  int   lastbar;
  int   fill1[100-29];
  
  float time;
  float atime;
  float htime;
  float dtime;
  float E_init;
  float E_kin;
  float E_ther;
  float E_pot;
  float Radiation;
  float Esum;
  float Rsum;
  float cpu;
  float time_end;
  float tout;
  float padding;
  float Tlost;
  float Qlost;
  float Ulost;
  float delta_min;
  float delta_max;
  float T_min;
  float avisc;
  float bvisc;
  float eta2;
  float rho_star;
  float c_star;
  float rmtot;
  float rmsep;
  float dnthres;
  float sft0;
  float sftmin;
  float sftmax;
  float h100;
  float box100;
  float rmgas;
  float rmdark;
  float omega0;
  float xlambda0;
  float h0t0;
  float omegab0;
  float sigma80;
  float ztime0;
  float e0;
  float fill2[100-43];
};

/** Convenient typedef */
typedef struct io_deva_header_struct io_deva_header_struct_t;

/** Convenient typedef */
typedef io_deva_header_struct_t *io_deva_header_t;


#endif /* IO_DEVA_HEADER_DEF_H */
