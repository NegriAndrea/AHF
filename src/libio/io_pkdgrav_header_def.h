#ifndef IO_PKDGRAV_HEADER_DEF_H 
#define IO_PKDGRAV_HEADER_DEF_H


/**
 * \file io_pkdgrav_header_def.h
 *
 * Provides the structure definition for the Pkdgrav header
 * structure. Including useful typedefs.
 */


/**********************************************************************\
 *    Includes                                                        * 
 \**********************************************************************/
#include <inttypes.h>
#include "io_pkdgrav_header_def.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
 \**********************************************************************/

struct io_pkdgrav_header_struct {
  unsigned int np[6];                 /*!< number of particles of each type in this file */
  double massarr[6];               /*!< mass of particles of each type. If 0, then the masses are explicitly
                                    stored in the mass-block of the snapshot file, otherwise they are omitted */
  double expansion;                  /*!< time of snapshot file */
  double redshift;              /*!< redshift of snapshot file */
  int flagsfr;                 /*!< flags whether the simulation was including star formation */
  int flagfeedback;            /*!< flags whether feedback was included (obsolete) */
  unsigned int nall[6];   /*!< total number of particles of each type in this snapshot. This can be
                           different from npart if one is dealing with a multi-file snapshot. */
  int flagcooling;             /*!< flags whether cooling was included  */
  int numfiles;                /*!< number of files in multi-file snapshot */
  double boxsize;               /*!< box-size of simulation in case periodic boundaries were used */
  double omega0;                /*!< matter density in units of critical density */
  double omegalambda;           /*!< cosmological constant parameter */
  double hubbleparameter;           /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flagstellarage;          /*!< flags whether the file contains formation times of star particles */
  int flagmetals;              /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int nallhighw[6];   /*!< High word of the total number of particles of each type (needed to combine with npartTotal to allow >2^31 particles of a given type) */
  int flagdoubleprecision;     /*!< flags that snapshot contains double-precision instead of single precision */
  
  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                 or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                 For snapshots files, the value informs whether the simulation was evolved from
                                 Zeldoch or 2lpt ICs. Encoding is as follows:
                                 FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                 FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                 FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                 FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                 FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                                 All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */
  
  char fill[18];                /*!< fills to 256 Bytes */
  char names[15][2];

  double MsolUnit;
  double KpcUnit;
  double ErgPerGmUnit;
  double KmPerSecUnit;
};

/** Convenient typedef */
typedef struct io_pkdgrav_header_struct io_pkdgrav_header_struct_t;

/** Convenient typedef */
typedef io_pkdgrav_header_struct_t *io_pkdgrav_header_t;


#endif /* IO_PKDGRAV_HEADER_DEF_H */

