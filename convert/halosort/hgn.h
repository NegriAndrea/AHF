#ifndef _HGN_HGN_H
#define _HGN_HGN_H

#include <stdio.h>
#include "rgadget.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief The file structure object
 */
typedef struct hgnfile {
	FILE *fp;		/*!< The FILE stream  */
	char *idline;		/*!< The id line */
	int level;		/*!< The level of the HGN processed */
	char *author;		/*!< The author comment in the file */
	char *hf;		/*!< The halo finder name  */
	char *date;		/*!< Yhe date */
    PartId count;		/*!< The count of particles */
} HGN;

typedef struct _HaloShell {
#define SHELL_RADIAL 1
#define SHELL_MASS 2
    int type;
    double metric;
    double lambda;
    double b, c;
} HaloShell;

// what an individual line looks like
typedef struct  {
    PartId id;
    double mass;
    double age;
    double z;
} PartDesc;


/** \brief Descriptor for an individual sub halo
 */
typedef struct HGNHalo {
  PartId haloid;		/*!< The haloes ID */
  PartId partcount;	/*!< The number of particles in the halo */
  PartDesc *partlist;	/*!< The allocated list of particles  */
   
  double Xc, Yc, Zc;
  double VXc, VYc, VZc;
  double Rad;
  double Mass;
  double Vmax, Rmax;
  double r2;
  double a,b,c;
  double lambda;
    double plummerSpin;
  long unsigned npart;
#define SS_SIZES 9
    HaloShell shells[SS_SIZES];
  double Lang[3];		/* angular momentum vector */
} HGNHalo;

/** The file opening interface to reading a HGN file */
HGN *openhaloidfile(const char *filename);
/** Read the HGN structure from a stream (e.g. STDIN) */
HGN *openhaloidstream(FILE *fp);
/** Fetch the next halo descriptor in the file 
 */

HGNHalo *getnexthalo(HGN *hgn);
/** Free an allocated halo */
void freeHGNHalo(HGNHalo *halo);
/** Finished with the file */
void closehgn(HGN *hgn);

#ifdef __cplusplus
}
#endif

#endif
