/******************************************************************************\
** File I/O (fio) Routines
**
** This module abstracts access to different format files such that a code can
** easily read all supported formats.
**
** To read a file, generally the generic "fioOpen" routine is called and the
** format is automatically detected by examining the file header.
**
** The following routines are most useful.
**
**   fioOpen     - Open the specified file (autodetects format)
**   fioClose    - Closes an open file.
**   fioGetN     - Total number of particles (or # of a specific species)
**   fioSeek     - Advance in the file to the specified particle
**   fioSpecies  - Returns the species of the next particle to read
**   fioReadDark - Reads a dark particle
**   fioReadSph  - Reads an SPH particle
**   fioReadStar - Reads a star particle
**
** Example (sequential read of all particles):
**
**   FIO fio;
**   uint64_t N, i;
**
**   fio = fioOpen("test.std");
**   if (fio==NULL) ...
**
**   N = fioGetN(fio,FIO_SPECIES_ALL);
**   for( i=0; i<N; i++ ) {
**       switch(fioSpecies(fio)) {
**       case FIO_SPECIES_DARK:
**           fioReadDark(fio,...);
**           break;
**       case FIO_SPECIES_SPH:
**           fioReadSph(fio,...);
**           break;
**       case FIO_SPECIES_STAR:
**           fioReadStar(fio,...);
**           break;
**           }
**       default:
**           perror("invalid/unknown particle type");
**           }
**       }
**
**  fioClose(fio);
**
** Example (read all dark matter particles)
**
**   FIO fio;
**   uint64_t nDark, i;
**
**   fio = fioOpen("test.std");
**   if (fio==NULL) ...
**
**   nDark = fioGetN(fio,FIO_SPECIES_DARK);
**   fioSeek(fio,0,FIO_SPECIES_DARK);
**   for( i=0; i<nDark; i++ ) {
**       fioReadDark(fio,...);
**       }
**
**   fioClose(fio);
**
\******************************************************************************/
#ifndef FIO_H
#define FIO_H

#include <stdint.h>
#include <assert.h>

/*
** These are the valid file formats
*/
typedef enum {
    FIO_FORMAT_TIPSY,
    FIO_FORMAT_HDF5,
    FIO_FORMAT_GRAFIC
    } FIO_FORMAT;

/*
** Here are the valid flags for Create
*/
#define FIO_FLAG_DOUBLE_POS    1
#define FIO_FLAG_DOUBLE_VEL    2
#define FIO_FLAG_COMPRESS_MASS 4
#define FIO_FLAG_COMPRESS_SOFT 8
#define FIO_FLAG_CHECKPOINT    16 /* Restart - normally double */
#define FIO_FLAG_POTENTIAL     32 /* Include the potential */
#define FIO_FLAG_DENSITY       64 /* Include the density */

typedef enum {
    FIO_MODE_READING,
    FIO_MODE_WRITING
    } FIO_MODE;

/*
** These are the valid data types for attributes.
*/
typedef enum {
    FIO_TYPE_FLOAT=0,
    FIO_TYPE_DOUBLE,
    FIO_TYPE_UINT32,
    FIO_TYPE_UINT64,
    FIO_TYPE_UINT8,
    FIO_TYPE_INT,
    FIO_TYPE_STRING,
    } FIO_TYPE;

/*
** Particle species are a special thing.  Codes must know explicitly what data
** each particle will contain, so we can enumerate them here.
*/
typedef enum {
    FIO_SPECIES_ALL=0,
    FIO_SPECIES_DARK,
    FIO_SPECIES_SPH,
    FIO_SPECIES_STAR,
    FIO_SPECIES_LAST /* Must be last */
    } FIO_SPECIES;

typedef uint64_t fioSpeciesList[FIO_SPECIES_LAST];

typedef struct {
    uint64_t       iFirst;      /* Starting particle index */
    char *         pszFilename; /* Filename of this file */
    fioSpeciesList nSpecies;    /* # of each species in this file */ 
    } fioFileInfo;

typedef struct {
    int iFile;                /* Current file */
    int nFiles;               /* Total number of files */
    fioFileInfo *fileInfo;    /* Array of information for each file */
    } fioFileList;

/* This structure should be treated as PRIVATE.  Call the "fio" routines. */
typedef struct fioInfo {
    FIO_FORMAT eFormat;
    FIO_MODE   eMode;
    int        mFlags;
    fioSpeciesList nSpecies;
    //uint64_t nSpecies[FIO_SPECIES_LAST];

    /* This is for multi-file support */
    fioFileList fileList;

    void (*fcnClose)(struct fioInfo *fio);
    int  (*fcnSeek) (struct fioInfo *fio,uint64_t iPart,FIO_SPECIES eSpecies);
    FIO_SPECIES (*fcnSpecies) (struct fioInfo *fio);

    int  (*fcnReadDark) (struct fioInfo *fio,
	uint64_t *piOrder,double *pdPos,double *pdVel,
			 float *pfMass,float *pfSoft,float *pfPot,float *pfDen);
    int  (*fcnReadSph) (
	struct fioInfo *fio,uint64_t *piOrder,double *pdPos,double *pdVel,
	float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
	float *pfTemp, float *pfMetals);
    int  (*fcnReadStar) (struct fioInfo *fio,
	uint64_t *piOrder,double *pdPos,double *pdVel,
	 float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
	float *pfMetals, float *pfTform);

    int  (*fcnWriteDark) (struct fioInfo *fio,
	uint64_t iOrder,const double *pdPos,const double *pdVel,
	float fMass,float fSoft,float fPot,float fDen);
    int  (*fcnWriteSph) (
	struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
	float fMass,float fSoft,float fPot,float fDen,
	float fTemp,float fMetals);
    int  (*fcnWriteStar) (struct fioInfo *fio,
	uint64_t iOrder,const double *pdPos,const double *pdVel,
	float fMass,float fSoft,float fPot,float fDen,
	float fMetals,float fTform);

    int  (*fcnGetAttr)(struct fioInfo *fio,
	const char *attr, FIO_TYPE dataType, void *data);
    int  (*fcnSetAttr)(struct fioInfo *fio,
	const char *attr, FIO_TYPE dataType, void *data);
    } *FIO;

/******************************************************************************\
** Generic Routines
\******************************************************************************/

/*
** Auto-detects the file format by looking at header information.
*/
FIO fioOpen(const char *fileName,double dOmega0,double dOmegab);
FIO fioOpenMany(int nFiles, const char * const *fileNames,double dOmega0,double dOmegab);

/*
** Close an open file of any format.
*/
static inline void fioClose(struct fioInfo *fio) {
    (*fio->fcnClose)(fio);
    }

/*
** Return the format of the currently open file.
*/
static inline FIO_FORMAT fioFormat(FIO fio) {
    return fio->eFormat;
    }

/*
** Return the mode of the currently open file.
*/
static inline FIO_MODE fioMode(FIO fio) {
    return fio->eMode;
    }

/*
** Returns the number of particles of a given species (or ALL).
*/
static inline uint64_t fioGetN(FIO fio,FIO_SPECIES eSpecies) {
    assert(eSpecies>=FIO_SPECIES_ALL && eSpecies<FIO_SPECIES_LAST);
    return fio->nSpecies[eSpecies];
    }

/*
** Seek to the N'th particle or the N'th particle of a given species.
*/
static inline int fioSeek(struct fioInfo *fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    return (*fio->fcnSeek)(fio,iPart,eSpecies);
    }

/*
** Returns the species at the current file position (what will next be read)
*/
static inline FIO_SPECIES fioSpecies(struct fioInfo *fio) {
    return (*fio->fcnSpecies)(fio);
    }

/*
** Read a particle.  Must already be positioned at the appropriate particle.
** Normally fioSpecies() is called to determine which type to read, or a
** specific seek is performed to start reading a particular type.
*/
static inline int fioReadDark(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen) {
    return (*fio->fcnReadDark)(fio,piOrder,pdPos,pdVel,pfMass,pfSoft,pfPot,pfDen);
    }
static inline int  fioReadSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfTemp,float *pfMetals) {
    return (*fio->fcnReadSph)(fio,piOrder,pdPos,pdVel,pfMass,pfSoft,pfPot,pfDen,
			      pfTemp,pfMetals);
    }
static inline int fioReadStar(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfMetals,float *pfTform) {
    return (*fio->fcnReadStar)(fio,piOrder,pdPos,pdVel,pfMass,pfSoft,pfPot,pfDen,
			       pfMetals,pfTform);
    }
/*
** Write a particle.  Must already be positioned at the appropriate particle.
*/
static inline int fioWriteDark(
    FIO fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen) {
    return (*fio->fcnWriteDark)(fio,iOrder,pdPos,pdVel,fMass,fSoft,fPot,fDen);
    }
static inline int  fioWriteSph(
    FIO fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float fTemp,float fMetals) {
    return (*fio->fcnWriteSph)(fio,iOrder,pdPos,pdVel,fMass,fSoft,fPot,fDen,
			      fTemp,fMetals);
    }
static inline int fioWriteStar(
    FIO fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float fMetals,float fTform) {
    return (*fio->fcnWriteStar)(fio,iOrder,pdPos,pdVel,fMass,fSoft,fPot,fDen,
			       fMetals,fTform);
    }
/*
** Returns the value of a given attribute.  Only "dTime" is available for
** Tipsy files, but HDF5 supports the inclusion of any arbitary attribute.
*/
static inline int fioGetAttr(FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    return (*fio->fcnGetAttr)(fio,attr,dataType,data);
    }

/*
** Sets an arbitrary attribute.  Only supported for HDF5; other formats
** return 0 indicating that it was not successful.
*/
static inline int fioSetAttr(FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    return (*fio->fcnSetAttr)(fio,attr,dataType,data);
    }

static inline int fioGetFlags(FIO fio) {
    return fio->mFlags;
    }

/******************************************************************************\
** TIPSY FORMAT
\******************************************************************************/

FIO fioTipsyCreate(const char *fileName,int mFlags,int bStandard,
		   double dTime,uint64_t nSph, uint64_t nDark, uint64_t nStar);
FIO fioTipsyAppend(const char *fileName,int mFlags,int bStandard);
FIO fioTipsyCreatePart(const char *fileName,int bAppend,int mFlags,int bStandard,
		       double dTime, uint64_t nSph, uint64_t nDark, uint64_t nStar,
		       uint64_t iStart);
int fioTipsyIsDouble(FIO fio);
int fioTipsyIsStandard(FIO fio);

/******************************************************************************\
** HDF5 FORMAT
\******************************************************************************/

/*
** Create an HDF5 file.
*/
FIO fioHDF5Create(const char *fileName,int mFlags);

/******************************************************************************\
** GRAFIC FORMAT
\******************************************************************************/

#endif
