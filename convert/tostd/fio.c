#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#ifdef HAVE_GSL
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#endif

#ifndef NO_XDR
#include <rpc/types.h>
#include <rpc/xdr.h>
#else
typedef struct {
    FILE *fp;
    int bEncoding;
} XDR;

#define XDR_ENCODE 1
#define XDR_DECODE 2

void xdrstdio_create(XDR *xdr,FILE *fp,int op) {
    assert(op==XDR_ENCODE || op==XDR_DECODE);
    xdr->bEncoding = (op==XDR_ENCODE);
    xdr->fp = fp;
}

void xdr_destroy(XDR *xdr) {
}

int xdr_double(XDR *xdr, double *d) {
    uint64_t v;
    unsigned char c[8];
    union {
	double   d;
	uint64_t v;
    } u;
    if (xdr->bEncoding) {
        u.d = *d;
	v = u.v;
	c[0] = (v>>56) & 255;
	c[1] = (v>>48) & 255;
	c[2] = (v>>40) & 255;
	c[3] = (v>>32) & 255;
	c[4] = (v>>24) & 255;
	c[5] = (v>>16) & 255;
	c[6] = (v>>8)  & 255;
	c[7] = (v>>0)  & 255;
	return fwrite(c,sizeof(c),1,xdr->fp);
    }
    else {
	if (fread(c,sizeof(c),1,xdr->fp)!=1) return 0;
	u.v = 0;   u.v |= c[0];
	u.v <<= 8; u.v |= c[1];
	u.v <<= 8; u.v |= c[2];
	u.v <<= 8; u.v |= c[3];
	u.v <<= 8; u.v |= c[4];
	u.v <<= 8; u.v |= c[5];
	u.v <<= 8; u.v |= c[6];
	u.v <<= 8; u.v |= c[7];
	*d = u.d;
	return 1;
    }
}

int xdr_float(XDR *xdr, float *f) {
    uint32_t v;
    unsigned char c[4];
    union {
	float    f;
	uint32_t v;
    } u;
    if (xdr->bEncoding) {
        u.f = *f;
	v = u.v;
	c[0] = (v>>24) & 255;
	c[1] = (v>>16) & 255;
	c[2] = (v>>8)  & 255;
	c[3] = (v>>0)  & 255;
	return fwrite(c,sizeof(c),1,xdr->fp);
    }
    else {
	if (fread(c,sizeof(c),1,xdr->fp)!=1) return 0;
	u.v = (c[0]<<24) | (c[1]<<16) | (c[2]<<8) | c[3];
	*f = u.f;
	return 1;
    }
}

int xdr_u_int(XDR *xdr, uint32_t *u) {
    uint32_t v;
    unsigned char c[4];
    if (xdr->bEncoding) {
        v = *u;
	c[0] = (v>>24) & 255;
	c[1] = (v>>16) & 255;
	c[2] = (v>>8)  & 255;
	c[3] = (v>>0)  & 255;
	return fwrite(c,sizeof(c),1,xdr->fp);
    }
    else {
	if (fread(c,sizeof(c),1,xdr->fp)!=1) return 0;
	v = (c[0]<<24) | (c[1]<<16) | (c[2]<<8) | c[3];
	*u = v;
	return 1;
    }
}
#endif

#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
#include <wordexp.h>
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
#include <glob.h>
#endif

#include "fio.h"

/*
** This uses the best available seek routine to move to the specified
** 64-bit offset into a file.  The new "fseeko" function is preferred,
** but we fallback to multiple fseek calls if required.
*/
static int safe_fseek(FILE *fp,uint64_t lStart) {
#ifdef HAVE_FSEEKO_ACTUALLY_OFTEN_BORKED_SO_DONT_USE_THIS
    uint64_t lSeeked;
    int iErr;
    errno = 0;
    iErr = fseeko(fp,lStart,SEEK_SET);
    if (iErr) {
	perror("safe_fseek failed");
	exit(errno);
	}
    lSeeked = ftello(fp);
    if(lSeeked != lStart) {
	fprintf(stderr,"fseeko() borked: %lu != %lu\n", lStart, lSeeked);
	exit(ESPIPE);
	}
    return 0;
#else
    /* fseek fails for offsets >= 2**31; this is an ugly workaround */
    static const off_t MAX_OFFSET = 2147483640;
    int iErr;

    if (lStart > MAX_OFFSET) {
	errno = 0;
	rewind(fp);
	if (errno) {
	    perror("safe_fseek failed");
	    exit(errno);
	    }
	while (lStart > MAX_OFFSET) {
	    fseek(fp,MAX_OFFSET,SEEK_CUR);
	    lStart -= MAX_OFFSET;
	    }
	iErr = fseek(fp,lStart,SEEK_CUR);
	if (iErr) {
	    perror("safe_fseek failed");
	    exit(errno);
	    }
	}
    else {
	iErr = fseek(fp,lStart,SEEK_SET);
	if (iErr) {
	    perror("safe_fseek failed");
	    exit(errno);
	    }
	}
    return 0;
#endif
    }

/*
** This is used for the HDF5 group name.
*/
const char *fioSpeciesName(FIO_SPECIES eSpecies) {
    switch(eSpecies) {
    case FIO_SPECIES_DARK: return "dark";
    case FIO_SPECIES_SPH:  return "sph";
    case FIO_SPECIES_STAR: return "star";
    default:
	break;
	}
    return NULL;
    }

static void fioFree(FIO fio) {
    if (fio->fileList.fileInfo) {
	if (fio->fileList.fileInfo[0].pszFilename) free(fio->fileList.fileInfo[0].pszFilename);
	free(fio->fileList.fileInfo);
	}
    }

static void fileScanFree(fioFileList *list) {
    if (list->fileInfo) {
	free(list->fileInfo[0].pszFilename);
	free(list->fileInfo);
	}
    }

/*
** Given a list of one or more files, this function will expand any wildcards
** present (if possible) and return a complete list of matching files.
*/
static int fileScan( fioFileList *list, int nFiles, const char * const *szFilenames) {
    int i, nScan, iIdx, nSize;
    char *pszFilename;
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
    wordexp_t files;
    int flag = 0;
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
    glob_t files;
    int flag = GLOB_ERR;
#endif

    /*
    ** Run through the list of files and expand them producing a total
    ** count of files, and the list in whatever format we can.
    */
    nSize = 0;
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
    for( iIdx = 0; iIdx<nFiles; iIdx++ ) {
	wordexp(szFilenames[iIdx], &files, flag);
	if ( files.we_wordc <= 0 ) {
	    return iIdx+1;
	    }
	flag |= WRDE_APPEND;
	}
    nScan = files.we_wordc;
    for(i=0; i<nScan; i++)
	nSize += strlen(files.we_wordv[i])+1;
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
    for( iIdx = 0; iIdx<nFiles; iIdx++ ) {
	if (glob(szFilenames[iIdx],flag,NULL,&files) || files.gl_pathc==0) {
	    return iIdx+1;
	    }
	flag |= GLOB_APPEND;
	}
    nScan = files.gl_pathc;
    for(i=0; i<nScan; i++)
	nSize += strlen(files.gl_pathv[i])+1;
#else
    nScan = nFiles;
    for( iIdx = 0; iIdx<nFiles; iIdx++ )
	nSize += strlen(szFilenames[iIdx]) + 1;
#endif

    /*
    ** Allocate space for the list of files.  We use a single buffer for the
    ** names -- we have already calculated how big this needs to be.
    */
    list->nFiles = nScan;
    list->fileInfo = malloc(sizeof(fioFileInfo)*(nScan+1));
    assert(list->fileInfo);
    pszFilename = malloc(nSize);
    assert(pszFilename);

    /*
    ** Great, now just copy the filenames into the new buffer
    */
    for( i=0; i<nScan; i++ ) {
	list->fileInfo[i].pszFilename = pszFilename;
	list->fileInfo[i].iFirst = 0;
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
	strcpy( pszFilename, files.we_wordv[i] );
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
	strcpy( pszFilename, files.gl_pathv[i] );
#else
	strcpy( pszFilename, szFilenames[i] );
#endif
	pszFilename += strlen(pszFilename) + 1;
	}

    /* Finished: free memory used by wildcard routines */
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
    wordfree(&files);
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
    globfree(&files);
#endif
    return 0;
    }

/******************************************************************************\
** Stubs that are called when a feature is not supported.
\******************************************************************************/

static int  fioNoReadDark(
    struct fioInfo *fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen) {
    fprintf(stderr,"Reading dark particles is not supported\n");
    abort();
    return 0;
    }

static int fioNoReadSph(
    struct fioInfo *fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft, float *pfPot,float *pfDen,
    float *pfTemp, float *pfMetals) {
    fprintf(stderr,"Reading SPH particles is not supported\n");
    abort();
    return 0;
    }

static int fioNoReadStar(
    struct fioInfo *fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,float *pfMetals, float *pfTform) {
    fprintf(stderr,"Reading star particles is not supported\n");
    abort();
    return 0;
    }

static int  fioNoWriteDark(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen) {
    fprintf(stderr,"Writing dark particles is not supported\n");
    abort();
    return 0;
    }

static int fioNoWriteSph(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float fTemp,float fMetals) {
    fprintf(stderr,"Writing SPH particles is not supported\n");
    abort();
    return 0;
    }

static int fioNoWriteStar(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,float fMetals,float fTform) {
    fprintf(stderr,"Writing star particles is not supported\n");
    abort();
    return 0;
    }

static int fioNoGetAttr(
    FIO fio,const char *attr, FIO_TYPE dataType, void *data) {
    return 0;
    }

static int fioNoSetAttr(
    FIO fio,const char *attr, FIO_TYPE dataType, void *data) {
    return 0;
    }

static int fioNoSeek(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fprintf(stderr,"Seeking is not supported\n");
    abort();
    return 0;
    }

static void fioNoClose(FIO fio) {
    fprintf(stderr,"Closing the file is not supported -- seriously?\n");
    abort();
    }

static FIO_SPECIES fioNoSpecies(FIO fio) {
    return FIO_SPECIES_LAST;
    }

/******************************************************************************\
** Generic Initialization - Provide default functions where possible
\******************************************************************************/

static void fioInitialize(FIO fio, FIO_FORMAT eFormat, FIO_MODE eMode, int mFlags) {
    int i;

    fio->eFormat = eFormat;
    fio->eMode   = eMode;
    fio->mFlags  = mFlags;

    for(i=0; i<FIO_SPECIES_LAST; i++) fio->nSpecies[i] = 0;

    fio->fileList.fileInfo= NULL;
    fio->fileList.iFile   = 0;
    fio->fileList.nFiles  = 0;

    /* At least some of these need to be overridden later */
    fio->fcnClose     = fioNoClose;
    fio->fcnSeek      = fioNoSeek;
    fio->fcnSpecies   = fioNoSpecies;
    fio->fcnReadDark  = fioNoReadDark;
    fio->fcnReadSph   = fioNoReadSph;
    fio->fcnReadStar  = fioNoReadStar;
    fio->fcnWriteDark = fioNoWriteDark;
    fio->fcnWriteSph  = fioNoWriteSph;
    fio->fcnWriteStar = fioNoWriteStar;
    fio->fcnGetAttr   = fioNoGetAttr;
    fio->fcnSetAttr   = fioNoSetAttr;
    }

/*
** Runs through the file specific species counts and accumulates them.
** The file specific total is set, as well as the global totals.
*/
static void fioTabulateSpecies(FIO fio) {
    int iFile, iSpec;

    assert(FIO_SPECIES_ALL==0);

    /* We will accumulate, so zero the global totals */
    for( iSpec=0; iSpec<FIO_SPECIES_LAST; iSpec++)
	fio->nSpecies[iSpec] = 0;

    /* Now accumulate for each file */
    for( iFile=0; iFile<fio->fileList.nFiles; iFile++ ) {
	fio->fileList.fileInfo[iFile].nSpecies[FIO_SPECIES_ALL] = 0;
	for( iSpec=1; iSpec<FIO_SPECIES_LAST; iSpec++) {
	    fio->nSpecies[iSpec] += fio->fileList.fileInfo[iFile].nSpecies[iSpec];
	    fio->fileList.fileInfo[iFile].nSpecies[FIO_SPECIES_ALL] += fio->fileList.fileInfo[iFile].nSpecies[iSpec];
	    }
	fio->nSpecies[FIO_SPECIES_ALL] += fio->fileList.fileInfo[iFile].nSpecies[FIO_SPECIES_ALL];
	}
    }


/******************************************************************************\
** TIPSY FORMAT
\******************************************************************************/

#define TIO_BUFFER_SIZE (1024*1024)

typedef struct {
    double dTime;
    unsigned nBodies;
    unsigned nDim;
    unsigned nSph;
    unsigned nDark;
    unsigned nStar;
    unsigned nPad;
    } tipsyHdr;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
    } tipsyDark;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
    } tipsySph;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
    } tipsyStar;

#define TIPSY_FD_SIZE(t,d) (sizeof(float)*((t)-(d))+sizeof(double)*(d))
#define TIPSY_DARK_SIZE(dark,dbl) ((dark)*TIPSY_FD_SIZE(9,dbl))
#define TIPSY_SPH_SIZE(sph,dbl)   ((sph)*TIPSY_FD_SIZE(12,dbl))
#define TIPSY_STAR_SIZE(star,dbl) ((star)*TIPSY_FD_SIZE(11,dbl))
#define TIPSY_TOTAL_SIZE(dark,gas,star,dbl) sizeof(tipsyHdr)+TIPSY_DARK_SIZE(dark,dbl)+TIPSY_SPH_SIZE(gas,dbl)+TIPSY_STAR_SIZE(star,dbl)

typedef struct {
    struct fioInfo fio; /* "base class" */
    FILE *fp;
    double dTime;
    char *fpBuffer;
    XDR xdr;
    off_t nHdrSize;     /* Size of header; 0 if second or subsequent file */
    uint64_t iOrder;    /* Current particle index */
    int bDoublePos;
    int bDoubleVel;
    } fioTipsy;

static int xdrHeader(XDR *pxdr,tipsyHdr *ph) {
    if (!xdr_double(pxdr,&ph->dTime)) return 0;
    if (!xdr_u_int(pxdr,&ph->nBodies)) return 0;
    if (!xdr_u_int(pxdr,&ph->nDim)) return 0;
    if (!xdr_u_int(pxdr,&ph->nSph)) return 0;
    if (!xdr_u_int(pxdr,&ph->nDark)) return 0;
    if (!xdr_u_int(pxdr,&ph->nStar)) return 0;
    if (!xdr_u_int(pxdr,&ph->nPad)) return 0;
    return 1;
    }

static int tipsyGetAttr(FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if ( strcmp(attr,"dTime")==0 ) {
	if (!tio->nHdrSize) return 0;
	switch(dataType) {
	case FIO_TYPE_FLOAT: *(float *)(data) = tio->dTime; break;
	case FIO_TYPE_DOUBLE:*(double *)(data) = tio->dTime; break;
	default: return 0;
	    }
	return 1;
	}

    return 0;
    }

static FIO_SPECIES tipsySpecies(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tio->iOrder<tio->fio.nSpecies[FIO_SPECIES_SPH]) return FIO_SPECIES_SPH;
    else if (tio->iOrder<tio->fio.nSpecies[FIO_SPECIES_SPH]+tio->fio.nSpecies[FIO_SPECIES_DARK])
	return FIO_SPECIES_DARK;
    else if (tio->iOrder<tio->fio.nSpecies[FIO_SPECIES_SPH]+tio->fio.nSpecies[FIO_SPECIES_DARK]+tio->fio.nSpecies[FIO_SPECIES_STAR])
	return FIO_SPECIES_STAR;
    else return FIO_SPECIES_LAST;
    }

/*
** Given an offset, calculate the particle number
*/
static uint64_t tipsyParticle(uint64_t iByte,uint64_t nHdrSize,
    uint64_t nSph, uint64_t nDark, uint64_t nStar,
    int nSphSize, int nDarkSize, int nStarSize ) {
    uint64_t iPart,n;

    iPart = 0;
    iByte -= nHdrSize;

    n = iByte / nSphSize;
    if ( n > nSph ) n = nSph;
    iPart += n;
    iByte -= nSphSize * n;

    n = iByte / nDarkSize;
    if ( n > nDark ) n = nDark;
    iPart += n;
    iByte -= nDarkSize * n;

    n = iByte / nStarSize;
    assert( n<=nStar );
    iPart += n;
    iByte -= nStarSize * n;

    assert(iByte==0);

    return iPart;
    }


/*
** Calculate the absolute offset of a given particle in a Tipsy file.
*/
static uint64_t tipsyOffset(
    uint64_t iPart,uint64_t nHdrSize,
    uint64_t nSph, uint64_t nDark, uint64_t nStar,
    int nSphSize, int nDarkSize, int nStarSize ) {
    uint64_t iByte;

    iByte = nHdrSize;
    if (iPart<nSph) {
	iByte += iPart * nSphSize;
	}
    else {
	iByte += nSph * nSphSize;
	iPart -= nSph;

	if (iPart<nDark) {
	    iByte += iPart * nDarkSize;
	    }
	else {
	    iByte += nDark * nDarkSize;
	    iPart -= nDark;

	    if (iPart<=nStar) {
		iByte += iPart * nStarSize;
		}
	    else {
		errno = ESPIPE;
		return 0;
		}
	    }
	}
    return iByte;
    }

/*
** Compares the current iOrder with the file list and switches to a different
** file if necessary.
*/
static int tipsySwitchFile(FIO fio, int bSeek) {
    fioTipsy *tio = (fioTipsy *)fio;

    assert(fio->fileList.iFile>=0 && fio->fileList.iFile<fio->fileList.nFiles);
    if (tio->iOrder<fio->fileList.fileInfo[fio->fileList.iFile].iFirst
	|| tio->iOrder>=fio->fileList.fileInfo[fio->fileList.iFile+1].iFirst ) {
	int bStandard = fioTipsyIsStandard(fio);
	if (bStandard) xdr_destroy(&tio->xdr);
	fclose(tio->fp);
	if (tio->iOrder < fio->fileList.fileInfo[fio->fileList.iFile].iFirst) {
	    while(fio->fileList.iFile!=0)
		if (fio->fileList.fileInfo[--fio->fileList.iFile].iFirst<tio->iOrder) break;
	    }
	else if (tio->iOrder >= fio->fileList.fileInfo[fio->fileList.iFile+1].iFirst) {
	    while(++fio->fileList.iFile<fio->fileList.nFiles)
		if (fio->fileList.fileInfo[fio->fileList.iFile+1].iFirst>tio->iOrder) break;
	    }
	assert(tio->iOrder >= fio->fileList.fileInfo[fio->fileList.iFile].iFirst
	       && tio->iOrder < fio->fileList.fileInfo[fio->fileList.iFile+1].iFirst);
	tio->fp = fopen(fio->fileList.fileInfo[fio->fileList.iFile].pszFilename,"rb");
	if (tio->fp == NULL) printf("Fail: %d %s\n",fio->fileList.iFile,fio->fileList.fileInfo[fio->fileList.iFile].pszFilename);
	if (tio->fp == NULL) return 1;
	if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);
	if (bStandard) xdrstdio_create(&tio->xdr,tio->fp,XDR_DECODE);
	bSeek = 1;
	}
    if (bSeek) {
	uint64_t iByte;
	off_t iOffset;
	int rc;

	/* For file fragments, the first particle is not 0.  Adjust accordingly. */
	iByte = tipsyOffset(tio->iOrder,0,
			    tio->fio.nSpecies[FIO_SPECIES_SPH],
			    tio->fio.nSpecies[FIO_SPECIES_DARK],
			    tio->fio.nSpecies[FIO_SPECIES_STAR],
			    TIPSY_SPH_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
			    TIPSY_DARK_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
			    TIPSY_STAR_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel))
	    - tipsyOffset(fio->fileList.fileInfo[fio->fileList.iFile].iFirst,0,
			  tio->fio.nSpecies[FIO_SPECIES_SPH],
			  tio->fio.nSpecies[FIO_SPECIES_DARK],
			  tio->fio.nSpecies[FIO_SPECIES_STAR],
			  TIPSY_SPH_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
			  TIPSY_DARK_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
			  TIPSY_STAR_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel))
	    + (fio->fileList.iFile?0:tio->nHdrSize);
	iOffset = iByte;
	assert(iOffset==iByte);
	rc = safe_fseek(tio->fp,iByte); assert(rc==0);
	}

    return 0;
    }

/* DARK PARTICLES */
static int tipsyReadNativeDark(FIO fio,
    uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];
    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);

    if (tipsySwitchFile(fio,0)) return 0;
    *piOrder = tio->iOrder++;

    rc = fread(pfMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fread(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdPos[d] = fTmp[d];
	}
    if (tio->bDoubleVel) {
	rc = fread(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdVel[d] = fTmp[d];
	}
    rc = fread(pfSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    *pfDen = 0.0f;
    return 1;
    }

static int tipsyWriteNativeDark(FIO fio,
    uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
    assert(iOrder == tio->iOrder++);
    rc = fwrite(&fMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fwrite(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdPos[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    if (tio->bDoubleVel) {
	rc = fwrite(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdVel[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    rc = fwrite(&fSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyReadStandardDark(FIO fio,
    uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tipsySwitchFile(fio,0)) return 0;
    *piOrder = tio->iOrder++;
    if (!xdr_float(&tio->xdr,pfMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    if (!xdr_double(&tio->xdr,pdPos+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdPos[d] = fTmp;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    if (!xdr_double(&tio->xdr,pdVel+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdVel[d] = fTmp;
	    }
	}
    if (!xdr_float(&tio->xdr,pfSoft)) return 0;
    if (!xdr_float(&tio->xdr,pfPot)) return 0;
    *pfDen = 0.0f;
    return 1;
    }

static int tipsyWriteStandardDark(FIO fio,
    uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;
    double dTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
//    assert(iOrder == tio->iOrder++ + tio->iStart);  // JW: needs to handle non-contiguous iOrders
    if (!xdr_float(&tio->xdr,&fMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    dTmp = pdPos[d]; if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdPos[d]; if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    dTmp = pdVel[d]; if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdVel[d]; if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    if (!xdr_float(&tio->xdr,&fSoft)) return 0;
    if (!xdr_float(&tio->xdr,&fPot)) return 0;
    return 1;
    }

/* SPH PARTICLES */
static int tipsyReadNativeSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft, float *pfPot,float *pfDen,
    float *pfTemp, float *pfMetals) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tipsySwitchFile(fio,0)) return 0;
    *piOrder = tio->iOrder++;
    rc = fread(pfMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fread(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdPos[d] = fTmp[d];
	}
    if (tio->bDoubleVel) {
	rc = fread(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdVel[d] = fTmp[d];
	}
    rc = fread(pfDen,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfTemp,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfMetals,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }


static int tipsyWriteNativeSph(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float fTemp,float fMetals) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
    assert(iOrder == tio->iOrder++);

    rc = fwrite(&fMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fwrite(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdPos[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    if (tio->bDoubleVel) {
	rc = fwrite(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdVel[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    rc = fwrite(&fDen,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fTemp,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fMetals,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyReadStandardSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft, float *pfPot,float *pfDen,
    float *pfTemp, float *pfMetals) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tipsySwitchFile(fio,0)) return 0;
    *piOrder = tio->iOrder++;
    if (!xdr_float(&tio->xdr,pfMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    if (!xdr_double(&tio->xdr,pdPos+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdPos[d] = fTmp;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    if (!xdr_double(&tio->xdr,pdVel+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdVel[d] = fTmp;
	    }
	}
    if (!xdr_float(&tio->xdr,pfDen)) return 0;
    if (!xdr_float(&tio->xdr,pfTemp)) return 0;
    if (!xdr_float(&tio->xdr,pfSoft)) return 0;
    if (!xdr_float(&tio->xdr,pfMetals)) return 0;
    if (!xdr_float(&tio->xdr,pfPot)) return 0;
    return 1;
    }

static int tipsyWriteStandardSph(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float fTemp,float fMetals) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;
    double dTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
//    assert(iOrder == tio->iOrder++ + tio->iStart); //JW -- non-contiguous iOrder fix needed

    if (!xdr_float(&tio->xdr,&fMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    dTmp = pdPos[d];
	    if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdPos[d];
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    dTmp = pdVel[d];
	    if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdVel[d];
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    if (!xdr_float(&tio->xdr,&fDen)) return 0;
    if (!xdr_float(&tio->xdr,&fTemp)) return 0;
    if (!xdr_float(&tio->xdr,&fSoft)) return 0;
    if (!xdr_float(&tio->xdr,&fMetals)) return 0;
    if (!xdr_float(&tio->xdr,&fPot)) return 0;
    return 1;
    }

/* STAR PARTICLES */
static int tipsyReadNativeStar(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfMetals, float *pfTform) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tipsySwitchFile(fio,0)) return 0;
    *piOrder = tio->iOrder++;
    rc = fread(pfMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fread(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdPos[d] = fTmp[d];
	}
    if (tio->bDoubleVel) {
	rc = fread(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdVel[d] = fTmp[d];
	}
    rc = fread(pfMetals,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfTform,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    *pfDen = 0.0f;
    return 1;
    }

static int tipsyWriteNativeStar(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,float fMetals,float fTform) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
    assert(iOrder == tio->iOrder++);
    rc = fwrite(&fMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fwrite(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdPos[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    if (tio->bDoubleVel) {
	rc = fwrite(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdVel[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    rc = fwrite(&fMetals,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fTform,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyReadStandardStar(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfMetals, float *pfTform) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tipsySwitchFile(fio,0)) return 0;
    *piOrder = tio->iOrder++;
    if (!xdr_float(&tio->xdr,pfMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    if (!xdr_double(&tio->xdr,pdPos+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdPos[d] = fTmp;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    if (!xdr_double(&tio->xdr,pdVel+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdVel[d] = fTmp;
	    }
	}
    if (!xdr_float(&tio->xdr,pfMetals)) return 0;
    if (!xdr_float(&tio->xdr,pfTform)) return 0;
    if (!xdr_float(&tio->xdr,pfSoft)) return 0;
    if (!xdr_float(&tio->xdr,pfPot)) return 0;
    *pfDen = 0.0f;
    return 1;
    }

static int tipsyWriteStandardStar(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,float fMetals,float fTform) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;
    double dTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
//    assert(iOrder == tio->iOrder++ + tio->iStart); // JW non-contiguous iOrder fix needed
    if (!xdr_float(&tio->xdr,&fMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    dTmp = pdPos[d];
	    if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdPos[d];
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    dTmp = pdVel[d];
	    if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdVel[d];
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    if (!xdr_float(&tio->xdr,&fMetals)) return 0;
    if (!xdr_float(&tio->xdr,&fTform)) return 0;
    if (!xdr_float(&tio->xdr,&fSoft)) return 0;
    if (!xdr_float(&tio->xdr,&fPot)) return 0;
    return 1;
    }

static void tipsyCloseNative(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY);
    fclose(tio->fp);
    if (tio->fpBuffer) free(tio->fpBuffer);
    fileScanFree(&fio->fileList);
    free(tio);
    }

static void tipsyCloseStandard(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY);
    xdr_destroy(&tio->xdr);
    tipsyCloseNative(fio);
    }

int fioTipsyIsDouble(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY);
    return tio->bDoublePos;
    }

int fioTipsyIsStandard(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY);
    return tio->fio.fcnReadDark == tipsyReadStandardDark;
    }

static int tipsySeek(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;

    switch(eSpecies) {
    case FIO_SPECIES_STAR: iPart += tio->fio.nSpecies[FIO_SPECIES_DARK];
    case FIO_SPECIES_DARK: iPart += tio->fio.nSpecies[FIO_SPECIES_SPH];
    case FIO_SPECIES_SPH:  break;
    case FIO_SPECIES_ALL:  break;
    default:
	fprintf(stderr,"Invalid particle species\n");
	exit(1);
	}
    tio->iOrder = iPart;

    /* Seek to a different file if necessary */
    rc = tipsySwitchFile(fio,1); assert(rc==0);
    if (rc) return rc;

    return 0;
    }

static int tipsySeekStandard(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    /*fioTipsy *tio = (fioTipsy *)fio;*/
    int rc;
    /*assert(fio->eFormat == FIO_FORMAT_TIPSY);*/
    /*xdr_destroy(&tio->xdr);*/
    rc = tipsySeek(fio,iPart,eSpecies);
    /*xdrstdio_create(&tio->xdr,tio->fp,fio->eMode==FIO_MODE_READING?XDR_DECODE:XDR_ENCODE);*/
    return rc;
    }

static uint64_t tipsyCountSpecies(uint64_t N, uint64_t *iStart, uint64_t *iLast) {
    uint64_t n;

    if ( *iStart >= N ) {
	n = 0;
	*iStart -= N; *iLast -= N;
	}
    else if ( *iLast > N) {
	n = N - *iStart;
	*iStart = 0; *iLast -= N;
	}
    else {
	n = *iLast - *iStart;
	*iStart = *iLast = 0;
	}

    return n;
    }

static void tipsySetFunctions(fioTipsy *tio, int mFlags, int bStandard) {
    tio->fio.fcnGetAttr = tipsyGetAttr;
    tio->fio.fcnSpecies = tipsySpecies;

    tio->bDoubleVel = (mFlags&FIO_FLAG_DOUBLE_VEL) != 0;
    tio->bDoublePos = (mFlags&FIO_FLAG_DOUBLE_POS) != 0 || tio->bDoubleVel;

    if ( bStandard ) {
	tio->fio.fcnClose    = tipsyCloseStandard;
	tio->fio.fcnSeek     = tipsySeekStandard;
	tio->fio.fcnReadDark = tipsyReadStandardDark;
	tio->fio.fcnReadSph  = tipsyReadStandardSph;
	tio->fio.fcnReadStar = tipsyReadStandardStar;
	tio->fio.fcnWriteDark= tipsyWriteStandardDark;
	tio->fio.fcnWriteSph = tipsyWriteStandardSph;
	tio->fio.fcnWriteStar= tipsyWriteStandardStar;
	}
    else {
	tio->fio.fcnClose    = tipsyCloseNative;
	tio->fio.fcnSeek     = tipsySeek;
	tio->fio.fcnReadDark = tipsyReadNativeDark;
	tio->fio.fcnReadSph  = tipsyReadNativeSph;
	tio->fio.fcnReadStar = tipsyReadNativeStar;
	tio->fio.fcnWriteDark= tipsyWriteNativeDark;
	tio->fio.fcnWriteSph = tipsyWriteNativeSph;
	tio->fio.fcnWriteStar= tipsyWriteNativeStar;
	}
    }

/*
** The pad field can extend the 32-bit integers to 40 bits, increasing
** the maximum number of particles from 4e9 to 1e12 (1 trillion).
** This checks if that is what has happened here.
*/
static void tipsySussHeader(tipsyHdr *h,uint64_t *pN, uint64_t *pDark,uint64_t *pSph, uint64_t *pStar) {
    *pN = h->nPad & 0x000000ff;
    *pN <<= 32;
    *pN += h->nBodies;

    *pSph = h->nPad & 0x0000ff00;
    *pSph <<= 24;
    *pSph += h->nSph;

    *pDark = h->nPad & 0x00ff0000;
    *pDark <<= 16;
    *pDark += h->nDark;

    *pStar = h->nPad & 0xff000000;
    *pStar <<= 8;
    *pStar += h->nStar;

    /* There may be junk in the pad field, in which case we ignore it */
    if ( *pN != *pDark + *pSph + *pStar ) {
	*pN = h->nBodies;
	*pDark = h->nDark;
	*pSph  = h->nSph;
	*pStar = h->nStar;
	}
    }


/*
** Read the header and detect native or standard
*/
static int tipsyDetectHeader(fioTipsy *tio, int mFlags) {
    tipsyHdr h;
    uint64_t N,nDark,nSph,nStar;
    int rc;

    assert(tio->fio.eFormat == FIO_FORMAT_TIPSY);
    rc = fread(&h,sizeof(h),1,tio->fp);
    if (rc!=1) {
	fprintf(stderr,"Error reading Tipsy header\n");
	return 0;
	}

    /* Check for native binary format */
    if (h.nDim>=1 && h.nDim<=3) {
	tipsySussHeader(&h,&N,&nDark,&nSph,&nStar);
	if ( nDark + nSph + nStar != N ) {
	    fprintf(stderr,"Tipsy Header mismatch:"
		    " nDim=%u nDark=%lu nSph=%lu nStar=%lu nBodies=%lu\n",
		    h.nDim, nDark, nSph, nStar, N);
	    return 0;
	    }
	tipsySetFunctions(tio,mFlags,0);
	}

    /* Okay, this must be "standard" format */
    else {
	rewind(tio->fp);
	xdrstdio_create(&tio->xdr,tio->fp,XDR_DECODE);
	rc = xdrHeader(&tio->xdr,&h);
	if (rc!=1) {
	    perror("Error reading tipsy header");
	    return 0;
	    }
	tipsySussHeader(&h,&N,&nDark,&nSph,&nStar);

	if (h.nDim<1 || h.nDim>3 || nDark + nSph + nStar != N ) {
	    fprintf(stderr,"Tipsy Header mismatch:"
		    " nDim=%u nDark=%lu nSph=%lu nStar=%lu nBodies=%lu\n",
		    h.nDim, nDark, nSph, nStar, N);
	    xdr_destroy(&tio->xdr);
	    return 0;
	    }
	tipsySetFunctions(tio,mFlags,1);
	}
    tio->nHdrSize = sizeof(h);
    tio->dTime = h.dTime;
    tio->fio.nSpecies[FIO_SPECIES_DARK] = nDark;
    tio->fio.nSpecies[FIO_SPECIES_SPH]  = nSph;
    tio->fio.nSpecies[FIO_SPECIES_STAR] = nStar;
    return 1;
    }

static int tipsyWriteHeader(fioTipsy *tio, int bStandard) {
    FIO fio = &tio->fio;
    tipsyHdr h;
    uint64_t N, nSph, nDark, nStar;
    int rc;

    nSph  = fioGetN(fio,FIO_SPECIES_SPH);
    nDark = fioGetN(fio,FIO_SPECIES_DARK);
    nStar = fioGetN(fio,FIO_SPECIES_STAR);
    N = nSph + nDark + nStar;

    h.dTime = tio->dTime;
    h.nBodies = N;
    h.nDim    = 3;
    h.nSph    = nSph;
    h.nDark   = nDark;
    h.nStar   = nStar;
    h.nPad    = ((N>>32)&0xff) + ((nSph>>24)&0xff00)
	+ ((nDark>>16)&0xff0000) + ((nStar>>8)&0xff000000);

    tio->nHdrSize = sizeof(h);
    if (bStandard) {
	rc = xdrHeader(&tio->xdr,&h);
	if (rc!=1) {
	    perror("Error writing tipsy header");
	    return 0;
	    }
	}
    else {
	rc = fwrite(&h,sizeof(h),1,tio->fp);
	if (rc!=1) {
	    fprintf(stderr,"Error writing Tipsy header\n");
	    return 0;
	    }
	}

    return 1;
    }

static FIO tipsyOpen(fioFileList *fileList) {
    fioTipsy *tio;
    int i;
    off_t *nSizes;
    uint64_t nSize, nDark, nSph, nStar;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    fioInitialize(&tio->fio,FIO_FORMAT_TIPSY,FIO_MODE_READING,0);

    tio->iOrder = 0;
    tio->fio.fileList = *fileList;

    tio->fp = fopen(tio->fio.fileList.fileInfo[0].pszFilename,"rb");
    tio->fio.fileList.iFile = 0;
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    tio->fpBuffer = malloc(TIO_BUFFER_SIZE);
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);

    if ( !tipsyDetectHeader(tio,/*set below*/0) ) {
	fclose(tio->fp);
	if (tio->fpBuffer) free(tio->fpBuffer);
	free(tio);
	return NULL;
	}

    tio->fio.nSpecies[FIO_SPECIES_ALL] = 0;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];

    /*
    ** We know how many particles we are supposed to have at this point, and if
    ** this is a native or standard tipsy binary.  What we don't know is how
    ** many particles are in each file.  We have to stat the files to find this.
    */
    nSizes = malloc(sizeof(off_t)*tio->fio.fileList.nFiles);
    nSize = 0;
    for( i=0; i<tio->fio.fileList.nFiles; i++) {
	struct stat s;
	/* The file/directory needs to exist */
	if ( stat(tio->fio.fileList.fileInfo[i].pszFilename,&s) != 0 ) return NULL;
#ifdef _MSC_VER
	if ( !(s.st_mode&_S_IFREG) ) return NULL;
#else
	if ( !S_ISREG(s.st_mode) ) return NULL;
#endif
	nSize += (nSizes[i] = s.st_size);
	}

    nDark = tio->fio.nSpecies[FIO_SPECIES_DARK];
    nSph  = tio->fio.nSpecies[FIO_SPECIES_SPH];
    nStar = tio->fio.nSpecies[FIO_SPECIES_STAR];

    if (TIPSY_TOTAL_SIZE(nDark,nSph,nStar,0) == nSize ) {
	tio->bDoublePos = 0;
	tio->bDoubleVel = 0;
	}
    else if (TIPSY_TOTAL_SIZE(nDark,nSph,nStar,3) == nSize ) {
	tio->bDoublePos = 1;
	tio->bDoubleVel = 0;
	}
    else if (TIPSY_TOTAL_SIZE(nDark,nSph,nStar,6) == nSize ) {
	tio->bDoublePos = 1;
	tio->bDoubleVel = 1;
	}
    else {
	fprintf(stderr,"Tipsy file size mismatch!\n");
	return NULL;
	}

    /*
    ** Okay, now we need to calculate the starting particle number for each
    ** file fragment.
    */
    tio->fio.fileList.fileInfo[0].iFirst = 0;
    nSize = 0;
    for(i=1; i<=tio->fio.fileList.nFiles; i++ ) {
	nSize += nSizes[i-1];
	tio->fio.fileList.fileInfo[i].iFirst = tipsyParticle(
	    nSize,tio->nHdrSize,nSph,nDark,nStar,
	    TIPSY_SPH_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
	    TIPSY_DARK_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
	    TIPSY_STAR_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel));
	}
    assert(tio->fio.fileList.fileInfo[tio->fio.fileList.nFiles].iFirst==tio->fio.nSpecies[FIO_SPECIES_ALL]);

    free(nSizes);

    return &tio->fio;
    }

FIO fioTipsyCreate(const char *fileName,int mFlags,int bStandard,
		   double dTime,uint64_t nSph, uint64_t nDark, uint64_t nStar) {
    fioTipsy *tio;
    int i;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    fioInitialize(&tio->fio,FIO_FORMAT_TIPSY,FIO_MODE_WRITING,mFlags);

    tio->fio.fileList.fileInfo = NULL;

    tio->fio.fileList.fileInfo = malloc(sizeof(fioFileInfo)*2);
    assert(tio->fio.fileList.fileInfo);
    tio->fio.fileList.fileInfo[0].pszFilename= strdup(fileName);
    tio->fio.fileList.nFiles  = 1;

    tio->fio.nSpecies[FIO_SPECIES_SPH]  = nSph;
    tio->fio.nSpecies[FIO_SPECIES_DARK] = nDark;
    tio->fio.nSpecies[FIO_SPECIES_STAR] = nStar;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];

    tio->fio.fileList.fileInfo[0].iFirst = 0;
    tio->fio.fileList.fileInfo[1].iFirst = tio->fio.nSpecies[FIO_SPECIES_ALL];

    tio->iOrder = 0;
    tio->dTime = dTime;

    tio->fp = fopen(fileName,"wb");
    tio->fio.fileList.iFile = 0;
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    tio->fpBuffer = malloc(TIO_BUFFER_SIZE);
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);
    if (bStandard) xdrstdio_create(&tio->xdr,tio->fp,XDR_ENCODE);

    tipsySetFunctions(tio,mFlags,bStandard);
    tipsyWriteHeader(tio,bStandard);

    return &tio->fio;
    }

FIO fioTipsyAppend(const char *fileName,int mFlags,int bStandard) {
    fioTipsy *tio;
    int i;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    fioInitialize(&tio->fio,FIO_FORMAT_TIPSY,FIO_MODE_WRITING,mFlags);

    tio->iOrder = 0;
    tio->fio.fileList.fileInfo = malloc(sizeof(fioFileInfo)*2);
    assert(tio->fio.fileList.fileInfo);
    tio->fio.fileList.fileInfo[0].pszFilename= strdup(fileName);
    tio->fio.fileList.nFiles  = 1;

    tio->fp = fopen(fileName,"rb+");
    tio->fio.fileList.iFile = 0;
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    tio->fpBuffer = malloc(TIO_BUFFER_SIZE);
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);

    if ( !tipsyDetectHeader(tio,mFlags) ) {
	fclose(tio->fp);
	if (tio->fpBuffer) free(tio->fpBuffer);
	free(tio);
	return NULL;
	}
    tio->fio.nSpecies[FIO_SPECIES_ALL] = 0;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];

    tio->fio.fileList.fileInfo[0].iFirst = 0;
    tio->fio.fileList.fileInfo[1].iFirst = tio->fio.nSpecies[FIO_SPECIES_ALL];

    if (bStandard) xdrstdio_create(&tio->xdr,tio->fp,XDR_ENCODE);
    tipsySetFunctions(tio,mFlags,bStandard);

    return &tio->fio;
    }

/*
** Create a partial Tipsy file.  A partial Tipsy file is one that has been split
** into multiple files.  Only the first part has the header.  To rebuild the
** original file, the parts need only be concatenated together.
*/
FIO fioTipsyCreatePart(const char *fileName,int bAppend,int mFlags,int bStandard,
		       double dTime, uint64_t nSph, uint64_t nDark, uint64_t nStar,
		       uint64_t iStart) {
    fioTipsy *tio;
    uint64_t iLast;
    int i;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    fioInitialize(&tio->fio,FIO_FORMAT_TIPSY,FIO_MODE_WRITING,mFlags);

    tio->fio.fileList.fileInfo = malloc(sizeof(fioFileInfo)*2);
    assert(tio->fio.fileList.fileInfo);
    tio->fio.fileList.fileInfo[0].pszFilename= strdup(fileName);
    tio->fio.fileList.nFiles  = 1;

    tio->iOrder = 0;
    tio->dTime = dTime;

    for( i=0; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[i] = 0;

    /* Adjust the local number of particles in this fragment */
    iLast = nSph + nDark + nStar;
    tio->fio.nSpecies[FIO_SPECIES_SPH]  = nSph;
    tio->fio.nSpecies[FIO_SPECIES_DARK] = nDark;
    tio->fio.nSpecies[FIO_SPECIES_STAR] = nStar;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];

    tio->fio.fileList.fileInfo[0].iFirst = iStart;
    tio->fio.fileList.fileInfo[1].iFirst = iStart + tio->fio.nSpecies[FIO_SPECIES_ALL];

    tio->fp = fopen(fileName,bAppend?"rb+":"wb");
    tio->fio.fileList.iFile = 0;
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    tio->fpBuffer = malloc(TIO_BUFFER_SIZE);
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);
    if (bStandard) xdrstdio_create(&tio->xdr,tio->fp,XDR_ENCODE);

    tipsySetFunctions(tio,mFlags,bStandard);
    if (iStart) tio->nHdrSize = 0;
    else tipsyWriteHeader(tio,bStandard);
    return &tio->fio;
    }

/******************************************************************************\
** General HDF5 support routines.  Not "fio" specific.
\******************************************************************************/
#ifdef USE_HDF5
#ifndef H5_USE_16_API
#define H5_USE_16_API 1
#endif

#include <hdf5.h>

/* Returns the number of records in a dataset */
static hsize_t getSetSize(hid_t setID) {
    hid_t spaceID;
    hsize_t dims[2], maxs[2];

    spaceID = H5Dget_space(setID); assert(spaceID>=0);
    assert( H5Sis_simple(spaceID) > 0 );
    assert( H5Sget_simple_extent_ndims(spaceID) <= 3 );
    H5Sget_simple_extent_dims(spaceID,dims,maxs);
    H5Sclose(spaceID);
    return dims[0];
    }

/* Create a dataset inside a group (positions, velocities, etc.) */
static hid_t newSet(hid_t locID, const char *name, uint64_t chunk,
		    uint64_t count, int nDims, hid_t dataType ) {
    hid_t dataProperties, dataSpace, dataSet;
    hsize_t /*@alt uint64_t[]@*/ iDims[2];
    hsize_t /*@alt uint64_t[]@*/ iMax[2];
    hsize_t rc;

    /* Create a dataset property so we can set the chunk size */
    dataProperties = H5Pcreate( H5P_DATASET_CREATE );
    assert(dataProperties >= 0);
    if (chunk) {
	iDims[0] = chunk;
	iDims[1] = 1;
	rc = H5Pset_chunk( dataProperties, nDims>1?2:1, iDims ); assert(rc>=0);
	/* Also request the FLETCHER checksum */
	rc = H5Pset_filter(dataProperties,H5Z_FILTER_FLETCHER32,0,0,NULL); assert(rc>=0);
	}


    /* And the dataspace */
    iDims[0] = count;
    iDims[1] = nDims;
    if (chunk) iMax[0] = H5S_UNLIMITED;
    else       iMax[0] = count;
    iMax[1] = nDims;

    dataSpace = H5Screate_simple( nDims>1?2:1, iDims, iMax );
    assert( dataSpace != H5I_INVALID_HID );

    /* Create the data set */
    dataSet = H5Dcreate( locID, name,
			 dataType, dataSpace, dataProperties );
    assert(dataSet != H5I_INVALID_HID);
    H5Pclose( dataProperties );
    H5Sclose( dataSpace );

    return dataSet;
    }

/* Read part of a set from disk into memory */
static void readSet(
    hid_t set_id,           /* set from which to read the data */
    void *pBuffer,          /* Buffer for the data */
    hid_t memType,          /* Data type in memory */
    hsize_t iOffset,        /* Offset into the set */
    hsize_t nRows,          /* Number of rows to read */
    hsize_t nCols ) {       /* Number of columns (normally 1 or 3) */
    hid_t memSpace, diskSpace;
    hsize_t dims[2], start[2];
    herr_t rc;

    dims[0] = nRows;
    dims[1] = nCols;
    memSpace = H5Screate_simple(nCols>1?2:1,dims,0);
    diskSpace = H5Dget_space(set_id);
    start[0] = iOffset;
    start[1] = 0;

    rc = H5Sselect_hyperslab(diskSpace,H5S_SELECT_SET,start,0,dims,0); assert(rc>=0);
    rc = H5Dread(set_id,memType,memSpace,diskSpace,H5P_DEFAULT,pBuffer); assert(rc>=0);
    rc = H5Sclose(memSpace); assert(rc>=0);
    rc = H5Sclose(diskSpace); assert(rc>=0);
    }


/* Add an attribute to a group */
static int writeAttribute( hid_t groupID, const char *name,
			    hid_t dataType, void *data ) {
    hid_t dataSpace, attrID;
    herr_t rc;

    dataSpace = H5Screate(H5S_SCALAR); assert(dataSpace!=H5I_INVALID_HID);
    attrID = H5Acreate( groupID,name,dataType,dataSpace,H5P_DEFAULT );
    assert(attrID!=H5I_INVALID_HID);
    rc = H5Awrite( attrID, dataType, data ); assert(rc>=0);
    rc = H5Aclose( attrID ); assert(rc>=0);
    rc = H5Sclose(dataSpace); assert(rc>=0);
    return 1;
    }

/* Read an attribute from a group */
static int readAttribute( hid_t groupID, const char *name,
			  hid_t dataType, void *data ) {
    hid_t attrID;
    herr_t rc;

    attrID = H5Aopen_name( groupID,name );
    if ( attrID == H5I_INVALID_HID ) return 0;
    rc = H5Aread( attrID, dataType, data ); assert(rc>=0);
    rc = H5Aclose( attrID ); assert(rc>=0);
    return 1;
    }

/* Write part of a set from memory */
static void writeSet(
    hid_t set_id,           /* set into which to write the data */
    const void *pBuffer,    /* Buffer containing the data */
    hid_t memType,          /* Data type in memory */
    hsize_t iOffset,        /* Offset into the set */
    hsize_t nRows,          /* Number of rows to write */
    hsize_t nCols ) {       /* Number of columns (normally 1 or 3) */
    hid_t memSpace, diskSpace;
    hsize_t dims[2], start[2], size[2];

    dims[0] = nRows;
    dims[1] = nCols;
    memSpace = H5Screate_simple(nCols>1?2:1,dims,0);
    size[0] = iOffset + nRows;
    size[1] = nCols;
    H5Dextend(set_id,size);
    diskSpace = H5Dget_space(set_id);
    start[0] = iOffset;
    start[1] = 0;
    H5Sselect_hyperslab(diskSpace,H5S_SELECT_SET,start,0,dims,0);
    H5Dwrite(set_id,memType,memSpace,diskSpace,H5P_DEFAULT,pBuffer);
    H5Sclose(memSpace);
    H5Sclose(diskSpace);
    }

/******************************************************************************\
** HDF5 FORMAT
** 
** The structure is as follows:
**
** fioHDF5         - The HDF5 file
**   - IOBASE[]    - A single species of particle -- an HDF5 group
**     - IOORDER   - Particle iOrder information
**     - IOCLASS   - particle class map
**     - IOFIELD[] - Each data field (positions,velocities,etc.)
\******************************************************************************/

#define CHUNK_SIZE (65536)

#define GROUP_PARAMETERS "parameters"

#define ATTR_IORDER      "iOrder"

#define FIELD_POSITION    "position"
#define FIELD_VELOCITY    "velocity"
#define FIELD_MASS        "mass"
#define FIELD_SOFTENING   "softening"
#define FIELD_POTENTIAL   "potential"
#define FIELD_DENSITY     "density"
#define FIELD_TEMPERATURE "temperature"
#define FIELD_METALS      "metals"

#define FIELD_ORDER      "order"
#define FIELD_CLASS      "class"
#define FIELD_CLASSES    "classes"

#define DARK_POSITION    0
#define DARK_VELOCITY    1
#define DARK_POTENTIAL   2
#define DARK_DENSITY     3
#define DARK_N           4

#define SPH_POSITION    0
#define SPH_VELOCITY    1
#define SPH_POTENTIAL   2
#define SPH_DENSITY     3
#define SPH_TEMPERATURE 4
#define SPH_METALS      5
#define SPH_N           6

typedef struct {
    double v[3];
    } ioHDF5V3;

typedef uint64_t PINDEX;

typedef struct {
    void *pBuffer;            /* Pointers to data */
    hid_t memType;            /* ... of this type */
    hid_t setId;              /* ... in this set */
    uint32_t nValues;         /* ... with n-dimensions */
    uint32_t nChunk;          /* ... by this size chunk */
    uint32_t iSize;           /* ... each this size */
    } IOFIELD;

typedef struct {
    uint8_t iClass;
    PINDEX  iOrderStart;      /* Start index of this class */
    double  fMass;            /* Particle mass */
    double  fSoft;            /* Softening */
    PINDEX  nCount;           /* Number of particles with this class */
    } classEntry;

typedef struct {
    classEntry    Class[256]; /* Array of class information */
    uint_fast32_t nClasses;   /* ... of this many values */
    IOFIELD       fldClasses; /* Corresponding particle class information */
    IOFIELD       fldMass;
    IOFIELD       fldSoft;
    } IOCLASS;

typedef struct {
    PINDEX iStart;            /* Start iOrder number */
    PINDEX iNext;             /* Next iOrder number (starts at iStart) */
    IOFIELD fldOrder;
    hid_t   groupID;
    } IOORDER;

typedef struct {
    PINDEX        iOffset;    /* Particle offset into the file */
    PINDEX        nTotal;     /* Total number of particles in the file */
    uint_fast32_t iIndex;     /* Index of next particle in memory */
    uint_fast32_t nBuffered;  /* Number of buffered particles */
    hid_t group_id;           /* HDF5 group: /dark, /sph, etc. */
    IOORDER ioOrder;          /* The iOrder for each particle */
    IOCLASS ioClass;          /* The class (mass/softening) of each particle */
    IOFIELD *fldFields;       /* Data fields (position,velocity,etc.) */
    int     nFields;
    } IOBASE;

typedef struct {
    struct fioInfo fio;
    hid_t fileID;
    hid_t parametersID;
    hid_t stringType;
    int mFlags;
    FIO_SPECIES eCurrent;
    IOBASE base[FIO_SPECIES_LAST];
    } fioHDF5;

static hid_t fio2hdf(FIO_TYPE dataType,fioHDF5 *hio) {
    switch(dataType) {
    case FIO_TYPE_FLOAT:  return H5T_NATIVE_FLOAT;
    case FIO_TYPE_DOUBLE: return H5T_NATIVE_DOUBLE;
    case FIO_TYPE_UINT8:  return H5T_NATIVE_UINT8;
    case FIO_TYPE_UINT32: return H5T_NATIVE_UINT32;
    case FIO_TYPE_UINT64: return H5T_NATIVE_UINT64;
    case FIO_TYPE_INT:    return H5T_NATIVE_INT;
    case FIO_TYPE_STRING: return hio->stringType;
    default:
	fprintf(stderr,"Invalid type specified: %d\n", dataType);
	abort();
	break;
	}
    return H5I_INVALID_HID;
    }

static void field_reset(IOFIELD *ioClass) {
    ioClass->pBuffer = NULL;
    ioClass->memType = H5I_INVALID_HID;
    ioClass->setId   = H5I_INVALID_HID;
    ioClass->nValues = 0;
    ioClass->iSize   = 0;
    }

static void alloc_fields(IOBASE *base,int nFields) {
    int i;

    base->nFields = nFields;
    base->fldFields = malloc(sizeof(IOFIELD)*base->nFields);
    for(i=0; i<base->nFields; i++) {
	field_reset(&base->fldFields[i]);
	}
    }

/*
** field "class"
*/

#define DECLARE_FIELD_TYPE(data_t,hdf5_t)				\
    static inline int field_get_ ## data_t(data_t *pData, IOFIELD *field, int iIndex) { \
	if (pData) {							\
	    int i;							\
	    if (field->setId != H5I_INVALID_HID) {			\
		data_t *pBuffer = (data_t *)field->pBuffer;		\
		assert(field->memType == hdf5_t);			\
		iIndex *= field->nValues;				\
		for(i=0; i<field->nValues; ++i) pData[i] = pBuffer[iIndex+i]; \
		return 1;						\
		}							\
	    }								\
	return 0;							\
	}\
    static inline void field_add_ ## data_t(const data_t *pData, IOFIELD *field, int iIndex) { \
	if (pData && field->setId != H5I_INVALID_HID) {			\
	    data_t *pBuffer = (data_t *)field->pBuffer;			\
	    int i;							\
	    iIndex *= field->nValues;					\
	    for(i=0; i<field->nValues; ++i) pBuffer[iIndex+i] = pData[i]; \
	    }								\
	}

DECLARE_FIELD_TYPE(float,H5T_NATIVE_FLOAT)
DECLARE_FIELD_TYPE(double,H5T_NATIVE_DOUBLE)
DECLARE_FIELD_TYPE(uint64_t,H5T_NATIVE_UINT64)
DECLARE_FIELD_TYPE(uint32_t,H5T_NATIVE_UINT32)
DECLARE_FIELD_TYPE(uint8_t,H5T_NATIVE_UINT8)

static void field_read(IOFIELD *field, PINDEX iOffset, uint_fast32_t nBuffered) {
    if (field->setId != H5I_INVALID_HID) {
	readSet( field->setId,field->pBuffer,
		 field->memType, iOffset, nBuffered, field->nValues );
	}
    }

static void field_write(IOFIELD *field, PINDEX iOffset, uint_fast32_t nBuffered) {
    if (field->setId != H5I_INVALID_HID && nBuffered>0) {
	writeSet( field->setId, field->pBuffer,
		  field->memType, iOffset, nBuffered, field->nValues );
	}
    }

static hsize_t field_size(IOFIELD *field) {
    return getSetSize(field->setId);
    }

static void field_open(IOFIELD *field, hid_t groupID, const char *name, hid_t memType, int nValues ) {
    size_t iSize;

    field->setId = H5Dopen(groupID,name);
    if (field->setId != H5I_INVALID_HID) {
	iSize = H5Tget_size(memType);
	field->nChunk = CHUNK_SIZE;
	field->pBuffer = malloc(iSize*nValues*field->nChunk);
	assert(field->pBuffer);
	field->memType = memType;
	field->iSize   = iSize;
	field->nValues = nValues;
	}
    else {
	field->pBuffer = NULL;
	field->memType = H5I_INVALID_HID;
	field->iSize   = 0;
	field->nValues = 0;
	}
    }

static void field_close(IOFIELD *field) {
    if (field->setId != H5I_INVALID_HID) {
	free(field->pBuffer);
	H5Dclose(field->setId);
	field->pBuffer = NULL;
	field->memType = H5I_INVALID_HID;
	field->iSize   = 0;
	field->nValues = 0;
	}
    }

static inline int field_isopen(IOFIELD *field) {
    return field->pBuffer != NULL;
    }

static void field_create(IOFIELD *field, hid_t groupID, const char *name, hid_t memType, hid_t diskType, int nValues ) {
    size_t iSize;

    field->nChunk = CHUNK_SIZE;
    field->setId = newSet(groupID,name,field->nChunk,0,nValues,diskType);
    assert(field->setId!=H5I_INVALID_HID);
    iSize = H5Tget_size(memType);
    field->pBuffer = malloc(iSize*nValues*field->nChunk);
    assert(field->pBuffer);
    field->memType = memType;
    field->iSize   = iSize;
    field->nValues = nValues;
    }

/*
** Particle class table
*/

/* Create the type for the class table */
static hid_t makeClassType(hid_t floatType, int bStart) {
    hid_t tid;
    hid_t dataType = sizeof(PINDEX)==4 ? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;

    tid = H5Tcreate (H5T_COMPOUND, sizeof(classEntry));
    assert(tid!=H5I_INVALID_HID);
    H5Tinsert(tid,"class",HOFFSET(classEntry,iClass), H5T_NATIVE_UINT8);
    H5Tinsert(tid,"mass", HOFFSET(classEntry,fMass), floatType);
    H5Tinsert(tid,"soft", HOFFSET(classEntry,fSoft), floatType);
    if ( bStart )
	H5Tinsert(tid,"start",HOFFSET(classEntry,iOrderStart), dataType);
    return tid;
    }


static void class_flush(IOCLASS *ioClass,hid_t group_id) {
    hid_t tid, set;

    if ( ioClass->nClasses > 0 ) {
	tid = makeClassType( H5T_NATIVE_DOUBLE, !field_isopen(&ioClass->fldClasses));
	set = newSet(group_id, FIELD_CLASSES,
		     256, ioClass->nClasses, 1, tid );

	writeSet( set, ioClass->Class, tid,
		  0, ioClass->nClasses, 1 );

	H5Dclose(set);
	H5Tclose(tid);
	}
    }

static void class_write(IOCLASS *ioClass, PINDEX iOffset, uint_fast32_t nBuffered) {
    field_write(&ioClass->fldMass,iOffset,nBuffered);
    field_write(&ioClass->fldSoft,iOffset,nBuffered);
    field_write(&ioClass->fldClasses,iOffset,nBuffered);
    }

static void class_read(IOCLASS *ioClass, PINDEX iOffset, uint_fast32_t nBuffered) {
    field_read(&ioClass->fldMass,iOffset,nBuffered);
    field_read(&ioClass->fldSoft,iOffset,nBuffered);
    field_read(&ioClass->fldClasses,iOffset,nBuffered);
    }

static void class_add( IOBASE *base, PINDEX iOrder, float fMass, float fSoft ) {
    IOCLASS *ioClass = &base->ioClass;
    uint_fast32_t i;
    uint8_t iClass, jClass;

    if (field_isopen(&ioClass->fldMass)) {
	field_add_float(&fMass,&ioClass->fldMass,base->iIndex);
	fMass = 0.0;
	}
    if (field_isopen(&ioClass->fldSoft)) {
	field_add_float(&fSoft,&ioClass->fldSoft,base->iIndex);
	fSoft = 0.0;
	}

    /* See if we already have this class: Mass/Softening pair */
    for ( i=0; i<ioClass->nClasses; i++ ) {
	if ( ioClass->Class[i].fMass == fMass && ioClass->Class[i].fSoft == fSoft )
	    break;
	}
    iClass = i;

    /* Case 1: This is a new class */
    if ( iClass == ioClass->nClasses ) {
	assert( ioClass->nClasses < 256 ); /*TODO: handle this case */
	ioClass->Class[iClass].iClass = iClass;
	ioClass->Class[iClass].iOrderStart = iOrder;
	ioClass->Class[iClass].fMass = fMass;
	ioClass->Class[iClass].fSoft = fSoft;
	ioClass->Class[iClass].nCount= 0;
	ioClass->nClasses++;
	if ( field_isopen(&ioClass->fldClasses) )
	    field_add_uint8_t(&iClass,&ioClass->fldClasses,base->iIndex);
	}

    /* Case 2: This was the last class, and we might be compressing */
    else if ( iClass == ioClass->nClasses - 1 && !field_isopen(&ioClass->fldClasses) ) {
	}

    /* Case 3: A prior match meaning we cannot compress */
    else {
	if ( !field_isopen(&ioClass->fldClasses) ) {
	    uint64_t iIndex, iOffset;
	    int n;
	    field_create(&ioClass->fldClasses,base->group_id,
			 FIELD_CLASS, H5T_NATIVE_UINT8, H5T_NATIVE_UINT8, 1 );
	    ioClass->Class[ioClass->nClasses].iOrderStart = iOrder;
	    n = 0;
	    iOffset = 0;
	    for( jClass=0; jClass<ioClass->nClasses; jClass++ ) {
		for( iIndex = ioClass->Class[jClass].iOrderStart; iIndex<ioClass->Class[jClass+1].iOrderStart; iIndex++ ) {
		    field_add_uint8_t(&jClass,&ioClass->fldClasses,n);
		    if (++n == ioClass->fldClasses.nChunk) {
			field_write(&ioClass->fldClasses, iOffset, n );
			iOffset += n;
			n = 0;
			}
		    }
		}
	    }
	field_add_uint8_t(&iClass,&ioClass->fldClasses,base->iIndex);
	}
    ioClass->Class[i].nCount++;
    }

static void class_reset(IOCLASS *ioClass) {
    field_reset(&ioClass->fldClasses);
    }

static void class_open(IOCLASS *ioClass, hid_t groupID) {
    hid_t tid, set;

    set = H5Dopen( groupID, FIELD_CLASSES );
    if ( set != H5I_INVALID_HID ) {
	/* Each particle can have their own class ID, bit is isn't required. */
	field_open(&ioClass->fldClasses,groupID, FIELD_CLASS, H5T_NATIVE_UINT8, 1 );
	tid = makeClassType( H5T_NATIVE_DOUBLE, !field_isopen(&ioClass->fldClasses) );
	ioClass->nClasses = getSetSize(set);
	readSet( set, ioClass->Class, tid, 0, ioClass->nClasses, 1 );
	H5Dclose(set);
	H5Tclose(tid);
	}
    else field_reset(&ioClass->fldClasses);
    field_open(&ioClass->fldMass,groupID, FIELD_MASS, H5T_NATIVE_FLOAT, 1 );
    field_open(&ioClass->fldSoft,groupID, FIELD_SOFTENING, H5T_NATIVE_FLOAT, 1 );
    }

static void class_create(IOCLASS *ioClass, hid_t groupID,int bMass, int bSoft) {
    field_reset(&ioClass->fldClasses);
    field_reset(&ioClass->fldMass);
    field_reset(&ioClass->fldSoft);
    ioClass->nClasses = 0;
    if (bMass) {
	field_create(&ioClass->fldMass,groupID, FIELD_MASS, H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT, 1 );
	}
    if (bSoft) {
	field_create(&ioClass->fldSoft,groupID, FIELD_SOFTENING, H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT, 1 );
	}
    }

static void class_close(IOCLASS *ioClass) {
    field_close(&ioClass->fldClasses);
    field_close(&ioClass->fldMass);
    field_close(&ioClass->fldSoft);
    }


static void class_get(float *pfMass,float *pfSoft,IOCLASS *ioClass,PINDEX iOrder,uint_fast32_t iIndex) {
    uint8_t iClass;

    /* The particles were sorted by class to save space */
    if (field_isopen(&ioClass->fldClasses)) {
	field_get_uint8_t(&iClass,&ioClass->fldClasses,iIndex);
	assert( iClass < ioClass->nClasses );
	}
    else {
	assert(ioClass->nClasses>=1);
	for ( iClass=0; iClass<ioClass->nClasses; iClass++ )
	    if ( ioClass->Class[iClass].iOrderStart > iOrder )
		break;
	assert( iClass>0 );
	--iClass;
	}
    *pfMass = ioClass->Class[iClass].fMass;
    *pfSoft = ioClass->Class[iClass].fSoft;

    if (field_isopen(&ioClass->fldMass)) {
	field_get_float(pfMass,&ioClass->fldMass,iIndex);
	}
    if (field_isopen(&ioClass->fldSoft)) {
	field_get_float(pfSoft,&ioClass->fldSoft,iIndex);
	}
    }

/*
** IOORDER "class"
*/

void ioorder_open(IOORDER *order,hid_t group_id) {
    hid_t dataType = sizeof(PINDEX)==4 ? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;
    order->groupID = group_id;
    field_open(&order->fldOrder,group_id,FIELD_ORDER,dataType,1);
    order->iStart = 0;
    if (!field_isopen(&order->fldOrder)) {
	readAttribute(group_id,ATTR_IORDER,dataType,&order->iStart);
	}
    order->iNext = order->iStart;
    }

void ioorder_create(IOORDER *order, hid_t group_id, PINDEX iOrder) {
    field_reset(&order->fldOrder);
    order->groupID = group_id;
    order->iStart = order->iNext = iOrder;
    }

void ioorder_close(IOORDER *order) {
    field_close(&order->fldOrder);
    }

static void ioorder_read(IOORDER *order, PINDEX iOffset, uint_fast32_t nBuffered) {
    field_read(&order->fldOrder,iOffset,nBuffered);
    }

static PINDEX ioorder_get(IOORDER *order, PINDEX iOffset, uint_fast32_t iIndex) {
    if (field_isopen(&order->fldOrder)) {
	PINDEX iOrder;
	field_get_uint64_t(&iOrder,&order->fldOrder,iIndex);
	return iOrder;
	}
    else {
	return order->iStart + iOffset + iIndex;
	}
    }

static void ioorder_add(IOORDER *order, PINDEX iOrder) {
    /* Make sure that we are still in order */
    if (order->iNext == iOrder) ++order->iNext;
    //else abort();
    }

static void ioorder_write(IOORDER *order,PINDEX iOffset, uint_fast32_t nBuffered) {
    }

static void ioorder_flush(IOORDER *order) {
    hid_t dataType = sizeof(PINDEX)==4 ? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;
    writeAttribute( order->groupID, ATTR_IORDER,
		    dataType, &order->iStart );
    }

/*
** hdf5 interface
*/

static int hdf5GetAttr(
    FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    fioHDF5 *hio = (fioHDF5 *)fio;
    H5E_auto_t save_func;
    void *     save_data;
    int rc;

    assert(fio->eFormat == FIO_FORMAT_HDF5);
    assert(fio->eMode == FIO_MODE_READING);

    H5Eget_auto(&save_func,&save_data);
    H5Eset_auto(0,0);
    rc = readAttribute(hio->parametersID,attr,fio2hdf(dataType,hio),data);
    H5Eset_auto(save_func,save_data);
    return rc;
    }

static int hdf5SetAttr(
    FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    fioHDF5 *hio = (fioHDF5 *)fio;

    assert(fio->eFormat == FIO_FORMAT_HDF5);
    assert(fio->eMode == FIO_MODE_WRITING);

    return writeAttribute(hio->parametersID,attr,fio2hdf(dataType,hio),data);
    }

static FIO_SPECIES hdf5Species (struct fioInfo *fio) {
    fioHDF5 *hio = (fioHDF5 *)fio;
    return hio->eCurrent;
    }

/* Open the i'th file */
static int hdf5OpenOne(fioHDF5 *hio, int iFile) {
    H5E_auto_t save_func;
    void *     save_data;
    int i;

    assert(iFile<hio->fio.fileList.nFiles);
    hio->fio.fileList.iFile = iFile;

    /* Open the HDF5 file. */
    hio->fileID = H5Fopen(hio->fio.fileList.fileInfo[iFile].pszFilename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( hio->fileID < 0 ) {
	perror(hio->fio.fileList.fileInfo[iFile].pszFilename);
	abort();
	return 0;
	}

    /* Global parameters (dTime,etc.) are stored here */
    hio->parametersID = H5Gopen( hio->fileID, GROUP_PARAMETERS );
    if ( hio->parametersID == H5I_INVALID_HID ) {
	perror(hio->fio.fileList.fileInfo[iFile].pszFilename);
	abort();
	return 0;
	}

    /* Now open all of the available groups.  It's okay if some aren't there. */
    H5Eget_auto(&save_func,&save_data);
    H5Eset_auto(0,0);
    for( i=1; i<FIO_SPECIES_LAST; i++) {
	IOBASE *base = hio->base+i;

	switch(i) {
	case FIO_SPECIES_DARK:
	    base->group_id = H5Gopen(hio->fileID,fioSpeciesName(i));
	    if (base->group_id!=H5I_INVALID_HID) {
		base->iOffset = 0;
		base->iIndex = base->nBuffered = 0;

		alloc_fields(base,DARK_N);
		field_open(&base->fldFields[DARK_POSITION],base->group_id,
			   FIELD_POSITION, H5T_NATIVE_DOUBLE,3 );
		field_open(&base->fldFields[DARK_VELOCITY],base->group_id,
			   FIELD_VELOCITY, H5T_NATIVE_DOUBLE,3 );
		field_open(&base->fldFields[DARK_POTENTIAL],base->group_id,
			   FIELD_POTENTIAL, H5T_NATIVE_FLOAT,1 );
		if (base->fldFields[DARK_POTENTIAL].setId == H5I_INVALID_HID)
		    hio->fio.mFlags &= ~FIO_FLAG_POTENTIAL;
		field_open(&base->fldFields[DARK_DENSITY],base->group_id,
			   FIELD_DENSITY, H5T_NATIVE_FLOAT,1 );
		if (base->fldFields[DARK_DENSITY].setId == H5I_INVALID_HID)
		    hio->fio.mFlags &= ~FIO_FLAG_DENSITY;
		base->nTotal = hio->fio.fileList.fileInfo[iFile].nSpecies[i] = field_size(&base->fldFields[DARK_POSITION]);
		hio->fio.nSpecies[i] += hio->fio.fileList.fileInfo[iFile].nSpecies[i];
		class_open(&base->ioClass,base->group_id);
		/* iOrder can have a starting value if they are sequential, or a list */
		ioorder_open(&base->ioOrder,base->group_id);
		}
	    else base->nTotal = 0;
	    break;
	default:
	    hio->fio.fileList.fileInfo[iFile].nSpecies[i] = 0;
	    base->group_id = H5I_INVALID_HID;
	    base->nTotal = 0;
	    break;
	    }
	}
    H5Eset_auto(save_func,save_data);
    return 1;
    }

/* Close the current open file; does not destroy the context */
static void hdf5CloseOne(fioHDF5 *hio) {
    int i, j;

    for( i=1; i<FIO_SPECIES_LAST; i++) {
	IOBASE *base = &hio->base[i];
	if (base->group_id!=H5I_INVALID_HID) {
	    if (hio->fio.eMode==FIO_MODE_WRITING) {
		ioorder_flush(&base->ioOrder);
		class_flush(&base->ioClass,base->group_id);
		class_write(&base->ioClass,base->iOffset,base->iIndex);
		for(j=0; j<base->nFields; j++) {
		    field_write(&base->fldFields[j],
				base->iOffset, base->iIndex);
		    }
		}
	    ioorder_close(&base->ioOrder);
	    for(j=0; j<base->nFields; j++)
		field_close(&base->fldFields[j]);
	    class_close(&base->ioClass);
	    H5Gclose(base->group_id);
	    }
	}
    H5Gclose(hio->parametersID);
    H5Fflush(hio->fileID,H5F_SCOPE_GLOBAL);
    assert(H5Fget_obj_count(hio->fileID,H5F_OBJ_ALL)==1);
    H5Fclose(hio->fileID);
    }

static int hdf5Seek(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioHDF5 *hio = (fioHDF5 *)fio;
    IOBASE *base;
    int i, iFile;

    assert(fio->eFormat == FIO_FORMAT_HDF5);

    /* First identify the correct file */
    for(iFile=0; iFile<fio->fileList.nFiles; ++iFile) {
	if (iPart>=fio->fileList.fileInfo[iFile].nSpecies[eSpecies])
	    iPart -= fio->fileList.fileInfo[iFile].nSpecies[eSpecies];
	else
	    break;
	}

    /* Open the proper file if necessary */
    if (iFile!=fio->fileList.iFile) {
	hdf5CloseOne(hio);
	hdf5OpenOne(hio,iFile);
	}

    if (eSpecies==FIO_SPECIES_ALL) {
	for( i=1; i<FIO_SPECIES_LAST; i++) {
	    base = &hio->base[i];
	    if (iPart<base->nTotal) {
		eSpecies = i;
		break;
		}
	    iPart -= base->nTotal;
	    }
	assert(eSpecies!=FIO_SPECIES_ALL);
	}
    hio->eCurrent = eSpecies;
    base = &hio->base[eSpecies];
    base->iOffset = iPart;
    base->iIndex = base->nBuffered = 0;
    return 0;
    }

static int base_read(fioHDF5 *hio,IOBASE *base) {
    hsize_t N;
    int i;

    base->iOffset += base->nBuffered;
    base->iIndex = base->nBuffered = 0;
    if ( base->iOffset >= base->nTotal ) {
	hdf5CloseOne(hio);
	hdf5OpenOne(hio,++hio->fio.fileList.iFile);
	}

    N = base->nTotal - base->iOffset;
    base->nBuffered = N > CHUNK_SIZE ? CHUNK_SIZE : N;

    for(i=0; i<base->nFields; i++) {
	field_read(&base->fldFields[i],
		   base->iOffset, base->nBuffered);
	}
    ioorder_read(&base->ioOrder,base->iOffset,base->nBuffered);
    class_read(&base->ioClass,base->iOffset,base->nBuffered);
    return 1;
    }

static int base_write(IOBASE *base) {
    int i;

    if ( base->iIndex) {
	for(i=0; i<base->nFields; i++) {
	    field_write(&base->fldFields[i],
			base->iOffset, base->iIndex);
	    }
	ioorder_write(&base->ioOrder,base->iOffset,base->iIndex);
	class_write(&base->ioClass,base->iOffset,base->iIndex);
	}
    base->iOffset += base->iIndex;
    base->iIndex = base->nBuffered = 0;
    return 1;
    }


static void base_create(fioHDF5 *hio,IOBASE *base,int iSpecies,int nFields,uint64_t iOrder) {
    hid_t posType, velType;

    posType = (hio->mFlags&FIO_FLAG_CHECKPOINT) || (hio->mFlags&FIO_FLAG_DOUBLE_POS)
	? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
    velType = (hio->mFlags&FIO_FLAG_CHECKPOINT) || (hio->mFlags&FIO_FLAG_DOUBLE_VEL)
	? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;

    base->group_id = H5Gcreate(hio->fileID,fioSpeciesName(iSpecies),0);
    base->iOffset = 0;
    base->iIndex = base->nBuffered = 0;
    ioorder_create(&base->ioOrder,base->group_id,iOrder);
    class_create(&base->ioClass,base->group_id,
		 !(hio->mFlags&FIO_FLAG_COMPRESS_MASS),
		 !(hio->mFlags&FIO_FLAG_COMPRESS_SOFT));
    alloc_fields(base,nFields);
    field_create(&base->fldFields[DARK_POSITION],base->group_id,
		 FIELD_POSITION, H5T_NATIVE_DOUBLE, posType, 3 );
    field_create(&base->fldFields[DARK_VELOCITY],base->group_id,
		 FIELD_VELOCITY, H5T_NATIVE_DOUBLE, velType,3 );
    if ( hio->mFlags&FIO_FLAG_POTENTIAL) {
	field_create(&base->fldFields[DARK_POTENTIAL],base->group_id,
		     FIELD_POTENTIAL, H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT,1 );
	}
    if ( hio->mFlags&FIO_FLAG_DENSITY) {
	field_create(&base->fldFields[DARK_DENSITY],base->group_id,
		     FIELD_DENSITY, H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT,1 );
	}
    }

static int hdf5ReadDark(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen) {
    fioHDF5 *hio = (fioHDF5 *)(fio);
    IOBASE *base = &hio->base[hio->eCurrent];
    int i;

    assert(fio->eFormat == FIO_FORMAT_HDF5);
    assert(hio->eCurrent == FIO_SPECIES_DARK);

    /* If we have exhausted our buffered data, read more */
    if (base->iIndex == base->nBuffered) {
	base_read(hio,base);
	}

    /* Position and Velocity are always present */
    field_get_double(pdPos,&base->fldFields[DARK_POSITION],base->iIndex);
    field_get_double(pdVel,&base->fldFields[DARK_VELOCITY],base->iIndex);

    /* Potential is optional */
    if ( !field_get_float(pfPot,&base->fldFields[DARK_POTENTIAL],base->iIndex) )
	*pfPot = 0.0f;

    /* Density is optional */
    if ( !field_get_float(pfDen,&base->fldFields[DARK_DENSITY],base->iIndex) )
	*pfDen = 0.0f;

    /* iOrder is either sequential, or is listed for each particle */
    *piOrder = ioorder_get(&base->ioOrder,base->iOffset,base->iIndex);

    /* If each particles has a unique class, use that */
    class_get(pfMass,pfSoft,&base->ioClass,*piOrder,base->iIndex);

    /*
    ** Next particle.  If we are at the end of this species,
    ** advance to the start of the next species.
    */
    base->iIndex++;
    if (base->iIndex==base->nTotal) {
	for( i=hio->eCurrent+1; i<FIO_SPECIES_LAST; i++) {
	    base = &hio->base[i];
	    if ( base->nTotal ) {
		hio->eCurrent = i;
		base->iOffset = base->iIndex = base->nBuffered = 0;
		break;
		}
	    }
	}
    return 1;
    }

static int hdf5ReadSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft, float *pfPot,float *pfDen,
    float *pfTemp, float *pfMetals) {
    return 0;
    }

static int hdf5ReadStar(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfMetals, float *pfTform) {
    return 0;
    }

static int  hdf5WriteDark(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen) {
    fioHDF5 *hio = (fioHDF5 *)(fio);
    IOBASE *base = &hio->base[FIO_SPECIES_DARK];

    assert(fio->eFormat == FIO_FORMAT_HDF5);
    assert(fio->eMode == FIO_MODE_WRITING);

    /* First time for this particle type? */
    if (base->group_id == H5I_INVALID_HID) {
	base_create(hio,base,FIO_SPECIES_DARK,DARK_N,iOrder);
	}

    ioorder_add(&base->ioOrder,iOrder);
    class_add(base,iOrder,fMass,fSoft);
    field_add_double(pdPos,&base->fldFields[DARK_POSITION],base->iIndex);
    field_add_double(pdVel,&base->fldFields[DARK_VELOCITY],base->iIndex);
    field_add_float(&fPot,&base->fldFields[DARK_POTENTIAL],base->iIndex);
    field_add_float(&fDen,&base->fldFields[DARK_DENSITY],base->iIndex);

    /* If we have exhausted our buffered data, read more */
    if (++base->iIndex == CHUNK_SIZE) {
	base_write(base);
	}
    return 1;
    }

static int hdf5WriteSph(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float fTemp,float fMetals) {
    fioHDF5 *hio = (fioHDF5 *)(fio);
    IOBASE *base = &hio->base[FIO_SPECIES_SPH];
    assert(fio->eFormat == FIO_FORMAT_HDF5);
    assert(fio->eMode == FIO_MODE_WRITING);

    /* First time for this particle type? */
    if (base->group_id == H5I_INVALID_HID) {
	base_create(hio,base,FIO_SPECIES_SPH,SPH_N,iOrder);
	field_create(&base->fldFields[SPH_TEMPERATURE],base->group_id,
		     FIELD_TEMPERATURE, H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT, 1 );
	field_create(&base->fldFields[SPH_METALS],base->group_id,
		     FIELD_METALS, H5T_NATIVE_FLOAT, H5T_NATIVE_FLOAT, 1 );
	}

    ioorder_add(&base->ioOrder,iOrder);
    class_add(base,iOrder,fMass,fSoft);
    field_add_double(pdPos,&base->fldFields[DARK_POSITION],base->iIndex);
    field_add_double(pdVel,&base->fldFields[DARK_VELOCITY],base->iIndex);
    field_add_float(&fPot,&base->fldFields[DARK_POTENTIAL],base->iIndex);
    field_add_float(&fDen,&base->fldFields[DARK_DENSITY],base->iIndex);
    field_add_float(&fTemp,&base->fldFields[SPH_TEMPERATURE],base->iIndex);
    field_add_float(&fMetals,&base->fldFields[SPH_METALS],base->iIndex);

    /* If we have exhausted our buffered data, read more */
    if (++base->iIndex == CHUNK_SIZE) {
	base_write(base);
	}
    return 1;
    }

static int hdf5WriteStar(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,float fMetals,float fTform) {
    fprintf(stderr,"Writing star particles is not supported\n");
    abort();
    }

static void hdf5Close(FIO fio) {
    fioHDF5 *hio = (fioHDF5 *)(fio);
    hdf5CloseOne(hio);
    H5Tclose(hio->stringType);
    fileScanFree(&fio->fileList);
    free(hio);
    }

static FIO hdf5Open(fioFileList *fileList) {
    fioHDF5 *hio;
    int i;

    hio = malloc(sizeof(fioHDF5));
    assert(hio!=NULL);
    fioInitialize(&hio->fio,FIO_FORMAT_HDF5,FIO_MODE_READING,0);

    hio->fio.fileList = *fileList;
    hio->fio.fileList.iFile = 0;

    hio->stringType = H5Tcopy(H5T_C_S1);
    H5Tset_size(hio->stringType, 256);

    /* Scan the files - all but the first one */
    hio->fio.mFlags |= FIO_FLAG_DENSITY | FIO_FLAG_POTENTIAL;
    for( i=hio->fio.fileList.nFiles-1; i>0; --i) {
	hdf5OpenOne(hio,i);
	hdf5CloseOne(hio);
	}
    /* Open the first HDF5 file. */
    hdf5OpenOne(hio,0);

    fioTabulateSpecies(&hio->fio);

    /* Set the current species (the first one in the first file) */
    hio->eCurrent = 0;
    for( i=1; i<FIO_SPECIES_LAST; i++) {
	if (hio->eCurrent==0 && hio->fio.fileList.fileInfo[0].nSpecies[i]) hio->eCurrent=i;
	}

    hio->fio.fcnClose    = hdf5Close;
    hio->fio.fcnSeek     = hdf5Seek;
    hio->fio.fcnReadDark = hdf5ReadDark;
    hio->fio.fcnReadSph  = hdf5ReadSph;
    hio->fio.fcnReadStar = hdf5ReadStar;
    hio->fio.fcnGetAttr  = hdf5GetAttr;
    hio->fio.fcnSetAttr  = hdf5SetAttr;
    hio->fio.fcnSpecies  = hdf5Species;

    return &hio->fio;
    }

FIO fioHDF5Create(const char *fileName, int mFlags) {
    fioHDF5 *hio;
    int i;

    hio = malloc(sizeof(fioHDF5));
    assert(hio!=NULL);
    fioInitialize(&hio->fio,FIO_FORMAT_HDF5,FIO_MODE_WRITING,mFlags);

    hio->mFlags = mFlags;

    /* Open the HDF5 file. */
    hio->fileID = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if ( hio->fileID < 0 ) {
	free(hio);
	return NULL;
	}
    assert(H5Fget_obj_count(hio->fileID,H5F_OBJ_ALL)==1);

    hio->stringType = H5Tcopy(H5T_C_S1);
    H5Tset_size(hio->stringType, 256);

    /* Global parameters (dTime,etc.) are stored here */
    hio->parametersID = H5Gcreate( hio->fileID, GROUP_PARAMETERS, 0 );
    if ( hio->parametersID == H5I_INVALID_HID ) {
	abort();
	}

    for( i=1; i<FIO_SPECIES_LAST; i++) {
	IOBASE *base = hio->base+i;
	base->group_id = H5I_INVALID_HID;
	base->nTotal = 0;
	}

    hio->fio.fcnClose    = hdf5Close;
    hio->fio.fcnSeek     = hdf5Seek;
    hio->fio.fcnWriteDark= hdf5WriteDark;
    hio->fio.fcnWriteSph = hdf5WriteSph;
    hio->fio.fcnWriteStar= hdf5WriteStar;
    hio->fio.fcnGetAttr  = hdf5GetAttr;
    hio->fio.fcnSetAttr  = hdf5SetAttr;
    hio->fio.fcnSpecies  = hdf5Species;

    return &hio->fio;
    }

#endif

/******************************************************************************\
** GRAFIC FORMAT
\******************************************************************************/

typedef struct {
    int32_t n[3];
    float dx, o[3];
    float astart, omegam, omegav, H0;
    } GraficHdr4f;

typedef struct {
    int64_t n[3];
    float dx, o[3];
    float astart, omegam, omegav, H0;
    } GraficHdr8f;

typedef struct {
    int32_t n[3];
    double dx, o[3];
    double astart, omegam, omegav, H0;
    } GraficHdr4d;

typedef struct {
    int64_t n[3];
    double dx, o[3];
    double astart, omegam, omegav, H0;
    } GraficHdr8d;

typedef struct {
    FILE *fp;
    union {
	float *pFloat;
	double *pDouble;
	} data;
    int iPosition[3];
    int iIndex;
    int nPerSlab;
    int nSlabSize;
    off_t nHdrSize;
    int bDouble;
    GraficHdr8d hdr;
    } graficFile;

typedef struct {
    graficFile fp_velcx;
    graficFile fp_velcy;
    graficFile fp_velcz;
    graficFile fp_velbx;
    graficFile fp_velby;
    graficFile fp_velbz;
    } graficLevel;

typedef struct {
    struct fioInfo fio;
    uint64_t iOrder;
    double   dTime;
    double   vFactor;
    double   pFactor1;
    double   pFactor2;
    double   mValueCDM;
    double   mValueBar;
    double   sValue;
    int nLevels;
    graficLevel *level;
    } fioGrafic;

/*
** Read the header information from a GRAFIC velocity file.
** The data types are inferred from the size of the record.
*/
static void graficReadHdr(graficFile *gf) {
    uint32_t w1,w2;
    GraficHdr4f hdr4f;
    GraficHdr8f hdr8f;
    GraficHdr4d hdr4d;
    GraficHdr8d hdr8d;
    int rc;

    rc = fread(&w1,sizeof(w1),1,gf->fp);
    assert(rc==1);
    switch(w1) {
    case sizeof(GraficHdr4f):
	rc = fread(&hdr4f,w1,1,gf->fp);
	assert(rc==1);
	gf->hdr.n[0] = hdr4f.n[0];
	gf->hdr.n[1] = hdr4f.n[1];
	gf->hdr.n[2] = hdr4f.n[2];
	gf->hdr.dx = hdr4f.dx;
	gf->hdr.o[0] = hdr4f.o[0];
	gf->hdr.o[1] = hdr4f.o[1];
	gf->hdr.o[2] = hdr4f.o[2];
	gf->hdr.astart = hdr4f.astart;
	gf->hdr.omegam = hdr4f.omegam;
	gf->hdr.omegav = hdr4f.omegav;
	gf->hdr.H0 = hdr4f.H0;
	gf->bDouble = 0;
	break;
    case sizeof(GraficHdr8f):
	rc = fread(&hdr8f,w1,1,gf->fp);
	assert(rc==1);
	gf->hdr.n[0] = hdr8f.n[0];
	gf->hdr.n[1] = hdr8f.n[1];
	gf->hdr.n[2] = hdr8f.n[2];
	gf->hdr.dx = hdr8f.dx;
	gf->hdr.o[0] = hdr8f.o[0];
	gf->hdr.o[1] = hdr8f.o[1];
	gf->hdr.o[2] = hdr8f.o[2];
	gf->hdr.astart = hdr8f.astart;
	gf->hdr.omegam = hdr8f.omegam;
	gf->hdr.omegav = hdr8f.omegav;
	gf->hdr.H0 = hdr8f.H0;
	gf->bDouble = 0;
	break;
    case sizeof(GraficHdr4d):
	rc = fread(&hdr4d,w1,1,gf->fp);
	assert(rc==1);
	gf->hdr.n[0] = hdr4d.n[0];
	gf->hdr.n[1] = hdr4d.n[1];
	gf->hdr.n[2] = hdr4d.n[2];
	gf->hdr.dx = hdr4d.dx;
	gf->hdr.o[0] = hdr4d.o[0];
	gf->hdr.o[1] = hdr4d.o[1];
	gf->hdr.o[2] = hdr4d.o[2];
	gf->hdr.astart = hdr4d.astart;
	gf->hdr.omegam = hdr4d.omegam;
	gf->hdr.omegav = hdr4d.omegav;
	gf->hdr.H0 = hdr4d.H0;
	gf->bDouble = 1;
	break;
    case sizeof(GraficHdr8d):
	rc = fread(&hdr8d,w1,1,gf->fp);
	assert(rc==1);
	gf->hdr.n[0] = hdr8d.n[0];
	gf->hdr.n[1] = hdr8d.n[1];
	gf->hdr.n[2] = hdr8d.n[2];
	gf->hdr.dx = hdr8d.dx;
	gf->hdr.o[0] = hdr8d.o[0];
	gf->hdr.o[1] = hdr8d.o[1];
	gf->hdr.o[2] = hdr8d.o[2];
	gf->hdr.astart = hdr8d.astart;
	gf->hdr.omegam = hdr8d.omegam;
	gf->hdr.omegav = hdr8d.omegav;
	gf->hdr.H0 = hdr8d.H0;
	gf->bDouble = 1;
	break;
    default:
	fprintf(stderr,"Invalid GRAFIC header size\n");
	abort();
	}
    rc = fread(&w2,sizeof(w2),1,gf->fp);
    assert(rc==1);
    assert(w1==w2);

    gf->nHdrSize = ftell(gf->fp);

    assert(sizeof(off_t)>=8);
    assert(gf->nHdrSize==w1+2*sizeof(w1));

    gf->nPerSlab = (uint64_t)gf->hdr.n[0] * (uint64_t)gf->hdr.n[1];
    gf->iIndex = gf->nPerSlab;
    gf->iPosition[0] = -1;
    gf->iPosition[1] = gf->iPosition[2] = 0;

    // Find out if we have single or double precision
    rc = fread(&w1,sizeof(w1),1,gf->fp);
    assert(rc==1);

    if (w1 == sizeof(double)*gf->nPerSlab) {
	gf->bDouble = 1;
	gf->nSlabSize = sizeof(double)*gf->nPerSlab;
	gf->data.pDouble = malloc(gf->nSlabSize);
	assert(gf->data.pDouble!=NULL);
	}
    else if (w1 == sizeof(float)*gf->nPerSlab) {
	gf->bDouble = 0;
	gf->nSlabSize = sizeof(float)*gf->nPerSlab;
	gf->data.pFloat = malloc(gf->nSlabSize);
	assert(gf->data.pFloat!=NULL);
	}
    else {
	fprintf(stderr,"Invalid GRAFIC slab size\n");
	abort();
	}
    rc = fseek(gf->fp,-sizeof(w1),SEEK_CUR);
    assert(rc==0);
    }

/*
** Open a single GRAFIC file and read in the header information
*/
static int graficOpen(graficFile *gf,const char *fileName) {
    gf->fp = fopen(fileName,"rb");
    if ( gf->fp != NULL ) graficReadHdr(gf);
    return gf->fp != NULL;
    }

/*
** Return the next velocity from the file, reading more if necessary.
*/
static double graficRead(graficFile *gf) {
    int rc;
    uint32_t w;
    assert( gf->iIndex <= gf->nPerSlab);
    if ( ++gf->iPosition[0] == gf->hdr.n[0] ) {
	gf->iPosition[0] = 0;
	if ( ++gf->iPosition[1] == gf->hdr.n[1] ) {
	    gf->iPosition[1] = 0;
	    ++gf->iPosition[2];
	    }
	}
    if ( gf->iIndex == gf->nPerSlab) {
	gf->iIndex = 0;

	/* Read an entire slab */
	rc = fread(&w,sizeof(w),1,gf->fp);
	assert(rc==1 && w==gf->nSlabSize);

	if ( gf->bDouble )
	    rc = fread(gf->data.pDouble,sizeof(double),gf->nPerSlab,gf->fp);
	else
	    rc = fread(gf->data.pFloat,sizeof(float),gf->nPerSlab,gf->fp);
	assert(rc==gf->nPerSlab);

	rc = fread(&w,sizeof(w),1,gf->fp);
	assert(rc==1 && w==gf->nSlabSize);
	}
    if ( gf->bDouble ) return gf->data.pDouble[gf->iIndex++];
    else return gf->data.pFloat[gf->iIndex++];
    }

static void graficSeekFile(graficFile *gf,uint64_t iPart) {
    uint64_t
	sFloat,     /* Bytes in each float (could be double) */
	iSlab;      /* Index of the slab */
    uint64_t iByte; /* Byte offset into the file */
    off_t iOffset;
    uint32_t w;
    int rc;

    gf->iPosition[2] = iPart / gf->nPerSlab;
    gf->iPosition[1] = (iPart-gf->iPosition[2]*gf->nPerSlab) / gf->hdr.n[0];
    gf->iPosition[0] = iPart % gf->hdr.n[0] - 1;

    /* Calculate the slab, particle in slab and byte offset */
    sFloat = gf->bDouble ? sizeof(double) : sizeof(float);
    iSlab = iPart / gf->nPerSlab;
    assert( iSlab < gf->hdr.n[2] );
    iPart -= iSlab*gf->nPerSlab;

    iByte = gf->nHdrSize
	+ iSlab * (sFloat*gf->nPerSlab + 2*sizeof(uint32_t))
	+ iPart * sFloat + sizeof(uint32_t);
    iOffset = iByte;
    assert(iOffset==iByte);

    rc = safe_fseek(gf->fp,iByte); assert(rc==0);

    /* Read the remainder of this slab */
    if ( gf->bDouble )
	rc = fread(gf->data.pDouble+iPart,sizeof(double),gf->nPerSlab-iPart,gf->fp);
    else
	rc = fread(gf->data.pFloat+iPart,sizeof(float),gf->nPerSlab-iPart,gf->fp);
    assert(rc==gf->nPerSlab-iPart);

    gf->iIndex = iPart;

    /* Also verify that the FORTRAN record length is correct */
    rc = fread(&w,sizeof(w),1,gf->fp);
    assert(rc==1 && w==gf->nSlabSize);
    }

/*
** Compare two GRAFIC headers for equality
*/
static int graficCompare(graficFile *a,graficFile *b) {
    return a->hdr.n[0] == b->hdr.n[0]
	&& a->hdr.n[1] == b->hdr.n[1]
	&& a->hdr.n[2] == b->hdr.n[2]
	&& a->hdr.dx == b->hdr.dx
	&& a->hdr.o[0] == b->hdr.o[0]
	&& a->hdr.o[1] == b->hdr.o[1]
	&& a->hdr.o[2] == b->hdr.o[2]
	&& a->hdr.astart == b->hdr.astart
	&& a->hdr.omegam == b->hdr.omegam
	&& a->hdr.omegav == b->hdr.omegav
	&& a->hdr.H0 == b->hdr.H0;
    }

static int graficGetAttr(FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    fioGrafic *gio = (fioGrafic *)fio;
    assert(fio->eFormat == FIO_FORMAT_GRAFIC);
    if ( strcmp(attr,"dTime")==0 ) {
	switch(dataType) {
	case FIO_TYPE_FLOAT: *(float *)(data) = gio->dTime; break;
	case FIO_TYPE_DOUBLE:*(double *)(data) = gio->dTime; break;
	default: return 0;
	    }
	return 1;
	}

    return 0;
    }

static int graficSeek(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioGrafic *gio = (fioGrafic *)fio;

    assert(fio->eFormat == FIO_FORMAT_GRAFIC);

    /* Turn Species/Particle into absolute particle index */
    switch(eSpecies) {
    case FIO_SPECIES_DARK:
	iPart += gio->fio.nSpecies[FIO_SPECIES_SPH];
	break;
    case FIO_SPECIES_SPH:
    case FIO_SPECIES_ALL:
	break;
    default:
	fprintf(stderr,"Seeking to an unsupported particle type: %d\n",eSpecies);
	abort();
	}

    gio->iOrder = iPart;

    /* Now seek within the appropriate files */
    if (iPart>=gio->fio.nSpecies[FIO_SPECIES_SPH]) {
	iPart -= gio->fio.nSpecies[FIO_SPECIES_SPH];
	graficSeekFile(&gio->level[0].fp_velcx,iPart);
	graficSeekFile(&gio->level[0].fp_velcy,iPart);
	graficSeekFile(&gio->level[0].fp_velcz,iPart);
	}
    else {
	graficSeekFile(&gio->level[0].fp_velbx,iPart);
	graficSeekFile(&gio->level[0].fp_velby,iPart);
	graficSeekFile(&gio->level[0].fp_velbz,iPart);
	}
    return 0;
    }

/* Return the cell center in GRAFIC coordinates (Mpc at Z=0) */
static double graficCellCenter(graficLevel *lvl,int iDim) {
    return (0.5+lvl->fp_velcx.iPosition[iDim]) * lvl->fp_velcx.hdr.dx + lvl->fp_velcx.hdr.o[iDim];
    }

static double wrap(double v) {
    if (v<-0.5) v += 1.0;
    else if (v>=0.5) v -= 1.0;
    return v;
    }

/* Apply scaling factors to position and velocity */
static void graficSetPV(FIO fio,double *r,double *v,double x,double y,double z) {
    fioGrafic *gio = (fioGrafic *)fio;
    r[0] = wrap((graficCellCenter(&gio->level[0],0) + x * gio->pFactor1) * gio->pFactor2 - 0.5);
    r[1] = wrap((graficCellCenter(&gio->level[0],1) + y * gio->pFactor1) * gio->pFactor2 - 0.5);
    r[2] = wrap((graficCellCenter(&gio->level[0],2) + z * gio->pFactor1) * gio->pFactor2 - 0.5);
    v[0] = x * gio->vFactor;
    v[1] = y * gio->vFactor;
    v[2] = z * gio->vFactor;
    }

static int graficReadDark(FIO fio,
			  uint64_t *piOrder,double *pdPos,double *pdVel,
			  float *pfMass,float *pfSoft,float *pfPot,float *pfDen) {
    fioGrafic *gio = (fioGrafic *)fio;
    assert(fio->eFormat == FIO_FORMAT_GRAFIC);
    assert(gio->iOrder >= gio->fio.nSpecies[FIO_SPECIES_SPH]);
    *piOrder = gio->iOrder++;
    graficSetPV(fio,pdPos,pdVel,
		graficRead(&gio->level[0].fp_velcx),
		graficRead(&gio->level[0].fp_velcy),
		graficRead(&gio->level[0].fp_velcz) );
    *pfMass = gio->mValueCDM;
    *pfSoft = gio->sValue;
    if ( pfPot) *pfPot = 0.0f;
    if ( pfDen) *pfDen = 0.0f;
    return 1;
    }

static int graficReadSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfTemp, float *pfMetals) {
    fioGrafic *gio = (fioGrafic *)fio;
    assert(fio->eFormat == FIO_FORMAT_GRAFIC);
    assert(gio->iOrder < gio->fio.nSpecies[FIO_SPECIES_SPH]);
    *piOrder = gio->iOrder++;
    graficSetPV(fio,pdPos,pdVel,
		graficRead(&gio->level[0].fp_velbx),
		graficRead(&gio->level[0].fp_velby),
		graficRead(&gio->level[0].fp_velbz) );
    *pfMass = gio->mValueBar;
    *pfSoft = gio->sValue;
    if (pfPot) *pfPot = 0.0;
    if (pfDen) *pfDen = 0.0;
    if (pfTemp) *pfTemp = 0.0;
    if (pfSoft) *pfSoft = 0.0;
    if (pfMetals) *pfMetals = 0.0;
    return 1;
    }

static void graficClose(FIO fio) {
    fioGrafic *gio = (fioGrafic *)fio;
    if ( gio->level[0].fp_velcx.fp!=NULL ) fclose(gio->level[0].fp_velcx.fp);
    if ( gio->level[0].fp_velcy.fp!=NULL ) fclose(gio->level[0].fp_velcy.fp);
    if ( gio->level[0].fp_velcz.fp!=NULL ) fclose(gio->level[0].fp_velcz.fp);
    if ( gio->level[0].fp_velbx.fp!=NULL ) fclose(gio->level[0].fp_velbx.fp);
    if ( gio->level[0].fp_velby.fp!=NULL ) fclose(gio->level[0].fp_velby.fp);
    if ( gio->level[0].fp_velbz.fp!=NULL ) fclose(gio->level[0].fp_velbz.fp);
    free(gio);
    }

typedef struct {
    double omegam;
    double omegav;
    } ddplus_ctx;

#ifdef HAVE_GSL
static double ddplus(double a,void *ctx) {
    ddplus_ctx *dc = ctx;
    double eta;
    if ( a == 0.0 ) return 0.0;
    eta = sqrt(dc->omegam/a + dc->omegav*a*a + 1.0 - dc->omegam - dc->omegav);
    return 2.5/(eta*eta*eta);
    }
static double dplus(double a,double omegam,double omegav) {
    double eta;
    ddplus_ctx ddctx;
    gsl_function F;
    double result, error;
    gsl_integration_workspace *W;

    ddctx.omegam = omegam;
    ddctx.omegav = omegav;
    eta = sqrt(omegam/a + omegav*a*a + 1.0 - omegam - omegav);
    W = gsl_integration_workspace_alloc (1000);
    F.function = &ddplus;
    F.params = &ddctx;
    gsl_integration_qag(&F, 0.0, a, 0.0, 1e-8, 1000,
        GSL_INTEG_GAUSS61, W, &result, &error); 
    gsl_integration_workspace_free(W);
    return eta/a * result;
    }
#else
static double ddplus(void *ctx,double a) {
    ddplus_ctx *dc = ctx;
    double eta;
    if ( a == 0.0 ) return 0.0;
    eta = sqrt(dc->omegam/a + dc->omegav*a*a + 1.0 - dc->omegam - dc->omegav);
    return 2.5/(eta*eta*eta);
    }
#ifdef HAVE_ROMBERG
extern double dRombergO(void *CTX, double (*func)(void *,double), double a,
                 double b, double eps);
#else
#define MAXLEV 13

/*
 ** Romberg integrator for an open interval.
 */

static double dRombergO(void *CTX,double (*func)(void *, double),double a,double b,double eps) {
    double tllnew;
    double tll;
    double tlk[MAXLEV+1];
    int n = 1;
    int nsamples = 1;

    tlk[0] = tllnew = (b-a)*(*func)(CTX, 0.5*(b+a));
    if (a == b) return tllnew;

    eps*=0.5;

    do {
        /*
 *          * midpoint rule.
 *                   */
        double deltax;
        double tlktmp;
        int i;

        nsamples *= 3;
        deltax = (b-a)/nsamples;
        tlktmp = tlk[0];
        tlk[0] = tlk[0]/3.0;

        for (i=0;i<nsamples/3;i++) {
            tlk[0] += deltax*(*func)(CTX,a + (3*i + 0.5)*deltax);
            tlk[0] += deltax*(*func)(CTX,a + (3*i + 2.5)*deltax);
            }

        /*
 *          * Romberg extrapolation.
 *                   */

        for (i=0;i<n;i++) {
            double tlknew = (pow(9.0, i+1.)*tlk[i] - tlktmp)
                            /(pow(9.0, i+1.) - 1.0);

            tlktmp = tlk[i+1];
            tlk[i+1] = tlknew;
            }
        tll = tllnew;
        tllnew = tlk[n];
        n++;
        } while ((fabs((tllnew-tll)/(tllnew+tll)) > eps) && (n < MAXLEV));

    assert((fabs((tllnew-tll)/(tllnew+tll)) < eps));

    return tllnew;
    }
#endif

static double dplus(double a,double omegam,double omegav) {
    double eta;
    ddplus_ctx ddctx;

    ddctx.omegam = omegam;
    ddctx.omegav = omegav;
    eta = sqrt(omegam/a + omegav*a*a + 1.0 - omegam - omegav);
    return eta/a * dRombergO(&ddctx,ddplus,0,a,1e-8);
}
#endif

static double fomega(double a,double omegam,double omegav) {
    double eta, omegak;
    if ( omegam == 1.0 && omegav == 0.0 ) return 1.0;
    omegak=1.0-omegam-omegav;
    eta=sqrt(omegam/a+omegav*a*a+omegak);
    return (2.5/dplus(a,omegam,omegav)-1.5*omegam/a-omegak)/(eta*eta);
}

static double dladt( double a, double omegam, double omegav ) {
    double eta;
    eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav);
    return a * eta;
    }

static FIO_SPECIES graficSpecies(FIO fio) {
    fioGrafic *gio = (fioGrafic *)fio;
    assert(fio->eFormat == FIO_FORMAT_GRAFIC && fio->eMode==FIO_MODE_READING);
    if (gio->iOrder<gio->fio.nSpecies[FIO_SPECIES_SPH]) return FIO_SPECIES_SPH;
    else if (gio->iOrder<gio->fio.nSpecies[FIO_SPECIES_SPH]+gio->fio.nSpecies[FIO_SPECIES_DARK])
	return FIO_SPECIES_DARK;
    else return FIO_SPECIES_LAST;
    }

static FIO graficOpenDirectory(const char *dirName,double dOmega0,double dOmegab) {
    fioGrafic *gio;
    struct stat s;
    size_t n;
    int i;
    char *fileName;

    /*
    ** GRAFIC files are found in a specific directory, so verify
    ** that a directory was given as input.
    */
    if ( stat(dirName,&s) == 0 ) {
#ifdef _MSC_VER
	if ( !(s.st_mode&_S_IFDIR) ) {
#else
	if ( !S_ISDIR(s.st_mode) ) {
#endif
	    errno = ENOTDIR;
            return NULL;
	    }
	}

    gio = malloc(sizeof(fioGrafic));
    assert(gio!=NULL);
    gio->fio.eFormat = FIO_FORMAT_GRAFIC;
    gio->fio.eMode   = FIO_MODE_READING;

    gio->fio.fileList.fileInfo = malloc(sizeof(fioFileInfo)*2);
    assert(gio->fio.fileList.fileInfo);
    gio->fio.fileList.fileInfo[0].pszFilename= NULL;
    gio->fio.fileList.nFiles  = 1;

    gio->fio.fcnClose    = graficClose;
    gio->fio.fcnSeek     = graficSeek;
    gio->fio.fcnReadDark = graficReadDark;
    gio->fio.fcnReadSph  = graficReadSph;
    gio->fio.fcnReadStar = fioNoReadStar;
    gio->fio.fcnGetAttr  = graficGetAttr;
    gio->fio.fcnSpecies  = graficSpecies;

    gio->nLevels = 1;
    gio->level = malloc(sizeof(graficLevel));
    gio->level[0].fp_velcx.fp = gio->level[0].fp_velcy.fp = gio->level[0].fp_velcz.fp = NULL;
    gio->level[0].fp_velbx.fp = gio->level[0].fp_velby.fp = gio->level[0].fp_velbz.fp = NULL;

    for( i=0; i<FIO_SPECIES_LAST; i++)
	gio->fio.nSpecies[i] = 0;

    n = strlen(dirName) + 1;
    fileName = malloc(n + 1 + strlen("ic_velcx"));
    assert(fileName!=NULL);
    strcpy(fileName,dirName);
    strcat(fileName,"/");

    strcpy(fileName+n,"ic_velcx");
    if ( !graficOpen(&gio->level[0].fp_velcx,fileName) ) {
	free(fileName);
	graficClose(&gio->fio);
	return NULL;
	}
    strcpy(fileName+n,"ic_velcy");
    if ( !graficOpen(&gio->level[0].fp_velcy,fileName) ) {
	free(fileName);
	graficClose(&gio->fio);
	return NULL;
	}
    strcpy(fileName+n,"ic_velcz");
    if ( !graficOpen(&gio->level[0].fp_velcz,fileName) ) {
	free(fileName);
	graficClose(&gio->fio);
	return NULL;
	}
    gio->iOrder = 0L;

    assert(graficCompare(&gio->level[0].fp_velcx,&gio->level[0].fp_velcy));
    assert(graficCompare(&gio->level[0].fp_velcx,&gio->level[0].fp_velcz));

    gio->fio.nSpecies[FIO_SPECIES_DARK] = (uint64_t)gio->level[0].fp_velcx.hdr.n[0]
	* (uint64_t)gio->level[0].fp_velcx.hdr.n[1]
	* (uint64_t)gio->level[0].fp_velcx.hdr.n[2];


    if (dOmegab>0.0) {
	strcpy(fileName+n,"ic_velbx");
	if ( !graficOpen(&gio->level[0].fp_velbx,fileName) ) {
	    free(fileName);
	    graficClose(&gio->fio);
	    return NULL;
	    }
	strcpy(fileName+n,"ic_velby");
	if ( !graficOpen(&gio->level[0].fp_velby,fileName) ) {
	    free(fileName);
	    graficClose(&gio->fio);
	    return NULL;
	    }
	strcpy(fileName+n,"ic_velbz");
	if ( !graficOpen(&gio->level[0].fp_velbz,fileName) ) {
	    free(fileName);
	    graficClose(&gio->fio);
	    return NULL;
	    }
	gio->fio.nSpecies[FIO_SPECIES_SPH] = gio->fio.nSpecies[FIO_SPECIES_DARK];

	assert(graficCompare(&gio->level[0].fp_velbx,&gio->level[0].fp_velby));
	assert(graficCompare(&gio->level[0].fp_velbx,&gio->level[0].fp_velbz));
	}

    for( i=1; i<FIO_SPECIES_LAST; i++)
	gio->fio.nSpecies[FIO_SPECIES_ALL] += gio->fio.nSpecies[i];
    gio->dTime = gio->level[0].fp_velcx.hdr.astart;

    assert(gio->level[0].fp_velcx.hdr.n[0]==gio->level[0].fp_velcx.hdr.n[1]&&gio->level[0].fp_velcx.hdr.n[1]==gio->level[0].fp_velcx.hdr.n[2]);

    /* Makes position dimensionless (i.e., be between 0 and 1) */
    gio->pFactor2 = 1.0 / (gio->level[0].fp_velcx.hdr.n[0]*gio->level[0].fp_velcx.hdr.dx);
    gio->pFactor1 = gio->level[0].fp_velcx.hdr.astart / (
        fomega(gio->level[0].fp_velcx.hdr.astart,gio->level[0].fp_velcx.hdr.omegam,gio->level[0].fp_velcx.hdr.omegav)
        * gio->level[0].fp_velcx.hdr.H0
        * dladt(gio->level[0].fp_velcx.hdr.astart,gio->level[0].fp_velcx.hdr.omegam,gio->level[0].fp_velcx.hdr.omegav) );

    gio->vFactor  = sqrt(8*M_PI/3) * gio->pFactor2 / (gio->level[0].fp_velcx.hdr.H0*gio->level[0].fp_velcx.hdr.astart);
    gio->mValueCDM  = (gio->level[0].fp_velcx.hdr.omegam-dOmegab) / gio->fio.nSpecies[FIO_SPECIES_DARK];
    gio->mValueBar  = dOmegab / gio->fio.nSpecies[FIO_SPECIES_DARK];
    gio->sValue  = 1.0 / (50.0 * gio->level[0].fp_velcx.hdr.n[0]);
    free(fileName);

    return &gio->fio;
    }


FIO fioGraficOpenMany(int nFiles, const char * const *dirNames,double dOmega0,double dOmegab) {
    fioGrafic *gio;
    int i;

    gio = malloc(sizeof(fioGrafic));
    assert(gio!=NULL);
    gio->fio.eFormat = FIO_FORMAT_GRAFIC;
    gio->fio.eMode   = FIO_MODE_READING;

    fileScan(&gio->fio.fileList,nFiles,dirNames);

    /* Verify that all "files" are really directories */
    for( i=0; i<gio->fio.fileList.nFiles; i++) {
	struct stat s;
	/* The file/directory needs to exist */
	if ( stat(gio->fio.fileList.fileInfo[i].pszFilename,&s) != 0 ) return NULL;
#ifdef _MSC_VER
	if ( !(s.st_mode&_S_IFDIR) ) return NULL;
#else
	if ( !S_ISDIR(s.st_mode) ) return NULL;
#endif
	}
    gio->nLevels = gio->fio.fileList.nFiles;
    gio->level = malloc(gio->nLevels*sizeof(graficLevel));
    assert(gio->level);
    return &gio->fio;
    }


/******************************************************************************\
** Generic Routines
\******************************************************************************/

/* Attempt to determinate the file type by examining it */
FIO fioOpenMany(int nFiles, const char * const *fileNames,
		double dOmega0,double dOmegab) {
    struct stat s;
    const char *fileName;
    fioFileList fileList;

    /* Turn any wildcard file names into real files before we do anything */
    fileScan(&fileList,nFiles,fileNames);

    /* Test only the first file; they need to be the same type anyway */
    fileName = fileList.fileInfo[0].pszFilename;

    /* The file/directory needs to exist */
    if ( stat(fileName,&s) != 0 ) return NULL;

    /* If given a directory, then it must be a GRAFIC file */
#ifdef _MSC_VER
	if ( s.st_mode&_S_IFDIR ) {
#else
	if ( S_ISDIR(s.st_mode) ) {
#endif
	return graficOpenDirectory(fileName,dOmega0,dOmegab);
	}

#ifdef USE_HDF5
    else if ( H5Fis_hdf5(fileName) ) {
	return hdf5Open(&fileList);
	}
#endif

    /* Try tipsy as a last resort */
    else {
	return tipsyOpen(&fileList);
	}

    return NULL;
    }

FIO fioOpen(const char *fileName,double dOmega0,double dOmegab) {
    return fioOpenMany(1,&fileName,dOmega0,dOmegab);
    }
