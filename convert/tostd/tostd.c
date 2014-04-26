#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include "fio.h"

#define BUFFER_SIZE (8*1024*1024)

#define OPT_DOUBLE 'd'
#define OPT_NATIVE 'n'
#define OPT_HDF5   '5'
#define OPT_POTENTIAL 'p'
#define OPT_DENSITY   'r'

int main( int argc, char *argv[] ) {
    int bError = 0;
    int bDouble = 0;
    int bNative = 0;
    int bHDF5 = 0;
    int bPotential = 0;
    int bDensity = 0;
    uint64_t N, nSph, nDark, nStar, i;
    double dTime;

    uint64_t iOrder;
    double r[3], v[3];
    float fMass, fSoft, fPot, fRho, u, fMetals, fTimer;

    FIO fioIn, fioOut;
    FIO_SPECIES eSpecies;
    int inNameIndex = 0;
    int inNameCount = 0;
    const char *outName = 0;

    //! Parse command line arguments (flags).
    for (;;) {
	int c, option_index=0;

	static struct option long_options[] = {
		{ "double",       0, 0, OPT_DOUBLE},
		{ "native",       0, 0, OPT_NATIVE},
		{ "hdf5",         0, 0, OPT_HDF5},
		{ "potential",    0, 0, OPT_POTENTIAL},
		{ "density"  ,    0, 0, OPT_DENSITY},
		{ NULL,   0, 0, 0 },
	    };

	c = getopt_long( argc, argv, "dn5pr",
			 long_options, &option_index );
	if ( c == -1 ) break;

	switch (c) {
	case OPT_DOUBLE:
	    bDouble = 1;
	    break;
	case OPT_NATIVE:
	    bNative = 1;
	    break;
	case OPT_HDF5:
	    bHDF5 = 1;
	    break;
	case OPT_POTENTIAL:
	    bPotential = 1;
	    break;
	case OPT_DENSITY:
	    bDensity = 1;
	    break;
	default:
	    bError = 1;
	    break;
	    }
	}

    if (bNative && bHDF5) {
	fprintf(stderr, "Specify only one of --hdf5 or --native\n" );
	bError = 1;
	}
#ifndef USE_HDF5
   if (bHDF5) {
        fprintf(stderr, "HDF5 support was not compiled in.\n" );
        bError = 1;
        }
#endif


    if ( optind < argc ) {
	inNameIndex = optind++;
	inNameCount = argc - optind;
	}
    else {
	fprintf(stderr, "Missing input file(s)\n" );
	bError = 1;
	}

    if ( optind < argc )
	outName = argv[argc-1];
    else {
	fprintf(stderr, "Missing Tipsy output file\n" );
	bError = 1;
	}

    if ( bError ) {
	fprintf(stderr, "Usage: %s [-p] <input...> <outtipsy>\n"
		"  -d,--double    Output double precision positions\n"
		"  -n,--native    Output a native tipsy binary\n"
#ifdef USE_HDF5
		"  -5,--hdf5      Output in HDF5 format\n"
		"  -p,--potential Included potentials in HDF5 output\n"
#endif
		,argv[0] );
	exit(1);
	}

    fioIn = fioOpenMany(inNameCount,(const char * const *)&argv[inNameIndex],0.0,0.0);
    if (fioIn==NULL) {
	perror(argv[inNameIndex]);
	exit(errno);
	}
    N     = fioGetN(fioIn,FIO_SPECIES_ALL);
    nSph  = fioGetN(fioIn,FIO_SPECIES_SPH);
    nDark = fioGetN(fioIn,FIO_SPECIES_DARK);
    nStar = fioGetN(fioIn,FIO_SPECIES_STAR);
    if (!fioGetAttr(fioIn,"dTime",FIO_TYPE_DOUBLE,&dTime)) dTime = 0.0;

#ifdef USE_HDF5
    if (bHDF5) {
	int mFlag = FIO_FLAG_COMPRESS_MASS | FIO_FLAG_COMPRESS_SOFT;
	if (bDouble) mFlag |= FIO_FLAG_DOUBLE_POS | FIO_FLAG_DOUBLE_VEL;
	if (bPotential) mFlag |= FIO_FLAG_POTENTIAL;
	if (bDensity) mFlag |= FIO_FLAG_DENSITY;
	fioOut = fioHDF5Create(outName,mFlag);
	}
#else
    if (0) {}
#endif
    else {
	fioOut = fioTipsyCreate(outName,bDouble,!bNative,dTime,nSph,nDark,nStar);
	}

    if (fioOut==NULL) {
	perror(outName);
	exit(errno);
	}
    fioSetAttr(fioOut,"dTime",FIO_TYPE_DOUBLE,&dTime);

    for( i=0; i<N; i++ ) {
        eSpecies = fioSpecies(fioIn);
        switch(eSpecies) {
        case FIO_SPECIES_SPH:
            fioReadSph(fioIn,&iOrder,r,v,&fMass,&fSoft,&fPot,&fRho,&u,&fMetals);
	    fioWriteSph(fioOut,iOrder,r,v,fMass,fSoft,fPot,fRho,u,fMetals);
            break;
        case FIO_SPECIES_DARK:
            fioReadDark(fioIn,&iOrder,r,v,&fMass,&fSoft,&fPot,&fRho);
            fioWriteDark(fioOut,iOrder,r,v,fMass,fSoft,fPot,fRho);
            break;
        case FIO_SPECIES_STAR:
            fioReadStar(fioIn,&iOrder,r,v,&fMass,&fSoft,&fPot,&fRho,&fMetals,&fTimer);
            fioWriteStar(fioOut,iOrder,r,v,fMass,fSoft,fPot,fRho,fMetals,fTimer);
            break;
        default:
            fprintf(stderr,"Unsupported particle type: %d\n",eSpecies);
            abort();
            }
	}

    fioClose(fioOut);
    fioClose(fioIn);

    return 0;
    }
