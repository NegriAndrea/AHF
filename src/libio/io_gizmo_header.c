/* $Id: io_gizmo_header.c,v 1.3 2007/02/05 16:18:01 knolli Exp $ */

/**
 * \file io_gizmo_header.c
 *
 * Provides functions for reading and writing the header of Gizmo
 * files.
 */

#ifdef WITH_HDF5

/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

#include "../define.h"

#include "io_gizmo_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/
#define SKIP(f) {fseek(f, 4L, SEEK_CUR);}


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_gizmo_header_t
io_gizmo_header_get(io_logging_t log, io_gizmo_t f)
{
	io_gizmo_header_t dummy;
	long skipsize;
        hid_t fpin, hdf5_headergrp, hdf5_attribute;;

#ifdef FOPENCLOSE
  //fprintf(stderr,"FOPENCLOSE: reading header information from %s ... ",f->fname);
  //fpin = fopen(f->fname,IO_FILE_MODE_READ);
    fpin = H5Fopen(f->fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fpin < 0) {
        io_logging_fatal(log,"io_gizmo_header_get(): could not open file %s for reading",f->fname);
        return NULL;
    }
  //fprintf(stderr,"file successfully opened ");
#else
  fpin = f->file;
#endif
  
	/* Some sanity checks */
	if ((f == NULL) || (f->file < 0))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

	/* Create the header structure array */
	dummy = (io_gizmo_header_t)malloc((size_t)GIZMO_HEADER_SIZE+1); // make +1 larger because of trailing '\0' inserted by io_util_readstring()
	if (dummy == NULL) {
		io_logging_memfatal(log, "Gizmo header structure %d",GIZMO_HEADER_SIZE+1);
		return NULL;
	}
    
        hdf5_headergrp = H5Gopen(fpin, "/Header");

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, dummy->np);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
        H5Aread(hdf5_attribute, H5T_NATIVE_UINT, dummy->nall);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
        H5Aread(hdf5_attribute, H5T_NATIVE_UINT, dummy->nallhighw);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
        H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, dummy->massarr);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
        H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &dummy->expansion);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &dummy->numfiles);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_DoublePrecision");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &dummy->flagdoubleprecision);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Redshift");
        H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &dummy->redshift);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "BoxSize");
        H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &dummy->boxsize);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Omega0");
        H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &dummy->omega0);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "OmegaLambda");
        H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &dummy->omegalambda);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "HubbleParam");
        H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &dummy->hubbleparameter);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Sfr");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &dummy->flagsfr);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Feedback");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &dummy->flagfeedback);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Cooling");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &dummy->flagcooling);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_StellarAge");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &dummy->flagstellarage);
        H5Aclose(hdf5_attribute);

        hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_Metals");
        H5Aread(hdf5_attribute, H5T_NATIVE_INT, &dummy->flagmetals);
        H5Aclose(hdf5_attribute);

        H5Gclose(hdf5_headergrp);

/* @TODO: Rennehan : BEGIN BLOCK READING CAN REMOVE BELOW */
        //dummy->flagentropyu = 0;

	f->header = dummy;
  
#ifdef FOPENCLOSE
        H5Fclose(fpin);
        //fprintf(stderr,"and closed (temporarily)\n");
#endif

	return dummy;
}

extern void
io_gizmo_header_del(io_logging_t log, io_gizmo_header_t *header)
{
	if ( (header == NULL) || (*header == NULL) )
		return;

	free(*header);

	*header = NULL;

	return;
}

extern void
io_gizmo_header_write(io_logging_t log,
                       io_gizmo_header_t header,
                       io_gizmo_t f)
{
	if ( (header == NULL) || (f == NULL))
		return;

	if (f->mode != IO_FILE_WRITE)
		return;

	if (f->file < 0)
		return;

	if (header != f->header) {
		io_logging_msg(log, INT32_C(1),
		               "Writing a different header than stored in "
		               "the file object to the file.");
	}

	/* TODO: Write the header */

	return;
}


extern void
io_gizmo_header_log(io_logging_t log, io_gizmo_header_t header)
{
	io_logging_msg(log, INT32_C(5),
	               "Headerobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  No of particles in file:       %" PRIi32,
	               header->np[0]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[1]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[2]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[3]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[4]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIi32,
	               header->np[5]);
	io_logging_msg(log, INT32_C(5),
	               "  Mass of particle species:      %e",
	               header->massarr[0]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[1]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[2]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[3]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[4]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %e",
	               header->massarr[5]);
	io_logging_msg(log, INT32_C(5),
	               "  Expansion:                     %e",
	               header->expansion);
	io_logging_msg(log, INT32_C(5),
	               "  Redshift:                      %e",
	               header->redshift);
	io_logging_msg(log, INT32_C(5),
	               "  Flagsfr:                       %" PRIi32,
	               header->flagsfr);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Feedback:                 %" PRIi32,
	               header->flagfeedback);
	io_logging_msg(log, INT32_C(5),
	               "  No of particles in total:      %" PRIu32,
	               header->nall[0]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[1]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[2]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[3]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[4]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nall[5]);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Cooling:                  %" PRIi32,
	               header->flagcooling);
	io_logging_msg(log, INT32_C(5),
	               "  Number of files:               %" PRIi32,
	               header->numfiles);
	io_logging_msg(log, INT32_C(5),
	               "  Boxsize:                       %e",
	               header->boxsize);
	io_logging_msg(log, INT32_C(5),
	               "  Omega0:                        %e",
	               header->omega0);
	io_logging_msg(log, INT32_C(5),
	               "  OmegaLambda:                   %e",
	               header->omegalambda);
	io_logging_msg(log, INT32_C(5),
	               "  Hubble parameter:              %e",
	               header->hubbleparameter);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Stellar Age:              %" PRIi32,
	               header->flagstellarage);
	io_logging_msg(log, INT32_C(5),
	               "  Flag Metals:                   %" PRIi32,
	               header->flagmetals);
	io_logging_msg(log, INT32_C(5),
	               "  No of particles in total (HW): %" PRIu32,
	               header->nallhighw[0]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[1]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[2]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[3]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[4]);
	io_logging_msg(log, INT32_C(5),
	               "                                 %" PRIu32,
	               header->nallhighw[5]);
	/*io_logging_msg(log, INT32_C(5),
	               "  Flag Entropy insted U:         %" PRIi32,
	               header->flagentropyu);*/

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/

#endif // WITH_HDF5