/**
 * \file io_util.c
 *
 * Implements functions to perfom input/output operations.
 */


/**********************************************************************\
 *    Includes                   * 
\**********************************************************************/
#include <stdlib.h>
#include <string.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif

#include "io_util.h"
#ifdef WITH_HDF
#include <hdf5.h>
#endif

/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/
#define LOCAL_MAX_NUMFILES_TO_TEST 10000


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern void
io_util_sexchange(void *p, size_t s)
{
	int n;
	unsigned char ptmp,*pc;

	pc = (unsigned char *)p;

	for(n=0; n < s/2; n++) {
		ptmp = pc[n];
		pc[n] = pc[s - n - 1];
		pc[s - n - 1] = ptmp;
	}

	return;
}

extern int
io_util_readbytes(FILE *f, void *n, size_t len, io_file_swap_t swap)
{
	if (fread(n, len, 1, f) != 1)
		return 0;

	if (swap == IO_FILE_IS_SWAPPED)
		io_util_sexchange(n, len);

	return 1;
}

extern size_t
io_util_readstring(FILE *f, char *s, size_t n)
{
	int c;
	size_t i;

	s[0] = '\0';
	for (i=0; i<n; i++) {
		c = fgetc(f);
		if (c == EOF)
			return i;
		s[i] = c;
		s[i+1] = '\0';
	}

	return i;
}

extern size_t
io_util_readline(FILE *f, char *s, size_t n)
{
	int c;
	size_t i;

	s[0] = '\0';
	for (i=0; i<n; i++) {
		c = fgetc(f);
		if ((c == '\n') || (c == EOF))
			return i;
		s[i] = c;
		s[i+1] = '\0';
	}
	return i;
}

extern char *
io_util_strdup(const char *str)
{
	char *dummy;

	dummy = (char *)malloc((strlen(str)+1)*sizeof(char));
	if (dummy == NULL)
		return NULL;

	return strcpy(dummy, str);
}

extern char *
io_util_split_pathfname(const char *str, char **path, char **stem)
{
	int64_t len, i;

	if ( (str == NULL) || (path == NULL) || (stem == NULL) )
		return NULL;

	len = strlen(str);

	/* Figure out where to cut */
	if (str[len-1] == '/') {
		*path = io_util_strdup(str);
		*stem = io_util_strdup("");
		return *path;
	}
	i = 1;
	while ( (i <= len-1) && (str[len-1-i] != '/'))
		i++;
	if ( (i == len) && (str[0] != '/') ) {
		*path = io_util_strdup("./");
	} else {
		*path = (char *)malloc(sizeof(char) * (len-i+1));
		if (*path == NULL)
			return NULL;
		strncpy(*path, str, len-i);
		(*path)[len-i] = '\0';
	}

	*stem = (char *)malloc(sizeof(char) * (i+1));
	if ( *stem == NULL) {
		free(*path);
		*path = NULL;
		return NULL;
	}
	strncpy(*stem, str+(len-i), i);
	(*stem)[i] = '\0';

	return *path;
}

extern int32_t
io_util_findfiles(const char *path,
                  const char *stem,
                  const char *format,
                  const char *suffix,
                  char ***fnames)
{
	int32_t numfiles, i;
	int32_t len;
	char *buff;
	FILE *f;
	char *full_format;

	if ( (path == NULL) || (stem == NULL) || (format == NULL) || (suffix == NULL) ||(fnames == NULL))
		return INT32_C(0);

	/* Adding additional 10 bytes for the number and the '\0' */
	len  = strlen(path);
	len += strlen(stem);
	len += strlen(suffix);
	buff = (char *)malloc(sizeof(char)*(len+10));
	if (buff == NULL)	
		return INT32_C(0);

	/* Construct the full format string used */
	len  = strlen("%%s%%s%%s") + 1;
	len += strlen(format);
	full_format = (char *)malloc(sizeof(char)*(len));
	if (full_format == NULL) {
		free(buff);
		return INT32_C(0);
	}
	/*
	 * Used format is: '%s%s<format>%s', i.e. for format = %i:
	 * fullformat = "%s%s%i%s"
	 */
	sprintf(full_format, "%%s%%s%s%%s", format);

	i = 0;
	f = NULL;

	do  {
		if (f != NULL)
     {
#ifdef DEBUG_MGADGET
      fprintf(stderr,"io_util_findfiles: now closing %s\n",buff);
#endif
			fclose(f);
     }
		sprintf(buff, full_format, path, stem, i, suffix);
		f=fopen(buff, "r");
#ifdef DEBUG_MGADGET
    if(f!=NULL)
      fprintf(stderr,"io_util_findfiles: successfully opened %s\n",buff);
#endif
		i++;
	} while ( (i<LOCAL_MAX_NUMFILES_TO_TEST) && (f != NULL) );
	i--;

	*fnames = (char **)malloc(sizeof(char *) * i);
	if (*fnames == NULL) {
		free(buff);
		return INT32_C(0);
	}
	
	for (numfiles=0; numfiles<i; numfiles++) {
		sprintf(buff, full_format, path, stem, numfiles, suffix);
		(*fnames)[numfiles] = io_util_strdup(buff);
		if ( (*fnames)[numfiles] == NULL ) {
			while (numfiles > 0) {
				numfiles--;
				free((*fnames)[numfiles]);
			}
			free(*fnames);
			free (buff);
			return INT32_C(0);
		}
	}

	/* Clean up */
	free(full_format);
	free(buff);
	
	return numfiles;
}

extern bool *
io_util_getminearr(int32_t numfiles)
{
	int32_t i = 0;
	bool *dummy = NULL;
	int rank = 0;
	int size = 1;

	/* First get the memory */
	dummy = (bool *)malloc(sizeof(bool)*numfiles);
	if (dummy == NULL)
		return NULL;

#ifdef WITH_MPI
	/* Get the MPI dimensions */
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/* Now figure out the things to read and not to read */
	for (i=0; i<numfiles; i++) {
		if (rank == i%size)
			dummy[i] = true;
		else
			dummy[i] = false;
	}

	return dummy;
}

#ifdef WITH_HDF5

void io_util_readhdf5(io_logging_t log,
                    uint8_t select,
                    void* f,
                    char dataset_name[],
                    int type,
                    int dataset_num_values,
                    hid_t hdf5_datatype,
                    hid_t hdf5_grp[],
                    void * CommBuffer)
{
    int rank;
    hid_t hdf5_dataset, hdf5_dataspace_in_memory, hdf5_dataspace_in_file;
    hsize_t dims[2], count[2], start[2];

    /* Read dataset from particle group */
    hdf5_dataset = H5Dopen(hdf5_grp[type], dataset_name);

    // cast to the right structure
    switch (select) {
        case 0:
            // gizmo file
            dims[0] = ((io_gizmo_t) f)->header->np[type];
            break;
        case 1:
            // pkdgrav file
            dims[0] = ((io_pkdgrav_t) f)->header->np[type];
            break;
        default:
            fprintf(stderr, "Error in casting void pointer, wrong flag %u\n", select);
            exit(0);
    }

    dims[1] = dataset_num_values;
    if (dims[1] == 1)
    {
        rank = 1;
    } else {
        rank = 2;
    }

    hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
    hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);

    start[0] = 0; // was pcsum;
    start[1] = 0;

    H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, dims, NULL);

    H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory, hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);

    H5Tclose(hdf5_datatype);
    H5Sclose(hdf5_dataspace_in_memory);
    H5Sclose(hdf5_dataspace_in_file);
    H5Dclose(hdf5_dataset);
}

void io_util_readhdf5_gizmo(io_logging_t log,
                    io_gizmo_t f,
                    char dataset_name[],
                    int type,
                    int dataset_num_values,
                    hid_t hdf5_datatype,
                    hid_t hdf5_grp[],
                    void * CommBuffer)
{

    io_util_readhdf5(log,
                    0,
                    (void*) f,
                    dataset_name,
                    type,
                    dataset_num_values,
                    hdf5_datatype,
                    hdf5_grp,
                    CommBuffer);

}

void io_util_readhdf5_pkdgrav(io_logging_t log,
                    io_pkdgrav_t f,
                    char dataset_name[],
                    int type,
                    int dataset_num_values,
                    hid_t hdf5_datatype,
                    hid_t hdf5_grp[],
                    void * CommBuffer)
{

    io_util_readhdf5(log,
                        1,
                        f,
                        dataset_name,
                        type,
                        dataset_num_values,
                        hdf5_datatype,
                        hdf5_grp,
                        CommBuffer);

}
#endif // WITH_HDF5

/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
