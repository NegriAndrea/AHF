#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

//*********************************************************************************************************
// DEFINES
//*********************************************************************************************************
#define END_OF_TRUNK                 -1
#define NO_SINGLE_TRUNK              -2

void   get_patch_level_range (patch_t **, int64_t *, int *, int *);
double patch_property        (patch_t);
void   write_patch_geom      (char *, patch_t );
void   get_patch_centre      (patch_t , double *, double *, double *);
#endif
