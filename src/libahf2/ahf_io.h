#ifndef AHF_IO_H
#define AHF_IO_H

/* Required for the type definitions */
#include "../tdef.h"

extern void
ahf_io_WriteProfiles(const char    *fprefix,
                     HALO          *halos,
                     unsigned long *idx,
                     int           numHalos);

extern void
ahf_io_WriteParticles(const char    *fprefix,
                      HALO          *halos,
                      unsigned long *idx,
                      int           numHalos);

extern void
ahf_io_WriteParticlesSTARDUST(const char    *fprefix,
                              HALO          *halos,
                              unsigned long *idx,
                              int           numHalos);

extern void
ahf_io_WriteHalos(const char    *fprefix,
                  HALO          *halos,
                  unsigned long *idx,
                  int           numHalos);
#ifdef AHFdisks
extern void
ahf_io_WriteDisks(const char    *fprefix,
                  HALO          *halos,
                  unsigned long *idx,
                  int           numHalos);
#endif

#ifdef AHFcentrefile
extern void
ahf_io_WriteCenterfile(const char *fprefix,
                       HALO       *halos,
                       int        numHalos);

#endif

#if ((defined AHFsubstructure) && !defined AHFrestart)
extern void
ahf_io_WriteSubstructure(const char    *fprefix,
                         HALO          *halos,
                         unsigned long *idx,
                         int           numHalos);

#endif

#ifdef AHFgeom
extern void
ahf_io_WriteHalosGeom(const char    *fprefix,
                      HALO          *halos,
                      unsigned long *idx,
                      int           numHalos);
#endif /* AHFgeom */

#ifdef AHFbinary
void ahf_binwrite_open_files (FILE **f, FILE **f_info, char *prefix, char *suffix);
void ahf_binwrite_profiles   (char *prefix, HALO *halos, unsigned long *idx, int numHalos);
void ahf_binwrite_particles  (char *prefix, HALO *halos, unsigned long *idx, int numHalos);
void ahf_binwrite_halos      (char *prefix, HALO *halos, unsigned long *idx, int numHalos);
#endif /* AHFbinary */

#endif   /* AHF_IO_H */
