#ifndef INCLUDE_LIBIO_H
#define INCLUDE_LIBIO_H

#include "include.h"

int      read_particles         (char filename[MAXSTRING], int);
int      read_particles_bin     (char filename[MAXSTRING], int);
int      write_mtree            (char OutFile[MAXSTRING]);

#endif
