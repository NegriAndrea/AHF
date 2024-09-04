#ifndef INCLUDE_LIBMTREE_H
#define INCLUDE_LIBMTREE_H

#include "include.h"
int      particle_halo_mapping  (int);
int      cross_correlation      (char OutFile[MAXSTRING]);
void     clean_connection       (uint64_t, int, int);
void     create_mtree           (uint64_t, int, int);
void     create_mtree_qsort     (uint64_t, int, int);

#endif
