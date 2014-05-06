#include "../param.h"
#include "../tdef.h"

#ifndef ALLOC_STRUCT_INCLUDED
#define ALLOC_STRUCT_INCLUDED

partptr c_part       (long block_size);
gasptr  c_gas        (long block_size);

nptr    r_node       (nptr  node_ptr,  long new_size);
pqptr   r_pquad      (pqptr pquad_ptr, long old_size, long new_size);
cqptr   r_cquad      (cqptr cquad_ptr, long old_size, long new_size);
nqptr   r_nquad      (nqptr nquad_ptr, long old_size, long new_size);

nptr    c_node       (long block_size);
pqptr   c_pquad      (long block_size);
cqptr   c_cquad      (long block_size);
nqptr   c_nquad      (long block_size);

void    dest_node    (nptr  *ptr);
void    dest_nquad   (nqptr *ptr);
void    dest_cquad   (cqptr *ptr);
void    dest_pquad   (pqptr *ptr);

void    alloc_quads  (gridls *cur_grid, long l1dim);

void    free_pquad   (pqptr cur_pquad);
void    free_nquad   (nqptr cur_nquad);
void    free_cquad   (cqptr cur_cquad);
void    free_grid    (gridls *cur_grid,   int *no_grids);


#if (defined AHF || defined AHF2)
void    c_profile(HALO *cur_halo, int nbins);
void    dest_profile(HALO *cur_halo);
#endif

#endif

