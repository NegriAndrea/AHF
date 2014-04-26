#include "../param.h"
#include "../tdef.h"

#ifndef GET_NODES_INCLUDED
#define GET_NODES_INCLUDED

void     get_TSCnodes     (gridls *curgrid, pqptr curpquad, cqptr curcquad, nqptr curnquad, nptr tsc_box[3][3][3], long *z, long *y, long *x);
void     get_MHDnodes     (gridls *cur_grid, pqptr cur_pquad, long z, long y, long x, nptr MHDnodes[5][5][5]);
void     get_SIXnodes     (gridls *cur_grid, pqptr cur_pquad, cqptr cur_cquad, nqptr cur_nquad, nptr tsc_nodes[3][3][3], int z, int y, int x);
boolean  test_mhd         (nptr mhd_nodes[5][5][5]);
boolean  test_tsc         (nptr tsc_nodes[3][3][3]);
int      count_tsc        (nptr tsc_nodes[3][3][3]);
void     check_node       (gridls *cur_grid, nptr node, long z_node, long y_node, long x_node);


# ifdef WITH_MPI
#  include "../libsfc/sfc.h"
/**
 * \brief This function searches in a given grid for a node at
 *        the position given by its sfc key.
 *
 * \param *curgrid  The grid we are searching on.
 * \param key       The SFC key describing the node we are
 *                  looking for.
 * \param ctype     Specifies the type of the SFC key.
 * \param bits      The number of bits used for one dimension of
 *                  the SFC key.
 *
 * \return Returns a pointer to the node in the grid, or NULL,
 *         if no node could be found.
 */
extern nptr
get_node_from_key(gridls *curgrid,
                  sfc_key_t key,
                  sfc_curve_t ctype,
                  uint32_t bits);

/**
 * TODO
 */
extern nptr
get_node_from_pos(gridls *curgrid,
                  uint32_t *pos);

# endif /* WITH_MPI */

#endif

