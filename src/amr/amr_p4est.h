/* C API to the p4est-based adaptive mesh layer.
 *
 * The coarse mesh is a brick of nbx*nby*nbz unit trees (nbz ignored when
 * dim == 2). Every leaf is a Nek element. Leaves are indexed 0..nloc-1 in
 * p4est order (tree, then z-order); ghost leaves are indexed nloc+g.
 *
 * Face numbering is p4est's: 0:-x 1:+x 2:-y 3:+y (4:-z 5:+z in 3D).
 * Nek's preprocessor order is 1:-y 2:+x 3:+y 4:-x (5:-z 6:+z): the
 * Fortran side must translate.
 *
 * Face types from amr_face_info:
 *   0 domain boundary
 *   1 same-size neighbor: nbr[0]
 *   2 neighbor is LARGER (this leaf is the small side): nbr[0], sub = which
 *     subface of the large neighbor this leaf touches
 *   3 neighbors are SMALLER (this leaf is the large side): nbr[0..nnbr-1],
 *     ordered by this leaf's face corners
 * nface/orient are the neighbor's face number and relative orientation.
 */
#ifndef AMR_P4EST_H
#define AMR_P4EST_H

#ifdef __cplusplus
extern "C" {
#endif

#include <sc.h>

#define AMR_MAXH 4

int  amr_init (sc_MPI_Comm comm, int dim, int nbx, int nby, int nbz,
               int px, int py, int pz, int minlevel, int maxlevel);
void amr_finalize (void);
int  amr_dim (void);
int  amr_num_local (void);
long amr_num_global (void);

/* tree id, level and integer quadrant coordinates of every local leaf */
void amr_leaves (int *tree, int *level, int *qx, int *qy, int *qz);
/* physical (brick) lower/upper corners of a leaf; xlo,xhi hold 3 doubles */
void amr_leaf_box (int leaf, double *xlo, double *xhi);

/* flag[i] = +1 refine leaf i, -1 coarsen (only if all siblings are -1),
 * 0 keep. Coarsens, refines, then 2:1 face-balances. Returns the new local
 * leaf count. Invalidates the mesh: call amr_build_mesh again.
 * Single-rank-consistent: no repartitioning, leaves never change rank. */
int  amr_adapt (const int *flag);
/* How new leaf i relates to the previous mesh (valid after amr_adapt):
 * type 0 = same leaf as old leaf `first`
 * type 1 = descendant of old leaf `first` (count = 1)
 * type 2 = ancestor of `count` old leaves starting at `first` */
void amr_map (int leaf, int *type, int *first, int *count);

/* ghost layer + face neighbor tables; call after init and after adapt */
int  amr_build_mesh (void);
int  amr_num_ghost (void);
int  amr_ghost_rank (int g);
void amr_face_info (int leaf, int face, int *type, int *nface, int *orient,
                    int *sub, int *nbr, int *nnbr);

/* global leaf count and an order-independent hash of (tree,level,coords),
 * for comparing forests across rank counts */
void amr_checksum (long *nleaf, long *hash);

#ifdef __cplusplus
}
#endif
#endif
