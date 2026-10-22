/* Dimension-specific entry points (implemented in amr_impl.h, compiled once
   as amr_p4est_2d.c and once as amr_p4est_3d.c). Not for use outside amr/. */
#ifndef AMR_INTERNAL_H
#define AMR_INTERNAL_H

#include <sc.h>

#define AMR_MAXH 4 /* max neighbors across one face: 2^(dim-1) */

#define AMR_DECLARE(P)                                                       \
  int P##init (sc_MPI_Comm comm, int nbx, int nby, int nbz, int px, int py, \
               int pz, int minlevel, int maxlevel);                          \
  void P##finalize (void);                                                   \
  int P##num_local (void);                                                   \
  long P##num_global (void);                                                 \
  void P##leaves (int *tree, int *level, int *qx, int *qy, int *qz);        \
  void P##leaf_box (int leaf, double *xlo, double *xhi);                     \
  int P##adapt (const int *flag);                                            \
  void P##map (int leaf, int *type, int *first, int *count);                 \
  int P##build_mesh (void);                                                  \
  int P##num_ghost (void);                                                   \
  int P##ghost_rank (int g);                                                 \
  void P##face_info (int leaf, int face, int *type, int *nface,              \
                     int *orient, int *sub, int *nbr, int *nnbr);            \
  void P##checksum (long *nleaf, long *hash);

AMR_DECLARE (amr2_)
AMR_DECLARE (amr3_)

#endif
