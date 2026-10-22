/* Dimension dispatcher for the C API in amr_p4est.h. */

#include <p4est.h>
#include <p8est.h>
#include "amr_p4est.h"
#include "amr_internal.h"

static int amr_ndim = 0;
static int amr_sc_ready = 0;

int
amr_init (sc_MPI_Comm comm, int dim, int nbx, int nby, int nbz, int px,
          int py, int pz, int minlevel, int maxlevel)
{
  if (!amr_sc_ready) {
    sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
    p4est_init (NULL, SC_LP_ESSENTIAL); /* shared by p4est and p8est */
    amr_sc_ready = 1;
  }
  amr_ndim = dim;
  return dim == 2 ?
    amr2_init (comm, nbx, nby, nbz, px, py, pz, minlevel, maxlevel) :
    amr3_init (comm, nbx, nby, nbz, px, py, pz, minlevel, maxlevel);
}

#define DISPATCH_VOID(name, args2, args3) \
  do { if (amr_ndim == 2) amr2_##name args2; else amr3_##name args3; } while (0)
#define DISPATCH_RET(name, args2, args3) \
  return amr_ndim == 2 ? amr2_##name args2 : amr3_##name args3

void
amr_finalize (void)
{
  DISPATCH_VOID (finalize, (), ());
  amr_ndim = 0;
}

int
amr_dim (void)
{
  return amr_ndim;
}

int
amr_num_local (void)
{
  DISPATCH_RET (num_local, (), ());
}

long
amr_num_global (void)
{
  DISPATCH_RET (num_global, (), ());
}

void
amr_leaves (int *tree, int *level, int *qx, int *qy, int *qz)
{
  DISPATCH_VOID (leaves, (tree, level, qx, qy, qz),
                 (tree, level, qx, qy, qz));
}

void
amr_leaf_box (int leaf, double *xlo, double *xhi)
{
  DISPATCH_VOID (leaf_box, (leaf, xlo, xhi), (leaf, xlo, xhi));
}

int
amr_adapt (const int *flag)
{
  DISPATCH_RET (adapt, (flag), (flag));
}

void
amr_map (int leaf, int *type, int *first, int *count)
{
  DISPATCH_VOID (map, (leaf, type, first, count),
                 (leaf, type, first, count));
}

int
amr_build_mesh (void)
{
  DISPATCH_RET (build_mesh, (), ());
}

int
amr_num_ghost (void)
{
  DISPATCH_RET (num_ghost, (), ());
}

int
amr_ghost_rank (int g)
{
  DISPATCH_RET (ghost_rank, (g), (g));
}

void
amr_face_info (int leaf, int face, int *type, int *nface, int *orient,
               int *sub, int *nbr, int *nnbr)
{
  DISPATCH_VOID (face_info,
                 (leaf, face, type, nface, orient, sub, nbr, nnbr),
                 (leaf, face, type, nface, orient, sub, nbr, nnbr));
}

void
amr_checksum (long *nleaf, long *hash)
{
  DISPATCH_VOID (checksum, (nleaf, hash), (nleaf, hash));
}
