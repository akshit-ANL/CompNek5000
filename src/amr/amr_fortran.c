/* Fortran-callable entry points (all arguments by reference).
   Names follow the repo convention in core/name.h; the amrf_ prefix keeps
   them distinct from the C API in amr_p4est.h. Indices are 0-based, as in
   the C API; convert at the Fortran call site. */

#include "name.h"
#include "amr_p4est.h"
#include <sc_config.h>

#define amrf_init FORTRAN_UNPREFIXED (amrf_init, AMRF_INIT)
#define amrf_finalize FORTRAN_UNPREFIXED (amrf_finalize, AMRF_FINALIZE)
#define amrf_num_local FORTRAN_UNPREFIXED (amrf_num_local, AMRF_NUM_LOCAL)
#define amrf_leaves FORTRAN_UNPREFIXED (amrf_leaves, AMRF_LEAVES)
#define amrf_leaf_box FORTRAN_UNPREFIXED (amrf_leaf_box, AMRF_LEAF_BOX)
#define amrf_adapt FORTRAN_UNPREFIXED (amrf_adapt, AMRF_ADAPT)
#define amrf_map FORTRAN_UNPREFIXED (amrf_map, AMRF_MAP)
#define amrf_build_mesh FORTRAN_UNPREFIXED (amrf_build_mesh, AMRF_BUILD_MESH)
#define amrf_face_info FORTRAN_UNPREFIXED (amrf_face_info, AMRF_FACE_INFO)
#define amrf_checksum FORTRAN_UNPREFIXED (amrf_checksum, AMRF_CHECKSUM)

/* fcomm: Fortran MPI communicator handle; returns the local leaf count */
int
amrf_init (int *fcomm, int *dim, int *nbx, int *nby, int *nbz, int *px,
           int *py, int *pz, int *minlevel, int *maxlevel)
{
#ifdef SC_ENABLE_MPI
  sc_MPI_Comm comm = MPI_Comm_f2c (*fcomm);
#else
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
#endif
  return amr_init (comm, *dim, *nbx, *nby, *nbz, *px, *py, *pz, *minlevel,
                   *maxlevel);
}

void
amrf_finalize (void)
{
  amr_finalize ();
}

int
amrf_num_local (void)
{
  return amr_num_local ();
}

void
amrf_leaves (int *tree, int *level, int *qx, int *qy, int *qz)
{
  amr_leaves (tree, level, qx, qy, qz);
}

void
amrf_leaf_box (int *leaf, double *xlo, double *xhi)
{
  amr_leaf_box (*leaf, xlo, xhi);
}

int
amrf_adapt (int *flag)
{
  return amr_adapt (flag);
}

void
amrf_map (int *leaf, int *type, int *first, int *count)
{
  amr_map (*leaf, type, first, count);
}

int
amrf_build_mesh (void)
{
  return amr_build_mesh ();
}

void
amrf_face_info (int *leaf, int *face, int *type, int *nface, int *orient,
                int *sub, int *nbr, int *nnbr)
{
  amr_face_info (*leaf, *face, type, nface, orient, sub, nbr, nnbr);
}

void
amrf_checksum (long *nleaf, long *hash)
{
  amr_checksum (nleaf, hash);
}
