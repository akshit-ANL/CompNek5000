/* Dimension-generic implementation. Included by amr_p4est_2d.c and, after
   p4est_to_p8est.h, by amr_p4est_3d.c. */

#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#include <p4est_bits.h>
#define AMRFN(n) amr2_##n
#else
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#include <p8est_bits.h>
#define AMRFN(n) amr3_##n
#endif

#include <stdlib.h>
#include <string.h>
#include "amr_internal.h"

typedef struct
{
  sc_MPI_Comm comm;
  p4est_connectivity_t *conn;
  p4est_t *p4est;
  p4est_ghost_t *ghost;
  p4est_mesh_t *mesh;
  int maxlevel;

  int nold;                      /* snapshot of leaves before the last adapt */
  int *old_tree, *old_level, *old_x, *old_y, *old_z;
  int *map_type, *map_first, *map_count; /* new leaf -> old leaves */
}
amr_state_t;

static amr_state_t S;

static void
free_mesh (void)
{
  if (S.mesh != NULL) {
    p4est_mesh_destroy (S.mesh);
    S.mesh = NULL;
  }
  if (S.ghost != NULL) {
    p4est_ghost_destroy (S.ghost);
    S.ghost = NULL;
  }
}

static void
free_snapshot (void)
{
  free (S.old_tree);
  free (S.old_level);
  free (S.old_x);
  free (S.old_y);
  free (S.old_z);
  free (S.map_type);
  free (S.map_first);
  free (S.map_count);
  S.old_tree = S.old_level = S.old_x = S.old_y = S.old_z = NULL;
  S.map_type = S.map_first = S.map_count = NULL;
  S.nold = 0;
}

static void
init_flag (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  *(int *) quadrant->p.user_data = 0;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  return *(int *) quadrant->p.user_data == 1 && quadrant->level < S.maxlevel;
}

static int
coarsen_fn (p4est_t * p4est, p4est_topidx_t which_tree,
            p4est_quadrant_t * quadrant[])
{
  int i;
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    if (*(int *) quadrant[i]->p.user_data != -1) {
      return 0;
    }
  }
  return 1;
}

int
AMRFN (init) (sc_MPI_Comm comm, int nbx, int nby, int nbz, int px, int py,
              int pz, int minlevel, int maxlevel)
{
  memset (&S, 0, sizeof (S));
  S.comm = comm;
  S.maxlevel = maxlevel;
#ifndef P4_TO_P8
  S.conn = p4est_connectivity_new_brick (nbx, nby, px, py);
#else
  S.conn = p8est_connectivity_new_brick (nbx, nby, nbz, px, py, pz);
#endif
  S.p4est = p4est_new_ext (comm, S.conn, 0, minlevel, 1, sizeof (int),
                           init_flag, NULL);
  return (int) S.p4est->local_num_quadrants;
}

void
AMRFN (finalize) (void)
{
  free_mesh ();
  free_snapshot ();
  if (S.p4est != NULL) {
    p4est_destroy (S.p4est);
  }
  if (S.conn != NULL) {
    p4est_connectivity_destroy (S.conn);
  }
  memset (&S, 0, sizeof (S));
}

int
AMRFN (num_local) (void)
{
  return (int) S.p4est->local_num_quadrants;
}

long
AMRFN (num_global) (void)
{
  return (long) S.p4est->global_num_quadrants;
}

/* visit all local leaves in order; fn(user, index, tree, quadrant) */
typedef void (*leaf_fn) (void *user, int idx, int tree,
                         p4est_quadrant_t * q);

static void
for_each_leaf (leaf_fn fn, void *user)
{
  p4est_topidx_t t;
  size_t i;
  int idx = 0;

  for (t = S.p4est->first_local_tree; t <= S.p4est->last_local_tree; ++t) {
    p4est_tree_t *tree = p4est_tree_array_index (S.p4est->trees, t);
    for (i = 0; i < tree->quadrants.elem_count; ++i) {
      fn (user, idx++, (int) t,
          p4est_quadrant_array_index (&tree->quadrants, i));
    }
  }
}

typedef struct
{
  int *tree, *level, *qx, *qy, *qz;
}
leaves_out_t;

static void
leaves_cb (void *user, int idx, int tree, p4est_quadrant_t * q)
{
  leaves_out_t *o = (leaves_out_t *) user;
  o->tree[idx] = tree;
  o->level[idx] = q->level;
  o->qx[idx] = q->x;
  o->qy[idx] = q->y;
#ifdef P4_TO_P8
  o->qz[idx] = q->z;
#else
  o->qz[idx] = 0;
#endif
}

void
AMRFN (leaves) (int *tree, int *level, int *qx, int *qy, int *qz)
{
  leaves_out_t o;
  o.tree = tree;
  o.level = level;
  o.qx = qx;
  o.qy = qy;
  o.qz = qz;
  for_each_leaf (leaves_cb, &o);
}

typedef struct
{
  int target;
  double *xlo, *xhi;
}
box_out_t;

static void
box_cb (void *user, int idx, int tree, p4est_quadrant_t * q)
{
  box_out_t *o = (box_out_t *) user;
  p4est_qcoord_t len;

  if (idx != o->target) {
    return;
  }
  len = P4EST_QUADRANT_LEN (q->level);
#ifndef P4_TO_P8
  p4est_qcoord_to_vertex (S.conn, tree, q->x, q->y, o->xlo);
  p4est_qcoord_to_vertex (S.conn, tree, q->x + len, q->y + len, o->xhi);
#else
  p8est_qcoord_to_vertex (S.conn, tree, q->x, q->y, q->z, o->xlo);
  p8est_qcoord_to_vertex (S.conn, tree, q->x + len, q->y + len, q->z + len,
                          o->xhi);
#endif
}

void
AMRFN (leaf_box) (int leaf, double *xlo, double *xhi)
{
  box_out_t o;
  o.target = leaf;
  o.xlo = xlo;
  o.xhi = xhi;
  for_each_leaf (box_cb, &o);
}

static void
set_flag_cb (void *user, int idx, int tree, p4est_quadrant_t * q)
{
  *(int *) q->p.user_data = ((const int *) user)[idx];
}

static void
snapshot (void)
{
  int n = (int) S.p4est->local_num_quadrants;
  leaves_out_t o;

  free_snapshot ();
  S.nold = n;
  S.old_tree = (int *) malloc (sizeof (int) * (n + 1));
  S.old_level = (int *) malloc (sizeof (int) * (n + 1));
  S.old_x = (int *) malloc (sizeof (int) * (n + 1));
  S.old_y = (int *) malloc (sizeof (int) * (n + 1));
  S.old_z = (int *) malloc (sizeof (int) * (n + 1));
  o.tree = S.old_tree;
  o.level = S.old_level;
  o.qx = S.old_x;
  o.qy = S.old_y;
  o.qz = S.old_z;
  for_each_leaf (leaves_cb, &o);
}

static void
quad_from (p4est_quadrant_t * q, int level, int x, int y, int z)
{
  memset (q, 0, sizeof (*q));
  q->level = (int8_t) level;
  q->x = x;
  q->y = y;
#ifdef P4_TO_P8
  q->z = z;
#endif
}

/* relation of old leaf (index o) to new quadrant n in the same tree:
   0 unrelated, 1 equal, 2 old is ancestor of new, 3 old is descendant */
static int
related (int o, int ntree, const p4est_quadrant_t * n)
{
  p4est_quadrant_t oq;

  if (S.old_tree[o] != ntree) {
    return 0;
  }
  quad_from (&oq, S.old_level[o], S.old_x[o], S.old_y[o], S.old_z[o]);
  if (p4est_quadrant_is_equal (&oq, n)) {
    return 1;
  }
  if (p4est_quadrant_is_ancestor (&oq, n)) {
    return 2;
  }
  if (p4est_quadrant_is_ancestor (n, &oq)) {
    return 3;
  }
  return 0;
}

typedef struct
{
  int o;                        /* running index into the old leaves */
}
map_state_t;

static void
map_cb (void *user, int idx, int ntree, p4est_quadrant_t * n)
{
  map_state_t *m = (map_state_t *) user;
  p4est_quadrant_t oq;
  int rel;

  /* skip old leaves that precede n in (tree, z-order) and are unrelated */
  while (m->o < S.nold) {
    rel = related (m->o, ntree, n);
    if (rel) {
      break;
    }
    quad_from (&oq, S.old_level[m->o], S.old_x[m->o], S.old_y[m->o],
               S.old_z[m->o]);
    if (S.old_tree[m->o] < ntree ||
        (S.old_tree[m->o] == ntree && p4est_quadrant_compare (&oq, n) < 0)) {
      ++m->o;
    }
    else {
      break;
    }
  }
  S.map_first[idx] = -1;
  S.map_count[idx] = 0;
  S.map_type[idx] = -1;
  if (m->o >= S.nold) {
    return;
  }
  rel = related (m->o, ntree, n);
  if (rel == 1) {
    S.map_type[idx] = 0;
    S.map_first[idx] = m->o;
    S.map_count[idx] = 1;
    ++m->o;
  }
  else if (rel == 2) {
    /* old leaf stays the parent for further new leaves inside it */
    S.map_type[idx] = 1;
    S.map_first[idx] = m->o;
    S.map_count[idx] = 1;
  }
  else if (rel == 3) {
    int k = 0;
    S.map_type[idx] = 2;
    S.map_first[idx] = m->o;
    while (m->o < S.nold && related (m->o, ntree, n) == 3) {
      ++k;
      ++m->o;
    }
    S.map_count[idx] = k;
  }
}

int
AMRFN (adapt) (const int *flag)
{
  int n;
  map_state_t m;

  free_mesh ();
  snapshot ();

  for_each_leaf (set_flag_cb, (void *) flag);

  p4est_coarsen_ext (S.p4est, 0, 0, coarsen_fn, init_flag, NULL);
  p4est_refine_ext (S.p4est, 0, S.maxlevel, refine_fn, init_flag, NULL);
  p4est_balance_ext (S.p4est, P4EST_CONNECT_FACE, init_flag, NULL);

  n = (int) S.p4est->local_num_quadrants;
  S.map_type = (int *) malloc (sizeof (int) * (n + 1));
  S.map_first = (int *) malloc (sizeof (int) * (n + 1));
  S.map_count = (int *) malloc (sizeof (int) * (n + 1));
  m.o = 0;
  for_each_leaf (map_cb, &m);
  return n;
}

void
AMRFN (map) (int leaf, int *type, int *first, int *count)
{
  *type = S.map_type[leaf];
  *first = S.map_first[leaf];
  *count = S.map_count[leaf];
}

int
AMRFN (build_mesh) (void)
{
  p4est_mesh_params_t params;

  free_mesh ();
  S.ghost = p4est_ghost_new (S.p4est, P4EST_CONNECT_FACE);
  p4est_mesh_params_init (&params);
  params.compute_tree_index = 1;
  params.btype = P4EST_CONNECT_FACE;
  S.mesh = p4est_mesh_new_params (S.p4est, S.ghost, &params);
  return (int) S.mesh->ghost_num_quadrants;
}

int
AMRFN (num_ghost) (void)
{
  return (int) S.mesh->ghost_num_quadrants;
}

int
AMRFN (ghost_rank) (int g)
{
  return S.mesh->ghost_to_proc[g];
}

void
AMRFN (face_info) (int leaf, int face, int *type, int *nface, int *orient,
                   int *sub, int *nbr, int *nnbr)
{
  const int F = P4EST_FACES;
  const int H = P4EST_HALF;
  const int FR = F * H;
  const p4est_locidx_t q2q = S.mesh->quad_to_quad[F * leaf + face];
  const int v = (int) S.mesh->quad_to_face[F * leaf + face];
  int k;

  *sub = 0;
  *nnbr = 1;
  if (v >= 0 && v < FR) {
    *nface = v % F;
    *orient = v / F;
    if (q2q == leaf && *nface == face) {
      *type = 0;
      *nnbr = 0;
      nbr[0] = -1;
    }
    else {
      *type = 1;
      nbr[0] = (int) q2q;
    }
  }
  else if (v >= FR) {
    const int vv = v - FR;
    *type = 2;
    *sub = vv / FR;
    *orient = (vv % FR) / F;
    *nface = vv % F;
    nbr[0] = (int) q2q;
  }
  else {
    const int vv = v + FR;
    const p4est_locidx_t *half =
      (const p4est_locidx_t *) sc_array_index (S.mesh->quad_to_half, q2q);
    *type = 3;
    *orient = vv / F;
    *nface = vv % F;
    *nnbr = H;
    for (k = 0; k < H; ++k) {
      nbr[k] = (int) half[k];
    }
  }
}

typedef struct
{
  long long sum;
}
hash_t;

static void
hash_cb (void *user, int idx, int tree, p4est_quadrant_t * q)
{
  hash_t *h = (hash_t *) user;
  unsigned long long x = 1469598103934665603ULL;
  unsigned long long v[5];
  int i;

  v[0] = (unsigned long long) tree;
  v[1] = (unsigned long long) q->level;
  v[2] = (unsigned long long) q->x;
  v[3] = (unsigned long long) q->y;
#ifdef P4_TO_P8
  v[4] = (unsigned long long) q->z;
#else
  v[4] = 0;
#endif
  for (i = 0; i < 5; ++i) {
    x ^= v[i] + 0x9e3779b97f4a7c15ULL + (x << 6) + (x >> 2);
    x *= 1099511628211ULL;
  }
  h->sum += (long long) (x & 0xFFFFFFFFFFULL);
}

void
AMRFN (checksum) (long *nleaf, long *hash)
{
  hash_t h;
  long long glob = 0;

  h.sum = 0;
  for_each_leaf (hash_cb, &h);
  sc_MPI_Allreduce (&h.sum, &glob, 1, sc_MPI_LONG_LONG_INT, sc_MPI_SUM,
                    S.comm);
  *nleaf = (long) S.p4est->global_num_quadrants;
  *hash = (long) glob;
}
