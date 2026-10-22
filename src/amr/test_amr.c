/* Standalone test of the AMR layer: refine/coarsen a moving front, then
   check mesh invariants. Usage: test_amr <dim 2|3> [nbx nby nbz]
   Exit status is the number of failed checks. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sc_config.h>
#include "amr_p4est.h"

static int nfail = 0;

#define CHECK(cond, ...)                                                     \
  do {                                                                       \
    if (!(cond)) {                                                           \
      if (nfail < 20) {                                                      \
        printf ("FAIL: ");                                                   \
        printf (__VA_ARGS__);                                                \
        printf ("\n");                                                       \
      }                                                                      \
      ++nfail;                                                               \
    }                                                                        \
  } while (0)

static int dim, nb[3];

static double
box_volume (const double *lo, const double *hi)
{
  double v = 1.;
  int d;
  for (d = 0; d < dim; ++d) {
    v *= hi[d] - lo[d];
  }
  return v;
}

static void
check_mesh (const char *tag)
{
  const int n = amr_num_local ();
  const int nfaces = 2 * dim;
  int i, f, k, nbr[AMR_MAXH], nnbr, type, nface, orient, sub;
  int *tree = malloc (sizeof (int) * n), *level = malloc (sizeof (int) * n);
  int *qx = malloc (sizeof (int) * n), *qy = malloc (sizeof (int) * n);
  int *qz = malloc (sizeof (int) * n);
  double vol = 0., dom = 1., lo[3], hi[3], nlo[3], nhi[3];
  int nbdry = 0, expect_bdry = 0, d;

  amr_build_mesh ();
  amr_leaves (tree, level, qx, qy, qz);

  for (d = 0; d < dim; ++d) {
    dom *= nb[d];
  }
  for (i = 0; i < n; ++i) {
    amr_leaf_box (i, lo, hi);
    vol += box_volume (lo, hi);
    for (f = 0; f < nfaces; ++f) {
      const int ax = f / 2;
      const double plane = (f % 2) ? hi[ax] : lo[ax];
      const double h = hi[ax] - lo[ax];
      amr_face_info (i, f, &type, &nface, &orient, &sub, nbr, &nnbr);
      if (plane == 0. || plane == (double) nb[ax]) {
        ++expect_bdry;
      }
      CHECK ((type == 0) == (plane == 0. || plane == (double) nb[ax]),
             "%s leaf %d face %d: boundary flag %d vs geometry", tag, i, f,
             type == 0);
      if (type == 0) {
        ++nbdry;
        continue;
      }
      if (type == 1) {
        const int j = nbr[0];
        int t2, nf2, or2, sub2, nb2[AMR_MAXH], nn2;
        amr_leaf_box (j, nlo, nhi);
        CHECK (level[j] == level[i], "%s same-size level mismatch", tag);
        CHECK (nface == (f ^ 1), "%s leaf %d face %d: nface %d", tag, i, f,
               nface);
        CHECK (fabs ((f % 2 ? nlo[ax] : nhi[ax]) - plane) < 1e-12,
               "%s leaf %d face %d: neighbor not adjacent", tag, i, f);
        amr_face_info (j, nface, &t2, &nf2, &or2, &sub2, nb2, &nn2);
        CHECK (t2 == 1 && nb2[0] == i && nf2 == f,
               "%s leaf %d face %d: asymmetric same-size link", tag, i, f);
      }
      else if (type == 2) {
        const int j = nbr[0];
        int t2, nf2, or2, sub2, nb2[AMR_MAXH], nn2, found = 0;
        amr_leaf_box (j, nlo, nhi);
        CHECK (level[j] == level[i] - 1, "%s 2:1 violated (%d,%d)", tag,
               level[i], level[j]);
        CHECK (fabs ((f % 2 ? nlo[ax] : nhi[ax]) - plane) < 1e-12,
               "%s leaf %d face %d: large neighbor not adjacent", tag, i, f);
        amr_face_info (j, nface, &t2, &nf2, &or2, &sub2, nb2, &nn2);
        CHECK (t2 == 3 && nf2 == f, "%s large neighbor not type 3", tag);
        for (k = 0; k < nn2; ++k) {
          if (nb2[k] == i) {
            ++found;
          }
        }
        CHECK (found == 1, "%s leaf %d not listed by its large neighbor", tag,
               i);
        (void) h;
      }
      else {
        CHECK (type == 3 && nnbr == (dim == 2 ? 2 : 4), "%s bad type 3", tag);
        for (k = 0; k < nnbr; ++k) {
          const int j = nbr[k];
          int t2, nf2, or2, sub2, nb2[AMR_MAXH], nn2;
          CHECK (level[j] == level[i] + 1, "%s 2:1 violated", tag);
          amr_face_info (j, nface, &t2, &nf2, &or2, &sub2, nb2, &nn2);
          CHECK (t2 == 2 && nb2[0] == i && nf2 == f,
                 "%s leaf %d face %d: small neighbor %d does not point back",
                 tag, i, f, j);
        }
      }
    }
  }
  CHECK (fabs (vol - dom) < 1e-10, "%s volume %g != domain %g", tag, vol,
         dom);
  CHECK (nbdry == expect_bdry, "%s boundary face count", tag);
  printf ("%-28s leaves=%6d  vol=%g  boundary faces=%d\n", tag, n, vol,
          nbdry);
  free (tree);
  free (level);
  free (qx);
  free (qy);
  free (qz);
}

/* old/new relationship from amr_map must be geometrically consistent */
static void
check_map (const double *oldlo, const double *oldhi, int nold,
           const char *tag)
{
  const int n = amr_num_local ();
  int i, d, type, first, count, o;
  double lo[3], hi[3], tot = 0.;

  for (i = 0; i < n; ++i) {
    amr_map (i, &type, &first, &count);
    CHECK (type >= 0, "%s new leaf %d has no old counterpart", tag, i);
    if (type < 0) {
      continue;
    }
    amr_leaf_box (i, lo, hi);
    if (type == 0) {
      for (d = 0; d < dim; ++d) {
        CHECK (fabs (lo[d] - oldlo[3 * first + d]) < 1e-12 &&
               fabs (hi[d] - oldhi[3 * first + d]) < 1e-12,
               "%s type0 box differs", tag);
      }
    }
    else if (type == 1) {
      for (d = 0; d < dim; ++d) {
        CHECK (lo[d] >= oldlo[3 * first + d] - 1e-12 &&
               hi[d] <= oldhi[3 * first + d] + 1e-12,
               "%s type1 not inside parent", tag);
      }
    }
    else {
      double ov = 0.;
      for (o = first; o < first + count; ++o) {
        ov += box_volume (oldlo + 3 * o, oldhi + 3 * o);
      }
      CHECK (fabs (ov - box_volume (lo, hi)) < 1e-12,
             "%s coarsened leaf volume != children volume", tag);
    }
    tot += box_volume (lo, hi);
  }
  (void) nold;
  (void) tot;
}

static void
snapshot_boxes (double **lo, double **hi, int *n)
{
  int i;
  *n = amr_num_local ();
  *lo = malloc (sizeof (double) * 3 * (*n + 1));
  *hi = malloc (sizeof (double) * 3 * (*n + 1));
  for (i = 0; i < *n; ++i) {
    amr_leaf_box (i, *lo + 3 * i, *hi + 3 * i);
  }
}

/* refine near x = xc (width w in x), coarsen far from it */
static int
front_cycle (double xc, double w, int minlevel, const char *tag)
{
  const int n = amr_num_local ();
  int *flag = calloc (n + 1, sizeof (int)), *level = malloc (sizeof (int) * n);
  int *t = malloc (sizeof (int) * n), *qx = malloc (sizeof (int) * n);
  int *qy = malloc (sizeof (int) * n), *qz = malloc (sizeof (int) * n);
  double *olo, *ohi, lo[3], hi[3];
  int i, nold, nnew;

  amr_leaves (t, level, qx, qy, qz);
  for (i = 0; i < n; ++i) {
    double c;
    amr_leaf_box (i, lo, hi);
    c = 0.5 * (lo[0] + hi[0]);
    if (fabs (c - xc) < w) {
      flag[i] = 1;
    }
    else if (level[i] > minlevel) {
      flag[i] = -1;
    }
  }
  snapshot_boxes (&olo, &ohi, &nold);
  nnew = amr_adapt (flag);
  check_map (olo, ohi, nold, tag);
  check_mesh (tag);
  free (olo);
  free (ohi);
  free (flag);
  free (level);
  free (t);
  free (qx);
  free (qy);
  free (qz);
  return nnew;
}

int
main (int argc, char **argv)
{
  int cycle;
  long nleaf, hash;
  char tag[64];
  const int minlevel = 1, maxlevel = 4;

  dim = argc > 1 ? atoi (argv[1]) : 2;
  nb[0] = argc > 2 ? atoi (argv[2]) : 8;
  nb[1] = argc > 3 ? atoi (argv[3]) : (dim == 2 ? 1 : 2);
  nb[2] = argc > 4 ? atoi (argv[4]) : (dim == 3 ? 2 : 1);

#ifdef SC_ENABLE_MPI
  MPI_Init (&argc, &argv);
#endif
  amr_init (sc_MPI_COMM_WORLD, dim, nb[0], nb[1], nb[2], 0, 0, 0, minlevel,
            maxlevel);
  check_mesh ("initial uniform");

  for (cycle = 0; cycle < 4; ++cycle) {
    sprintf (tag, "refine front x=%.1f #%d", 0.5 * nb[0], cycle);
    front_cycle (0.5 * nb[0], 0.4, minlevel, tag);
  }
  for (cycle = 0; cycle < 4; ++cycle) {
    sprintf (tag, "move front x=%.1f #%d", 0.25 * nb[0], cycle);
    front_cycle (0.25 * nb[0], 0.4, minlevel, tag);
  }
  amr_checksum (&nleaf, &hash);
  printf ("checksum: leaves=%ld hash=%ld\n", nleaf, hash);

  amr_finalize ();
#ifdef SC_ENABLE_MPI
  MPI_Finalize ();
#endif
  printf (nfail ? "RESULT: %d FAILED CHECKS\n" : "RESULT: PASS (%d)\n", nfail);
  return nfail != 0;
}
