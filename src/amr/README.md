# AMR layer (p4est)

Status: the forest/neighbor layer is written and tested standalone. It is not
yet connected to the Nek data structures, so no solver run uses it.

## Quick start on the HPC (gcc + OpenMPI)

Nothing here has been run on a cluster yet; this is the intended first session.

    git clone <fork-url> && cd CompNek5000 && git checkout amr
    export NEK_SOURCE_ROOT=$PWD
    # 1. build p4est (downloads p4est 2.8.6; compute nodes may lack internet,
    #    so run on a login node, or use: cd 3rd_party/p4est && ./install fetch)
    cd 3rd_party/p4est && CC=mpicc ./install && cd ../..
    # 2. standalone AMR test (see "Build the test" below), 2D and 3D
    # 3. baseline solver cases, no AMR yet:
    #    dg_examples/sod2d  (Cantera-free; README there)
    #    dg_examples/dmr2d  (Cantera-free double Mach reflection)
    #    userparam05 in base.par selects LLF / HLL / HLLC / blend

Other branches: `riemann-solvers` (HLL/HLLC/blend, already in `amr`) and
`master` (adds `dmr2d`).

## Files

| file | role |
|---|---|
| `amr_p4est.h` | C API (face types, indexing, adapt/map contract) |
| `amr_impl.h`, `amr_p4est_2d.c`, `amr_p4est_3d.c` | one implementation compiled twice (p4est / p8est) |
| `amr_p4est.c` | dispatches on dimension at run time |
| `amr_fortran.c` | `amrf_*` entry points for the F77 solver (0-based indices) |
| `test_amr.c` | standalone invariant test (2D and 3D) |

## What is verified (test_amr, serial p4est 2.8.6)

Refine/coarsen of a moving front plus 2:1 face balance, in 2D and 3D, checking:
leaf volumes sum to the domain; boundary faces match geometry; same-size,
large-neighbor and 4-small-neighbor (3D) links are symmetric; levels differ by
at most 1 across faces; the old-to-new map is geometrically consistent
(equal / inside parent / union of children). Not yet exercised: MPI runs with
ghost leaves, periodic bricks.

## Build the test

    P=$NEK_SOURCE_ROOT/3rd_party/p4est   # after: PPLIST="... P4EST" ./makenek
    mpicc -I. -I$P/include $(ls *.c | grep -v amr_fortran) -L$P/lib -lp4est -lsc -lm -o test_amr
    mpirun -np 4 ./test_amr 2 8 1 1 && mpirun -np 4 ./test_amr 3 4 2 2

(`amr_fortran.c` also needs `core/name.h` on the include path.) On a build
where p4est was configured with MPI, `test_amr` calls `MPI_Init` itself.

## Design decisions made

- Coarse mesh = a brick of trees, one tree per initial Nek element; every leaf
  is a Nek element. General (curved/unstructured) coarse meshes need a
  connectivity built from Nek's vertex data and are not handled yet.
- Only face 2:1 balance is enforced; that is all the face mortars need.
- No repartitioning inside `amr_adapt`: leaves stay on their rank, so the
  old-to-new map is local. Rebalancing needs element data migration and is a
  separate step.
- p4est face numbers (-x,+x,-y,+y,-z,+z) differ from Nek's (-y,+x,+y,-x,-z,+z);
  the Fortran side must translate.

## Next steps

1. Fortran module: build Nek element geometry (xm1/ym1, metrics, `dg_face`,
   normals, areas) from the leaves, starting from the initial mesh.
2. Replace the face `gs_op` exchange with: conforming faces as now, hanging
   faces through an L2-projection mortar (project, flux, project back).
3. Solution transfer from `amr_map` (copy, prolongation, restriction).
4. Indicator and adapt loop; then 2D Sod (`dg_examples/sod2d`), double Mach
   reflection, 3D.
5. MPI: ghost data exchange for off-rank neighbors, then repartitioning.
