# dmr2d: double Mach reflection

Woodward & Colella (1984), J. Comput. Phys. 54, 115. A Mach 10 shock hits a
30-degree ramp. The frame is wedge-aligned: the domain is [0,4] x [0,1], the
bottom is a reflecting wall for x > 1/6 and post-shock inflow for x < 1/6, the
top boundary follows the exact moving shock, the left is post-shock inflow and
the right is outflow. Compare density at t = 0.2.

- Gas: calorically perfect, gamma = 1.4, R = 1, one species, no chemistry.
- Pre-shock: rho = 1.4, p = 1, u = 0. Post-shock: rho = 8, p = 116.5,
  |u| = 8.25 at 30 degrees below the x axis.
- Cantera is not needed. `src/cantera_stubs.f` provides the link symbols and
  aborts if a Cantera routine is ever called.
- Shock capturing: entropy viscosity (`ifavisc`) plus the positivity limiter,
  with the current LLF interface flux. This is the baseline for later work
  (HLL/HLLC blending, AMR).

## Status

Written without a compiler or MPI, so it has not been built or run. Expect to
fix small things on the first build.

## Mesh

96 x 24 elements (h = 1/24, so x = 1/6 falls on an element edge). Build the
tools once, then generate `base.re2` and `base.ma2`:

    cd $NEK_SOURCE_ROOT/tools && ./maketools genbox reatore2 genmap
    # in this directory (each tool prompts for a name; adjust to its prompt):
    genbox      # input.box -> box.rea
    reatore2    # box.rea   -> base.re2
    genmap      # base.re2  -> base.ma2

## Build and run

    ./makenek base
    mpirun -np 4 ./nek5000 base

`makenek` follows `bin/makenek` but uses the stock `nekconfig` (no Cantera).
`PPLIST` is kept the same as det2d; trim it if you do not need PARRSB/HYPRE.

## Tuning

`dt` in `base.par` is fixed (5e-5). Check the printed C/aC/dC from
`compress_cfl` and adjust. SIZE uses lx1 = 4 (p = 3); `lpmin` limits the
smallest MPI rank count (lelt = lelg/lpmin + 3).

## Interface flux selection (riemann-solvers branch)

`userparam05` in `base.par` selects the interface flux:

| value | flux |
|---|---|
| 0 | original LLF with the KEP central part (default) |
| 1 | HLL |
| 2 | HLLC |
| 3 | HLLC blended into HLL by a face pressure-jump sensor |

For 3, `alpha = clamp((s - userparam06)/(userparam07 - userparam06), 0, 1)`
with `s = |pR - pL|/(pR + pL)`, and `F = alpha*F_HLL + (1-alpha)*F_HLLC`.
Run the same mesh with 0, 2 and 3 and compare density at t = 0.2.
