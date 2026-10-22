# sod2d: Sod shock tube in a thin 2D strip

Quasi-1D test, x in [0,1], diaphragm at 0.5, gamma = 1.4, t = 0.2. Left state
(rho, u, p) = (1, 0, 1), right (0.125, 0, 0.1). All walls are slip.
Perfect gas, one species, no Cantera (same setup as `dmr2d`).

This is the first AMR test problem. Today it runs on a uniform 64 x 8 mesh and
serves as the reference; the AMR run will start from an 8 x 1 element coarse
mesh and refine near the waves to the same finest size (3 levels).

## Status

Written without a compiler or MPI for the solver, so not yet built or run.
`exact_sod.py` was checked against the known star state
(p* = 0.30313, u* = 0.92745, rho* = 0.42632 / 0.26557).

## Run

    cd $NEK_SOURCE_ROOT/tools && ./maketools genbox reatore2 genmap
    # here: genbox (input.box) -> reatore2 -> genmap, giving base.re2 / base.ma2
    ./makenek base
    mpirun -np 4 ./nek5000 base
    python exact_sod.py            # merges prof_*.dat, prints L1/Linf density error

Set `userparam05` in `base.par` to 0, 2 or 3 to compare LLF, HLLC and the
HLLC/HLL blend against the exact solution.
