# Nek5000 + ParaView 5.13 Catalyst 2

This is a new Catalyst 2 adaptor. It uses the current Nek5000 in-situ hooks and
does not reuse the legacy Catalyst 1 API from `AtzoriMarco/InSituPackage`.

## Current data model

Each MPI rank contributes its local spectral elements as an unstructured
Conduit Blueprint mesh. Adjacent GLL points form linear quads in 2D or hexes in
3D. Point fields are:

- `velocity` (three components)
- `temperature` when `ifheat` is enabled

Pressure is deliberately not exposed yet. In the usual Pn-Pn-2 formulation it
must first be mapped from the pressure grid to the velocity grid; the old
adaptor incorrectly treats both arrays as if they had identical extents.

## Build on the HPC system

Load the same MPI/compiler stack used to build ParaView 5.13. The Catalyst API,
Conduit headers, and libraries must come from that installation. In the case's
`makenek`, set variables similar to the following (adjust include/library paths
to the installed ParaView/Catalyst layout):

```sh
FC=mpifort
CC=mpicc
CXX=mpicxx
LD=mpifort

CATALYST2=1
CATALYST2_CXXFLAGS="-I${CATALYST_ROOT}/include/catalyst-2.0 -I${CATALYST_ROOT}/include/conduit"
CATALYST2_LIBS="-L${CATALYST_ROOT}/lib -Wl,-rpath,${CATALYST_ROOT}/lib -lcatalyst -lconduit -lstdc++"
```

Keep the Fortran MPI wrapper as the linker so it supplies the Fortran runtime;
`-lstdc++` supplies the C++ runtime needed by the adaptor and Catalyst.
Some installations place Conduit in versioned include directories or require
additional transitive libraries. Prefer the flags exported by the site module
or Catalyst installation when available.

## Eddy smoke test

The `short_tests/eddy` case is small and already exercises a 2D unsteady flow.
Prepare it in a scratch/run directory using the repository's normal short-test
procedure, then copy the pipeline and build:

```sh
cp "$NEK_SOURCE_ROOT/docs/in_situ/catalyst_pipeline.py" .
mkdir -p catalyst-output
export NEK_CATALYST_SCRIPT="$PWD/catalyst_pipeline.py"
export NEK_CATALYST_OUTPUT="$PWD/catalyst-output"
export NEK_CATALYST_FREQUENCY=10
export CATALYST_IMPLEMENTATION_PATHS="$PARAVIEW_ROOT/lib/catalyst"

./makenek clean
./makenek eddy_uv
mpiexec -n 2 ./nek5000
```

The precise `CATALYST_IMPLEMENTATION_PATHS` directory is installation-specific.
It must contain ParaView's Catalyst implementation library. A successful test
has all of these signals:

1. The log prints `Catalyst 2 initialized` and `Catalyst 2 finalized` once.
2. The run completes with two or more MPI ranks and no Blueprint verification
   or missing-implementation errors.
3. `catalyst-output` contains time-indexed `.pvtu` files and their rank-local
   pieces.
4. Opening a `.pvtu` in ParaView 5.13 shows the eddy mesh and a three-component
   `velocity` point array.

For an adaptor-only diagnostic, temporarily set
`CATALYST_IMPLEMENTATION_PATHS` to the Catalyst `stub` implementation supplied
by the installation. That validates API calls but will not execute the Python
pipeline or write visualization files.

## Runtime controls

- `NEK_CATALYST_SCRIPT`: pipeline path; defaults to `catalyst_pipeline.py`.
- `NEK_CATALYST_OUTPUT`: output directory used by the supplied pipeline.
- `NEK_CATALYST_FREQUENCY`: extraction interval; defaults to 10 steps.
- `CATALYST_IMPLEMENTATION_PATHS`: ParaView Catalyst implementation location.

The output directory must exist before launching the MPI job.
