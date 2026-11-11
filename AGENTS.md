# BISICLES (bisicles-uob)

BISICLES is an adaptive mesh ice sheet model. The repository is C++ built on
Chombo 3.2, with some Fortran/Python components. The main simulation executable
is `driver2d` (`code/exec2D/driver2d...ex`).

## Repository layout

- `code/` — C++/Fortran source, build system, tests, tools.
  - `code/src/` — core model implementation (`AmrIce`, velocity solvers, calving,
    thermodynamics, friction, I/O).
  - `code/exec2D/` — build/run directory for the main `driver2d` executable.
  - `code/lib/` — `libBisicles2d` and Chombo helper libraries.
  - `code/test/` — small unit/integration tests.
  - `code/regression/` — end-to-end regression problems with `inputs.regression`.
  - `code/filetools/` — pre/post-processing utilities (`amrtocf`, `diagnostics`, ...).
  - `code/mk/` — BISICLES make configuration and host-specific Make.defs.
- `docs/` — user/developer documentation (Markdown, published via mkdocs).
- `examples/` — example setups and Python analysis scripts.
- `applications/` — large application workflows (Greenland, Antarctica, ISMIP).
- `code/regression/` results and `*.hdf5` are generated outputs; do not hand-edit.

## Key documentation

- Build/install: [docs/readme.md](docs/readme.md)
- Site-specific notes: [docs/sites.md](docs/sites.md)
- Common build/runtime issues: [docs/trouble.md](docs/trouble.md)
- Regression tests: [docs/regression.md](docs/regression.md)
- Driver options: [docs/driver.md](docs/driver.md)
- Python interface: [docs/pythoninterface.md](docs/pythoninterface.md)

## Prerequisites

- GNU/Linux, GNU-compatible `make`, `csh`, C++ and Fortran compilers.
- MPI environment for parallel builds (`mpiCC`, `mpicc`).
- Chombo 3.2 source tree expected beside this repository
  (`$BISICLES_HOME/Chombo`), or set `CHOMBO_HOME`.
- HDF5; optional FFTW, netCDF, PETSc, and Python (strongly recommended).
- Optional dependencies require matching `Make.defs` settings.

## Build system

BISICLES uses Chombo-style GNU make. Two host-specific definitions are needed:

1. Chombo: `$CHOMBO_HOME/mk/Make.defs.local` (from `docs/Make.defs.local.*`).
2. BISICLES: `code/mk/Make.defs.$(uname -n)`, or `code/mk/Make.defs.none`.
   The repository tracks several `Make.defs.ubuntu_*` examples and a template at
   `code/mk/Make.defs.template`.

If the host file is missing, copy a template/example:
```
cd code/mk
cp Make.defs.template Make.defs.$(uname -n)
```

Common builds (run from repository root):
```
cd code
make all OPT=TRUE MPI=TRUE
```
Smaller scoped builds:
```
cd code/exec2D && make all OPT=TRUE MPI=TRUE
cd code/test     && make all OPT=TRUE MPI=TRUE
cd code/filetools && make all OPT=TRUE MPI=TRUE
```
Optimization/MPI flags are selected by `OPT=TRUE` and `MPI=TRUE`. Without them
you get unoptimized, serial executables. PETSc builds use `USE_PETSC=TRUE`.

Clean:
```
cd code
make clean OPT=TRUE MPI=TRUE
make realclean OPT=TRUE MPI=TRUE   # also removes generated hdf5 outputs
```

The current environment provides `code/mk/Make.defs.it038128` (symlink to
`Make.defs.ubuntu_24.4`) and a prebuilt sibling Chombo tree at
`/mnt/raid_10_0/stephen/Development/Chombo`.

## Running

The executable takes an inputs file:
```
cd code/exec2D
./driver2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex ../path/to/inputs.problem
```
or under MPI:
```
mpirun -np 4 ./driver2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex inputs.problem
```

Many examples and regression tests are driven by Python input-generation scripts
and Python extensions (`PythonIBC`, `PythonBasalFriction`, `PythonVelocitySolver`).
Run these from their own directory so Python modules can be imported, and ensure
`HAVE_PYTHON` was enabled at build time.

Default outputs are HDF5 (`*.hdf5`) plus optional CF-format files and `pout.*`.

## Tests

Unit/integration tests live in `code/test`. End-to-end regression problems live
under `code/regression/*` and are run with their local `inputs.regression`
(or problem-specific inputs):
```
cd code/regression/conservation
mpirun ../../exec2D/driver2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex inputs.regression
```
Regression suites exercised in CI: `twistyStream`, `conservation`, `plot-cf`,
and `ASE-control`. `ASE-control` also has sequence/resume/continuation input
files that must be run in order.

Regeneration of large `*.hdf5` outputs is expected during development; keep
generated artifacts out of commits. When changing model code, rebuild the
`libBisicles2d` library and rerun the regression test covering the affected area.

## Coding conventions

- Header files use `.H`; implementations use `.cpp`; Chombo macros and Chombo
  includes dominate the codebase.
- Follow existing Chombo style: class members, `ParmParse` input parameters,
  and `CH_assert`/message macros are standard.
- C++ is compiled as C++11; generated filenames and object directories encode
  the configuration (compiler, optimization, MPI, PETSc).
- Fortran components use Chombo Fortran (`.ChF`) plus
  `code/src/ColumnThermodynamics.f90`.
- Python extension modules are imported by name from the run directory; avoid
  renaming without updating input files.
- Do not edit generated files in `d/`, `f/`, `o/`, `p/`, compiled `*.ex`, or
  host-specific symlinks unless intentionally changing the environment.

## Change workflow for agents

1. Read the relevant docs before editing (`docs/readme.md`, `docs/trouble.md`).
2. Make focused changes in `code/src` or relevant supporting code.
3. Rebuild with `cd code && make all OPT=TRUE MPI=TRUE`.
4. Run the closest regression test and inspect `pout.*` / output success.
5. Do not revert or bulk-clean unrelated user changes; this repository commonly
   has many generated and in-progress files in the working tree.
6. Leave generated outputs uncommitted unless the task explicitly requires them.
