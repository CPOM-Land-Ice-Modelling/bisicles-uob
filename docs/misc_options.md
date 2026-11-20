# Miscellaneous run-time options

This page was rated by opencode/Deep Seek Flash 4.1 Free. YMMV

This page collects run-time options that are not covered by the more
focused pages: [input/output](io.md), [geometry](geometry.md),
[stresses](stress.md), [thermodynamics](thermodynamics.md), [ice
velocity](velocity.md), [surface fluxes](surfaceflux.md),
[LevelData](leveldatainterface.md), [Python](pythoninterface.md), and
[calving](calvingfront.md). It is mainly intended as a reference for
less commonly used, advanced, or experimental settings.

## How options are read

Almost every BISICLES option is read from the run-time configuration
file (usually called `inputs.*`) using a Chombo `ParmParse` object. A
`ParmParse` object is created with a *prefix*, and then queries keys
within that prefix. For example

    ParmParse pp("foo");
    Real bar = 0.0;             // default value
    pp.query("bar", bar);       // reads foo.bar if it is present

reads the entry

    foo.bar = 1.0

from the configuration file and stores `1.0` in `bar`. If `foo.bar` is
not present, `query` leaves `bar` at its default. By contrast, `get`
raises an error when the entry is missing.

The most important conventions are

-   A prefix and key are joined with a dot, so `ParmParse pp("amr")` plus
    `pp.query("cfl", cfl)` reads `amr.cfl`.
-   Comments begin with `#`. Later entries override earlier entries.
-   Values can be scalars, strings, or space-separated lists. Lists are
    read with `queryarr`/`getarr`.
-   Boolean values accept `true`/`false` (and `1`/`0`).
-   Options can also be supplied as extra command-line arguments after
    the inputs file, e.g.
    `driver2d...ex inputs.run amr.maxLevel=4`. Chombo parses arguments
    of the form `prefix.key=value` and they override the file.

The tables below list the option name (including its prefix), the
default value in the current code, and a short description. Defaults
are only a guide: some options are only used when another option enables
them.

## Physical constants

These are read by the `constants` prefix. They are used to convert
between model units and SI units, and to compute buoyancy and driving
stress.

| Option | Default | Meaning |
| --- | --- | --- |
| `constants.ice_density` | `910.0` | Density of ice, kg m^-3. |
| `constants.sea_water_density` | `1028.0` | Density of sea water, kg m^-3. Used for flotation and sea level. |
| `constants.mantle_density` | `3313.0` | Mantle density, kg m^-3. Used by GIA and some surface fluxes. |
| `constants.gravity` | `9.81` | Gravitational acceleration, m s^-2. |
| `constants.seconds_per_unit_time` | seconds in a tropical year | Number of seconds in one model time unit. Set to `1.0` to run in seconds (see [input/output](io.md#run-time-configuration-inputs-files)). |

Note that the Bueler GIA topography flux reads the sea water density as
`constants.seaWater_density` (camel case), unlike the rest of the code,
which uses `constants.sea_water_density`.

## Main program controls

These options are read by the `main` prefix, or by the [driver](driver.md)
directly.

| Option | Default | Meaning |
| --- | --- | --- |
| `main.domain_size` | required | Size of the physical domain, one value per horizontal dimension, in metres. |
| `main.maxTime` | required | Simulation end time, in model time units. |
| `main.maxStep` | required | Maximum number of time steps to take. |
| `main.newTimeStep` | `false` | Use the alternative fixed-step predictor-corrector driver (`AmrIce::timeStep2`) instead of the standard adaptive time stepping. Requires `amr.fixed_dt` to be positive. |
| `main.damage_model` | `false` | Enable the continuum damage model. This wraps the selected constitutive relation in a damage relation and replaces the calving model with one that can calve on damage. See [damage model](#continuum-damage-model). |
| `main.melange_model` | `false` | Enable the melange transport model. See [melange model](#melange-model). |
| `main.poutBaseName` | `pout` | Base name for the per-process log files (`pout.0`, `pout.1`, ...). Also described in [input/output](io.md#log-output-pout-files). |

`main.newTimeStep` is experimental. It requires a fixed time step and is
normally only used for particular coupled or idealised configurations.

## Time stepping and ice evolution

These options are read from the `amr` prefix.

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.evolve_thickness` | `true` | Evolve the ice thickness. Setting this to `false` holds the geometry fixed. |
| `amr.evolve_velocity` | `true` | Recompute the velocity field as the simulation advances. |
| `amr.velocity_solve_interval` | `1` | Recompute the velocity every n time steps. (Some older documentation calls this `amr.velocity_solver_interval`; the code reads `amr.velocity_solve_interval`.) |
| `amr.time_step_n_iter` | `1` | Number of corrector iterations in the alternative time stepper: `0` is forward Euler, `1` is Heun, larger values approach Crank-Nicolson. Only used with `main.newTimeStep`. |
| `amr.velocity_exit` | `1.0e+6` | If the maximum velocity exceeds this value (in model velocity units), write an error plot file and abort. |
| `amr.max_dt` | effectively unlimited | Upper limit on the adaptive time step. |
| `amr.fixed_dt` | `0.0` | If positive, use this fixed time step instead of the CFL-limited step. Also required by `main.newTimeStep`. |
| `amr.diffusion_treatment` | `none` | Thickness diffusion treatment: `none` or `implicit`. The `explicit` option is not currently supported. |
| `amr.additional_diffusivity` | `0.0` | Constant added to the thickness diffusivity. |
| `amr.finest_timestep_level` | `-1` | If non-negative, evolve thickness and internal energy only on levels up to this level, then interpolate the result to finer levels. `-1` evolves all levels. |
| `amr.temporal_accuracy` | `1` | Time integration order selector. Currently `1` is the only working option. |
| `amr.additional_velocity` | `false` | Include an additional (shallow-ice) velocity contribution in the L1L2 constitutive relation. |
| `amr.domainLoIndex` | `0` (per direction) | Lower index of the base-level domain. Useful for domains whose index space does not start at zero. |
| `amr.limit_margin_velocity` | `false` | Limit excessive velocities at ice margins after the velocity solve. |
| `amr.reset_floating_friction` | `true` | Set the basal friction coefficient to zero on floating ice before the velocity solve, using a thickness-above-flotation interpolation at the grounding line. |

### Thickness advection

When `amr.temporal_accuracy` is less than `3`, thickness is advected with
the Chombo `PatchGodunov` scheme. Its options are read from the
`amr.thickness_advection` prefix.

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.thickness_advection.normal_pred_order` | `2` | Order of the normal predictor. |
| `amr.thickness_advection.fourth_order_slopes` | `true` | Use fourth-order limited slopes rather than second-order slopes. |
| `amr.thickness_advection.artificial_viscosity` | `0.0` | Artificial viscosity coefficient. Zero disables artificial viscosity. |

### Holding ice fixed and prescribing mass balance

These options are used to keep parts of the ice sheet from evolving by
modifying the basal mass balance, or to hold the grounding line.

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.grounded_ice_stable` | `false` | Modify the basal mass balance so that grounded ice does not change thickness. |
| `amr.floating_ice_stable` | `false` | Modify the basal mass balance so that floating ice does not change thickness. |
| `amr.grounded_ice_basal_flux_is_dhdt` | `false` | Set the grounded-ice basal flux to balance surface flux and divergence, i.e. impose `dh/dt = 0`. |
| `amr.floating_ice_basal_flux_is_dhdt` | `false` | Set the floating-ice basal flux to balance surface flux and divergence, i.e. impose `dh/dt = 0`. |
| `amr.floating_ice_basal_flux_is_min_dhdt` | `false` | Modify the floating-ice basal flux so that floating ice cannot thin, but is allowed to thicken. |
| `amr.grounding_line_stable` | `false` | Modify the surface mass balance near the grounding line to keep grounded ice just above flotation, discouraging grounding-line migration. |
| `amr.mask_stable_sources` | `false` | Apply the stable-source modifications only where a supplied mask is zero. |
| `amr.stable_sources_mask_file` | empty | hdf5 file containing the stable-sources mask. |
| `amr.stable_sources_mask_name` | `stable_sources_mask` | Variable name of the stable-sources mask in that file. |
| `amr.grounding_line_proximity_scale` | `1.0e+4` | Length scale (metres) for the exponentially decaying grounding-line proximity field used by `groundingLineLocalizedFlux` and the Python interface. |
| `amr.grounding_line_proximity_calc_type` | `0` | Selects the coefficient used when solving for the grounding-line proximity field. `0` is the historical behaviour; values greater than zero use a squared coefficient. |

### Ice fraction and source masking

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.mask_sources` | `true` | Multiply surface and basal thickness sources by the ice area fraction. |
| `amr.frac_sources` | `true` | Synonym for `amr.mask_sources`. |
| `amr.evolve_ice_frac` | derived | Advect the ice area fraction. By default this is true unless thickness evolution is inhibited by one of the stable-ice options. |
| `amr.evolve_ice_frac_2` | `false` | Use the alternative ice-fraction advection scheme. Setting this true also enables `amr.evolve_ice_frac`. |
| `amr.floating_check_ocean_connected` | `false` | When holding floating ice stable, only apply the adjustment to floating ice that is connected to the open ocean. |

### Removing remote ice

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.eliminate_remote_ice` | `false` | Identify floating ice that is not connected to grounded ice and remove it. Very thin isolated grounded cells are removed at the same time. |
| `amr.eliminate_remote_ice_max_iter` | `10` | Maximum number of iterations in the remote-ice elimination. |
| `amr.eliminate_remote_ice_tol` | `1.0` | Thickness tolerance used when eliminating remote ice. |
| `amr.eliminate_remote_ice_after_regrid` | `false` | Also run remote-ice elimination immediately after regridding. |

### Initial velocity guess

By default BISICLES can construct an inexpensive initial guess before
the first nonlinear velocity solve. These options control that guess.

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.do_initial_velocity_solve` | `true` | Perform the velocity solve when the model is initialized. |
| `amr.do_initial_velocity_guess` | `false` | Construct a semi-informed initial guess before the first solve. |
| `amr.initial_velocity_guess_type` | `0` | `0`: sliding-law guess; `1`: constant-viscosity solve; `2`: a supplied function. |
| `amr.initial_velocity_guess_const_mu` | `1.0e+7` | Constant viscosity used when `amr.initial_velocity_guess_type = 1`. |
| `amr.initial_velocity_guess_solver_type` | `1` | Solver for the constant-viscosity guess: `0` Picard, `1` JFNK. |
| `amr.initial_velocity_guess_const_vel` | `0 0` | Constant initial velocity used together with the constant-viscosity guess and a nonlinear basal friction law. |
| `amr.initial_velocity_function_type` | `constant` | When `amr.initial_velocity_guess_type = 2`, the function type. `flowline` reads a piecewise-linear flowline profile. |
| `amr.initial_velocity_function_flowline_dx` | required for flowline | Spacing of the flowline data. |
| `amr.initial_velocity_function_flowline_file` | required for flowline | File containing the flowline data. |
| `amr.initial_velocity_function_flowline_set` | required for flowline | Name of the data set in that file. |

## Thermodynamics smoothing

The `Thermodynamics` prefix controls optional horizontal smoothing of the
internal energy and till water depth. Smoothing is applied with a
Helmholtz solve every `smooth_interval` steps.

| Option | Default | Meaning |
| --- | --- | --- |
| `Thermodynamics.smooth_interval` | `-1` (disabled) | Smooth internal energy every n steps. |
| `Thermodynamics.length_scale` | `sqrt(IceConductivity()/IceHeatCapacity())` | Smoothing length scale, in metres. |
| `Thermodynamics.till_water_smooth_interval` | `-1` (disabled) | Smooth the till water depth every n steps. |
| `Thermodynamics.till_water_length_scale` | `sqrt(IceConductivity()/IceHeatCapacity())` | Smoothing length scale for till water. |

## Mesh refinement and tagging

### Grid generation

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.tagCap` | `100` | Maximum level index considered when tagging cells for refinement. |
| `amr.max_box_size` | `64` | Maximum size of a patch on the finest level. |
| `amr.max_base_grid_size` | same as `max_box_size` | Maximum size of a base-level patch. |
| `amr.blockFactor` | `-1` | Chombo block factor for grid generation. `-1` uses the Chombo default. |
| `amr.fill_ratio` | `-1` | Chombo fill ratio for grid generation. `-1` uses the Chombo default. |
| `amr.nestingRadius` | `1` (`3` with PETSc) | Required nesting radius between levels. PETSc builds force at least `3`. |
| `amr.interpolate_zb` | `true` | Interpolate the bedrock topography when refining, rather than re-reading it. |
| `amr.regrid_thickness_interpolation_method` | `0` | Thickness interpolation on regrid: `0` interpolate thickness, `1` interpolate upper and lower surfaces and use fine bedrock on grounded cells, `2` take the maximum thickness. |
| `amr.tags_grow` | `1` | Grow the tagged region by this many cells in every direction. |
| `amr.tags_grow_dir` | `0 0` | Grow tags further in each direction; the actual growth in direction `d` is `max(0, tags_grow_dir[d] - tags_grow)`. |
| `amr.tagSubsetBoxesFile` | empty | File containing per-level boxes; tagging is restricted to these boxes. |
| `amr.gridsFile` | empty | Use predefined grids read from this file instead of generating them at start-up. |

The predefined grids file used by `amr.gridsFile` has the format

    numberOfLevels
    numGridsOnLevel1
    lo hi
    lo hi
    ...
    numGridsOnLevel2
    lo hi
    ...

where `lo hi` are Chombo index-space boxes. Level 0 is generated
automatically from the domain and `amr.max_base_grid_size`. The file
must contain at most `amr.maxLevel + 1` levels.

### Tagging criteria

If no tagging criterion is enabled, BISICLES falls back to
`amr.tag_on_grad_velocity = true`. Otherwise the criteria below are
combined (a cell is tagged if any enabled criterion tags it).

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.tag_on_grad_velocity` | `false` | Tag where the undivided gradient of the horizontal velocity exceeds `amr.tagging_val`. |
| `amr.tagging_val` | `0.1` | Threshold for `amr.tag_on_grad_velocity`. |
| `amr.tag_on_laplacian_velocity` | `false` | Tag where the undivided Laplacian of the velocity exceeds `amr.lap_vel_tagging_val`. |
| `amr.tag_on_grounded_laplacian_velocity` | `false` | As above, but only on grounded ice. |
| `amr.lap_vel_tagging_val` | `1.0` (required when enabled) | Threshold for either Laplacian criterion. |
| `amr.lap_vel_tagging_max_basal_friction_coef` | `1.2e+300` | Only apply the Laplacian criterion where the basal friction coefficient is below this value. |
| `amr.tag_on_strain_rate_invariant` | `false` | Tag where the second invariant of the strain-rate tensor, scaled by `dx^2`, exceeds `amr.strain_rate_invariant_tagging_val`. |
| `amr.strain_rate_invariant_tagging_val` | `0.1` (required when enabled) | Threshold for the strain-rate-invariant criterion. |
| `amr.tag_on_velocity_rhs` | `false` | Tag where the magnitude of the velocity right-hand side multiplied by `dx` exceeds `amr.velocity_rhs_tagging_val`. |
| `amr.velocity_rhs_tagging_val` | `1.0` (required when enabled) | Threshold for the velocity-RHS criterion. |
| `amr.tag_grounding_line` | `false` | Tag cells adjacent to the grounding line where the velocity exceeds a threshold. |
| `amr.grounding_line_tagging_min_vel` | `200.0` (required when enabled) | Minimum speed for grounding-line tagging. |
| `amr.grounding_line_tagging_max_basal_friction_coef` | `1.2e+300` | Only apply grounding-line tagging where the basal friction coefficient is below this value. |
| `amr.tag_vel_dx` | `false` | Tag where `|u| dx` exceeds `amr.vel_dx_tagging_val`, up to a specified finest level. This mimics the meshes produced by observed-velocity-based controllers. |
| `amr.vel_dx_tagging_val` | `0.0` (required when enabled) | Threshold for `amr.tag_vel_dx`. |
| `amr.vel_dx_finest_level_grounded` | `-1` | Refine on `|u| dx` below this level on grounded ice. |
| `amr.vel_dx_finest_level_floating` | same as grounded | As above, on floating ice. |
| `amr.tag_thin_cavity` | `false` | Tag floating ice thinner than `amr.tag_thin_cavity_thickness`. |
| `amr.tag_thin_cavity_thickness` | tiny thickness | Thickness threshold for `amr.tag_thin_cavity`. |
| `amr.tag_ice_margin` | `false` | Tag where the ice area fraction is close to `0.5`, on levels below `amr.margin_finest_level`. |
| `amr.margin_finest_level` | `maxLevel + 1` | Highest level on which ice-margin tagging is applied. |
| `amr.tag_all_ice` | `false` | Tag every cell containing ice. |
| `amr.tag_all_ice_on_level_0` | `false` | Tag every cell containing ice on level 0 only. |
| `amr.tag_entire_domain` | `false` | Tag the whole domain (useful for uniform refinement tests). |
| `amr.tag_on_div_H_grad_vel` | `false` | Tag where `div(H grad u)` exceeds `amr.div_H_grad_vel_tagging_val`. |
| `amr.div_H_grad_vel_tagging_val` | `1.0` (required when enabled) | Threshold for the `div(H grad u)` criterion. |
| `amr.tag_in_regions` | `false` | Tag inside rectangular regions read from `amr.tag_regions_file`. |

#### Tagging in regions

When `amr.tag_in_regions = true`, the file named by
`amr.tag_regions_file` must contain

    numberOfRegions
    xlo ylo xhi yhi resolution
    xlo ylo xhi yhi resolution
    ...

with one line per region. Coordinates are in physical space. A region is
tagged on a level while its requested resolution is finer than that
level's mesh spacing.

#### Tagging inside a subset of boxes

`amr.tagSubsetBoxesFile` names a file that masks the tagged cells: after
the usual criteria have run on a level, any tags outside the boxes listed
for that level are discarded. The file lists levels in order, each as
either `level N` or `domain box`, followed by the number of boxes and
then the boxes in Chombo index-space syntax:

    level 0
    1
    0 0 99 99
    level 1
    2
    0 0 49 49
    50 50 99 99

#### Python tagging

When BISICLES is built with Python support,

    amr.tag_python = true
    amr.tag_python_module = mymodule
    amr.tag_python_function = mytag

calls `mytag(x, y, dx, thickness, topography)` in `mymodule.py` for each
cell. A positive return value tags the cell. This is a convenient way to
refine around arbitrary geometric features without a LevelData file.

## Velocity solver controls

The choice of velocity solver is described in [ice
velocity](velocity.md). The options below tune the solvers themselves.

### `amr` solver controls

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.max_solver_iterations` | `-1` | Override the solver's own maximum iteration count. `-1` leaves the solver default. |
| `amr.velocity_solver_tolerance` | `1.0e-10` | Override the solver tolerance. |
| `amr.do_initial_velocity_solve` | `true` | See [initial velocity guess](#initial-velocity-guess). |
| `amr.do_initial_velocity_guess` | `false` | See [initial velocity guess](#initial-velocity-guess). |

### JFNK solver

The JFNK solver reads the `JFNKSolver` prefix. The most commonly used
options are demonstrated in [ice velocity](velocity.md) and
[Pine Island Glacier](pineisland.md); the remainder are listed here for
reference.

| Option | Default | Meaning |
| --- | --- | --- |
| `JFNKSolver.maxIter` | `15` | Maximum number of Newton iterations. |
| `JFNKSolver.maxIter_init` | same as `maxIter` | Maximum number of Newton iterations for the initial solve. |
| `JFNKSolver.absTol` | `1.0e-10` | Absolute convergence tolerance. |
| `JFNKSolver.relTol` | `1.0e-10` | Relative convergence tolerance. |
| `JFNKSolver.normType` | `0` | Norm used for the convergence test. |
| `JFNKSolver.verbosity` | `5` | Solver logging level. |
| `JFNKSolver.solverType` | `0` (Relax) | Linear solver: `0` Relax, `1` BiCGStab, `2` GMRES (disabled), `3` CG, `4` PETSc. |
| `JFNKSolver.bottom_solver_type` | `0` (BiCGStab) | Multigrid bottom solver: `0` BiCGStab, `1` PETSc, `2` Relax. |
| `JFNKSolver.RelaxRelTol` | `1.0e-3` | Relative tolerance for the Relax linear solver. `RelaxTol` is a synonym. |
| `JFNKSolver.maxRelaxIter` | `10` | Maximum iterations for the Relax solver. |
| `JFNKSolver.RelaxHang` | `0.25` | Hang parameter for the Relax solver. |
| `JFNKSolver.BiCGStabRelTol` | `1.0e-3` | Relative tolerance for BiCGStab. |
| `JFNKSolver.maxBiCGStabIter` | `10` | Maximum iterations for BiCGStab. |
| `JFNKSolver.CGRelTol` | `1.0e-3` | Relative tolerance for CG. |
| `JFNKSolver.maxCGIter` | `10` | Maximum iterations for CG. |
| `JFNKSolver.numMGSmooth` | `8` | Number of smoothing sweeps in multigrid. |
| `JFNKSolver.numMGIter` | `1` | Number of multigrid V-cycles per Krylov iteration. |
| `JFNKSolver.vtopSafety` | `0.5` | Safety factor for the viscous-tensor operator. |
| `JFNKSolver.vtopRelaxTol` | `1.0e-2` | Relaxation tolerance for the viscous-tensor operator. |
| `JFNKSolver.vtopRelaxMinIter` | `8` | Minimum number of viscous-tensor smoothing sweeps. |
| `JFNKSolver.mg_solver_depth` | `-1` | Maximum multigrid depth. `-1` uses all available levels. |
| `JFNKSolver.mgCoefficientAverageType` | `0` (arithmetic) | Chombo coefficient averaging type for multigrid. |
| `JFNKSolver.mgProlongType` | linear interpolation | Chombo prolongation type for multigrid. |
| `JFNKSolver.lazyGSRB` | Chombo default | Enable or disable lazy Gauss-Seidel red-black smoothing (Chombo trunk builds only). |
| `JFNKSolver.h` | `0.025` | Step size for the finite-difference Jacobian approximation. |
| `JFNKSolver.hAdaptive` | `false` | Choose the finite-difference step adaptively. |
| `JFNKSolver.err` | `1.0e-6` | Error parameter used by adaptive `h`. |
| `JFNKSolver.umin` | `1.0e-6` | Minimum speed used by adaptive `h`. |
| `JFNKSolver.switchRate` | `1.5` | Rate at which the solver switches from Picard to Newton iterations. |
| `JFNKSolver.minPicardIterations` | `0` | Minimum number of Picard iterations before switching to Newton. |
| `JFNKSolver.minStepFactor` | `1.0` | Smallest acceptable Newton step reduction factor. |
| `JFNKSolver.muMin` | `0.0` | Minimum effective viscosity. |
| `JFNKSolver.muMax` | `1.2e+300` | Maximum effective viscosity. |
| `JFNKSolver.scale` | model time units per tropical year | Solve internally for `scale * u`, where `scale` converts the velocity to the solver's preferred units; the returned velocity is unscaled. This is `1` when running in years. |
| `JFNKSolver.writeResiduals` | `false` | Write residual information at each iteration, useful for diagnosing failed solves. |
| `JFNKSolver.artificialDragCoef` | `0.0` | Coefficient of an artificial drag applied everywhere, intended to improve conditioning. |
| `JFNKSolver.artificialDragPower` | `8.0` | Exponent of the artificial drag. |
| `JFNKSolver.eliminateFastIce` | `false` | Remove fast ice at the end of each Newton iteration. |
| `JFNKSolver.eliminateFastIceSpeed` | `1.0e+5` | Speed above which ice is considered "fast". |
| `JFNKSolver.eliminateFastIceEdgeOnly` | `false` | Only remove fast ice adjacent to ice-free regions. |
| `JFNKSolver.eliminateRemoteIceTol` | `1.0` | Thickness tolerance for remote-ice removal during the solve. Inherits `amr.eliminate_remote_ice_tol` if set. |
| `JFNKSolver.eliminateRemoteIceMaxIter` | `10` | Maximum iterations for remote-ice removal during the solve. Inherits `amr.eliminate_remote_ice_max_iter` if set. |

### PETSc nonlinear solver

These options are read by the `petsc` prefix when using the experimental
PETSc SNES velocity solver (`amr.velocity_solver_type = 2`). General
PETSc linear-solver options are set in a `.petsrc` file, as described in
[ice velocity](velocity.md#option-petsc-linear-solvers).

| Option | Default | Meaning |
| --- | --- | --- |
| `petsc.maxIter` | `20` | Maximum number of nonlinear iterations. |
| `petsc.absNLTol` | `1.0e-30` | Absolute nonlinear tolerance. |
| `petsc.relNLTol` | `1.0e-6` | Relative nonlinear tolerance. |
| `petsc.minPicardIterations` | `3` | Minimum number of Picard iterations before Newton. |
| `petsc.plotResidual` | `false` | Write the residual history for plotting. |

### Picard solver initial guess

When the constant-viscosity initial guess uses the Picard solver
(`amr.initial_velocity_guess_type = 1` and
`amr.initial_velocity_guess_solver_type = 0`), these options apply.

| Option | Default | Meaning |
| --- | --- | --- |
| `picardSolver.linearsolver_tolerance` | `1.0e-4` | Tolerance for the linear solve. |
| `picardSolver.max_picard_iterations` | `1` | Maximum number of Picard iterations. |

## Basal friction tuning

The basal friction law is selected with `main.basalFrictionRelation` and
described in [stresses](stress.md#basal-friction-relation). The following
options tune the built-in laws.

### Power law

These are read from the `BasalFrictionPowerLaw` prefix.

| Option | Default | Meaning |
| --- | --- | --- |
| `BasalFrictionPowerLaw.fastSlidingSpeed` | negative (disabled) | Speed above which an additional Joughin-style sliding term is applied. Negative disables the term. |
| `BasalFrictionPowerLaw.highThickness` | negative (effectively infinite) | Effective thickness used by that fast-sliding term. |
| `BasalFrictionPowerLaw.oceanConnectivityCoef` | `1.0` | Weight between thickness above flotation and full thickness in the effective-pressure factor, from `0` to `1`. |

The `m` and `includeEffectivePressure` options are described in
[stresses](stress.md#power-law-basal-friction-relation).

A separate option, `InitialThickness.criticalHeightAboveFlotation`
(default `0.0`, disabled), reduces the basal friction coefficient near
the grounding line. When it is positive, cells whose thickness above
flotation is below that value have their friction coefficient scaled
down, using the initial thickness supplied by the `InitialThickness`
surface flux. This requires a valid `InitialThickness` flux, e.g.

    InitialThickness.type = LevelData
    InitialThickness.fileFormat = initial%4d.2d.hdf5
    InitialThickness.n = 1
    InitialThickness.name = thk

### Pressure-limited law

These are read from the `BasalFrictionPressureLimitedLaw` prefix. The
`model` and `coefficient` options are described in
[stresses](stress.md#pressure-limited-basal-friction-relation).

| Option | Default | Meaning |
| --- | --- | --- |
| `BasalFrictionPressureLimitedLaw.power` | `1.0` | Exponent of the pressure-limited law. |
| `BasalFrictionPressureLimitedLaw.maxTillWaterDepth` | `-1.0` (no limit) | Maximum till water depth used when forming the effective pressure. |
| `BasalFrictionPressureLimitedLaw.tillPressureFactor` | `0.97` | Till pressure factor (the `alpha` of Van Pelt and Oerlemans 2012). |

## Additional plot-file contents

The main plot-file options are described in
[input/output](io.md#plot-files-plot-2dhdf5-files). The following options
are also available.

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.write_preSolve_plotfiles` | `false` | Write a plot file immediately before each velocity solve, with `preSolve` in the name. |
| `amr.write_vonmises_stress` | `false` | Write the von Mises stress. |
| `amr.write_topography_sources` | `false` | Write the bedrock/uplift source term from the topography flux. |
| `amr.write_ismip6` | `false` | Write the extra fields required by ISMIP6: surface and basal temperature, and the drag coefficient. Forces `write_base_velocities` and `write_mask` on. |
| `amr.write_map_file` | `false` | Write a sigma-mapped AMR hierarchy rather than the standard flattened hierarchy. |

Older input files and documentation use `amr.write_fluxVel`,
`amr.write_baseVel`, and `amr.write_viscousTensor`. The current code
reads only the underscore names `amr.write_flux_velocities`,
`amr.write_base_velocities`, and `amr.write_viscous_tensor`.

## Restart controls

Most restart options are described in
[input/output](io.md#checkpoints-chk-2dhdf5-files-and-restarts). A few
further options are useful when restarting.

| Option | Default | Meaning |
| --- | --- | --- |
| `amr.restart_step` | checkpoint value | Override the step number stored in the checkpoint file. Also sets the step from which checkpoint overwriting is measured. |
| `amr.restart_read_ice_frac` | `true` if the checkpoint contains it | Set to `false` to ignore the ice area fraction stored in a checkpoint, for backward compatibility with older files. |
| `amr.check_exit` | `false` | Exit immediately after writing a checkpoint (except at step zero). Useful for testing checkpoint/restart or for splitting a run into short jobs. |
| `amr.cf_restart_file` | empty | CF-format file to read CF restart data from after loading the AMR checkpoint. If empty, `CFIO.imply_restart` decides whether the default file name is used. |
| `amr.force_gia_init` | `false` | Re-initialise the topography flux (GIA model) after restarting, for example when the checkpoint was written without GIA. |

## Diagnostics and conservation

| Option | Default | Meaning |
| --- | --- | --- |
| `DomainDiagnosticData.max_err_mass_con` | `-1.0` (disabled) | If positive, abort when the domain-integrated mass-conservation error exceeds this value. |

Several `amr.report_*` options (`report_discharge`,
`report_time_interval`, `report_sum_grounded_ice`, `report_ice_area`,
`report_total_flux`, `report_calving`) are accepted but are currently
placeholders and have no effect. The domain diagnostics written to CF
plot files are controlled through the `CFIO` options described in
[input/output](io.md#cf-plot-files-plot-cf2dhdf5-files).

## Optional model components

### Continuum damage model

Setting `main.damage_model = true` enables the continuum damage model.
Damage is transported with the ice and reduces the effective viscosity.
The calving model is wrapped so that heavily damaged floating ice can
calve.

| Option | Default | Meaning |
| --- | --- | --- |
| `damage.shear_alpha` | `-1.0` (disabled) | Coefficient for damage production by shear. Values greater than zero enable the shear source. |
| `damage.shear_critical` | `1.0e-10` | Shear threshold below which no shear damage is produced. |
| `damage.water_source.type` | `zeroFlux` | A [surface flux](surfaceflux.md) providing water to the damage model, e.g. from surface melt. Uses the usual `damage.water_source.*` options. |
| `damage.critical_damage` | `1.0` | Damage threshold used by the damage calving model. |
| `damage.critical_thickness` | `0.0` | Additional thickness threshold used by the damage calving model. |

### Melange model

Setting `main.melange_model = true` enables a melange thickness field
that receives calved ice and can spread diffusively.

| Option | Default | Meaning |
| --- | --- | --- |
| `melange_model.external_source.type` | `zeroFlux` | A [surface flux](surfaceflux.md) providing the external melange source (or sink), e.g. from an ocean or iceberg model. Uses the usual `melange_model.external_source.*` options. |
| `melange_model.diffusion_factor` | `0.0` | Diffusivity factor for melange transport. Zero disables diffusion. |
| `melange_model.increment_method` | `1` | Source method: `0` increments melange by ice thickness through explicit calls; `1` uses the thickness removed by calving during the time step. |

## Additional surface flux types

These flux types are documented here because they are not covered by
[surface fluxes](surfaceflux.md). They can be used anywhere a surface
flux is expected, e.g. `surfaceFlux.type`, `basalFlux.type`, or
`topographyFlux.type`.

### targetThicknessFlux

Applies a flux that relaxes the ice thickness toward a target field on a
given timescale. It is useful for idealised spin-up experiments.

    surfaceFlux.type = targetThicknessFlux
    surfaceFlux.timescale = 100.0
    surfaceFlux.target.type = constantFlux
    surfaceFlux.target.flux_value = 500.0

`timescale` is required. The `target` subsection is itself a surface
flux.

### FloatingDivUHLocalizedFlux

Confines a supplied flux to the floating ice, and scales it by the local
value of `div(uH)`. The supplied flux is given in a `flux` subsection and
the mesh spacing must be provided.

    basalFlux.type = FloatingDivUHLocalizedFlux
    basalFlux.mesh_spacing = 1000.0
    basalFlux.flux.type = constantFlux
    basalFlux.flux.flux_value = -10.0

### ISMIP6 ocean forcing

`ISMIP6OceanForcing` computes a basal melt rate from ISMIP6-style thermal
forcing and ocean temperature anomaly data. In the current code the type
name is spelled `IMSIP6OceanForcing` (without the second `I`).

    basalFlux.type = IMSIP6OceanForcing
    basalFlux.file_format = thermal_forcing_%4d.nc
    basalFlux.start_year = 1995
    basalFlux.end_year = 2100

| Option | Default | Meaning |
| --- | --- | --- |
| `file_format` | required | `printf`-style format used to build the yearly thermal-forcing file names. |
| `start_year` | `1995` | First year of forcing. |
| `end_year` | `2100` | Last year of forcing. |
| `name` | `thermal_forcing_0000` | Variable name for the thermal forcing. |
| `basin_file` | empty | Optional file containing a basin-number field. |
| `basin_var_name` | `basinNumber` | Variable name for the basin numbers. |
| `n_basin` | `16` | Number of basins (the IMBIE2 default). |
| `deltaT_specified` | `true` | Whether an ocean temperature anomaly is supplied. |
| `deltaT_file` | required if `deltaT_specified` | File containing the temperature anomaly. |
| `deltaT_var_name` | `deltaT_basin` | Variable name for the anomaly. |
| `anomaly` | `false` | Treat the supplied field as an anomaly. |
| `n_layer` | `30` | Number of vertical forcing layers (ISMIP6 default). |
| `dz` | `60.0` | Layer spacing in metres (ISMIP6 default). |
| `basin_mean_min_thickness` | `0.0` | Minimum thickness used when forming basin means. |
| `local` | `true` | Use the local (rather than non-local) melt-rate formula. |
| `local_slope` | `false` | Use the local-slope variant. |
| `gamma0` | formula default | Heat-transfer coefficient. Defaults depend on `local`. |
| `factor` | formula default | Directly set the melt-rate factor, overriding `gamma0`. Should be negative. |
| `tf_is_source` | `false` | Interpret the supplied thermal forcing as a source rather than as ISMIP6 thermal forcing. |

All of these options live under the flux prefix, e.g.
`basalFlux.file_format`.

### Bueler GIA topography flux

When built with FFTW support, `topographyFlux.type = buelerGIA` enables
the Bueler GIA model. It requires a rectangular, uniformly spaced domain
covering the whole ice sheet.

    topographyFlux.type = buelerGIA
    topographyFlux.pad = 2
    topographyFlux.nlayers = 1
    topographyFlux.visc = 1.0e+21
    topographyFlux.flex = 1.0e+5

| Option | Default | Meaning |
| --- | --- | --- |
| `pad` | required | Number of grid cells added around the domain for the FFT. |
| `inside_box` | `true` | Treat the problem domain as a box inside a larger GIA domain. |
| `box_lox`, `box_hix`, `box_loy`, `box_hiy` | large values | Extent of the GIA box in index space. |
| `nlayers` | required | Number of mantle layers; currently only `1` and `2` are implemented. |
| `visc` | required | Mantle viscosity (one value for one layer, two for two layers). |
| `thk` | required for two layers | Thickness of the upper layer. |
| `flex` | required | Flexural rigidity. |
| `includeElas` | `false` | Include elastic lithosphere effects. |
| `lame1`, `lame2` | `3.42666e+10`, `2.66e+10` | Lamé parameters used when `includeElas` is true. |
| `oceanLoad` | `false` | Include ocean loading. |
| `ELRA` | `false` | Use the ELRA relaxation model. |
| `ELRAtau` | `1.0e+3` | ELRA relaxation time, in years. |
| `dt` | `0` | GIA time step. |
| `init` | `false` | Initialise the GIA model from a supplied uplift field. |
| `init.file`, `init.name` | required if `init` | File and variable name for the initial uplift field. |

The Bueler GIA flux uses `constants.ice_density`, `constants.gravity`,
`constants.mantle_density`, and `constants.seaWater_density` (note the
camel-case spelling described above).

## Experimental and developer options

These options are accepted by the code but are either experimental or
intended for development and debugging.

| Option | Default | Meaning |
| --- | --- | --- |
| `bodges.frac_rhs_bodge` | `false` | Divide the velocity right-hand side by the ice area fraction. This is an experimental hack to speed up ice in calving zones. |
| `amr.basal_length_scale` | `0.0` | Accepted but currently unused. |
| `CalvingModel.preserveSea` | `false` | Accepted but currently unused. |
| `CalvingModel.preserveLand` | `false` | Accepted but currently unused. |
| `geometry.basalSlope` | `0 0` | Read and stored but currently unused. |
