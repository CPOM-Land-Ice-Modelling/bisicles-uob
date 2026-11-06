# BISICLES calving front

BISICLES represents the loss of ice at the ice sheet margin -- the
calving front -- through a *calving model*, which is selected with the
run-time option

    CalvingModel.type = ...

If you have built the doxygen code documentation, the C++ class
hierarchy underlying all this is described at
[classCalvingModel.html](../code/doc/doxygen/html/classCalvingModel.html)

## Calving models and run-time options

The entry point for all calving is the static factory method
`CalvingModel::parseCalvingModel`, which reads the `CalvingModel.type`
entry and dispatches to the appropriate model class. The run-time options
for each of the calving models are described below.

### All models

To impose a calving front at a domain boundary, set the `CalvingModel.front_hi/lo` options.
For example, if the domain is `0 < x < W, 0 < y < H` with a calving front at
` x = W ` and ice divides elsewhere, set 

    CalvingModel.front_hi = 1 0    #marine boundary at x = W, 
    CalvingModel.front_lo = 0 0    

Calving models by default reduce mean ice thickess and ice area fraction until both are zero.
It is possible (but a bodge!) to limit reducution of mean ice thickness in previously ice covered regions 

   CalvingModel.min_thickness = 10 

### NoCalvingModel

    CalvingModel.type = NoCalvingModel

This is the simplest calving model: no ice is removed and the usual result
will be calving front advance at the ice velocity.

### DomainEdgeCalvingModel

CalvingModel.type = DomainEdgeCalvingModel

This model is deprecated, but probably everything will work as before.
It is equivalent to the `NoCalvingModel`, with the front imposed at the
domain edge.

### FixedFrontCalvingModel

    CalvingModel.type = FixedFrontCalvingModel

This model fixes the front, or rather does not allow advance: the calving
criterion marks the open ocean (and open land) cells at the start of
each timestep and prevents ice from entering them

No run-time options are read by this model.

### FlotationCalvingModel

    CalvingModel.type = FlotationCalvingModel

This model removes all floating ice: the calving criterion is satisfied
everywhere the ice is floating.

No run-time options are read by this model.

### BennCalvingModel

    CalvingModel.type = BennCalvingModel
    BennCalvingModel.includeBasalCrevasses = false
    BennCalvingModel.criticalStress = 0.0
    BennCalvingModel.stressMeasure = FirstPrincipalStress

This is a crevasse calving model, following Benn (2007) and Nick et al
(2011). The `crevasseWaterDepth` (or the older-style `waterDepth`) sets the
water depth used in the crevasse depth calculation, and the `factor` is a
dimensionless multiplier on the crevasse depth. Basal crevasses can be
included with `includeBasalCrevasses`. Calving occurs only where `R > criticalStress`
for a range of `stressMeasure` choices: `FirstPrincipalStress`,
`AlongFlowNormalStress`, `Trace`, or `MAX_STRESS_MEASURE`.

### ThicknessCalvingModel

    CalvingModel.type = ThicknessCalvingModel
    CalvingModel.calving_thickness = 0.0
    CalvingModel.calving_depth = 0.0
    CalvingModel.factor_mu_coef = false

This model calves floating ice whose effective thickness is below the
calving thickness. The effective thickness is `h/f' , or 'h/ f *mu_coef'
`factor_mu_coef` is set. The `min_thickness`, `calving_depth`, and the
`start_time` and `end_time` window control the behaviour of the model.

### MaskedCalvingModel

    CalvingModel.type = MaskedCalvingModel

The masked calving model uses a mask -- specified as a surface flux, and
keyed under the `CalvingModel.mask` prefix -- to decide where calving
occurs. For example

    CalvingModel.mask.type = LevelData
    CalvingModel.mask.n = 1
    CalvingModel.mask.timeStep = 1.0e+10
    CalvingModel.mask.startTime = 0.0
    CalvingModel.mask.fileFormat = calvingmask.2d.hdf5
    CalvingModel.mask.name = calvmask

### CompositeCalvingModel

    CalvingModel.type = CompositeCalvingModel
    CalvingModel.nElements = 2
    CalvingModel.element0.type = FixedFrontCalvingModel
    CalvingModel.element1.type = ThicknessCalvingModel

Composite calving models allow several calving models to be applied
serially; the ordering is important.

### VariableRateCalvingModel

    CalvingModel.type = VariableRateCalvingModel

The variable-rate calving model allows a calving rate to be specified
through any of the surfaceFlux classes.

    CalvingModel.CalvingRate.type = LevelData
    CalvingModel.CalvingRate.n = 1
    CalvingModel.CalvingRate.timeStep = 1.0e+10
    CalvingModel.CalvingRate.startTime = 0.0
    CalvingModel.CalvingRate.fileFormat = calvingrate.2d.hdf5
    CalvingModel.CalvingRate.name = calrate

### RateAuBuhatCalvingModel

    CalvingModel.type = RateAuBuhatCalvingModel
    CalvingModel.proportion.type = constantFlux
    CalvingModel.proportion.flux_value = 0.5
    CalvingModel.independent.type = constantFlux
    CalvingModel.independent.flux_value = 0.1
    CalvingModel.independent_normal = false
    CalvingModel.vector = true

This model is also available under the alias
`RateProportionalToSpeedCalvingModel`. The rate is
`r = a * u + b * u / |u|`, combining a part proportional to the ice flow
speed with an independent part. The constants are set by surface fluxes
(`proportion` and `independent`).

### VonMisesCalvingModel

    CalvingModel.type = VonMisesCalvingModel
    CalvingModel.start_time = -1.2345678e+300
    CalvingModel.end_time = 1.2345678e+300
    CalvingModel.scale.type = constantFlux
    CalvingModel.scale.flux_value = 1.0
    CalvingModel.independent.type = constantFlux
    CalvingModel.independent.flux_value = 0.1
    CalvingModel.vector = true

The rate is `r = a * u / |u| + (sigma_vm / sigma_max) * u`, where `a` is an
independent part and the second term is proportional to the Von Mises
stress. The constants are again set by surface fluxes (`independent` and
`scale`).

### CriterionToRateCalvingModel

    CalvingModel.type = CriterionToRateCalvingModel
    CalvingModel.CriterionCalvingModel.type = FixedFrontCalvingModel
    CalvingModel.calving_rate = 1.0e+4

This model converts a calving criterion (another calving model, here
keyed under the `CalvingModel.CriterionCalvingModel` prefix) into a
calving rate: the scalar `calving_rate` is applied where the criterion is
satisfied.

## Specifying calving rates

Most calving models provide a calving rate or a calving velocity rather
than a calving criterion. The rate-based models combine a rate with the
ice flow:

- The `VariableRateCalvingModel` takes a cell-centred calving rate from a
  surface flux.

- The `RateAuBuhatCalvingModel` combines a velocity-proportional calving
  rate `r = a * u + b * u / |u|` with an independent calving rate `b`.

- The `VonMisesCalvingModel` combines a speed-proportional calving rate
  with a Von Mises-stress-proportional calving rate.

- The `CriterionToRateCalvingModel` converts a calving criterion into a
  scalar calving rate.

## Calving in the time step

Calving is applied to the ice thickness and ice fraction at the
`PostThicknessAdvection` stage of the time step, after thickness
advection and after the calving criterion has been evaluated. The ice
fraction is evolved through advection when `amr.evolve_ice_frac` is set
(and through the newer scheme selected by `amr.evolve_ice_frac_2`).

In this way, calving and ice fraction advection are tightly coupled: the
calving rate acts on the ice fraction field as it is being advected with
the ice.
