# Restart goal and acceptance criteria

Prepared 2026-09-18. Design only; the user will move to a new cluster.

## Scientific question

Can a thermally self-consistent disk-plus-envelope model, compatible with independently measured geometry and heating, reproduce TMC1A's aperture-matched continuum, spatial structure, H2O band and silicate region? If not, establish which assumption is inadequate before reporting an abundance posterior.

## What the previous search contributes

- Keep the exact nine-anchor continuum data/measurement contract, the excluded first point, H2O observation v3, Method-2 tests and signed-aperture kernel as benchmarks.
- Retain the nominal mass2.25/q2.75/cavity17.5 solution, cavity20 near-tie and mass2.5/q3.25 family as controls.
- Treat the latest DHS vmax0.1/amax0.4/xH2O0.04 model as a useful conditional dust variant. It is not a final solution: its H2O red wing is weak, its silicate dip is too deep, and temperature was reused.
- Reuse compact predictions for comparisons on the same observable; regenerate images when changing apertures or spatial geometry.
- Preserve old-estimator experiments as history. Do not combine their scores into a new likelihood.

## Stage 0: freeze data and reproduce controls

Record the actual cube versions/CRDS contexts, extraction centre, aperture radius, background treatment, masks, overlap handling and PSF. Resolve the approximately 0.122-arcsec centre difference between the old continuum and new H2O extraction. Keep 1.0-arcsec and 0.35-arcsec measurements distinct. Do not assume the older MRT product equals either one.

Rebuild or validate MCFOST spatial Method 2 on the new cluster. Reproduce a few archived wavelengths for the old continuum model, its cavity20 near-tie and the latest ice variant using their stated frozen inputs. Include a short-wave scattering point and a silicate point. Verify numerical closure where coeval products are requested and separately test aperture prediction stability with pixel scale and seed. The archived images cannot be presumed recoverable from the retired cluster.

Exit: the new code reproduces the archived controls within a declared Monte Carlo/numerical tolerance small compared with observational errors; all scientific parameter and opacity definitions have a portable manifest. A provisional target is less than 2% in representative aperture fluxes, tested rather than assumed for each relevant regime.

## Stage 1: test thermal consistency and literature geometry

Proposed 12 physical models = inclinations [50,60,70] degrees × source-heating targets [1.9,2.7] Lsun × [original coated dust, latest DHS ice dust]. Keep nominal mass2.25, q2.75 and cavity17.5 initially to isolate responses. Recompute radiative-equilibrium temperature for each physical model and document stellar/accretion heating separately.

Evaluate nine continuum and four silicate anchors for each (156 image calculations, plus twelve temperature solves). Use the 85-node H2O operator only on an initial shortlist of three (255 additional wavelength nodes). Counts exclude numerical control repeats. Start from 128k packets but change numerical settings where the convergence evidence requires it.

Exit: determine whether the near-IR/H2O/silicate conflict survives fresh temperatures, and how inclination and heating change it. A poor fixed-mass pilot at a literature inclination does not rule out that inclination after mass/dust refitting. Choose the next search axes from these response patterns.

## Stage 2: fit continuum and spatial constraints jointly

Use 48–96 initial space-filling/adaptive physical settings rather than the entire Cartesian product. Initial candidate bounds are in config/restart_proposal.json. They are exploratory design bounds, not literature confidence limits. Give inclination 50–65 degrees priority and retain 70 as a control. Condition dimensions using morphology/kinematics to avoid inferring six loosely constrained parameters from only nine continuum anchors.

Vary column/mass, grain size and heating together with viewing geometry, and test whether q/cavity/disk parameters need to be released. Include the old minima and boundary directions. Predict held-out continuum wavelengths and spatial/radial profiles at a few diagnostic bands. Match position angle only when comparing spatial data; the old PA=0 was not a literature morphology fit.

Exit: improved out-of-sample predictions without large structured residuals, demonstrated numerical stability, and a relevant minimum bracketed within the tested bounds. If a best model sits at a boundary, extend or physically justify that boundary before inference.

## Stage 3: ice and silicate with physical dust structure

Evaluate dust representations and laboratory optical constants using new temperatures and explicit ice survival/distribution assumptions. Preserve separate disk and envelope grain populations. Compare the original mantle-volume fraction with disjoint species mass fractions only after converting using the correct densities/definitions. Do not interpret DHS vmax as measured bulk porosity.

Use a common forward operator for both model and observation. Do not count the same flux spectrum and its derived optical depths as independent data. Include calibration and continuum covariance. The different 1-arcsec and 0.35-arcsec observables need separate aperture predictions and an appropriate joint error model. Add CO/CO2 after the H2O/silicate model passes; gas-phase co.dat is not a CO-ice opacity component.

Exit: continuum and band profiles can be explained within declared observational plus model/numerical uncertainty, or a clear model failure is quantified. A broad red wing is not evidence that a laboratory material is wrong by itself.

## Stage 4: production inference

Only after validating the likelihood, run MCMC/nested sampling or an emulator validated with withheld direct MCFOST calculations. Preserve a finite budget and monitor sampler/seed/resolution convergence. Report posterior predictive spectra/images and dependence on dust family, distance and geometry assumptions. If the model remains inadequate, report the failure instead of a narrow conditional abundance error bar.

## Deliverables and stopping rule

Deliver a versioned data contract, portable environment/build record, a parameter catalogue with physical units, a complete run manifest, numerical validation plots, posterior/score-landscape diagnostics and grey-JWST-plus-anchor spectral plots. Success is reliable prediction and an honest account of identifiable parameters; it is not merely a lower value of the old screening score.

No cluster address, reservation, account, executable path or wall-time is assumed here. The next implementation should obtain those from configuration and use independent OpenMP model tasks only after benchmarking the new hardware.
