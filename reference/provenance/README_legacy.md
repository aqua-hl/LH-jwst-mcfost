# JWST–MCFOST workflow for TMC1A

For the current state and the planned move to a new cluster, start with
[START_HERE.md](START_HERE.md) and the [18 September retrospective](reports/restart_review_2026-09-18/REPORT.md).
The original workflow notes below are preserved as history.

This folder keeps the real JWST spectrum separate from the synthetic MCFOST
tutorial. Step 1 is implemented by `extract_tmc1a_sed.py`.

The script takes the on-source coordinate currently used in
`JWST_notebook.ipynb`, extracts a fixed 1 arcsec radius aperture from all three
NIRSpec IFU cubes and all twelve MIRI MRS cubes, and converts the cube values
from MJy/sr to integrated Jy.

Run it from the project directory in the `lesschool` environment:

```bash
conda activate lesschool
python jwst_mcfost/extract_tmc1a_sed.py
```

It creates:

- `output/tmc1a_sed_unstitched.ecsv`: wavelength, flux, formal error,
  lambda-F-lambda, valid aperture coverage, instrument, and segment.
- `output/tmc1a_sed_unstitched.png`: a quick-look plot with every cube segment
  shown separately.

“Unstitched” is intentional. At this point there is no new background
subtraction, overlap rescaling, spectral binning, gas-line masking, or point
source aperture correction. The next step is to inspect the source position,
aperture growth, background, and the NIRSpec/MIRI offset before altering the
spectrum for an MCFOST fit.

## Step 2: inspect the aperture and background

Run:

```bash
python jwst_mcfost/inspect_aperture_background.py
```

This does not modify the Step 1 spectrum. It makes continuum maps for all 15
segments, measures the peak and flux-weighted centroid, compares 0.5, 1.0, and
1.5 arcsec aperture radii, and tests a provisional 1.5--2.0 arcsec background
annulus. The annulus may contain real envelope or cavity emission and is partly
outside some cube fields, so its subtraction is diagnostic rather than final.

## MIRI registration test

After Step 2, run:

```bash
python jwst_mcfost/test_miri_registration.py
```

This uses the median MIRI centroid offset measured in Step 2 to re-extract the
MIRI cubes with a shifted 1 arcsec aperture. It compares that experiment with
the unchanged Step 1 spectrum, including the NIRSpec/MIRI overlap and the
fraction of the shifted aperture that remains inside each cube. It does not
replace the adopted spectrum automatically.

## Prepare the continuum SED

The fixed-coordinate extraction remains the adopted input. Prepare an
unscaled, continuum-focused spectrum at resolving power 100 with:

```bash
python jwst_mcfost/prepare_continuum_sed.py
```

The script masks low-coverage samples, NIRSpec detector gaps, and known narrow
gas lines, then robustly clips remaining positive spikes within logarithmic
bins. Every detector segment is binned independently on the same wavelength
grid; overlapping segments are retained and are never forced to agree. Broad
ice and silicate features remain in the data. For experimentation, use for
example `--resolution 200` to retain twice as much spectral detail.

The main `continuum_sed_R100.png` is intentionally clean and shows only native
samples used for binning plus the final bins. Rejected gas-line and positive
spike samples are shown separately in `continuum_mask_diagnostic_R100.png`.

## First MCFOST baseline comparison

Before fitting any parameters, compare the prepared JWST continuum to the
existing disk+envelope+cavity model:

```bash
python jwst_mcfost/compare_baseline_mcfost.py
```

This reads the exact saved baseline run in
`mcfost-tuto/runs/mock_truth_seeded/seed=24680`, whose parameter snapshot records
the command `mcfost disk_envelope_1.para -cavity 274.7 100 1 -seed 24680
-root_dir runs/mock_truth_seeded`. It applies only the known inverse-square
distance correction from the model's 140 pc to the adopted 147 pc. It does not
fit a scale or alter a physical parameter.

The products are `output/baseline_mcfost_comparison.png`, a point-by-point
`output/baseline_mcfost_comparison.ecsv`, and wavelength-region summaries in
`output/baseline_mcfost_regions.ecsv`.

Interpret this as a shape diagnostic only. The saved MCFOST SED samples the
spectrum at about R=20 and is the total spatially integrated model flux; the
JWST continuum is R=100 and came from a 1 arcsec radius aperture. A later model
run must improve the wavelength sampling and match the synthetic aperture
before a formal chi-squared fit is meaningful.

## Repeat the unchanged model on the JWST R=100 grid

The continuum table has 358 detector-segment rows but only 312 unique retained
wavelength bins. Generate a sorted MCFOST wavelength file, copy the existing
parameter file with only that filename changed, run the same 20-degree cavity
model, and make the updated comparison with:

```bash
bash jwst_mcfost/run_baseline_R100.sh
```

The generated inputs and their provenance live in
`model_inputs/baseline_R100`. The model output goes to
`runs/baseline_R100/data_th`, so the original tutorial runs are never
overwritten. No random seed is specified. The run retains the original photon
counts and uses 312 output wavelengths instead of 80, so the SED portion is
expected to take roughly four times longer. The completed baseline took 26m40s
on 14 processors.

When MCFOST finishes, the same script creates
`output/baseline_R100_mcfost_comparison.png` and its two ECSV tables. This fixes
the wavelength-sampling mismatch, but it is still a total-model versus 1 arcsec
aperture diagnostic rather than a formal fit.

## Representative 1 arcsec synthetic apertures

Use the completed R=100 temperature structure to ray-trace full-field images
at 1.0, 1.5, 2.2, 3, 6, 10, 20, and 28 microns, then measure the same fixed
1 arcsec aperture used for the JWST extraction:

```bash
bash jwst_mcfost/run_aperture_images.sh
```

The runner repeats the command-line cavity, reuses
`runs/baseline_R100/data_th/Temperature.fits.gz`, and makes 601 by 601 images
without changing the physical model. Completed wavelengths are safely reused
if the script is restarted; partial output directories cause it to stop.

`measure_aperture_images.py` centers the circle on the MCFOST FITS reference
pixel and uses a radius of 147 AU, corresponding to 1 arcsec at the adopted
147 pc distance. It reports both the intrinsic fraction and a version after
convolution with approximate wavelength-dependent Gaussian JWST PSFs. No
recentring or annulus subtraction is applied, matching the adopted spectrum.
The approximation is not a replacement for STPSF/WebbPSF and the eight anchor
fractions should not be interpolated through solid-state features.

Outputs are `output/aperture_fractions_1arcsec.ecsv`,
`output/aperture_fractions_1arcsec.png`, and
`output/aperture_images_1arcsec.png`.

## Near-IR viewing-inclination test

Hold the star, disk, envelope, grains, cavity, photon counts, and temperature
structure fixed while comparing 65, 75, and 82 degree viewing inclinations at
1.0, 1.5, 2.2, and 3.0 microns:

```bash
bash jwst_mcfost/run_inclination_test.sh
```

The generated working directories live under `runs/inclination_test`. Each
parameter copy differs only in the map-inclination line and reuses the verified
baseline temperature file. The analysis measures the same 1 arcsec aperture
and approximate PSF, then calibrates the relative image changes to the
completed 75-degree R=100 aperture SED. Results are
`output/inclination_nearir_test.ecsv` and
`output/inclination_nearir_test.png`.

## Near-IR envelope-mass test

After the inclination test, vary only the Zone 2 envelope dust mass by factors
of 0.5, 1, 1.5, 2, and 3 while keeping the 75-degree view, cavity, star, disk,
grains, and photon counts fixed:

```bash
bash jwst_mcfost/run_envelope_mass_test.sh
```

Every mass gets a fresh radiative-equilibrium temperature calculation. The
full 312-point SED is disabled during this diagnostic, avoiding its long
runtime. Images at 1.0, 1.5, 2.2, 2.7, 3.0, and 3.3 microns test both the
near-IR slope and the local shape of the 3 micron trough. The 2.7--3.3 micron
index is not a formal ice optical depth because both endpoints can remain in
the broad feature. The same 1 arcsec aperture and approximate PSF are then
applied. Results are
`output/envelope_mass_nearir_test.ecsv` and
`output/envelope_mass_nearir_test.png`.

The completed grid shows that envelope mass is not a sufficient one-parameter
solution. Around 2 times the baseline mass, the 1.0--1.5 micron fluxes approach
the data, but the model is already about 6 times too faint at 2.2 microns and
86 times too faint at 3.0 microns. Increasing mass also makes the local 3
micron trough much deeper than observed. This points to the wavelength-
dependent envelope dust opacity/scattering as the next controlled test.

## Near-IR cavity-opening test

Test whether the scattered-light escape geometry can improve the near-IR shape
by varying only the bipolar cavity half-opening angle through 10, 15, 20, 25,
and 30 degrees:

```bash
bash jwst_mcfost/run_cavity_opening_test.sh
```

The envelope mass remains at the baseline value and the viewing inclination
remains 75 degrees. Each cavity receives a fresh temperature calculation,
followed by images at 1.0, 1.5, 2.2, 2.7, 3.0, and 3.3 microns. The 20-degree
model is the existing geometry and provides the relative-flux reference.
Results are `output/cavity_opening_nearir_test.ecsv` and
`output/cavity_opening_nearir_test.png`.

This is a coarse sensitivity pre-test. Widening the cavity also removes more
envelope material, and the 70-cell polar grid makes the cavity edge move in
discrete angular steps. Do not interpolate these five results as a precise
best-fitting angle; a later fitting grid should refine the spatial grid and
rerun its own baseline.

The completed sweep confirms that cavity opening affects the flux but does not
repair the near-IR spectral shape. The six-anchor log-RMS changes only from
0.864 dex at 10 degrees to 0.799 dex at 30 degrees. Even at 30 degrees the
model is about 29 times too bright at 1 micron, while it is about 3.6, 5.0,
and 2.5 times too faint at 2.7, 3.0, and 3.3 microns, respectively. Therefore
30 degrees is only the lowest diagnostic score among these five trials, not a
physical best fit. The next controlled pre-test should change the envelope
dust opacity/scattering properties rather than refine the cavity angle.

## Envelope dust screening grid

Combine the three parameters that have useful near-IR leverage in a 27-model
screen:

- Zone 2 envelope dust mass: 0.5, 1.0, and 1.5 times the baseline value.
- Zone 2 maximum grain size: 0.3, 1.0, and 3.0 microns.
- Ice-mantle volume fraction: 0.10, 0.30, and 0.50, with the silicate-core
  fraction set to one minus the ice fraction.

Run it with:

```bash
bash jwst_mcfost/run_dust_grid_screen.sh
```

Inclination remains 75 degrees and the cavity remains
`-cavity 274.7 100 1.0`. Every model receives a fresh temperature calculation
and six 601 by 601 images at 1.0, 1.5, 2.2, 2.7, 3.0, and 3.3 microns. The
runner is restart-safe for completed products and refuses partial output
directories. Inputs, hashes, and model metadata live under
`runs/dust_grid_screen`.

After completion, `validate_dust_grid_screen.py` checks all 27 temperatures,
162 images, FITS shapes and units, logs, saved commands, parameter prefixes,
and component-flux closure before the scientific ranking is accepted.

The analysis uses the exact baseline grid member (mass factor 1, maximum grain
size 1 micron, ice fraction 0.30) as the dynamic aperture-flux reference. It
ranks models without fitting a free scale using the six-anchor logarithmic
RMS, worst mismatch, 2.2-micron residual, continuum shape, and local 3-micron
curvature. Products are `output/dust_grid_screen.png`, the per-anchor and
per-model ECSV tables, and an automatically selected diverse shortlist.

This remains a screening experiment, not the final JWST fit. Only shortlisted
models should receive independent full R=100 SED calculations and denser
mid-infrared aperture checks.

The completed run produced and validated all 27 temperature structures and all
162 images, with no structural, provenance, or Monte Carlo quality warning.
The lowest six-anchor RMS is the model with 1.5 times the baseline envelope
mass, 1-micron maximum grains, and a 0.10 ice fraction. Its RMS is 0.531 dex,
compared with 0.840 dex for the baseline member, but it is not yet a good fit:
the model remains about 6.0 times too bright at 1 micron and about 2.7, 5.1,
1.6, and 2.8 times too faint at 2.2, 2.7, 3.0, and 3.3 microns. The preferred
mass and ice coordinates are both grid edges, so these values are directions
for a refinement test rather than measured optima.

Before running full R=100 SEDs for every shortlisted model, use a lean local
screen at fixed envelope `amax = 1` micron: mass factors 1.25, 1.50, and 1.75
combined with ice fractions 0.10, 0.15, 0.20, and 0.25. This is 12 models. The
observed local 3-micron curvature lies between the current 0.10 and 0.30 ice
predictions, so pushing the ice fraction below 0.10 merely because it is the
lowest-RMS grid edge would be misleading. Once mass and ice are bracketed,
refine `amax` locally and calculate full R=100 SEDs for the best two or three
models plus the baseline.

## Local envelope mass/ice refinement

Run the 12-model refinement with:

```bash
bash jwst_mcfost/run_dust_refinement_screen.sh
```

This holds the envelope maximum grain size at 1 micron and combines mass
factors 1.25, 1.50, and 1.75 with ice fractions 0.10, 0.15, 0.20, and 0.25.
Every point receives a fresh temperature calculation and the same six
diagnostic images as the coarse screen. The physical point at mass factor 1.5
and ice fraction 0.10 is deliberately repeated; its new-to-old aperture flux
ratios measure Monte Carlo drift and are never used to renormalize the grid.

The completed run validated all 12 temperatures and all 72 images. The repeat
control has an RMS drift of 0.025 dex and a maximum drift of 0.041 dex (a
factor of 1.10). The lowest absolute RMS is again mass factor 1.5 and ice
fraction 0.10, at 0.523 dex, but this is statistically indistinguishable from
its old 0.531-dex realization and remains on the lower ice boundary. It still
overpredicts the 1-micron flux by a factor of 6.17 and underpredicts the
2.7-micron flux by a factor of 5.19, so it is not a satisfactory SED fit.

Two complementary models are retained for full R=100 calculations:

- Mass factor 1.50, ice fraction 0.10: lowest six-anchor RMS and smallest worst
  anchor residual, but its local 3-micron trough is too shallow.
- Mass factor 1.75, ice fraction 0.15: closest local 3-micron curvature, but it
  lies on the upper mass boundary and has a worse absolute RMS of 0.581 dex.

Because the top score gaps are comparable to the repeat-control variation,
their exact ranking is unresolved. These are complementary test cases, not
parameter estimates. The next defensible step is to compute an independent
full R=100 SED for both, compare the full continuum without a free scale, and
retain the unchanged baseline as a control.

Results are `output/dust_refinement_screen.png`, the per-anchor and per-model
ECSV tables, `output/dust_refinement_screen_finalists.ecsv`, and the strict
post-run report `output/dust_refinement_screen_validation.json`. Generated
inputs and MCFOST products live under `runs/dust_refinement_screen`.

## Full R=100 finalist SEDs

Compute independent temperature structures and spatially integrated R=100
SEDs for the two complementary finalists with:

```bash
bash jwst_mcfost/run_full_r100_finalists.sh
```

The runner freezes the exact refinement inputs, applies the same
`-cavity 274.7 100 1.0` geometry, and validates both complete products before
running `compare_full_r100_finalists.py`. No random seed or fitted flux scale is
used. The output is a useful full-wavelength check, but it is not aperture
matched: it compares the total model flux with the fixed 1 arcsec JWST
extraction.

That total-SED diagnostic marginally preferred mass factor 1.75 and ice
fraction 0.15 (0.553 dex logarithmic RMS) over mass factor 1.50 and ice
fraction 0.10 (0.563 dex). The 0.010-dex gap is too small to interpret as a
parameter measurement, especially before matching the spatial aperture.

Results are `output/full_R100_finalists_comparison.png`, the bin, regional, and
summary ECSV tables with the same prefix, and
`output/full_R100_finalists_validation.json`.

## Full R=100 synthetic 1 arcsec apertures

Ray-trace a wavelength-specific 601 by 601 image at every one of the 312 R=100
wavelengths for both finalists, reuse each finalist's already validated
temperature structure, and perform the matched-aperture comparison with:

```bash
bash jwst_mcfost/run_full_r100_aperture_images.sh
```

This is 624 monochromatic image calculations. The runner is restart-safe and
validates every image before reuse. After all images pass, the analyzer
measures an exact fractional 1 arcsec circular aperture centered on MCFOST
`CRPIX`, applies a detector-row-specific approximate Gaussian PSF, collapses
the 358 detector rows to 312 unique wavelength bins, and compares with no free
normalization, background subtraction, recentering, or aperture correction.
Each prediction is the model's own total R=100 SED multiplied by its measured
image aperture fraction and the known 140-to-147-pc inverse-square distance
factor. Raw monochromatic image totals are retained only as a ray-tracing
consistency diagnostic.

The completed run validated all 624 images and all 624 hashes are unique. The
aperture-matched result reverses the very small total-SED ordering: mass factor
1.50 and ice fraction 0.10 has a logarithmic RMS of 0.505 dex, versus 0.528 dex
for mass factor 1.75 and ice fraction 0.15. The corresponding typical mismatch
factors are 3.20 and 3.38. The first model is therefore the better starting
point, not a satisfactory fit or a measured optimum. It remains too bright
below about 1.5 microns and too faint by factors of several through the
2.5--13-micron ice and silicate regions. The ranking advantage is driven by
those broad dust-feature regions; the 0.024-dex overall gap is not large.

The synthetic aperture contains only about 34--40 percent of the model flux at
1 micron, about 80 percent at 3 microns, and more than 93 percent near 5
microns. This wavelength-dependent extended scattered light is why the matched
aperture materially changes the near-infrared residuals and why a total-SED
comparison was insufficient.

The Gaussian PSFs are an approximation, not STPSF/WebbPSF or a reconstructed
IFU-cube response. Cube covariance, MCFOST Monte Carlo uncertainty, the exact
per-cube valid-pixel footprint, and integration through finite spectral bins
are also omitted. Therefore this is an aperture-matched finalist ranking, not
a likelihood or a formal full-grid fit.

Results are `output/full_R100_aperture_finalists_comparison.png`, detector-row,
aperture-fraction, bin, regional, and summary ECSV tables with the same prefix,
and the strict report
`output/full_R100_aperture_images_validation.json`. The monochromatic products
and frozen provenance occupy about 3 GB under
`runs/full_R100_aperture_images`.

## Aperture-matched 36-model feature grid

The next screen varies the three envelope-dust coordinates jointly over a
wider, feature-focused grid:

- Envelope dust-mass factor: 1.25, 1.50, and 1.75.
- Envelope maximum grain size: 0.5, 1.0, and 2.0 microns.
- Ice volume fraction: 0.05, 0.10, 0.15, and 0.20, with the core fraction set
  to one minus the ice fraction.

Run or safely resume all 36 models with:

```bash
bash jwst_mcfost/run_feature_grid_screen.sh
```

Each model receives an independent temperature calculation and a reduced
20-point SED spanning the broad JWST continuum and solid-state features. The
runner also ray-traces a 601 by 601 image at each wavelength, for 720 images in
total. The analyzer measures each model's own 1 arcsec aperture fractions,
applies the same approximate Gaussian PSFs and 140-to-147-pc distance
correction, and uses no fitted normalization, background, recentering, or
aperture correction. Seven wavelength regions receive equal weight so that a
dense part of the spectrum cannot dominate merely by containing more samples.

The completed run validated all 36 temperature structures, all 36 SEDs, and
all 720 unique images. The nominal lowest score is
`m1p75_a0p5_i0p10`: mass factor 1.75, maximum grain size 0.5 microns, and ice
fraction 0.10, with a region-balanced logarithmic RMS of 0.475 dex (a typical
factor of 2.99). This is not a good fit or a parameter measurement. It is still
about five times too bright in the shortest near-IR region, roughly three
times too faint around the ice/short-MIRI region, and about five times too
faint across the silicate region; the largest individual discrepancy is a
factor of 9.8 near 9.7 microns. It agrees much better only in the mid- and
long-MIRI regions, at typical factors of about 1.3--1.5.

The preferred mass is the upper grid boundary and the preferred grain size is
the lower boundary. Nearly every model is region-Pareto-optimal, demonstrating
strong trade-offs between wavelength regions rather than a single isolated
minimum. More importantly, an independent repeat shows 0.066 dex RMS drift in
the reduced total SED and changes the control score by 0.024 dex, larger than
the score gaps among ranks 1--8. The image aperture fractions themselves are
repeat-stable. Therefore the exact top ordering is not numerically converged at
the current `nbr_photons_lambda = 5.12e3`.

Before extending the grid or starting MCMC, rerun only the reduced SEDs for
the top eight candidates and the repeat control with at least
`nbr_photons_lambda = 1.28e5`. Reuse the validated temperatures and aperture
fractions; there is no reason to repeat the 720 images yet. If the scores then
stabilize, the next boundary grid should test mass factors 1.75, 2.0, and 2.25;
maximum grain sizes 0.25, 0.5, 0.75, and 1.0 microns; and ice fractions 0.05,
0.10, and 0.15. If the sign-changing residual pattern persists, the next model
family must vary geometry, stellar luminosity/extinction, or dust composition
instead of asking a three-parameter MCMC to refine a structurally inadequate
model.

Results are `output/feature_grid_screen.png`, the row, aperture-fraction, bin,
region, summary, and shortlist ECSV tables with the same prefix, and
`output/feature_grid_screen_validation.json`. Frozen inputs and MCFOST products
occupy about 5.3 GB under `runs/feature_grid_screen`.

### MCFOST SED/image sampler diagnosis

Exact MCFOST 4.1.13 source inspection and the controlled plot
`output/mapsize_vs_pixelscale_4p5um.png` separate three numerical choices that
were previously conflated.  `-rt2` selects the RT2 source-function treatment;
it does not change the SED spatial quadrature.  SEDs use the method-1
log-polar/no-subpixel integral out to `2*Rmax`, whereas images use method-2
Cartesian pixels with adaptive subpixels.  At 4.5 microns, changing map extent
at fixed approximately 10-AU pixels leaves `image_sum/SED_total` flat near
0.62, so the method-1 outer annulus is not leaking flux.  Refining pixel scale
at fixed 6000-AU field instead drives the ratio toward one (non-monotonically),
identifying unconverged Cartesian sampling as the raw-image-total error.

The final policy is to replace the log-polar predictor with a converged
method-2 calculation.  A separate SED2 and `-img` launch is still not an exact
pair: image mode counts emitted packets, whereas SED mode targets accepted
detector packets, so equal numeric photon settings and seed do not give the
same Monte Carlo realization.  The corrected one-wavelength SED2 run therefore
writes both `sed_rt.fits.gz` and a coeval `RT.fits.gz` from its existing
`Stokes_ray_tracing` map.  Its method-2 pixel integral is accumulated in
float64; the stock float32 loop loses up to 1.313% on the existing 2401-squared
test.  Log-polar products and separately launched images are retained only as
controls.

A coarse method-2 SED closing against its coeval coarse map proves internal
consistency, not physical convergence.  Per-component same-realization closure
and stability between the two finest pixel grids are both required before
fitting.  The exact held-job state, custom build provenance, and replacement
grid gates are recorded in `CONTINUE_HERE.md`.

## Published central-NIRSpec spectrum test

As a separate sanity check, compare the current nominal grid winner directly
with `DATA/JWST_ver196_center.mrt` using:

```bash
python jwst_mcfost/test_current_model_against_mrt.py
```

This MRT table is not identical to the cube spectrum used above. It is a
stitched NIRSpec-only spectrum from the central 0.3 by 0.3 arcsec extraction,
whereas `continuum_sed_R100.ecsv` uses a fixed 1 arcsec-radius circle and also
contains MIRI. The MRT is typically about twice as bright as the local cube
extraction, demonstrating an additional reduction/calibration difference; the
two spectra must not be combined as duplicate measurements.

The test reuses the winner's existing SED and ten NIRSpec images. It robustly
bins the 9,358 valid MRT samples at the ten precomputed model anchors and
remeasures an approximate PSF-convolved 0.3 by 0.3 arcsec square centered on
MCFOST `CRPIX`. No MCFOST rerun or physical-parameter adjustment occurs. The
published brightest-spaxel centering, exact reconstructed-cube response, and
calibration/aperture correction are unavailable, so this remains a diagnostic
rather than an exact spatial match.

The absolute three-region RMS is 0.607 dex, a typical mismatch factor of 4.04.
A single diagnostic multiplier of 3.27 reduces the remaining shape RMS to
0.322 dex, still a factor of 2.10. The model is about 9.5 times too faint near
2.7 microns and 11.6 times too faint near 4.5 microns before scaling. This
three-region NIRSpec score must not be compared numerically with the main
seven-region NIRSpec+MIRI grid score.

Outputs are `output/current_model_vs_JWST_ver196_center.png`, the matching
`_anchors.ecsv` table, and `_summary.ecsv`.

## Cluster arrays: convergence first, literature screen second

The long nine-model high-photon convergence calculation has a Slurm-array
version for the ENS-Lyon reservation. Each array element owns one model
directory and runs one MCFOST process with 24 OpenMP threads; the default
`0-8%9` throttle permits all nine models (216 allocated CPUs) at once when
the shared reservation has room. A separate dependent job validates all nine
products before making the refined
ranking and plot. Cancellation, non-destructive partial-output recovery,
submission, monitoring, and retry commands are in
`HIGH_PHOTON_ARRAY_RECOVERY.md`.

After that convergence test identifies a numerically separated winner, the
literature-informed pilot can be submitted with
`submit_literature_screen.sh`. It tests inclinations 47, 50, 53, and 60 degrees
against central luminosities 2.3, 2.5, and 2.7 solar luminosities. Each of its
12 array elements computes an independent temperature structure, a 20-anchor
SED at 32,000 photon packets per wavelength, and 20 wavelength-specific
601-by-601 images. The analysis uses each model's own PSF-convolved 1-arcsec
aperture fraction at 140 pc, with no fitted flux scale or distance correction.
See `LITERATURE_SCREEN_CLUSTER.md` for its scientific caveats and cluster
commands. This remains a screening experiment; full production modeling uses
all 312 R=100 wavelengths only for the resulting shortlist.
