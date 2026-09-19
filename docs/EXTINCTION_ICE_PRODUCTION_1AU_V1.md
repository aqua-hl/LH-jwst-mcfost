# Production at approximately 1 AU/pixel

This run tests how much the spectrum constrains envelope dust and ice **within
the existing generic-ice coated-Mie prescription**. It measures conditional
responses and degeneracies; it need not find uniquely constrained parameters
to return a useful result. Geometry, heating, mass and grain-size alternatives
remain in the catalogue. A separate final numerical comparison ends the
pixel-scale investigation after recording its outcome.

## Frozen production design

- Retain all **96 physical parameter combinations** from v1.1, in their original
  order. They receive fresh temperatures and new images, not cached old fluxes.
- Add **48 controlled cases**: inclinations 45/50/55/60/65/70°, envelope dust
  masses 0.00015/0.000225 M☉, and ice **mass fractions** 0/0.02/0.04/0.08.
  Other parameters match the original nominal model. Selection does not use
  fit rankings. Each inclination has a matched mass contrast.
- Zero ice is a genuine single-component bare-silicate Mie species. Positive
  fractions use the existing generic ice at density 1.2 g cm⁻³ and silicate at
  3.5 g cm⁻³. The code converts mass fractions to mantle volume fractions.
  Total envelope dust mass includes ice; the outer grain-size distribution is
  fixed within a matched ladder. The retained 96 use their original 5% ice
  **volume** fraction, approximately 1.77% by mass.
- **2,048,000 packets** for temperature, SED and image settings; seed 43001.
  There is one fresh temperature solve per physical model, including every
  change in ice fraction. Production uses one seed per model.
- All images use **6001 × 6001 pixels over 6000 AU**, or 0.999833 AU/pixel.
  The odd grid keeps the source centred. Source-centred 1″ aperture, Gaussian
  PSF, 140-pc model distance and 147-pc measurement distance are preserved.
- `aperture_v2` keeps whole-field edge/convolution losses as recorded warnings.
  Finite positive flux, valid image geometry and aperture/PSF support remain
  mandatory. No numerical outliers or inconvenient bands are silently removed.

The H₂O 30 K constants are **not** used: their full-wavelength extension and
thermal application have not been validated. This run does not infer ice
temperature, laboratory composition, grain-shape preference, or a physical
abundance posterior. Ice survival is not temperature dependent in the adopted
prescription. True bare controls do not themselves certify an infinitesimal
mantle opacity-limit test.

## Observables and weights

| Block | Primary bands | Declared comparison weight |
|---|---:|---:|
| Near-IR continuum, approximately 1.1–4.55 µm | 9 | 45% |
| Broad 3-µm ice absorption windows | 4 | 30% |
| MIRI broad bands, including 9.7-µm silicate | 6 | 25% |

The near-IR receives **75% in total**. Each block is normalized by its number
of bands, so extra image wavelengths do not become extra independent scores.
MIRI w007/w008/w009 are check-only with **zero primary weight**. The six broad
MIRI windows include 8.2–8.6, 9.3–10.1, 10.7–11.3, 12.4–13.0, 18.1–18.8 and
19.3–20.3 µm. The 18.1–18.8-µm window uses CH4-short only.

Each observed band is a wavelength-width mean Fν from explicitly recorded
R100 source-bin overlaps within one detector segment. NIR continuum windows
respect the old feature mask; H₂O and silicate feature windows are explicitly
unmasked. The four H₂O windows use the **same 1″ spectrum** as this continuum
analysis, not the separate 0.35″ R2400 optical-depth contract. Do not combine
their results as if they were the same aperture or spectral operator.

Five model wavelengths per band give a composite Simpson mean; a nested
three-node estimate records integration sensitivity. A difference exceeding
1% is flagged, not hidden or converted into a convergence claim. These are
**95 primary image wavelengths plus three MIRI checks = 98 images/model**:
**144 temperature solves and 14,112 images** overall. The embedded difference
does not bound all unresolved spectral structure.

The baseline covariance uses 3% correlated gains per NIR segment and 2% per
MIRI sub-band, with 1%/5% MIRI alternatives. Conservative formal errors and
measured within-bin structure are retained without a square-root-of-count
reduction. A separate 3% common multiplicative-scale nuisance is declared.
These gain/structure choices are sensitivity assumptions, not measured
absolute-calibration or artifact uncertainties. The analysis also reports
NIR-only and equal-NIR/MIRI weighting scenarios.

Matched ice/bare and mass contrasts are assessed with covariance and a free
common-amplitude nuisance. The output identifies weak or degenerate responses
as well as detectable conditional contrasts. Rankings are exploratory
screens, not posteriors or an assertion of independent information per band.
The **75% near-IR weight applies to screening/ranking**. The separate response
norms and singular values use the full 19-band conditional covariance without
regional weights, with the three MIRI calibration scenarios reported explicitly.

## Resources

Each array task requests **64 CPUs and 160 GB**, runs one physical model, and
uses MCFOST `-max_mem 112`. At most **16 production tasks** run concurrently:
1024 CPUs and 2560 GB of requested aggregate memory. A 24-hour model-task limit
and four-hour timeout per MCFOST invocation are initial resource bounds;
completed stages resume if a task needs another allocation. Analysis requests
two CPUs, 8 GB and two hours.

The final diagnostic requests five tasks at the same per-task resources.
Submitting both arrays simultaneously permits **21 model tasks / 1344 CPUs**
before scheduler limits. To hold the combined ceiling at 16, set production's
`slurm.max_parallel` to 11 before configuration, or submit its array after
the final array completes. Jobs may remain pending for memory or CPU capacity.
These are resource requests, not measured runtime promises. Keep raw products
on cluster storage; 14,112 large images can require substantial space.

## Commit, pull, prepare and submit

Commit and push the staged code from the laptop. On the cluster, load the
working MCFOST environment and use the existing shared Python environment:

```bash
cd ~/LH-jwst-mcfost
git pull --ff-only

# Final numerical comparison: use the actual original numerical-run name.
~/.venvs/mcfost-v11/bin/python -B scripts/build_extinction_ice_crossover_v2.py \
  --source runs/extinction_ice_numerics_v2_git \
  --final-resolution --output runs/final_resolution_1au_v1

# Production is independent of that diagnostic's result.
~/.venvs/mcfost-v11/bin/python -B scripts/build_extinction_ice_production.py
```

If the original numerical run is named `runs/extinction_ice_numerics_v2`, use
that as `--source` instead. Only the final diagnostic needs its saved
temperature and receipts. Production builds entirely from committed reference
inputs; no downloaded results, package archive or laptop-only path is needed.
Builders refuse to overwrite an existing destination.

Copy the successful cluster `partition`, `account`, `reservation` and module
`setup_lines`, if applicable, into **both** new `machine.template.json` files.
Export the working `MCFOST_UTILS` if it is not already set. For each run:

```bash
final_run=runs/final_resolution_1au_v1
production_run=runs/extinction_ice_production_1au_v1

~/.venvs/mcfost-v11/bin/python -B "$final_run/code/configure_cluster.py" \
  "$final_run" --machine "$final_run/machine.template.json"
~/.venvs/mcfost-v11/bin/python -B "$production_run/code/configure_cluster.py" \
  "$production_run" --machine "$production_run/machine.template.json"

bash "$final_run/submit.sh" "$PWD/$final_run/machine.template.cluster.json"
bash "$production_run/submit.sh" "$PWD/$production_run/machine.template.cluster.json"
```

The configure step validates imports, immutable inputs, shared Python,
MCFOST executable and utilities without invoking MCFOST. Every worker repeats
the checks, avoiding the earlier node-dependent SciPy/Astropy environment.
Both submit scripts automatically submit a dependent **afterany** analysis
and preserve job-ID receipts under `submissions/`. An incomplete run therefore
still produces an explicit incomplete report. No analysis emails are sent.

## Check status and collect results

Use each printed array ID with:

```bash
squeue -r -j JOBID -o '%.22i %.10T %.10M %.30R'
sacct -j JOBID -X --format=JobID%22,State%20,NodeList%20,Elapsed,ExitCode

~/.venvs/mcfost-v11/bin/python -B "$production_run/code/production_task.py" \
  "$production_run" --status

# Manual analysis, if the dependent analysis needs rerunning:
~/.venvs/mcfost-v11/bin/python -B \
  "$production_run/code/analyze_extinction_ice_production.py" "$production_run"
```

Production output includes `results/summary.json`, `REVIEW.md`, rankings,
band predictions, excluded-model reasons, ice responses, controlled ice/mass
response diagnostics, and PNG/PDF figures. The final diagnostic records its
1% reference outcome and `investigation_closed`; closure is an operational
decision and never changes a failed/inconclusive result into a pass.

For an initial review download each run's **`results/` and `logs/`**, plus
`experiment.json`, `experiment.sha256`, `runtime_binding.json` where present,
`production_runtime_binding.json` for production, and
`machine.template.cluster.json`. For reproducible compact production
re-analysis also include `manifest.json`, `manifest.sha256`,
`production_experiment.json`, `inputs/`, `code/`, `README.md`,
`requirements-cluster.txt`, and per-model `measurements.json`, `status.json`
and frozen `*.para`/`wavelength.lambda` files. Preserve paths; omit raw FITS
only from the transfer, not from the cluster run. The analyzer verifies frozen
inputs and compact measurement identities without loading the raw images.

No pass by this final one-wavelength diagnostic certifies all 19 production
bands, all physical models, or the full optical-material assumptions. The
production report carries those limits while answering the practical question
of which conditional responses the spectrum can separate.
