# Running silicate search v2

This package implements `docs/SILICATE_SEARCH_V2.md`: 34 unique physical models,
48 predeclared comparisons, 98 images per model, and 512,000 photon packets for
temperature, SED and image calculations. Each model gets fresh temperatures.
Production uses 64 CPUs and 160 GB per task, with at most 16 simultaneous tasks.
No new optical constants are needed.

## What is frozen, and when

The initial package freezes the optical constants, observations, workflow and
34 composition/geometry cases. It deliberately contains no production masses
or rendered model files. On the cluster, Stage 0 measures five pure-species
DHS opacities (Draine, pyroxene, olivine, carbon and supplied H2O 30 K), on the
98 production wavelengths plus 2.2, 9.7 and 18 microns: 99 distinct wavelengths.
The wavelength file has no count header. No zero-k stress samples are added.

All five calculations must produce valid finite opacities and albedos. Their
receipts, FITS products and simulator identity are checked before the builder
uses their mass-weighted extinction at 2.2 microns to calculate matched envelope
masses. It then freezes `production/manifest.json`, `production/experiment.json`
and the 34 models. Compact-Mie estimates in the design are never run inputs.

The production array depends on successful Stage 0 completion. Analysis runs
automatically after the array ends, including after a partial failure.
Submission IDs are saved under `submissions/`.

## Launch on the cluster

After committing and pushing the staged files locally, run:

```bash
cd "$HOME/LH-jwst-mcfost"
git pull
source "$HOME/.venvs/mcfost-v11/bin/activate"
export MCFOST_UTILS="$HOME/leshouches/mcfost/utils"
run="$HOME/LH-jwst-mcfost/runs/silicate_search_v2_512k"

python -B "$run/code/configure_silicate_search_v2.py" "$run" \
  --machine "$run/machine.template.json"
bash "$run/submit.sh" "$run/machine.template.cluster.json"
```

Configuration checks and pins the active shared virtual environment, MCFOST
binary, utilities and frozen package; it performs no simulation and submits
no jobs. If needed, set the Slurm account/partition in `machine.template.json`
before configuration. Do not change frozen inputs or an in-use cluster machine
file. The package includes every required scientific table and Python source;
the cluster still provides MCFOST, its utilities, and the existing Python venv.

## Check progress and failures

Use the production array job ID printed at submission:

```bash
squeue -r -j ARRAY_JOB_ID -o "%.22i %.12T %.40R"
sacct -j ARRAY_JOB_ID -X --format=JobID%22,State%20,Elapsed,ExitCode
python -B "$run/code/silicate_search_v2_task.py" "$run" --status
```

An array waiting on `(Dependency)` while Stage 0 runs is expected. If Stage 0
fails, inspect `logs/stage0-JOB_ID.err` and its attempt's `mcfost.log` under
`stage0/`. Production will not start. For a model failure, inspect
`logs/array-ARRAY_JOB_ID_INDEX.err` and
`production/models/MODEL_ID/status.json`.

After resolving a failure and ensuring no prior tasks are still running,
the same `bash "$run/submit.sh" "$run/machine.template.cluster.json"`
command safely resubmits the workflow. Valid Stage 0 outputs and completed
models are cached; unfinished work resumes using preserved attempts. This
does not discard or overwrite the previous campaign.

The inherited `aperture_v2` policy records image-boundary diagnostics without
rejecting a model solely for those diagnostics. Finite positive signed image
and aperture fluxes and adequate aperture support remain required. Negative
total-I images are failures, not clipped into acceptance.

## Results and interpretation

Download `results/`, `stage0/`, `production/manifest.json`,
`production/manifest.sha256`, `production/experiment.json`,
`production/experiment.sha256`, `production/production_experiment.json`, and
`production/inputs/stage0_*.json` for the compact review. Keeping the entire
package preserves raw images and complete provenance. Generated outputs stay
on the cluster across `git pull` because they are ignored by Git.

Key products in `results/`:

- `CAMPAIGN_REVIEW.md` and `campaign_summary.json`: completeness, all declared
  comparisons, and the joint success result.
- `model_success.csv`: every model and every success criterion.
- `opacity_matched_catalogue.csv`: actual frozen masses, composition, opacity
  ratios and analytic mass factors.
- `matched_band_responses.csv` and `matched_feature_shapes.csv`: all 48
  predeclared comparisons, with missing pairs explicitly reported.
- `matched_response_map` and `all_model_criteria` figures (PDF and PNG).
- `summary.json` and `ranking.csv`: the inherited nine weighting/calibration
  scenarios as descriptive diagnostics.

Success means **one model simultaneously** meets: absolute 9.7-micron
continuum-normalized log residual at most 0.10 dex; raw s05/s06 residuals each
at most 0.15 dex; raw h02/h03 residuals each at most 0.15 dex; and RMS of the nine
near-IR continuum log residuals below 0.20 dex. The observed normalized feature
is calculated from the frozen bands, rather than hard-coded from rounded text.
Normalized ice depths are also reported, while h04 is outside the declared
success test. No fitted amplitude or segment gain is applied to these criteria.

All 34 cases and 48 comparisons are declared before simulation; weighted
rankings never select which comparisons are reported. This is an exploratory
composition/structure test, not a unique abundance solution, a posterior, or
proof of a unique physical cause. Mass matching is analytic, not a guarantee
of identical optical depths through MCFOST's discrete geometry. The supplied
H2O table and its accepted native out-of-range treatment remain unchanged.
One seed does not estimate numerical covariance; MIRI calibration and coarse
band integration limitations remain visible in the analysis.

To repeat analysis after additional tasks finish:

```bash
python -B "$run/code/silicate_search_v2_task.py" "$run" --analyze
```
