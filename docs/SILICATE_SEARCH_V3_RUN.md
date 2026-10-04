# Running silicate search v3 (512k)

This package implements `docs/SILICATE_SEARCH_V3.md`: 24 physical models and
two independent seed repeats, with 54 comparisons fixed before simulation.
It tests whether shallower envelope profiles and a near-infrared column
rematch reduce the 12–20 µm excess while preserving the near-infrared
continuum, ice band and 9.7 µm shape.

## Frozen design

| Setting | Values |
|---|---|
| Profile exponent, column factor | (−1.5, 1.00), (−1.25, 0.86), (−1.0, 0.83), (−1.0, 0.66) |
| Envelope silicate | Draine; laboratory pyroxene Mg50 |
| Carbon mass fraction | 15%, 20%, 25% |
| Inclination | 70° |
| Ice | Supplied H₂O 30 K; fraction derived from the fixed, material-specific ice-column target |
| Envelope grains | Separate DHS populations; 0.03–0.4 µm, exponent 2.75, 50 bins, vmax = 0.1 |
| Disk grains | Existing 70/30 mixture and Mie treatment, unchanged |
| Numerics | 512,000 packets for each temperature, SED and image setting; 1 AU pixels |
| Image geometry | 6001 × 6001 pixels across 6000 AU |
| Probes | 98 images per run, measuring 19 primary bands and three check bands |
| Seeds | 43001 for indices 0–23; 43002 for replicas 24 and 25 of indices 7 and 10 |
| Resources | 64 CPUs and 160 GB per task; at most 16 concurrent tasks; 24-hour limit |

Every run, including each seed repeat, solves fresh temperatures. `aperture_v2`
is retained, including the finite-positive signed-flux requirement. No images
are zero-clipped for scoring. The supplied ice table, density 0.94 g/cm³ and
MCFOST's existing treatment outside its wavelength coverage remain unchanged;
this campaign does not supply newly measured UV or far-infrared ice constants.

There is **no corrected MIRI reduction or collaborator update yet**, as
confirmed on 2026-10-04. The package retains the v2 frozen spectrum and
observation contract. Interpretations remain conditional on that reduction.
If a corrected reduction changes the 12–20 µm level by more than 0.1 dex,
recompute the residuals and projections and prepare a new frozen package
before launching it. The projections in the design document are motivations
for this experiment, not simulated v3 results.

## What happens on submission

1. Stage 0 uses two CPUs to measure four pure-species DHS opacity tables on the
   same 99 wavelengths as v2. It checks κ_ext(2.2), κ_abs(9.7) and κ_abs(18)
   against the recorded v2 values to relative tolerance 10⁻⁶; a zero reference
   must remain exactly zero. A mismatch stops the dependent array.
2. Only after Stage 0 passes does the builder derive and freeze all masses and
   ice fractions under `production/`. It uses the full-precision median ice
   columns from seven Draine and two pyroxene bracketed v2 crossings. Planning
   masses and rounded fractions are never substituted for measured opacities.
3. The 26-task array runs with `afterok` on Stage 0. Up to 16 tasks can run
   together, subject to cluster availability. This is 1024 CPUs at full use.
4. Analysis runs with `afterany` on the array and writes `results/`. Incomplete
   models and comparisons stay explicitly marked as missing.

The initial package intentionally has no production manifest or rendered model
parameters. All portable inputs, source code, reference values and the source
CSV hashes for the column/ice derivation are frozen in `bootstrap.json`.

## Launch on the cluster

After committing and pushing the prepared files locally, run:

```bash
cd "$HOME/LH-jwst-mcfost"
git pull
source "$HOME/.venvs/mcfost-v11/bin/activate"
export MCFOST_UTILS="$HOME/leshouches/mcfost/utils"
run="$HOME/LH-jwst-mcfost/runs/silicate_search_v3_512k"

python -B "$run/code/configure_silicate_search_v3.py" "$run" \
  --machine "$run/machine.template.json"
bash "$run/submit.sh" "$run/machine.template.cluster.json"
```

Configuration checks the shared virtual environment and pins its Python path,
MCFOST executable and utilities. It creates the submission scripts without
running MCFOST. `submit.sh` records all three job IDs in a new directory under
`submissions/`. Add any required Slurm account/partition to the machine template
before configuring. The absolute `run` path works even if the previous shell
was inside another run directory.

## Monitor and resume

Use the printed array ID in place of `ARRAY_ID`:

```bash
squeue -r -j ARRAY_ID -o "%.22i %.12T %.40R"
sacct -j ARRAY_ID -X --format=JobID%22,State%20,Elapsed,ExitCode
python -B "$run/code/silicate_search_v3_task.py" "$run" --status
```

File status is not a live scheduler query. Error logs are
`logs/stage0-JOB_ID.err` and `logs/array-ARRAY_ID_INDEX.err`. A Stage 0 failure
must be understood before retrying; it must not be bypassed with planning
opacities. For a stopped production array, fix the underlying error first.
Once its tasks have finished, resubmitting the complete chain is safe:

```bash
bash "$run/submit.sh" "$run/machine.template.cluster.json"
```

The pipeline revalidates and reuses completed work with matching receipts.
Retrying does not alter a seed, relax photometry checks, or repair a reproducible
numerical failure. To rerun analysis alone:

```bash
python -B "$run/code/silicate_search_v3_task.py" "$run" --analyze
```

## Results and interpretation

Download `results/` for the review, together with `stage0/receipt.json`, the
four compact Stage 0 opacity JSONs, `production/manifest.json`,
`production/production_experiment.json` and the Slurm logs. Keep raw images
and temperature products on the cluster until the review is complete.

`results/campaign_summary.json`, `model_success.csv`, `bridge_to_v2.csv`,
`seed_scatter.csv`, `ice_column_prediction.csv`, `nir_rematch.csv` and
`profile_outcomes.csv` expose the declared tests. The results also contain all
54 pair responses, nine diagnostic score scenarios and figures. Rankings do
not select the comparisons or replace the raw success criteria.

The Draine 20%-carbon anchor is compared with the v2 8%/12%-ice controls first.
Failure or absence of this bridge, or missing required runs, suppresses the
interpretive outcome conclusions. The seed checks report absolute paired
log-flux differences, not a covariance estimate or statistical error bar;
transferring those two checks to other configurations is an explicit
assumption. Any apparent p = −1 success still needs comparison with independent
envelope-mass constraints. The inner-flattened two-zone contingency is a
discussion outcome only and is not submitted by this package.

## Local preparation and verification

The package has already been prepared for this campaign. To build an additional
copy locally, with the returned v1/v2 evidence available:

```bash
python -B scripts/build_silicate_search_v3.py --output /path/to/new/package
python -B /path/to/new/package/code/silicate_search_v3_task.py \
  /path/to/new/package --validate
```

The builder refuses to replace an existing package. Synthetic simulator outputs
are used only in temporary test directories; scientific masses are frozen on
the cluster from the real Stage 0 results.
