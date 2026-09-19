# Fixed-temperature pixel-scale diagnostic: Git and Slurm

This experiment investigates the anomalously low 2.546353-µm image in the
50° nominal model. It reuses **one existing temperature file, seed 42004**, at
**2,048,000 image packets**. Each image seed **42001–42005** produces a matched
pair whose only changed image setting is pixel count. The temperature,
physical model, field width, photon budget and image seed are fixed within
each pair:

| Image setting | Map width | Pixels | Linear pixel size |
| --- | ---: | ---: | ---: |
| Original | 6,000 AU | 2,401 × 2,401 | 2.499 AU |
| Finer | 6,000 AU | 4,801 × 4,801 | 1.250 AU |

There are **5 array tasks, 10 images and no new temperature solves**. Each task
runs its original and finer images sequentially and requests **64 CPUs and
96 GB**, with a 24-hour Slurm limit and four-hour timeout per MCFOST call.
The MCFOST memory limit is **64 GB for both images**. The larger grid has about
four times as many image pixels. Its allocation accounts for the image arrays
replicated across the 64 threads; this is not the earlier 16-GB job setting.
The allowed concurrency ceiling remains 16, while this five-task array runs
at most **5 tasks, using 320 CPUs**. A separate two-CPU, 8-GB analysis job runs
after the array, including when a task fails.

The field remains 6,000 AU at both pixel scales. The temperature-comparison
arm has been removed. This is a test of pixel-scale sensitivity at one fixed
temperature; five paired seeds measure its repeatability.

## Prepare and launch

Commit and push the staged implementation from the laptop. On the cluster:

```bash
cd ~/LH-jwst-mcfost
git pull --ff-only

~/.venvs/mcfost-v11/bin/python -B scripts/build_extinction_ice_crossover_v2.py \
  --source runs/extinction_ice_numerics_v2_git \
  --output runs/extinction_ice_pixel_scale_v2
```

`--source` must identify the completed numerical experiment on the cluster.
If that experiment is named `runs/extinction_ice_numerics_v2`, use that path
instead. The builder reads the existing 2.048M-packet, 50° model at seed 42004,
validates its provenance, and copies its frozen inputs and temperature file
into the new bundle. It needs that model's raw output, not just the numerical
CSV summaries. It does not require the seed-42001 temperature model or modify
the source run.

Prepare once: the builder refuses to replace an existing output directory.
Use a fresh `--output` path for a new experiment and use it consistently below.
The script names retain `crossover_v2` for compatibility, but the new experiment
has `schema_version: 2` and `diagnostic_id: fixed_temperature_pixel_scale_v2`.
The updated cluster helper rejects old 20-image temperature-crossover bundles.
Do not use a previously generated `submit.sh` from that superseded design.

Preparation requires the shared Python dependencies but launches no MCFOST or
Slurm job. It reads the source runtime's MCFOST paths into the machine template
and sets the larger memory allocation needed by this experiment.

If needed, set the working Slurm `partition`, `account`, `reservation` or module
`setup_lines` under `slurm` in the new run's `machine.template.json`. Then:

```bash
pixel_run=runs/extinction_ice_pixel_scale_v2

~/.venvs/mcfost-v11/bin/python -B "$pixel_run/code/configure_cluster.py" \
  "$pixel_run" --machine "$pixel_run/machine.template.json"

bash "$pixel_run/submit.sh" \
  "$PWD/$pixel_run/machine.template.cluster.json"
```

The first command checks imports and frozen inputs, pins the shared Python,
MCFOST and utilities, and generates the launch files. The second submits the
array and dependent analysis. The cluster machine argument must be the
absolute path printed by configuration. Every worker repeats the package and
environment checks before running. Each submission preserves both job IDs
under `submissions/`; if analysis submission alone fails, follow the printed
instruction to submit analysis separately instead of resubmitting the array.

## Status, validation and outputs

Replace `JOBID` with the array job ID:

```bash
squeue -r -j JOBID -o '%.22i %.10T %.10M %.25R'
sacct -j JOBID -X --format=JobID%22,State%20,NodeList%20,Elapsed,ExitCode

~/.venvs/mcfost-v11/bin/python -B "$pixel_run/code/crossover_task.py" \
  "$pixel_run" --status

~/.venvs/mcfost-v11/bin/python -B "$pixel_run/code/crossover_task.py" \
  "$pixel_run" --validate
```

File status is not a live scheduler query. Task logs are
`logs/array-JOBID_INDEX.out` and `.err`. The dependent analysis has separate
`logs/analysis-JOBID.out` and `.err` files. Inspect an error before retrying;
completed images are reused only after their cached products and inputs have
been checked. A resumed image is not an additional independent seed.

The analysis can also be rerun manually after all tasks finish:

```bash
~/.venvs/mcfost-v11/bin/python -B \
  "$pixel_run/code/analyze_extinction_ice_crossover_v2.py" "$pixel_run"
```

Download the new run's **`results/` and `logs/`**, plus `experiment.json`,
`experiment.sha256` and `runtime_binding.json`. Keep the complete run on the
cluster so raw image and temperature checks remain possible. The report
compares paired pixel-scale changes at fixed temperature and image seed,
including signed component aperture fluxes when available. Missing or failed
images remain explicit; the analyzer does not drop the anomalous seed to
improve a scatter estimate.

The results are `summary.json`, `seed_fluxes.csv`, `scatter.csv`,
`spatial_changes.csv`, `crossover_diagnostics.png` and `.pdf`, plus `REVIEW.md`.
There is no temperature-contrast report. Summary status is
`complete_diagnostic` when all images are complete and validated, or
`incomplete`/`integrity_error` when they are not.
`production_convergence_certified` remains false in every case. Read the
summary as well as the scheduler state.

## Scope and the earlier near-IR anchors

The frozen temperature is reused byte for byte. The material, inclination,
source-centred 1″ aperture, Gaussian PSF convention, model distance of 140 pc
and measurement distance of 147 pc remain those of the numerical test. The
existing aperture-v2 policy keeps whole-image edge and convolution losses as
diagnostic warnings while retaining strict finite flux, image geometry and
aperture/PSF support checks.

The `continuum_production_v1.1` manifest uses the same **6,000-AU,
2,401 × 2,401** image geometry for all five near-IR anchors:

| Anchor | Wavelength (µm) |
| --- | ---: |
| w001 | 1.202812 |
| w002 | 1.498796 |
| w003 | 1.758853 |
| w004 | 2.191666 |
| w005 | 2.546353 |

The w005 wavelength is exactly the failing numerical probe. Pixel-scale
sensitivity could therefore affect earlier near-IR parameter responses and
should be checked in representative production models. This one-model,
one-wavelength experiment cannot certify all five anchors or the earlier
parameter-response comparisons.

The numerical issue also has a bounded scientific interpretation. In the
audited 50° numerical comparison, the normal seed's direct thermal component
is 0.00492331 Jy out of 0.09627613 Jy, or **5.11% of the model flux**. The
earlier v1.1 nominal 70° model predicts 0.09009940 Jy against 0.12420141 Jy
observed at w005: **27.457% below the observation**. These fractions have
different denominators and refer to different inclinations. The measured
direct-thermal component alone is insufficient to bridge that earlier
0.03410201-Jy deficit. Its fraction is not a universal upper bound on all
image-integration errors or other models; scattered thermal emission also
changed in the audited image pair. Correcting numerical sampling therefore
does not establish a resolution of the model–data tension.

This diagnostic performs no observed-spectrum fitting, ice-material search,
broad MIRI-band scoring or automatic production certification. Completing the
10 images does not clear the earlier numerical accuracy gate.

The source, configuration, tests and this guide belong in Git. The existing
cluster experiment supplies the single raw temperature; do not add `runs/`,
raw FITS products or downloaded results to the commit.
