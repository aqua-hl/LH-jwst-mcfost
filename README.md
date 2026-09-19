# TMC1A modelling workspace

This is the consolidated local workspace. The former loose `jwst_mcfost` tree is preserved in a verified archive; its removal status is recorded in [the migration receipt](archive/migration_receipt.json). New work should begin here, not in the archived ENS-Lyon workflow.

## New local / cluster workflow

[Workflow guide](docs/WORKFLOW.md) · [Editable parameter grid](config/grid.example.json) · [Laptop smoke configuration](config/grid.smoke.json) · [Cluster machine template](config/machine.cluster.example.json)

[Local test plot](runs/laptop_smoke_v4/results/spectrum_comparison.png) · [Validation and timings](validation/WORKFLOW_VALIDATION.md)

`workflow.py` prepares a frozen model catalogue, runs independent physical models locally or as a Slurm array, measures Method-2 image aperture fluxes once, and generates compact rankings and plots. The cluster smoke test, full-resolution continuum control and initial twelve-model search have completed. A broader 96-model search is prepared below. The small laptop demonstration is not a scientific fit. The current physical template is the original coated-Mie disk/envelope model; the later DHS H₂O model needs a separate extension.

```bash
python -B workflow.py --help
python -B workflow.py status runs/laptop_smoke_v4
python -B workflow.py analyze runs/laptop_smoke_v4 --workers 2
```

New source is in `src/mcfost_grid/`, tests in `tests/`, and each prepared experiment owns its manifest, code snapshot, model outputs and `results/` directory under `runs/`. Active code does not modify the historical reference bundle.

### Quick Slurm array test

To check scheduling and log creation on the cluster, submit the small Bash-only
[array test](scripts/slurm_smoke.sh) from this workspace:

```bash
sbatch scripts/slurm_smoke.sh
# If required by your cluster:
# sbatch --partition=YOUR_PARTITION --account=YOUR_ACCOUNT scripts/slurm_smoke.sh
```

It runs four tasks, at most two concurrently, each requesting one CPU, 128 MB and
a two-minute time limit. Each `slurm-smoke-JOBID_TASKID.out` log in the submission
directory should end with `SUCCESS: task N completed`. This checks basic Slurm
execution; MCFOST and Python are not required.

To test real MCFOST calculations, activate a Python environment with
`requirements-workflow.txt` installed, put `mcfost` on `PATH`, and export
`MCFOST_UTILS` pointing to the utilities directory containing `Dust`, `Lambda`
and `Stellar_Spectra`. Then run:

```bash
python -B workflow.py prepare config/grid.smoke.json
python -B workflow.py slurm runs/laptop_smoke_v4 --machine config/machine.slurm-smoke.json
bash runs/laptop_smoke_v4/submit.sh "$PWD/config/machine.slurm-smoke.json"
```

Run preparation only once for a given run name. This submits two coarse models
with two CPUs and 4 GB each, followed by an analysis job. The smoke machine
configuration uses your default partition/account; add these under `slurm` if
your cluster requires them. Check progress with
`python -B workflow.py status runs/laptop_smoke_v4` and inspect the run's `logs/`
and `results/` directories.

If the smoke array fails with `MCFOST startup failed` or `SIGILL`, compare startup
on the login and compute nodes using the same environment and partition/account:

```bash
bash scripts/mcfost_startup_check.sh > mcfost-login.out 2>&1
sbatch scripts/mcfost_startup_check.sh
```

Compare `mcfost-login.out` with `mcfost-startup-JOBID.out`. The diagnostic records
CPU capabilities, executable checksum, linked libraries, utilities visibility
and the result of a bounded `mcfost -help` call with automatic updates disabled.
It performs no model calculation. An executable built for unsupported CPU
instructions is one possible cause of `SIGILL`; the crash alone does not
identify the instruction or prove that cause. Select a compatible build or
rebuild for the compute-node CPUs after checking the diagnostics. The
[MCFOST build instructions](https://mcfost.readthedocs.io/en/latest/installation.html)
describe source installation; [GCC's CPU-target documentation](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html)
explains why builds using `-march=native` may not run on other CPUs.
The [upstream MCFOST Makefile](https://github.com/cpinte/mcfost/blob/master/src/Makefile)
uses `-march=native` for GNU and `-xHOST` for Intel in ordinary non-release
builds. Check your own Makefile and compiler before selecting portable rebuild
flags; a source build on a different CPU can therefore require rebuilding.

### One full-resolution continuum control

After the smoke test passes, [grid.control.json](config/grid.control.json)
defines one nominal coated-Mie model: inclination 70 degrees, envelope dust mass
2.25e-4 solar masses, maximum envelope grain radius 0.4 micrometres, grain-size
exponent 2.75 and cavity half-opening 17.5 degrees. It computes a fresh
temperature and all nine continuum images with 128,000 photon packets, the
original density/grain grid and the wavelength-dependent 2401/1201-pixel image
presets. The configuration limits the catalogue to one model.

Use the same environment as the successful smoke test. Copy any required
partition, account, reservation or module settings from your working smoke
machine configuration into
[machine.slurm-control.json](config/machine.slurm-control.json), then run on the
cluster:

```bash
python -B workflow.py prepare config/grid.control.json
python -B workflow.py slurm runs/continuum_control_v1 --machine config/machine.slurm-control.json
bash runs/continuum_control_v1/submit.sh "$PWD/config/machine.slurm-control.json"
```

Preparation is needed only once. The array contains one model task requesting
64 CPUs, 16 GB and a 12-hour time limit; each simulator command has a one-hour
limit. These are initial resource bounds, not a measured runtime estimate. A
dependent analysis job generates the spectrum plot and compact predictions.

Check `python -B workflow.py status runs/continuum_control_v1` for one complete
model with nine completed anchors. After analysis,
`runs/continuum_control_v1/results/summary.json` should report
`ranked_model_count: 1` and `excluded_model_count: 0`. The run retains simulator
timings and logs for sizing subsequent work. Compare its predictions with the
archived nominal control before treating the new environment as scientifically
validated; a single run does not establish convergence.

### Completed twelve-model continuum search

The [downloaded control review](validation/continuum_control_v1/REVIEW.md)
confirms nine passing anchors and a score of 0.08687054 dex. Nominal fluxes differ
from the archive by up to 14.46%, so reproduction and convergence remain open.

The [run plan and cluster commands](docs/CONTINUUM_PRODUCTION_V1.md) describe
the prepared `continuum_production_v1` search: 12 models, nine anchors each,
**64 CPUs per task and up to 12 concurrent tasks** (768 CPUs at peak). It starts
from archived continuum controls, adds 50°/60° inclinations and lower mass/q
directions, and retains the original dust with 5% ice mantle volume. Detailed
ice/silicate features and heating changes are deferred.

Use [grid.continuum-production.json](config/grid.continuum-production.json) and
[machine.slurm-production.json](config/machine.slurm-production.json). The
downloaded cluster results now include all 12 models and 108 passing predictions,
including the three previously failed tasks. The nominal model scores best
among these twelve settings at 0.088998 dex; this does not rule out literature
geometry with other physical parameters varied. Its repeated fluxes differ from
the earlier control by up to 9.05%, so close rankings still need numerical checks. See the
[production results review](validation/continuum_production_v1/RESULTS_REVIEW.md)
and [spectrum plot](runs/continuum_production_v1/results/spectrum_comparison.png).

### Next batch: continuum production v1.1

[continuum_production_v1.1](docs/CONTINUUM_PRODUCTION_V1_1.md) is prepared with
**96 models, 512,000 photons, nine anchors per model, 64 CPUs per task and
16 concurrent tasks** (up to 1,024 CPUs). Six inclinations (45°, 50°, 55°, 60°,
65°, 70°) each sample the same 16 combinations of envelope mass, maximum grain
size, grain exponent, cavity opening and source heating. Three settings retain
historical controls; thirteen explore the broader joint parameter space.
The [catalogue](validation/continuum_production_v1.1/model_catalogue.csv) and
[coverage plot](validation/continuum_production_v1.1/coverageplot.png) show the
design. Dust/ice composition and the continuum-only anchors stay fixed.

More photons address possible sampling noise; neither this increase nor a
fixed seed alone establishes convergence. The earlier winner is conditional on
the small grid. This run tests other physical settings at every inclination.

Transfer the updated workspace and prepared run to the cluster. Activate the
working shared Python environment, set MCFOST on PATH and export MCFOST_UTILS.
Copy required partition/account/module settings into the
[machine template](config/machine.slurm-production-v1.1.json), then run:

```bash
python -B scripts/configure_slurm_environment.py runs/continuum_production_v1.1 --machine config/machine.slurm-production-v1.1.json
bash runs/continuum_production_v1.1/submit.sh "$PWD/config/machine.slurm-production-v1.1.cluster.json"
```

The helper pins Python, MCFOST and utilities paths in a separate cluster
configuration, checks imports and frozen inputs, and adds checks to both Slurm
jobs. The submission script schedules the array and dependent analysis. No jobs
have been submitted from this local workspace. See the
[run plan](docs/CONTINUUM_PRODUCTION_V1_1.md) for completion checks and numerical
limitations.

For the missing-SciPy failures on worker097 in array `10514475`, follow the
[shared-environment recovery commands](docs/RECOVER_CONTINUUM_V1_1.md). They keep
the running models and retry indices 8–10 and 16–95 after the original array ends.

For the later **image-boundary failures in models 48–53**, the returned aperture
audit passed at floating-point precision. The [aperture_v2 recovery](docs/APERTURE_V2_RECOVERY.md)
creates a separate run, preserves boundary warnings and strict aperture checks,
and reuses the original temperatures/images. Expected remaining work is
16 images and no new temperature solves once the other 90 models are complete.

## Scientific report

- [Paper-style report — PDF](reports/restart_review_2026-09-18/REPORT.pdf)
- [Paper-style report — HTML](reports/restart_review_2026-09-18/REPORT.html)
- [Editable manuscript — Markdown](reports/restart_review_2026-09-18/REPORT.md)
- [Figures — vector PDF collection](reports/restart_review_2026-09-18/FIGURES.pdf)
- [Restart objective and decision criteria](GOAL.md)

The report contains six scientific figures and a separate administrative appendix. Individual plots are supplied as vector PDF and 300-dpi PNG files. All figures derive from preserved measurements and calculations; no new radiative-transfer models were run during consolidation.

## Contents

| Directory | Purpose |
|---|---|
| `reports/` | Manuscript, figures, scientific audits and the self-contained report generator |
| `reference/observations/` | Native/R100 spectra and the corrected H₂O observation contract |
| `reference/predictions/` | Canonical continuum rankings, ice and silicate diagnostics |
| `reference/parameters/` | Representative parameter/temperature snapshots and wavelength files |
| `reference/source/` | Frozen measurement kernels, dependencies, regression tests and Method-2 patches |
| `reference/dust/`, `reference/runtime_assets/` | Selected production-matching dust and stellar-spectrum resources |
| `reference/provenance/` | Model/task manifests and historical working notes |
| `config/` | Restart proposal, editable grid, laptop test and new-cluster machine templates |
| `src/`, `tests/`, `docs/` | New portable simulation/analysis workflow and quick regression checks |
| `runs/` | Explicitly prepared experiments, attempts, compact measurements and generated plots |
| `archive/` | Complete legacy snapshot, file-level integrity ledger and restoration instructions |
| `validation/` | Migration checks, dependency inventory and test results |

## Reproduce the report locally

With NumPy, SciPy, Astropy, Matplotlib and Pillow installed:

```bash
python -B prepare_reference_bundle.py
python -B reports/restart_review_2026-09-18/build_review.py
python -B -m unittest discover -s reference/source -p 'test_*.py' -v
```

The report does not require `jwst_mcfost`, raw cubes, a simulator or cluster access. [requirements-report.txt](requirements-report.txt) records the tested package versions. Reference bytes are deliberately unchanged; historical extraction and simulation scripts are not a newly ported production package. The copied analyser is included for kernel parity tests, not for execution against an invented directory layout.

## Preservation and transfer limitations

[MIGRATION.md](MIGRATION.md) explains verification, restoration and deletion. [MIGRATION_CODE_AUDIT.md](MIGRATION_CODE_AUDIT.md) records code/dependency findings. The archive includes every file that remained in the old tree after moving this report, including hidden scientific files, logs, failed-attempt evidence and the earlier cache archive.

The raw JWST cubes and supplied ice packages remain unchanged in the parent project's `cubes/` and `DATA/` directories. They are catalogued in [the external-input manifest](validation/external_inputs.json) and must be transferred separately for re-extraction. This folder alone is sufficient for report regeneration, not for a complete observational re-reduction.

Most production images and temperature solutions were never retained locally. Their absence predates consolidation and is documented in the report; a backup cannot recover them. A new simulation environment still requires a compatible MCFOST build and complete utility installation. The selected copied assets are not the full utility distribution. Historical manifests preserve original absolute paths as provenance; they are not active path settings for the new cluster.

No old job script should be submitted unchanged. Use the new workflow's generated array scripts only after editing the machine configuration for the new cluster. Results now live inside each run rather than in one shared output directory.
