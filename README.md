# TMC1A modelling workspace

This is the consolidated local workspace. The former loose `jwst_mcfost` tree is preserved in a verified archive; its removal status is recorded in [the migration receipt](archive/migration_receipt.json). New work should begin here, not in the archived ENS-Lyon workflow.

## New local / cluster workflow

[Workflow guide](docs/WORKFLOW.md) · [Editable parameter grid](config/grid.example.json) · [Laptop smoke configuration](config/grid.smoke.json) · [Cluster machine template](config/machine.cluster.example.json)

[Local test plot](runs/laptop_smoke_v4/results/spectrum_comparison.png) · [Validation and timings](validation/WORKFLOW_VALIDATION.md)

`workflow.py` prepares a frozen model catalogue, runs independent physical models locally or as a Slurm array, measures Method-2 image aperture fluxes once, and generates compact rankings and plots. No jobs have been submitted. The small laptop demonstration is not a scientific fit; the production example has not been run. The current physical template is the original coated-Mie disk/envelope model; the later DHS H₂O model needs a separate extension.

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
