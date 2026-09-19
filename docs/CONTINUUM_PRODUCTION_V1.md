# Twelve-model continuum search

Prepared after reviewing the downloaded `continuum_control_v1` results. This is
the bounded continuum search selected by the user, with **64 CPUs per model and
up to 12 models running concurrently**. It is prepared locally; submission takes
place on the cluster.

The [downloaded results review](../validation/continuum_production_v1/RESULTS_REVIEW.md)
now verifies all 12 ranked models and 108 passing predictions, including indices
8–10. The nominal control remains best at 0.088998 dex; repeated nominal fluxes
differ by up to 9.05%, so numerical repeatability remains open. The design and
submission instructions below are retained as the run record.

## Starting point and scope

The starting point is the archived nominal continuum solution. The 50° and 60°
inclinations follow the literature-motivated geometry range discussed in the
[restart report](../reports/restart_review_2026-09-18/REPORT.md); 70° retains the
old control. The other settings below are controls and exploratory directions,
not a collection of literature best estimates or confidence limits.

Dust and ice are included. All models keep the original envelope prescription:
coated Mie silicate grains with a **5% generic ice mantle by volume**, maximum
grain radius 0.4 µm, and the existing minimum grain radius. The grain size
exponent q is varied in two families. The disk dust and disk structure remain
those of the template. The 5% mantle volume is not the 4% H₂O species mass
fraction used in the later DHS experiments.

Stellar parameters remain 4000 K, 2.5 R☉ and 0.5 M☉, with stellar accretion
7.3011756855e-8 M☉/yr. Every model receives a fresh equilibrium temperature
calculation. The calculation distance stays 140 pc and the aperture measurement
is transformed to 147 pc, preserving the existing observation contract.

Only the nine strict continuum anchors in a 1-arcsec-radius aperture are scored.
The proposed addition of two points immediately before 10 µm was reviewed and
the user chose to retain nine anchors: the existing mask excludes 5.5–8.05 µm
for ice and 8.0–13.3 µm for silicate, covering the entire intervening range.
The [anchor review and spectrum plot](../validation/continuum_anchor_extension/REVIEW.md)
record this decision and the evaluated alternatives.
Ice abundance, dust shape/composition, the H₂O profile and the silicate trough
are not fitted in this batch. Heating and the DHS dust alternative are deferred.
This is a limited continuum search, not the full geometry/heating/dust pilot in
[GOAL.md Stage 1](../GOAL.md), and it cannot yield an ice abundance inference.

## Model catalogue

Each family has two cavity half-opening angles, measured from the polar axis.
Array indices are zero based; mass is **envelope dust mass**, not gas-plus-dust
mass. The q exponent describes the grain size distribution.

| Array indices | Family | Inclination | Dust mass (M☉) | q | Cavity half-openings |
|---|---|---:|---:|---:|---|
| 0, 1 | Nominal control and cavity near-tie | 70° | 2.25e-4 | 2.75 | 17.5°, 20° |
| 2, 3 | Lower-inclination response | 50° | 2.25e-4 | 2.75 | 17.5°, 20° |
| 4, 5 | Intermediate-inclination response | 60° | 2.25e-4 | 2.75 | 17.5°, 20° |
| 6, 7 | Archived alternative mass/q family | 70° | 2.50e-4 | 3.25 | 17.5°, 20° |
| 8, 9 | Lower-mass boundary direction | 70° | 2.00e-4 | 2.75 | 17.5°, 20° |
| 10, 11 | Lower-q boundary direction | 70° | 2.25e-4 | 2.50 | 17.5°, 20° |

The [grid configuration](../config/grid.continuum-production.json) is the editable
design; the [prepared manifest](../runs/continuum_production_v1/manifest.json)
freezes the actual parameters, inputs and code. Index 0 has the same physical
parameters as the completed control. Do not edit frozen scientific inputs;
use a new run name for any scientific change.

The batch contains **12 temperature solves and 108 image calculations**, with
128,000 photon packets and the same density/grain grids as the control. Each
model has seven 2401-pixel images and two 1201-pixel images. These calculations
use the image Method-2 backend, not the optional coeval SED backend.

The [preparation checks](../validation/continuum_production_v1/preparation_checks.json)
passed: 351 frozen input hashes, all 120 temperature/image parameter files,
nominal-control input equality, model/anchor counts, Slurm resources and shell
syntax. These are preparation checks; MCFOST execution occurs on the cluster.

## What the control established

The [control review](../validation/continuum_control_v1/REVIEW.md) verified the
downloaded manifest identity, all nine reported quality passes, prediction
metadata, residuals and score. One model ranked, zero were excluded. Its regional
log RMS is **0.08687054 dex**. Largest fractional residuals are +30.38% at
13.80 µm and −28.02% at 5.39 µm.

Nominal fluxes differ from the archived nominal by up to **14.46%**, with seven
of nine anchors outside the provisional 2% reproduction target. Archive parity
and numerical convergence remain unresolved. The lower score is not evidence
of improved physics. Only the results directory was downloaded, so runtime,
memory, binary/utilities hashes and raw-image photometry could not be checked.

Rank this new batch on its common measurement contract. Compare its index-0
control with `continuum_control_v1` separately to assess repeatability; do not
merge archived and new scores into one ranking. Small score differences must be
interpreted against numerical uncertainty. A lower-inclination model failing at
fixed mass/dust/heating does not rule out that inclination after refitting.

## Resources and submission

[machine.slurm-production.json](../config/machine.slurm-production.json) requests:

- 12 array tasks, `--array=0-11%12`;
- one OpenMP process per task, `--cpus-per-task=64`;
- 16 GB and 12 hours per task; each MCFOST command has a one-hour timeout;
- at most 768 allocated CPUs and 192 GB for the model array;
- one dependent analysis job using 1 CPU, 2 GB and 10 minutes.

The scheduler may run fewer than twelve tasks at once. Resource limits inherit
the successful control setup and are not measured runtime estimates. Keep the
working portable MCFOST build and Python environment. Export the same
`MCFOST_UTILS`, and copy any required partition/account/reservation/module
settings from the working cluster machine configuration into the new one.

Transfer the updated workspace, including `runs/continuum_production_v1`, to the
cluster. From the workspace root on the cluster, activate the working environment
and run:

```bash
python -B workflow.py slurm runs/continuum_production_v1 --machine config/machine.slurm-production.json
bash runs/continuum_production_v1/submit.sh "$PWD/config/machine.slurm-production.json"
```

If you transferred only the source and configurations, first prepare the run
with `python -B workflow.py prepare config/grid.continuum-production.json`.
Preparation is needed only once and refuses to overwrite an existing run.
Regenerate Slurm scripts after changing machine resources. Use `submit.sh` so
the array and dependent analysis both receive the correct absolute paths.

## Completion and review

```bash
python -B workflow.py status runs/continuum_production_v1
cat runs/continuum_production_v1/results/summary.json
```

Expect 12 complete models with nine completed anchors each. After the analysis
job finishes, expect `ranked_model_count: 12`, `excluded_model_count: 0`, and 108
prediction rows. Check reported warnings and the spectrum/parameter-space plots.
The analysis runs after all array tasks end, including failures; a finished
analysis job alone does not establish a successful model batch.

The analysis automatically writes rankings, predictions, summary and plots into
the run's `results/`; it does not transfer them to the laptop. Download those
plus the run-level `runtime_binding.json`, each model's `runtime_binding.json`,
`status.json`, `measurements.json`, `temperature_complete.json`, command and
execution records, and simulator/Slurm logs. Keep the FITS products on the
cluster for independent photometry and convergence checks.

Do not rerun the ordinary analyzer locally on a results-only download: absent
per-model measurements would replace imported results with exclusions. The
dedicated control review script checks the downloaded tables without rewriting
them.

## Recovering the missing-Astropy failures

In array `10514458`, indices 8, 9 and 10 failed with
`ModuleNotFoundError: No module named 'astropy'`. The import is the first action
inside `run_model`, before MCFOST execution or creation of model status records.
This explains why those tasks may still appear as `not_started` in file-based
status. The error establishes that their Python could not import Astropy; it
does not establish which environment difference caused this on the three nodes.

Copy [retry_slurm_models.py](../scripts/retry_slurm_models.py) to the cluster.
Activate the shared Python environment used for the successful control. From the
workspace root, first verify its dependencies and interpreter:

```bash
python -c 'import sys, numpy, scipy, astropy, matplotlib; print(sys.executable); print(astropy.__version__)'
```

The interpreter and its packages must be accessible on compute nodes. If this
check fails, select the working environment before submitting. The helper does
not install or modify packages in an environment used by the nine running jobs.

Then submit only the three failed indices:

```bash
python -B scripts/retry_slurm_models.py runs/continuum_production_v1 --machine config/machine.slurm-production.json --indices 8,9,10 --submit
```

Without `--submit`, this command only prepares recovery files and checks local
imports; it does not query Slurm or submit anything. With `--submit`, it snapshots
this user's existing jobs named `continuum_production_v1`, submits three 64-CPU
tasks, and submits a final analysis that waits for both the retries and every
captured original job, including the original analysis. This prevents an older
analysis from overwriting the final report. The dependency uses Slurm's
[`afterany`](https://slurm.schedmd.com/sbatch.html) semantics so failures are
reported too.

The helper pins `sys.executable` without resolving a virtual-environment symlink,
checks frozen inputs, and checks the scientific imports again on each compute
node. New logs include the hostname, actual interpreter, package versions and
package paths. It writes separate scripts and a `recovery.json` receipt below
`runs/continuum_production_v1/recovery/`; the frozen model inputs, original
scripts and machine configuration stay unchanged. Regenerating a script alone
does not change a job already submitted to Slurm.

The login-node check cannot prove packages are visible on every compute node.
If a retry still fails at the import check, its new log gives the environment
details needed to fix that visibility. If final-analysis submission fails after
the retry array is accepted, the helper records the retry job ID and tells you
not to submit it again. Run analysis after all jobs for the run finish using the
same working Python. A fully recovered batch should rank 12 models and exclude
none.
