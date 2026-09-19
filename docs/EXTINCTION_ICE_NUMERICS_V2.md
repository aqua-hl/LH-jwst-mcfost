# V2 numerical test: Git and Slurm

This package runs the first numerical experiment in the V2 design. It is not
the 32-model extinction search or the 24-case material pilot.

The array has **20 tasks**: two nominal inclinations (50°, 70°), five seeds
(42001–42005), and two packet budgets (512,000 and 2,048,000). Each task computes
a fresh temperature followed by six Method-2 images at 1.202812, 2.546353,
3.076480, 3.296169, 9.7 and 18.5 µm: **20 temperatures and 120 images** in total.
Temperature and image packet counts change together. Ten independent frozen
two-model subruns keep seeds/budgets separate without altering existing runs.

It requests **64 CPUs and 16 GB per task**, at most **16 tasks** (1,024 CPUs),
with a 24-hour Slurm limit and a four-hour timeout per MCFOST invocation. These
are initial bounds, not measured runtimes. Keep 64 threads across replicas;
lowering array concurrency to suit the allocation is allowed. A dependent
analysis job requests two CPUs, 4 GB and one hour, and runs after the array even
when some tasks fail.

## Pull, prepare and launch

Commit and push the staged numerical-test files from the laptop. On the
cluster, pull them into the existing project and prepare the run directly:

```bash
cd ~/LH-jwst-mcfost
git pull --ff-only
source ~/.venvs/mcfost-v11/bin/activate
export MCFOST_UTILS=/home/hlo/leshouches/mcfost/utils

python -B scripts/build_extinction_ice_numerical_package.py --no-archive
```

Preparation creates `runs/extinction_ice_numerics_v2/` from the versioned source
and reference inputs. It invokes neither MCFOST nor Slurm and requires no tar
file or downloaded results. Keep the known-working `mcfost` on PATH. The shared
Python environment, MCFOST and its full utilities must already be installed.

Prepare only once. If this run directory already exists from an earlier
upload, use a fresh path, for example
`--output runs/extinction_ice_numerics_v2_git --no-archive`, and use that path
in every following command. Do not overwrite an existing experiment.

If required, add your working Slurm `partition`, `account`, `reservation` or
module `setup_lines` under `slurm` in
`runs/extinction_ice_numerics_v2/machine.template.json`. Then configure and submit:

```bash
numerical_run=runs/extinction_ice_numerics_v2

python -B "$numerical_run/code/configure_cluster.py" \
  "$numerical_run" \
  --machine "$numerical_run/machine.template.json"

bash "$numerical_run/submit.sh" \
  "$PWD/$numerical_run/machine.template.cluster.json"
```

The configure command checks dependencies and frozen files, pins Python/MCFOST/
utility paths, and writes the cluster JSON and Slurm scripts. It submits
nothing. The last command submits the array and dependent numerical analysis.
Compute nodes repeat environment checks before execution. The selected Python
must be the shared virtual-environment interpreter.

Do not run `workflow.py prepare` on the V2 design JSON; use the builder above.
Do not use ordinary `workflow.py analyze` on these subruns; observational
ranking is explicitly disabled for numerical-only experiments.

### Files to commit for this numerical stage

Required new or modified execution files:

- `config/extinction_ice_identifiability_v2.design.json`
- `scripts/build_extinction_ice_numerical_package.py`
- `scripts/extinction_ice_numerical_task.py`
- `scripts/configure_extinction_ice_numerics_v2.py`
- `scripts/analyze_extinction_ice_numerics_v2.py`
- `src/mcfost_grid/configuration.py`
- `src/mcfost_grid/analysis.py`
- `docs/EXTINCTION_ICE_NUMERICS_V2.md` (also copied into the prepared run)

The three `tests/test_extinction_ice_numerical_package.py`,
`tests/test_extinction_ice_numerics_analysis.py` and
`tests/test_extinction_ice_numerics_cluster.py` files accompany these changes.

The existing tracked `workflow.py`, the rest of `src/mcfost_grid/`,
`requirements-cluster-recovery.txt`, `reference/parameters/continuum_original_dust.para`,
`reference/observations/continuum_sed_R100.ecsv`, `reference/dust/Draine_Si_sUV.dat`,
`reference/runtime_assets/Dust/{ac_opct,ice_opct}.dat` and
`reference/runtime_assets/Stellar_Spectra/lte4000-3.5.NextGen.fits.gz` supply the
remaining dependencies. No new `runs/`, `transfer/`, downloaded results or V2
broad-band audit products need to be committed for this numerical test. The
wider design JSON contains planning references; this builder reads its embedded
model parameters, seed values, photon budgets and wavelengths directly.

## Status and completion

Replace JOBID with the array ID printed at submission:

```bash
squeue -r -j JOBID -o '%.22i %.10T %.10M %.25R'
sacct -j JOBID -X --format=JobID%22,State%20,NodeList%20,Elapsed,ExitCode

python -B runs/extinction_ice_numerics_v2/code/numerical_task.py \
  runs/extinction_ice_numerics_v2 --status
```

Task logs are `logs/array-JOBID_INDEX.out` and `.err`; analysis has separate
logs. Submission IDs are retained under `submissions/`. File-based status is
not a live scheduler query. If a task fails, inspect its log before resubmitting.
Submitting the package again with the same pinned environment safely reuses
validated completed tasks and resumes incomplete ones; it does not create new
independent samples or relax quality checks.

The numerical analyzer can also be rerun manually after the array finishes:

```bash
python -B runs/extinction_ice_numerics_v2/code/analyze_extinction_ice_numerics_v2.py \
  runs/extinction_ice_numerics_v2
```

Download **`runs/extinction_ice_numerics_v2/results/`**, together with `logs/`,
`experiment.json` and `runtime_binding.json`, for the initial review. Keep the
full run directory on the cluster for raw-image validation and possible follow-up.
Outputs include:

- `numerical_summary.json`: completeness, provenance and numerical decisions.
- `seed_fluxes.csv`: every valid replicate's six aperture fluxes.
- `numerical_scatter.csv`: five-seed scatter and uncertainty by wavelength/case/budget.
- `budget_changes.csv`: paired photon-budget changes and uncertainty.
- `numerical_convergence.png` and `.pdf`: comparison figure.

A completed array is not automatically a numerical pass. At 2.048M packets,
the upper end of the two-sided 95% Gaussian standard-deviation interval,
divided by mean flux, must be ≤1% at every probe. For each paired fractional
change `(high − low)/low`, require `|mean| + 95% t-interval half-width ≤1%`.
Missing replicas, mixed runtimes, invalid products and unresolved bounds cannot
pass. Five identical seed fluxes are reported as degenerate and require review.
The analyzer writes reports even when the experiment is incomplete; an exit
code of 2 means it did not establish a pass. Check the JSON for the reason.

## Scope and preserved assumptions

All calculations use the original **generic coated-Mie dust**. H₂O30K constants
are not needed here. The measurement is a source-centred 1″ signed aperture,
with the inherited Gaussian PSF convention. Model distance remains 140 pc;
aperture/flux evaluation uses 147 pc. Aperture-v2 quality policy keeps harmless
whole-image boundary losses as warnings while enforcing finite flux and
discrete aperture/PSF support. It does not certify image-total/SED closure.

There are no observed flux/error placeholders and no fit ranking. These six
monochromatic probes test sampling noise; they do not implement integrated
MIRI bands, the offset continuum extraction, the 0.35″ H₂O/LSF operator,
material changes, dense SEDs or spatial-resolution convergence. The physical-
contrast criterion (numerical uncertainty below one-third of a claimed
response) remains for later matched physical tests. A pass is conditional on
these two models, six probes and fixed spatial/grain resolution.

Frozen source and product hashes, fresh-temperature receipts, model locks and
the original resumable runner are retained. The package adds a common runtime
identity across all ten subruns, including binary, utilities, thread count,
Python and scientific-package versions. Previously prepared experiments keep
their original frozen code and results.

## Optional archive

The builder can still create a transfer archive when explicitly wanted:

```bash
python -B scripts/build_extinction_ice_numerical_package.py
```

The Git workflow above uses `--no-archive`. Both modes refuse to replace an
existing prepared run. No MCFOST calculation or Slurm submission is performed
during preparation or package validation.
