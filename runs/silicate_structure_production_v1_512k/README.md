# Silicate mass, material and envelope-structure production search

Run: `runs/silicate_structure_production_v1_512k`.

This freezes the approved 60-case search following the completed four-model
size pilot. It tests mass–size compensation, alternative silicate constants,
and redistribution of envelope dust. All cases are declared before their
results. No pilot measurements or old temperature solutions are imported.

| Arm | Settings | Cases |
|---|---|---:|
| Mass and size | Draine; total envelope dust mass [0.6, 1.0, 1.5, 2.25] × 10⁻⁴ M☉; silicate amax [0.4, 1.0] µm; inclination [50, 60, 70] degrees | 24 |
| Material | Olivine and pyroxene; mass [1.0, 1.5, 2.25] × 10⁻⁴ M☉; amax 0.4 µm; the same inclinations | 18 |
| Density profile | Draine, amax 0.4 µm; volume-density exponent −1.25 or −1.75; three reference columns and the same inclinations | 18 |

The mass–size arm includes four physical counterparts of the pilot, recomputed
within this campaign's manifest and runtime binding. The comparison baseline
for each material or profile case is a declared mass–size model, not a fitted
winner. Material and profile alternatives are not crossed with each other;
this experiment does not exhaust their interactions.

## Column matching

The baseline envelope has rho(r) proportional to r⁻¹·⁵ between 1 and 3000 AU.
For exponent p, let J0(p) be the integral of r^p dr over those radii and J2(p)
the integral of r^(p+2) dr. The central radial envelope column is proportional
to M J0/J2. A changed profile therefore has

    M(p) = Mref × [J2(p)/J0(p)] / [J2(-1.5)/J0(-1.5)].

The conical cavity and its angular factor are identical in every comparison;
all three viewing angles lie outside the cavity. Reference masses are
[1.0, 1.5, 2.25] × 10⁻⁴ M☉ in the baseline −1.5 profile. The actual input masses
equal to reference mass multiplied by 3.60037382695 at p=−1.25 or by
0.2393567089536 at p=−1.75.
These are analytic smooth-envelope column matches, not a certification of
identical columns on MCFOST's finite density grid, total source optical depths,
scattered light paths or aperture-integrated extinction. The disk is unchanged.
Both reference and actual masses are recorded in the catalogue.

## Fixed inputs and numerical settings

- 512,000 packets for temperature, SED and each image; seed 43001.
- A fresh temperature calculation for every model.
- 6001 × 6001 pixels over 6000 AU (0.999833 AU/pixel).
- 98 image wavelengths integrated into 19 bands: nine near-IR continuum,
  four near-IR ice and six MIRI bands. The three legacy MIRI points are checks.
- A 1-arcsec aperture; model distance 140 pc, comparison distance 147 pc.
- Envelope grains remain separate, nonporous DHS species with vmax=0.1:
  96% silicate and 4% supplied H2O30K by mass. H2O amax remains 0.4 µm.
- Grain exponent 2.75, amin 0.03 µm, 50 bins, cavity half-opening 17.5 degrees,
  stellar/accretion heating and all disk parameters remain fixed.
- The supplied optical constants are copied unchanged. Olivine density is
  3.71 g/cm³ and pyroxene density 3.20 g/cm³; no density rescaling to Draine.
- Native treatment outside the material tables is an explicit modelling
  assumption. The short H2O table retains its original 206 zero-k rows.
  Its 30 K laboratory label is not the computed dust temperature. Ice survival
  is not temperature dependent in this template.
- `aperture_v2` quality policy: full-field boundary diagnostics are recorded;
  finite positive signed image/aperture flux and aperture support remain
  mandatory. Fluxes are never clipped to make a model pass.

There are 5,880 image requests plus 60 fresh-temperature solves. Production
requests 64 CPUs and 160 GB per task, at most 16 tasks (1,024 CPUs), with the
MCFOST memory setting 112 GB. A 24-hour task limit and 4-hour subprocess timeout
are retained. These are limits, not runtime predictions.

## Run on the cluster

The complete frozen run directory is versioned with its source changes. After
committing and pushing the staged files on the laptop, run on the cluster:

```bash
cd ~/LH-jwst-mcfost
git pull --ff-only
source ~/.venvs/mcfost-v11/bin/activate
export MCFOST_UTILS="$HOME/leshouches/mcfost/utils"
run="$PWD/runs/silicate_structure_production_v1_512k"

python -B "$run/code/configure_cluster.py" "$run" \
  --machine "$run/machine.template.json"
bash "$run/submit.sh" "$run/machine.template.cluster.json"
```

Configuration pins the Python environment, executable and utilities paths. It
does not simulate or submit. Submission first runs a bounded 2-CPU dust
initialization check on a compute node for the four distinct envelope
material/size prescriptions. Only successful checks release the 60-model array.
This check is not a source temperature solution or proof of physical accuracy
of extrapolated optical constants. If it fails, inspect its error/log files;
the array's failed dependency is not evidence of 60 independent model failures.

The material check runs in a fresh attempt directory with MCFOST's default
output root. In MCFOST 4.1.14, the dust-property writer uses `data_dust/`
relative to the working directory even when initialization was given a
different `-root_dir`. The preflight therefore does not override that root.
It also rejects FITSIO errors even if MCFOST exits with status zero. The
optional `g.fits.gz` file is checked when present; this version writes it only
for the HG scattering option, whereas this campaign uses the exact phase
function. The wavelength, opacity, albedo, grain-opacity, phase-function and
polarizability files remain required. These details were verified in the
[writer for the cluster's reported source revision](https://github.com/cpinte/mcfost/blob/af0dec17cc305df4cc77f033418fef00f508c395/src/dust_prop.f90#L1384).

The custom `check.lambda` input contains one wavelength per numeric line and
**no count header**. The pinned MCFOST
[wavelength reader](https://github.com/cpinte/mcfost/blob/af0dec17cc305df4cc77f033418fef00f508c395/src/input.f90#L555)
counts every numeric row as a sample. The previous preflight incorrectly
prepended the sample count, making MCFOST return one extra wavelength (for
example, a leading `866` became an unintended 866-µm sample). The writer and
regression-test simulator now follow this format. Output dimensions and
ordered wavelength values must still match the requested grid; mismatches
report both expected and actual dimensions. Earlier attempts remain intact.

The revised wavelength gate requires valid products at every production image
probe, all temperature bin centres read from the frozen parameter files, and
100 broad-coverage samples from 0.1 to 3000 µm. Both ideal and source-rounded
thermal centres are included; the custom dust-property wavelength table is
single precision, so this tests thermal coverage at table precision, not a
native temperature solve. Nonfinite/negative opacity and invalid albedo at any
required sample still stop the array. A required wavelength always takes
priority if it coincides with a diagnostic sample.

Exact zero-k knots from the supplied H2O table and their immediate neighbours
are additional stress diagnostics. Their failures outside the required grid
are recorded in `material_preflight/receipt.json` rather than blocking this
production grid. The reported Draine/0.4-µm preflight had 97 nonfinite integrated
opacity samples and 4,850 nonfinite grain-opacity entries, with no affected
production image probes. This is consistent with the unguarded `log(k)`
[interpolation in the pinned source](https://github.com/cpinte/mcfost/blob/af0dec17cc305df4cc77f033418fef00f508c395/src/dust_prop.f90#L349);
it does not certify the temperature grid or the other material prescriptions.
All four prescriptions must pass the revised required-wavelength checks on
the cluster. The optical constants are unchanged, no NaNs are replaced, and
the stress-test FITS/logs remain available. Passing this gate cannot certify
arbitrary new wavelengths or the physical accuracy of the native table tails.

The analysis job runs automatically after the array, including partial/failed
arrays. Submission prints the job IDs and saves them under `submissions/`.

```bash
squeue -u "$USER" -r -o "%.22i %.12T %.12M %.40R"
python -B "$run/code/production_task.py" "$run" --status
```

To rerun analysis without submitting simulations:

```bash
python -B "$run/code/analyze_silicate_structure_production.py" "$run"
```

For recovery, rerunning the same `submit.sh` command keeps valid completed
measurements and resumes missing work under the existing per-model lock and
runtime checks. Do not rebuild or edit a prepared run containing results.

The builder is available for a deliberately new destination:

```bash
python -B scripts/build_silicate_structure_production.py --plan
python -B scripts/build_silicate_structure_production.py --output runs/ANOTHER_NAME
```

## Analysis and results to return

Raw matched band responses and regional residuals are primary. The silicate
core is also compared to its shoulders; the two halves of the 3-µm ice band and
the 18–20-µm bands remain visible. The existing 45% near-IR continuum / 30% ice /
25% MIRI score and all nine weighting/calibration scenarios are retained as
descriptive screens. A better scalar score alone does not resolve the joint
near-IR/ice/silicate tension.

The analysis writes tables, figures and its review to `results/`. Return that
folder together with `manifest.json`, `manifest.sha256`, `experiment.json`,
`production_experiment.json`, the observation contract, runtime bindings,
material-check receipts and per-model `measurements.json`/`status.json`. Raw
FITS are only needed for a subsequent image-level audit, not routine analysis.

All integration warnings and missing comparisons are retained. A single seed
does not provide Monte Carlo covariance. Four-percent ice is fixed: varying
total mass changes absolute ice mass, so this campaign cannot measure ice
abundance. The proposed density and mass levels are exploratory, not literature
confidence intervals. No likelihood, abundance posterior or unique
opacity-versus-geometry verdict is claimed.
