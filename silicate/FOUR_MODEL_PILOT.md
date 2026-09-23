# Four fresh-temperature silicate-size models

The requested diagnostic uses the **unchanged archived H2O30K table**, with
MCFOST's native treatment outside its tabulated wavelength range explicitly
accepted for this bounded comparison. A new full-range table is not a launch
prerequisite. Preparation freezes inputs and code; it does not simulate or
submit models. The 512k configuration uses a new run folder so any previously
prepared 2.048M bundle is preserved.

## Frozen design

| Array index | Inclination | Silicate amax | Total envelope dust mass |
|---|---:|---:|---:|
| 0 | 50 degrees | 0.4 um | 2.25e-4 solar masses |
| 1 | 50 degrees | 1.0 um | 2.25e-4 solar masses |
| 2 | 70 degrees | 0.4 um | 2.25e-4 solar masses |
| 3 | 70 degrees | 1.0 um | 2.25e-4 solar masses |

The archived near-IR-preferred v02 prescription is used: separate DHS species,
96% Draine silicate and 4% supplied H2O by mass, vmax=0.1, zero material
porosity, amin=0.03 um, q=2.75 and 50 size bins. Only **silicate** amax changes;
H2O amax stays 0.4 um. This deliberately differs from varying both grain
populations together in the old dust-axis search. The ice density assumption
is 0.94 g/cm3. The laboratory temperature of the ice constants is 30 K; this
does not impose a 30-K radiative-equilibrium temperature.

The total mass includes both species (2.16e-4 solar masses silicate and 9e-6
solar masses ice). This parameter is a dust mass, not a directly measured
line-of-sight column. The cavity half-angle stays 17.5 degrees; the star stays
4000 K, 2.5 solar radii and 0.5 solar masses, with stellar accretion rate
7.3011756855e-8 solar masses/year. Disk structure and its 70/30 silicate/carbon
Mie species are retained from `reference/parameters/ice_v02_dust.para`.

Each case computes a fresh equilibrium temperature at 512,000 packets,
then 98 images at the same packet count. No archived temperature is reused.
The density grid stays 100 x 70, with 20 inner radial cells. Images use
6001 x 6001 pixels over 6000 AU, seed 43001, image Method 2 and signed total-I
aperture photometry. Model distance is 140 pc and comparison distance 147 pc;
the aperture radius is 1 arcsec. This is not a new convergence experiment or
a certification that every image is converged.

The 98 wavelengths integrate the same 19 broad bands as the last production:
nine near-IR continuum, four ice and six MIRI bands, with three additional
checks. The primary results are **unscaled flux and residual changes for each
predeclared inclination pair**. Inherited 45/30/25 continuum/ice/MIRI weighted
screening scores and nuisance-gain fits are secondary diagnostics. They do not
select which models are reported. Quadrature flags and incomplete pairs remain
visible. Edge-flux and whole-image convolution losses remain warnings under
`aperture_v2`; finite positive signed image and aperture fluxes and aperture
support remain mandatory.

## H2O input policy for this pilot

The user chose to retain the earlier H2O prescription rather than delay this
four-model test for new constants. The input is
`reference/dust/H2O_30K_Leiden_mcfost.dat` (SHA256
`29b572f8907993870d0b15823e1182a9c4813b8cdb8f4fbdd40ef39cf28fb00c`).
All 3000 wavelength/n/k rows, the density/header assumptions and 206 zero-k
samples are preserved exactly. No floor, splice or replacement material is
introduced. The acceptance is recorded in
`silicate/constants/H2O_30K_historical.provenance.json` and copied to the run as
`inputs/ice_input_provenance.json`.

The table covers 0.25004000002–19.940540101 um; the thermal grid is 0.1–3000 um.
[MCFOST documents native extrapolation of finite-range optical tables](https://mcfost.readthedocs.io/en/latest/input.html).
This pilot accepts the running executable's treatment outside that range and
binds the executable hash at runtime. Its exact version-specific extrapolation
law is not independently certified here. The input receipt truthfully records
`coverage_passed=false`, `policy_accepted=true` and
`native_extrapolation_required=true`. These describe the accepted approximation,
not a completed numerical or physical validation.

The earlier successful H2O grids reused original-dust temperatures, as recorded
in `reports/restart_review_2026-09-18/ice_audit.md`. This pilot still computes
fresh temperatures. It tests silicate grain-size response conditional on the
shared ice treatment; that assumption need not cancel between different grain
sizes. Simulator failures and mandatory photometry failures remain failures.
The unrelated olivine/H2O coated-Mie crash is not a demonstrated failure of
the selected Draine/separate-DHS prescription.

An explicitly supplied full-range table remains an optional alternative via
`--ice-constants` and `--ice-provenance`. That policy retains its checks for
0.1–3000 um coverage, unchanged laboratory samples and documented tails. The
extension provenance example is for that optional mode, not a required input
for the default pilot. No thermal-tail investigation is added to this task.

No new silicate constants are requested. The existing Draine table ends at
1000 um: its file and the simulator's treatment outside that range through
3000 um are retained as an explicit baseline assumption. This pilot adds no
silicate tail and does not certify tail sensitivity. Ice survival/sublimation
is not newly modeled; the existing fixed-composition prescription is retained.

## Preparation and cluster commands

The design can be inspected now without any supplied file:

```bash
python -B scripts/build_silicate_size_pilot.py --plan
```

Commit/pull the staged source files and historical provenance, then prepare on
the cluster with the existing table:

```bash
cd ~/LH-jwst-mcfost
~/.venvs/mcfost-v11/bin/python -B scripts/build_silicate_size_pilot.py
```

Preparation refuses a missing or changed historical H2O table and missing or
inconsistent provenance before creating a run. The accepted short wavelength
coverage is recorded in the receipt rather than treated as a launch failure.
It also refuses to overwrite an existing run. It freezes all required code,
optical tables, parameter files and observations under
`runs/silicate_size_pilot_v1_512k`. No previous production results are required.

Set `MCFOST_UTILS` and ensure `mcfost` is the working cluster build. If needed,
edit partition/account/reservation in the new run's `machine.template.json`
before configuration. Then:

```bash
~/.venvs/mcfost-v11/bin/python -B \
  runs/silicate_size_pilot_v1_512k/code/configure_cluster.py \
  runs/silicate_size_pilot_v1_512k \
  --machine runs/silicate_size_pilot_v1_512k/machine.template.json

bash runs/silicate_size_pilot_v1_512k/submit.sh \
  "$PWD/runs/silicate_size_pilot_v1_512k/machine.template.cluster.json"
```

The array is `0-3%4`: 64 CPUs and 160 GB per model, up to 256 CPUs total.
The MCFOST memory setting is 112 GB. Each model has a 24-hour job limit and
each simulator command a four-hour timeout. Existing identical completed
stages can resume after interruption. Four models at 512k use
(4 x 512k)/(144 x 2048k), about 0.7% of the last production's nominal
temperature/image packet workload; this is not a
wall-clock prediction for DHS grains or difficult ray tracing.

An analysis job is submitted with `afterany`, so failures appear as incomplete
pairs rather than disappearing. Re-run the analysis after any recovery:

```bash
~/.venvs/mcfost-v11/bin/python -B \
  runs/silicate_size_pilot_v1_512k/code/production_task.py \
  runs/silicate_size_pilot_v1_512k --status

~/.venvs/mcfost-v11/bin/python -B \
  runs/silicate_size_pilot_v1_512k/code/analyze_silicate_size_pilot.py \
  runs/silicate_size_pilot_v1_512k
```

Download `results/`, `experiment.json`, `manifest.json`, input provenance and
the logs. The paired analyzer writes `pilot_summary.json`,
`paired_band_responses.csv`, `pilot_plot.png`, `pilot_plot.pdf`, and
`PILOT_REVIEW.md`, alongside validated per-model screening products.

## Interpretation of the opacity screen

The existing screen describes **bare Mie spheres**; these runs use a disjoint
DHS mixture. Moreover, its S = kappa_abs(9.7)/kappa_ext(2.2) ratio compares
9.7-um absorption relative to 2.2-um extinction at **fixed near-IR optical depth**, requiring a change in
dust mass. This pilot holds dust mass fixed.

For the screen's 0.4 to 1.0-um change, kappa_abs(9.7) rises 0.96%, while
kappa_ext(2.2) rises by a factor 3.916. S falls to 0.258 of its initial value.
Keeping near-IR extinction fixed would require lowering mass from 2.25e-4 to
5.745e-5 solar masses, which is not this experiment. The screen's advertised
0.32–0.97-dex gains therefore are not predictions for these four models.

Report the measured NIR, ice, silicate-core and shoulder responses, together
with the fresh-temperature effects. A response establishes sensitivity to this
specific silicate size change. A null response does not exclude different
silicate compositions, and neither outcome uniquely identifies opacity versus
envelope structure as the cause of the spectral residual.
