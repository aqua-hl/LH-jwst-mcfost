# Continuum production v1.1

The user selected **96 models, 512,000 photon packets, and 16 concurrent model
tasks** for `continuum_production_v1.1`. Each task uses 64 CPUs. This broadens the
joint search; the best score among the earlier 12 settings is not evidence that
70° is preferred over literature geometry when other parameters can change.

## Search design

Six inclinations share each of 16 dust/heating settings, so every viewing angle can
be compared after changing mass, grains, cavity and heating together. There are
16 models at each inclination and 64 models in the 50–65° interval emphasized
by the [literature review](../reports/restart_review_2026-09-18/literature_audit.md).

| Searched quantity | Values or bounds | Sampling |
|---|---|---|
| Inclination from pole | 45°, 50°, 55°, 60°, 65°, 70° | Every physical setting at all six angles |
| Envelope dust mass | 7.5e-5–4.0e-4 M☉ | Logarithmic |
| Maximum envelope grain radius | 0.1–3.0 µm | Logarithmic |
| Grain size exponent q | 2.5–4.0 | Linear |
| Cavity half-opening from pole | 10–30° | Linear |
| Source heating target | 1.5–3.0 L☉ | Linear; changes stellar accretion |

Three settings preserve the exact nominal, cavity20, and mass2.5/q3.25/cavity17.5
controls. The other 13 settings form a deterministic Latin-hypercube design:
each searched axis uses all 13 evenly spaced levels in its linear/logarithmic
coordinate, including both endpoints. Of 2,048 candidate arrangements, the
generator selects the one with the largest minimum pairwise distance in the
unit cube. These are exploratory design bounds, not literature confidence
intervals or posterior samples.

Each six-model block uses inclination order **70°, 45°, 50°, 55°, 60°, 65°**.
Indices 0–5 are the nominal setting, 6–11 the cavity20 setting, 12–17 the
archived alternative, and 18–95 the broader search. Index 0 is the old nominal
physical model at the higher photon count. This design samples 16 distinct
combinations of the other five axes; it is not a dense six-dimensional grid.

The [explicit grid](../config/grid.continuum-production-v1.1.json),
[model catalogue](../validation/continuum_production_v1.1/model_catalogue.csv),
[design provenance](../validation/continuum_production_v1.1/design.json), and
[prepared manifest](../runs/continuum_production_v1.1/manifest.json) record the
settings. The generator is [build_continuum_grid_v11.py](../scripts/build_continuum_grid_v11.py).
Its design seed is separate from the simulator seed.
The [coverage plot](../validation/continuum_production_v1.1/coverageplot.png)
shows the sampled settings. The
[preparation checks](../validation/continuum_production_v1.1/preparation_checks.json)
record verification of all 2,703 frozen input hashes, exact historical controls,
parameter files and Slurm resources. All 71 workflow tests pass; they include
mocked simulator runs, not a real local execution of this production batch.

## Heating, dust and measurements

The photosphere remains 4000 K and 2.5 R☉, with stellar mass 0.5 M☉. Its
luminosity is approximately 1.4415 L☉. The generator converts each total source
heating target H into accretion rate using
`Mdot = (H - Lstar) * Rstar / (G * Mstar)`, with consistent physical units.
MCFOST adds accretion luminosity as `G * Mstar * Mdot / Rstar` according to its
[stellar-source implementation](https://github.com/cpinte/mcfost/blob/main/src/stars.f90).
The exact old accretion rate is retained for the three controls. H is a model
source-heating target; it is not an assertion that observed bolometric or
aperture luminosity equals that input. The current upstream source was checked;
the exact cluster build's reported luminosities should also be retained.

The disk structure/dust, envelope density law, stellar photosphere, and dust
composition remain fixed. Envelope grains use the original coated-Mie
silicate/generic-ice mixture, **5% ice mantle by volume**, and minimum radius
0.03 µm. Changing grain size and heating requires a fresh equilibrium temperature
for every model; all 96 models receive one. DHS/composition changes and ice
abundance fitting remain deferred.

The nine existing continuum anchors remain unchanged, including the 1-arcsec
aperture, 140-pc simulation distance, 147-pc measurement distance and existing
PSF/aperture operator. No anchors enter the excluded ice/silicate ranges. Nine
fluxes alone do not remove degeneracies among all searched and fixed parameters.

## Photon counts and numerical controls

Temperature and image photon settings both increase from 128,000 to **512,000**;
the saved SED photon setting is also 512,000, though this backend measures images.
Density resolution, grain bins and wavelength-dependent image sizes are unchanged.
The full batch has **96 fresh temperature solves and 864 image calculations**.

A fixed `numerics.random_seed: 41001` is now passed to temperature and image
commands. MCFOST's [`-seed` option](https://mcfost.readthedocs.io/en/latest/options.html)
also nests output in `seed=41001`; the workflow handles that layout. The option
and photon settings are frozen in the manifest/configuration, and invocation
records preserve the actual command arguments. Old runs without this option
retain their existing behavior.

More packets address a plausible source of the previous 9.05% repeat difference.
MCFOST estimates temperature and scattering source functions through photon
sampling, then ray-traces images from those estimates
([MCFOST overview](https://mcfost.readthedocs.io/en/latest/overview.html)). They do
not establish convergence by themselves. The preceding cluster runs did not
explicitly fix a seed, so comparing them with this new control does not isolate
the photon-count effect. One shared seed is not a repeat ensemble or a guarantee
of identical results across builds/parallel executions. Repeat seeds and further
photon/resolution checks remain appropriate before interpreting small score gaps.

## Cluster resources and launch

[machine.slurm-production-v1.1.json](../config/machine.slurm-production-v1.1.json)
requests `--array=0-95%16`, one process and 64 CPUs per task, 16 GB per task,
and a 12-hour task limit. The command timeout is increased to four hours. At most
**1,024 CPUs and 256 GB** are requested concurrently by the model array.
The dependent analysis uses 2 CPUs, 4 GB and 30 minutes. These are resource bounds;
the downloaded previous results contained no timing or peak-memory records.
The photon workload is about 32 times the previous 12-model/128k batch, not a
measured wall-time estimate.

Transfer the updated workspace, including the prepared run and new helper, to
the cluster. Copy any required partition/account/reservation/module settings
from the successful cluster configuration into the new machine template.
Activate the working shared Python environment, put the portable `mcfost` build
on PATH, and export the working `MCFOST_UTILS`. From the workspace root **on the
cluster**, run:

```bash
python -B scripts/configure_slurm_environment.py runs/continuum_production_v1.1 --machine config/machine.slurm-production-v1.1.json
bash runs/continuum_production_v1.1/submit.sh "$PWD/config/machine.slurm-production-v1.1.cluster.json"
```

The first command validates imports/frozen inputs and writes a separate
`.cluster.json` with absolute Python, MCFOST and utilities paths. It regenerates
both Slurm scripts with compute-node dependency checks and environment logging.
This addresses the previous bare-Python/Astropy problem without changing the
working environment's packages. Interpreter and package files must be shared
with compute nodes. The template's placeholder Python deliberately requires
configuration before submission. The second command submits the array and its
dependent analysis. No scheduler is available in this local workspace.

If transferring only source/configurations, prepare the run first with
`python -B workflow.py prepare config/grid.continuum-production-v1.1.json`.
Do this only when the run directory does not already exist. Do not reprepare or
edit the frozen inputs of an existing run. A different scientific design needs
a new name; the environment helper also preserves an existing differing
cluster configuration because queued jobs may reference it.

## Completion and interpretation

Array `10514475` exposed differing system-Python packages between compute nodes.
The [recovery instructions](RECOVER_CONTINUUM_V1_1.md) set up a shared environment
and retry the failed tasks while preserving the running models and 16-task cap.

The later array `10514573` passed the shared-environment checks. Six models
(48–53) failed whole-image boundary checks; their aperture-invariance audit
passed. Use the separate [aperture_v2 recovery](APERTURE_V2_RECOVERY.md) to keep
the boundary warnings while resuming only the missing image calculations.

Use `python -B workflow.py status runs/continuum_production_v1.1`. After all model
tasks and the analysis finish, expect 96 complete/ranked models, nine completed
anchors each, 864 prediction rows and zero exclusions. Analysis runs after the
array ends even if some models fail, so inspect counts and warnings.

Compare matched inclination blocks and the best attainable score at each
inclination, retaining the parameters that changed with it. Report the best
tested settings and residual patterns; do not interpret a grid winner as a
unique solution or a rejection of literature geometry. Check boundaries and
numerical stability before adaptive refinement. No ice abundance or band-profile
constraint follows from these continuum anchors.

Download `results/` plus run/model runtime bindings, measurements, status and
temperature receipts, command/execution records, environment and Slurm logs.
Keep FITS products for photometry/convergence checks. Results-only downloads
support table review; ordinary local reanalysis requires the model measurements.
