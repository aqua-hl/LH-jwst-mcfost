# Silicate opacity workspace

Opened 2026-09-23 to do for silicate what the self-supplied H₂O 30 K table did
for ice: replace a default prescription with independently characterised optical
constants, and test whether that relieves the 9.7-µm deficit.

**Optical constants are never generated here.** Every table must be supplied by
a named source and recorded in [constants/MANIFEST.json](constants/MANIFEST.json)
with its hash, exactly as the H₂O table was. Alongside the historical Draine
baseline, `laboratory_silicates_v1/` now supplies two Dorschner et al. (1995)
amorphous-glass alternatives: olivine and pyroxene compositions, each tabulated
over 0.2–500 µm. Their raw/converted hashes and numerical rows were checked
against the supplied manifest and raw tables. They are candidate inputs for a
later composition comparison; the four-model size pilot uses Draine.

The package's H₂O file is byte-identical to our original short table, ending at
19.940540101 µm. It does not provide the missing UV/far-IR extension. The
included olivine/H₂O coated-Mie initialization failed with SIGSEGV; that is a
different grain prescription from the selected separate-DHS pilot, so it does
not establish that the pilot fails or diagnose the cause. Its 206 zero-k rows
remain unchanged.

The incoming README describes a `mcfost_grid.laboratory_silicates` module,
`grid.lab_*.example.json` configs and material-selection flags from another
checkout. Those source/config files are absent here; the package's reproduction
commands cannot be used in this checkout as written. The local commands are in
[FOUR_MODEL_PILOT.md](FOUR_MODEL_PILOT.md).

What can be done without a new table is to screen the levers in opacity space.
That screen has now run, and it changes what the next step should be.

The four-case implementation and cluster commands are in
[FOUR_MODEL_PILOT.md](FOUR_MODEL_PILOT.md). The selected inclinations are 50°
and 70°, with the near-IR-preferred v02 supplied-H2O DHS prescription. Code is
ready for preparation using the unchanged archived H2O table. The user has
accepted MCFOST's native out-of-range handling for this bounded pilot; the
assumption is frozen in `constants/H2O_30K_historical.provenance.json`.

## The target

For a simple screen normalized to a fixed near-infrared optical depth, the
relative silicate optical depth scales with

    S = κ_abs(9.7 µm) / κ_ext(2.2 µm)

the 9.7-µm opacity per unit near-infrared extinction. Baseline value, at the
production size distribution (0.03–0.4 µm, q = 2.75, compact):

    S = 1.967        feature contrast κ(9.7)/κ(shoulders) = 2.27

Reproduce with `python -B scripts/silicate_opacity.py`. The calculator is Mie
(Bohren & Huffman) with optional Bruggeman vacuum inclusion; it self-tests
against the Rayleigh and geometric-optics limits on every run.

## What the screen found

Grain growth is a large lever, and it saturates almost immediately:

| amax (µm) | 0.4 | 1.0 | 3.0 | 10 | 30 |
|---|---:|---:|---:|---:|---:|
| S / S(baseline) | 1.00 | **0.26** | 0.27 | 0.28 | 0.27 |

Going from 0.4 to 1 µm drops S by a factor 3.9. This is mostly a larger
2.2-µm extinction opacity: it rises from 1518.84 to 5948.34 cm²/g, while
9.7-µm absorption rises from 2987.79 to 3016.54 cm²/g (+0.96%). At fixed
dust mass this does **not** predict a factor-3.9 reduction in feature optical
depth. Keeping NIR extinction fixed instead would require reducing the mass
by that factor. Beyond 1 µm S changes little. At the baseline size, porosity makes S
**worse** (1.26× at 30% vacuum, 1.74× at 60%), because vacuum inclusion lowers
the near-infrared extinction faster than it lowers the 9.7-µm absorption.

## Why this argues against launching a production run yet

A factor-3.9 lever exists in opacity space. The completed catalogue says it does
almost nothing to the observed residual: regressing the 9.7-µm residual across
all 142 models gives **+0.03 dex per dex of amax**, against −2.69 for column.

The screen and regression use different comparisons: S uses NIR-normalized
column, while the regression conditions on mass and other sampled parameters.
They are not contradictory evidence for a unique controlling mechanism.
Temperature stratification and scattering can affect the emergent feature,
but those effects cannot be isolated from these summaries alone. The bare-Mie
screen also differs from the ice-bearing grain prescription in the runs.

There is a specific coverage gap: the catalogue pairs one amax with each
column, and **no model combines
amax ≥ 1 µm with the near-infrared-preferred column of 2.25 × 10⁻⁴ M☉**. The
+0.03 coefficient is therefore an inference from a design that never tested the
combination the opacity screen points at. (The design is well decorrelated —
correlation between log amax and log mass is 0.07, variance inflation 1.01 — so
the coefficient is estimable, but it rests on no direct contrast.)

## The next step, before any production run

Run **four models**: the near-infrared-preferred column, 2.25 × 10⁻⁴ M☉, at
amax = 0.4 and 1.0 µm, at 50° and 70°, with fresh temperatures and the
existing Draine constants. The user-selected v02 H2O prescription is retained
with native extrapolation explicitly accepted. Only silicate grain size
changes in each pair; each case recomputes temperature under that assumption.

- Measure the changes in NIR, ice, silicate-core and shoulder band residuals
  together. Improvement demonstrates sensitivity to this grain-size change
  under the fixed background assumptions.
- A null response does not rule out a different silicate composition. The
  screen's tabulated flux gains assume fixed NIR extinction and are not a
  quantitative prediction for these fixed-mass DHS models.

This supplies a direct size contrast missing from the larger catalogue, using
four models rather than 144. It informs the next experiment without declaring
material properties or envelope structure uniquely responsible.

## Optional full-range material replacement

For a future replacement that avoids native out-of-range extrapolation:

- Three columns (wavelength, n, k), with a bulk-density/sublimation-temperature
  header, in the same layout as `H2O_30K_Leiden_mcfost.dat`.
- Coverage across **0.1–3000 µm**, with the source of any added tails identified.
  This is an optional stricter input policy. The current pilot accepts the
  simulator's handling outside the unchanged H2O table's tabulated range.
- A stated provenance: the measurement or model the constants come from, the
  temperature they apply to, and any processing already applied.

## Files

- `constants/MANIFEST.json` — what is present, what is awaited, and the target.
- `screen_draine.json` — the full opacity screen, every variation.
- `../scripts/silicate_opacity.py` — the calculator.
- `FOUR_MODEL_PILOT.md` — selected design, supplied-file requirements and commands.
- `four_model_plan.json` — inspectable design without generated optical constants.

## Limits

- S is an opacity ratio. It does not predict emergent flux, which depends on
  geometry, temperature and radiative transfer. The four-model test measures
  one conditional response, not all of those effects separately.
- Mie theory assumes compact or homogeneously porous spheres. Distributions of
  hollow spheres and aggregates behave differently, particularly in feature
  shape, and are not covered by this screen.
- The screen varies one material. It says nothing about silicate composition —
  olivine against pyroxene, crystallinity, iron content — which is what a
  supplied table would actually change.
