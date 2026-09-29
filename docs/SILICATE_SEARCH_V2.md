# Silicate search v2: envelope composition at matched near-infrared extinction

Design prepared 2026-09-29, following the completed
[silicate structure campaign](SILICATE_STRUCTURE_PRODUCTION_V1.md)
(`runs/silicate_structure_production_v1_512k`, 60/60 models). This is a
design only. No run directory has been prepared and no jobs have been submitted.

The question is narrow: **what change to the envelope dust lets one model
reproduce the 9.7-µm silicate feature while keeping the near-infrared
continuum, the 3-µm ice band and the 18–20-µm bands?** The previous campaign
located the problem precisely. This one tests a specific, physically standard
answer to it, and separates that answer from the envelope-structure alternative.

## What the previous campaign established

The self-supplied H₂O 30 K ice in separate DHS grains largely fixed the 3-µm
band shape: 9 of 60 models fit both the 2.95- and 3.11-µm windows within
0.15 dex, against 0 of 142 with the generic MCFOST ice. That ice treatment is
retained here unchanged.

For the silicate, raw fluxes proved misleading, so the analysis used
continuum-normalised feature depths, which do not depend on the absolute flux
level. With Draine silicate, the baseline density profile (p = −1.5) and 4% ice,
each feature's depth scales close to linearly with envelope column:

| Feature | Observed core/continuum | Envelope dust mass matching that depth |
|---|---:|---:|
| Ice, 2.95 µm | 0.205 | ≈ 2.55 × 10⁻⁴ M☉ |
| Ice, 3.11 µm | 0.145 | ≈ 2.65 × 10⁻⁴ M☉ |
| Near-IR continuum level | — | ≈ 2.3 × 10⁻⁴ M☉ |
| Silicate, 9.7 µm | 0.463 | **≈ 0.93 × 10⁻⁴ M☉** |

The 9.7-µm feature wants about **2.8 times less silicate absorption per unit
near-infrared extinction** than the model supplies. How each previous lever
fared against that:

| Lever | Effect on the 9.7-µm feature | Cost |
|---|---|---|
| Silicate grain size at fixed mass | +5.6% in flux; κ_abs(9.7) barely changes | none, but no gain |
| Mass reduced with grain growth | ~10× brighter at 9.7 µm | near-IR and ice broken in all 9 pairs |
| Olivine or pyroxene (Dorschner, Mg50) | feature **deeper** at matched column | near-IR 10× off at 2.25 × 10⁻⁴ M☉ |
| Steeper density profile, p = −1.75 | shallower (0.170 → 0.264); ice depth unchanged | 18–20 µm overshoots ~1.5× |

Olivine models 27–29 matched the 9.7-µm *flux* to within ±0.1 dex only because
olivine raised the surrounding continuum; their feature was deeper than
Draine's. The steeper profile is the one lever that moved the feature in the
right direction without touching the ice, and it got about halfway.

## The hypothesis: carbon in the envelope, with the ice fraction free

The envelope is currently 96% silicate and 4% ice by mass, with **no carbon**.
The disk in the same model already uses 70/30 silicate/amorphous carbon, and
interstellar dust carries a comparable carbon fraction. Carbon is the natural
candidate because it adds near-infrared extinction with almost no 9.7-µm
feature. Compact-Mie mass opacities at the production size distribution
(0.03–0.4 µm, q = 2.75):

| Per gram | κ_ext(2.2 µm) | κ_abs(9.7 µm) | κ_abs(18 µm) / κ_abs(9.7 µm) |
|---|---:|---:|---:|
| Draine silicate | 1,519 | 2,988 | 0.357 |
| Amorphous carbon | 10,767 | 635 | — |
| H₂O 30 K ice | 198 | — | — |

The screening quantity is S = κ_abs(9.7 µm) / κ_ext(2.2 µm), the 9.7-µm
absorption per unit near-infrared extinction. With 20% carbon by mass, S falls to
**0.34–0.37 of the reference**, on the 0.36–0.38 target implied by the table
above. Solving for all three depth targets at once predicts an envelope of
approximately:

| Silicate | Carbon | Ice | Total envelope dust |
|---:|---:|---:|---:|
| 76% | 16% | 9% | 1.2 × 10⁻⁴ M☉ |

Carbon must be accompanied by a higher ice fraction. At matched near-infrared
extinction the total dust mass falls, so a fixed 4% ice fraction would thin the
ice band. The ice column the 3-µm band wants corresponds to roughly
1.0 × 10⁻⁵ M☉ of ice, which at 16–20% carbon is 9–10% of the envelope dust.

These are compact-Mie estimates. The production grains are DHS with
vmax = 0.1, which changes the absolute opacities, particularly for carbon;
Stage 0 below replaces them before any model is built.

### A risk the design is built to test

Carbon dilutes the 18-µm silicate feature by the same factor as the 9.7-µm one:
κ_abs(18)/κ_abs(9.7) moves only from 0.363 to 0.389 across 0–30% carbon.
Cutting the silicate absorption 2.6-fold at fixed temperature structure would,
crudely, brighten 18–20 µm by about 0.17 dex — the same limit the steeper
profile reached. The silicate material sets that ratio (Draine 0.357, pyroxene
0.493, olivine 0.581, compact Mie), so Stage 2 tests whether a material with a
relatively stronger 18-µm band can hold 18–20 µm down while carbon handles the
9.7-µm depth.

## Comparison rule: match near-infrared extinction, not dust mass

Every composition change is compared at **matched near-infrared extinction**.
The near-infrared continuum fixes the column; comparing compositions at equal
dust mass changes composition and column together. That is what confounded the
previous material arm. Iron-bearing olivine has far more near-infrared opacity
per gram than Draine silicate, so at equal mass its near-infrared came out 10×
too faint and its feature appeared deeper. Its compact-Mie S is in fact lower
than Draine's (1.45 against 1.97), so at matched extinction it may behave
quite differently.

For separate grain populations the mass opacity of a mixture is the
mass-weighted sum of the species opacities, so

    κ_ext,mix(2.2) = f_sil κ_sil + f_C κ_C + f_ice κ_ice
    M(composition) = M_ref × κ_ext,ref(2.2) / κ_ext,mix(2.2)

with the reference composition 96% Draine / 0% carbon / 4% ice at
M_ref = 2.25 × 10⁻⁴ M☉ and p = −1.5 — the configuration that fits the
near-infrared continuum and both ice windows in the previous campaign. For a
steeper profile the composition-matched mass is further multiplied by the
analytic radial-column factor of the previous campaign, 0.2393567089536 at
p = −1.75. Both matchings are analytic; neither certifies identical optical
depths on MCFOST's discrete, cavity-bearing grid.

## Design

### Stage 0 — DHS opacities (no radiative transfer)

Compute DHS mass opacities for each of the five species on its own — Draine,
pyroxene and olivine silicate, amorphous carbon, H₂O 30 K ice — at the 98
production image wavelengths plus 2.2, 9.7 and 18 µm, using the same grain
parameters as the models. Because the populations are separate, every mixture
follows from these five by mass weighting, so five dust-property calculations
cover the whole design.

Stage 0 produces the matched masses for Stages 1 and 2 and freezes them into
the catalogue before any model is built. It selects nothing by fit. It uses the
dust-initialisation check already in the production pipeline and inherits its
gate: finite, non-negative opacities and valid albedos at every production
probe.

Stage 0 also reports S and κ_abs(18)/κ_abs(9.7) for every Stage 1 and 2
composition, so the compact-Mie predictions above can be checked against the
actual grain treatment before the grid runs.

### Stage 1 — composition grid

Draine silicate, p = −1.5, mass matched to near-infrared extinction.

| Axis | Levels |
|---|---|
| Carbon mass fraction | 0, 10, 20, 30% |
| Ice mass fraction | 4, 8, 12% |
| Inclination | 60°, 70° |

**24 models.** Silicate takes the remaining mass fraction. Indicative
compact-Mie matched masses, which Stage 0 supersedes:

| Carbon | Ice | Silicate | Total mass (M☉) | Ice mass (M☉) | S / S_ref |
|---:|---:|---:|---:|---:|---:|
| 0% | 4% | 96% | 2.25 × 10⁻⁴ | 9.0 × 10⁻⁶ | 1.00 |
| 0% | 8% | 92% | 2.33 × 10⁻⁴ | 1.9 × 10⁻⁵ | 0.99 |
| 0% | 12% | 88% | 2.43 × 10⁻⁴ | 2.9 × 10⁻⁵ | 0.99 |
| 10% | 4% | 86% | 1.38 × 10⁻⁴ | 5.5 × 10⁻⁶ | 0.56 |
| 10% | 8% | 82% | 1.41 × 10⁻⁴ | 1.1 × 10⁻⁵ | 0.55 |
| 10% | 12% | 78% | 1.44 × 10⁻⁴ | 1.7 × 10⁻⁵ | 0.54 |
| 20% | 4% | 76% | 0.99 × 10⁻⁴ | 4.0 × 10⁻⁶ | 0.37 |
| 20% | 8% | 72% | 1.01 × 10⁻⁴ | 8.1 × 10⁻⁶ | 0.36 |
| 20% | 12% | 68% | 1.03 × 10⁻⁴ | 1.2 × 10⁻⁵ | 0.34 |
| 30% | 4% | 66% | 0.78 × 10⁻⁴ | 3.1 × 10⁻⁶ | 0.26 |
| 30% | 8% | 62% | 0.79 × 10⁻⁴ | 6.3 × 10⁻⁶ | 0.25 |
| 30% | 12% | 58% | 0.80 × 10⁻⁴ | 9.6 × 10⁻⁶ | 0.24 |

The carbon axis brackets the S target, with 10% short of it and 30% beyond.
The ice levels bracket the ice-column target at 0–20% carbon, which includes the
predicted optimum near 16%. At 30% carbon the 12% level falls just below it,
since that level bounds the carbon axis from above rather than being a
candidate.

### Stage 2 — 18–20-µm discriminator

Carbon fixed at the predicted 20%, inclination 60°, mass matched as above.

| Axis | Levels |
|---|---|
| Silicate material | Draine, pyroxene (Mg50), olivine (Mg50) |
| Density profile exponent | −1.5, −1.75 |
| Ice mass fraction | 8, 12% |

**12 cases, of which 2 are shared with Stage 1** (Draine, p = −1.5, 20% carbon,
8% and 12% ice, 60°): **10 new models**. The carbon fraction is fixed from the
opacity screen, not from Stage 1 results, so the whole catalogue is declared
before any model runs.

### Declared comparisons

Every comparison is fixed in advance and reported regardless of fit.

| Comparison | Holds fixed | Count |
|---|---|---:|
| Carbon: 10, 20, 30% against 0% | ice, inclination | 18 |
| Ice: 8, 12% against 4% | carbon, inclination | 16 |
| Material: pyroxene, olivine against Draine | carbon, ice, profile | 8 |
| Profile: −1.75 against −1.5 | material, carbon, ice | 6 |

**34 models and 48 declared comparisons in total.**

## Fixed inputs

Everything outside the axes above is carried over unchanged from the previous
campaign:

- Separate, non-porous DHS grain species with vmax = 0.1; amin 0.03 µm,
  amax 0.4 µm, q = 2.75, 50 size bins for silicate, carbon and ice alike. Using
  the silicate size distribution for carbon is an assumption of this design.
- Self-supplied H₂O 30 K ice at density 0.94 g cm⁻³, table unchanged, with the
  previous campaign's accepted native treatment outside its 0.25–19.94-µm range.
- Cavity half-opening 17.5°; star 4000 K, 2.5 R☉, 0.5 M☉; accretion rate
  7.3011756855 × 10⁻⁸ M☉ yr⁻¹; the disk and its 70/30 silicate/carbon Mie
  species unchanged.
- 512,000 photon packets for temperature, SED and image; seed 43001; a fresh
  temperature solution for every model; density grid 100 × 70 with 20 inner
  radial cells.
- 6001 × 6001 pixels over 6000 AU; image Method 2; signed total-I photometry in
  a source-centred 1″ aperture; model and comparison distances 140 and 147 pc.
- 98 image wavelengths integrated into the same 19 bands: nine near-infrared
  continuum, four ice and six MIRI, plus three unweighted checks.
- `aperture_v2` quality policy: boundary diagnostics recorded; finite positive
  signed image and aperture flux and aperture support mandatory.

**Implementation prerequisite:** the envelope must accept a third separate grain
species. The previous campaign's template carried two (silicate and ice), so the
builder and its template checks need extending before Stage 1 can be prepared.

## How results will be judged

Continuum-normalised feature depths come first and raw fluxes second, because
the previous campaign showed a raw flux can match while the feature is wrong.

- **Silicate 9.7 µm:** core over the power-law continuum between the 8.4- and
  11.0-µm bands.
- **Ice:** the 2.95- and 3.11-µm bands over the continuum between 2.55 and
  3.91 µm, and their raw residuals.
- **18–20 µm:** raw residuals of the 18.45- and 19.8-µm bands. No continuum band
  lies beyond 20.3 µm in the observation contract, so this band pair cannot be
  continuum-normalised.
- **Near-infrared continuum:** RMS of the nine band residuals.

The existing 45/30/25 weighted screening score and its nine weighting and
calibration scenarios are reported as descriptive diagnostics only. They do not
select which models are discussed.

### Success, declared in advance

**One model** meeting all four at once:

| Diagnostic | Tolerance |
|---|---|
| 9.7-µm core/continuum against observed 0.463 | within 0.10 dex |
| 18.45- and 19.8-µm bands | each within 0.15 dex |
| 2.95- and 3.11-µm bands | each within 0.15 dex |
| Near-infrared continuum | RMS below 0.20 dex |

### What each outcome would mean

| Outcome | Reading |
|---|---|
| Carbon brings 9.7 µm onto the data and one silicate material holds 18–20 µm | The spectrum constrains envelope dust composition: silicate, carbon and ice mass fractions. This answers the collaborators' question directly. |
| Carbon brings 9.7 µm onto the data but 18–20 µm overshoots for every material | Composition is not sufficient. The problem lies in the envelope's temperature structure, and the density profile becomes the primary axis. |
| Carbon brings 9.7 µm onto the data but the near-infrared colour breaks | The carbon grain size or optical constants need revisiting, not the fraction. |
| The 9.7-µm feature does not respond as Stage 0 predicts | Opacity ratios do not carry over into full radiative transfer here; envelope structure remains the lever. |

A null or partial result is reported as such. No abundance is fitted to force
agreement.

## Cost

| | Previous campaign | This design |
|---|---:|---:|
| Models | 60 | 34 |
| Images (98 per model) | 5,880 | 3,332 |
| Fresh temperature solves | 60 | 34 |
| Photon packets | 512,000 | 512,000 |

About 57% of the previous campaign's nominal workload, plus five Stage 0
dust-property calculations of a few minutes each. Resources as before: 64 CPUs
and 160 GB per task, at most 16 concurrent tasks, MCFOST memory setting 112 GB,
a 24-hour task limit and a 4-hour subprocess timeout. These are limits, not
runtime predictions; three-species DHS grains may run more slowly.

## Input tables

| File | Source | Density (g cm⁻³) | SHA-256 |
|---|---|---:|---|
| `Draine_Si_sUV.dat` | Draine astronomical silicate, MCFOST distribution | 3.5 | `3a6af751…1dea80` |
| `ac_opct.dat` | amorphous carbon, Rouleau & Martin (1991) | 1.8 | `ef8f0984…a44de2` |
| `H2O_30K_Leiden_mcfost.dat` | self-supplied H₂O ice, 30 K | 0.94 | `29b572f8…8fb00c` |
| `Pyroxene_Mg05Fe05SiO3_Dorschner1995_mcfost.dat` | pyroxene Mg50, Dorschner et al. (1995) | 3.20 | `05c1ebf1…8b44a2` |
| `Olivine_MgFeSiO4_Dorschner1995_mcfost.dat` | olivine Mg50, Dorschner et al. (1995) | 3.71 | `41e34acb…5f3965` |

No new optical constants are required. The carbon table is the one the disk
already uses; the other four are the previous campaign's frozen inputs. Full
hashes are recorded by the builder at preparation.

## Limits

- The depth targets rest on optical depth scaling linearly with column. That
  held across the four sampled columns, but at the previous campaign's fixed
  profile and composition. Carbon changes the heating and will shift the
  temperature structure; fresh temperatures account for that, which is why the
  targets are tested rather than assumed.
- The composition prediction and the 18–20-µm estimate are compact-Mie
  calculations. Stage 0 replaces the first with DHS values; only the radiative
  transfer settles the second.
- The matched masses are analytic. They do not certify identical sightline
  optical depths, scattered-light paths or aperture-integrated extinction on the
  discrete MCFOST grid.
- Carbon shares the silicate size distribution by assumption. A different
  carbon size distribution changes its near-infrared colour.
- Inclination is not constrained by this spectrum and is sampled only at 60°
  and 70°.
- Out of scope: the 3.36-µm red wing of the ice band, which none of these
  changes is expected to fix; an ice-abundance posterior; and any claim beyond
  this single source.
- One image seed per model supplies no numerical covariance. The success
  tolerances are descriptive, not statistical acceptance regions.

## Reproducing the numbers in this document

The compact-Mie opacities, the S ratios and the indicative matched masses come
from `scripts/silicate_opacity.py`, which self-tests against the Rayleigh and
geometric-optics limits on every run. The feature-depth targets come from the
continuum-normalised shapes in
`runs/silicate_structure_production_v1_512k/results/matched_feature_shapes.csv`.
