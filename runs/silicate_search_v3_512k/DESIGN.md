# Silicate search v3: envelope density profile with carbon

Design prepared 2026-10-04, following [silicate search v2](SILICATE_SEARCH_V2.md)
(`runs/silicate_search_v2_512k`, 34/34 models, 48/48 comparisons) and its review,
[validation/silicate_search_v2_512k/FINDINGS.md](../validation/silicate_search_v2_512k/FINDINGS.md).
This is a design only. No run directory has been prepared and no jobs have been
submitted.

The question is narrow: **can a shallower envelope density profile, with the
carbon fraction raised to keep the 9.7-µm feature, bring the 12–20-µm bands
onto the data while the near-infrared continuum and the 3-µm ice band stay
there?** v2 produced models that fit everything from 1 to 11 µm and fail only
at 12–20 µm. This run tests the one standard lever not yet tried against that
failure. Its outcome decides whether the silicate search continues.

## What v2 established

**Carbon worked as predicted.** v2 predicted that about 16% carbon by mass
would bring the 9.7-µm feature onto the data. The normalised 9.7-µm residual
crosses zero at about 15% (Draine, 8% ice, 60°).

**For the first time, one model fits the near-infrared continuum, the ice band
and the 9.7-µm feature together.** Twelve of the 34 models meet the 9.7-µm
shape and near-infrared criteria together, and two of them (models 8 and 29)
also meet the raw ice criterion. None of the 206 models in the three earlier
campaigns met even the first two together: every earlier model that
reproduced the 9.7-µm shape had a near-infrared RMS of at least 0.71 dex.

**What still fails: 12–20 µm is too bright in all 30 carbon models**, by +0.17
to +0.38 dex in the 18.45- and 19.8-µm bands. Along the matched carbon ladder
(Draine, 8% ice, p = −1.5, 60°; log₁₀ model/observed):

| Carbon | 9.7 shape | 8.4 µm | 11.0 µm | 12.7 µm | 13.8 µm* | 18.45 µm | 19.8 µm | 27.5 µm* |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0% | −0.41 | −0.47 | −0.48 | −0.12 | +0.10 | +0.09 | +0.11 | +0.11 |
| 10% | −0.09 | −0.02 | +0.09 | +0.20 | +0.28 | +0.27 | +0.25 | +0.11 |
| 20% | +0.09 | +0.17 | +0.35 | +0.32 | +0.34 | +0.32 | +0.29 | +0.08 |
| 30% | +0.20 | +0.28 | +0.50 | +0.39 | +0.37 | +0.34 | +0.30 | +0.06 |

\* Unweighted check bands.

Carbon raises the 9.7-µm shape faster than it raises 18–20 µm (+0.32 against
+0.18 dex from 0 to 10%), but it raises both. At matched 2.2-µm extinction the
mixtures' κ_abs(18)/κ_abs(9.7) stays at 0.373–0.399 from 0% to 30% carbon while
S falls fourfold: carbon lowers absorption at 9.7 and 18 µm per unit
near-infrared extinction by nearly the same factor. Correcting 9.7 µm and
brightening 18–20 µm are one effect, not two. The excess is confined to roughly
12–20 µm; the 27.5-µm check does not brighten with carbon.

**How each available lever moves the two targets** (median change, dex; from
`lever_responses.csv`):

| Lever | Campaign | 9.7 shape | 18.45 µm | 19.8 µm | NIR RMS |
|---|---|---:|---:|---:|---:|
| Draine → pyroxene, 20% carbon | v2 | −0.061 | −0.036 | −0.025 | −0.010 |
| Draine → olivine, 20% carbon | v2 | +0.014 | −0.039 | −0.031 | +0.057 |
| p −1.5 → −1.75, 20% carbon | v2 | +0.035 | +0.048 | +0.038 | +0.059 |
| Inclination 60° → 70°, all carbon levels | v2 | −0.002 | −0.069 | −0.069 | +0.018 |
| Silicate amax 0.4 → 1.0 µm, 0% carbon | v1 | +0.013 | +0.027 | +0.027 | −0.179 |
| **p −1.5 → −1.25**, lowest column (1.0 × 10⁻⁴ M☉), 0% carbon | v1 | −0.071 | **−0.122** | **−0.106** | −0.124 |
| **p −1.5 → −1.25**, matched column (2.25 × 10⁻⁴ M☉), 0% carbon | v1 | −0.251 | **−0.236** | **−0.211** | +0.162 |

The silicate material is a weak lever: raising κ_abs(18)/κ_abs(9.7) by 31–59%
moves 18–20 µm by only −0.02 to −0.04 dex. The steeper profile tested in v2's
Stage 2 and larger grains both move 18–20 µm the wrong way. Inclination lowers
12–20 µm by 0.06–0.08 dex at every carbon level without touching the 9.7-µm
shape; this design fixes it at 70°. **Apart from inclination, a shallower
profile is the only lever that lowers 18–20 µm by roughly as much as it deepens
the 9.7-µm feature, or more**, and it has never been run with carbon.

**The 3-µm band fixes the ice column, not the ice fraction.** On the nine ice
ladders where the crossing is bracketed, the mean raw 2.95/3.11-µm residual
crosses zero at a column-equivalent envelope ice mass of 0.86–1.14 × 10⁻⁵ M☉
(0.99–1.30 × 10⁻⁵ M☉ on four extrapolated ladders). Over the same ladders the
total dust mass changes threefold and the ice fraction at the crossing runs
4.1–12.8%. The ice fraction that fits rises with carbon only because carbon
lowers the total dust mass at matched near-infrared extinction.

## The hypothesis: too much warm dust close to the protostar

A shallower profile moves envelope mass outward, away from the warm inner region
that emits at 12–20 µm. Carbon pushes the other way on both targets, and it
moves the 9.7-µm shape further than it moves 18–20 µm. Together the two levers
can in principle reach both targets, which neither reaches alone.

The physical picture is standard: rotating infall flattens the envelope density
inside the centrifugal radius relative to r⁻¹·⁵. A single global power law is the
crudest version of that redistribution. It is the one the current template
already supports, since v1 and v2 varied the exponent without code changes.

### Match the near-infrared, not the central column

v1 and v2 compared profiles at a matched central radial column. For a shallower
profile that is the wrong constraint. At the matched column, v1's step to
p = −1.25 dimmed the 1.1- and 1.25-µm bands by 0.72 and 0.61 dex: the extra mass
at large radii extincts the scattered near-infrared light. The near-infrared
continuum is what fixes the column, so the column must be re-matched to it.

v1 shows where that column lies. Interpolating its column ladders (Draine, no
carbon, 70°), the near-infrared RMS is lowest at 1.93 × 10⁻⁴ M☉ for p = −1.25,
**0.86 of the analytic match**, against 2.25 × 10⁻⁴ M☉ for p = −1.5 (the edge
of the sampled range). The best achievable RMS rises from 0.106 to 0.167 dex,
because the profile also changes the near-infrared colour. Between the two
optima the 9.7-µm shape deepens by 0.132 dex while 18.45 and 19.8 µm fall by
0.168 and 0.154 dex.

The two halves of the spectrum respond to different columns. On the steeper
step, the only profile step run both with and without carbon, the 20%-carbon
models moved the 1.1/1.25-µm bands by 0.79 of v1's response at the matched
column. Their 9.7-µm and 18–20-µm responses were only 0.19–0.34 of it. The
near-infrared follows the near-infrared extinction column, which composition
matching holds fixed. The MIRI bands follow silicate absorption, which carbon
dilutes. So the column re-match carries over to carbon models, but the MIRI
gain may shrink.

### What the arithmetic suggests

Adding v1's near-infrared-matched step to the measured 20%-carbon, 8%-ice, 70°
models (FINDINGS F6). The strong variant keeps the v1 step as measured. The weak
variant scales its 9.7-µm and 18–20-µm parts by the carbon weakening measured on
the steeper step. Ice is set by column in this design and is not projected.

| Profile response | Silicate | 9.7 shape | 18.45 µm | 19.8 µm | NIR RMS | Met (9.7, 18–20, NIR) |
|---|---|---:|---:|---:|---:|---:|
| Measured, p = −1.5 | Draine | +0.086 | +0.254 | +0.216 | 0.102 | 2 |
| Measured, p = −1.5 | Pyroxene | +0.020 | +0.211 | +0.186 | 0.059 | 2 |
| Strong, p = −1.25 re-matched | Draine | −0.046 | +0.086 | +0.062 | 0.133 | 3 |
| Strong, p = −1.25 re-matched | Pyroxene | −0.112 | +0.044 | +0.032 | 0.129 | 2 |
| Weak, p = −1.25 re-matched | Draine | +0.062 | +0.197 | +0.168 | 0.133 | 2 |
| Weak, p = −1.25 re-matched | Pyroxene | −0.005 | +0.155 | +0.138 | 0.129 | 2 |

With the column re-matched, the projected near-infrared RMS is about 0.13 dex
in every row. In the strong case Draine meets all three projected criteria, and
pyroxene's feature turns too deep, which more carbon corrects. In the weak case
18.45 µm stays 0.005–0.05 dex beyond its tolerance, which is why the design
extends to p = −1.0. No projection is made for p = −1.0, which no campaign has
sampled. These rows motivate the test; they do not predict its outcome.

## Comparison rule

Every model is matched to the reference near-infrared extinction, scaled to the
reference central radial column, and then re-matched to the near-infrared by a
column factor c:

    κ_ext,mix(2.2) = f_sil κ_sil + f_C κ_C + f_ice κ_ice
    M = M_ref × c × κ_ext,ref(2.2) / κ_ext,mix(2.2) × F(p)

The reference is 96% Draine / 0% carbon / 4% ice at M_ref = 2.25 × 10⁻⁴ M☉ and
p = −1.5. F(p) is the analytic radial-column factor of the previous campaigns,
M(p)/M(−1.5) = [J(p+2)/J(p)] / [J(0.5)/J(−1.5)] with J(q) the integral of r^q
from 1 to 3000 AU.

| p | F(p) | Column factor c | Source of c |
|---:|---:|---:|---|
| −1.5 | 1 | 1.00 | reference |
| −1.25 | 3.60037382695 | 0.86 | v1 near-infrared optimum at 70°, rounded |
| −1.0 | 10.07433956241 | 0.83 and 0.66 | 0.86² = 0.74, bracketed by ±0.05 dex |

The p = −1.0 factor extends the measured log-factor linearly in p and is
uncertain, so two levels bracket it. The near-infrared bands change steeply with
column, so the bracket spans a wide range of near-infrared fits.

**The ice fraction is set by the ice column.** The radial ice column
f_ice × M_ref × c × κ_ext,ref/κ_ext,mix is held at a per-silicate target.
Because κ_ext,mix is linear in f_ice, the fraction follows in closed form:

    f_ice = T A / (M_ref c κ_ext,ref − T (κ_ice − κ_sil)),   A = (1 − f_C) κ_sil + f_C κ_C

with all κ at 2.2 µm. Each target T is the median column-equivalent crossing on
that silicate's bracketed v2 ladders:

| Silicate | Target ice column (M☉) | v2 crossings used |
|---|---:|---:|
| Draine | 9.14 × 10⁻⁶ | 7 |
| Pyroxene Mg50 | 1.11 × 10⁻⁵ | 2 |

Pyroxene needed about 20% more ice column than Draine in v2. A single target
would leave its ice residual roughly 0.1 dex from zero for a reason unrelated to
12–20 µm.

Holding the column rather than the fraction is supported by v1. At p = −1.25 and
its near-infrared-matched column, a fixed ice fraction left the mean
2.95/3.11-µm residual 0.117 dex higher than at p = −1.5. Raising the fraction by
1/c offsets 0.105 dex of that, at the median ice-ladder slope in F5. A lower
column factor therefore raises the ice fraction, up to 17.4% (pyroxene,
p = −1.0, c = 0.66, 25% carbon).

All matchings are analytic or interpolated. None certifies identical optical
depths on MCFOST's discrete, cavity-bearing grid.

## Design

### Stage 0: DHS opacities (no radiative transfer)

Recompute the pure-species DHS opacities for the four species used here
(Draine, pyroxene, amorphous carbon and H₂O 30 K) exactly as v2 did: same
tables, grain parameters and 99 wavelengths. The builder derives every matched
mass and ice fraction from these and freezes them before any model is built.

Stage 0 also checks reproducibility. Its κ_ext(2.2), κ_abs(9.7) and κ_abs(18)
must agree with v2's Stage 0 to 1 part in 10⁶. The v2 values are recoverable
exactly from v2's frozen catalogue, as `scripts/silicate_search_v3_plan.py`
does. A disagreement stops preparation until it is explained.

### Grid

| Axis | Levels |
|---|---|
| Profile and column (p, c) | (−1.5, 1.00) anchor; (−1.25, 0.86); (−1.0, 0.83); (−1.0, 0.66) |
| Carbon mass fraction | 15, 20, 25% |
| Silicate | Draine, pyroxene (Mg50) |
| Ice mass fraction | derived: holds the silicate's target ice column |
| Inclination | 70° |

**24 physical models**, plus two seed replicas (Draine and pyroxene, 20% carbon,
p = −1.25, c = 0.86, photon seed 43002): **26 runs.** Planning values, which
Stage 0 supersedes.

Ice mass fraction:

| Silicate | Carbon | (−1.5, 1.00) | (−1.25, 0.86) | (−1.0, 0.83) | (−1.0, 0.66) |
|---|---:|---:|---:|---:|---:|
| Draine | 15% | 7.9% | 9.1% | 9.4% | 11.7% |
| Draine | 20% | 9.2% | 10.6% | 11.0% | 13.6% |
| Draine | 25% | 10.4% | 12.1% | 12.5% | 15.5% |
| Pyroxene | 15% | 8.4% | 9.7% | 10.1% | 12.5% |
| Pyroxene | 20% | 10.0% | 11.6% | 12.0% | 15.0% |
| Pyroxene | 25% | 11.7% | 13.5% | 14.0% | 17.4% |

Total envelope dust mass (10⁻⁴ M☉):

| Silicate | Carbon | (−1.5, 1.00) | (−1.25, 0.86) | (−1.0, 0.83) | (−1.0, 0.66) |
|---|---:|---:|---:|---:|---:|
| Draine | 15% | 1.16 | 3.61 | 9.77 | 7.85 |
| Draine | 20% | 1.00 | 3.11 | 8.41 | 6.76 |
| Draine | 25% | 0.88 | 2.73 | 7.38 | 5.93 |
| Pyroxene | 15% | 1.32 | 4.11 | 11.10 | 8.91 |
| Pyroxene | 20% | 1.11 | 3.44 | 9.30 | 7.46 |
| Pyroxene | 25% | 0.95 | 2.96 | 8.00 | 6.42 |

At c = 1.00, S/S_ref runs from 0.432 to 0.285 (Draine) and 0.457 to 0.285
(pyroxene) across 15–25% carbon. κ_abs(18)/κ_abs(9.7) is 0.38–0.41 for Draine
and 0.51–0.53 for pyroxene across the whole grid.

Why these levels:

- **Profile and column.** −1.25 is v1's shallower level, now with carbon and at
  its near-infrared-matched column. −1.0 is one step further, sampled for the
  first time, for the weak case. −1.75 is excluded: it moved 18–20 µm the wrong
  way in v1 and v2. No configuration keeps the analytic column at a shallower
  profile, since v1 shows its near-infrared fails there.
- **Carbon.** At p = −1.5 the 9.7-µm zero crossing lies near 15%. A shallower
  profile deepens the feature, so the crossing moves to higher carbon, while the
  lower column partly offsets that. 15–25% brackets it. 18–20 µm rises only
  slowly above 20% (+0.02 dex from 20% to 30% in v2), so the upper level costs
  little there.
- **Silicate.** Draine is the reference. Pyroxene gave v2's best near-infrared
  RMS (0.038) and 9.7-µm shape (+0.022), and lowered 18–20 µm by 0.03–0.04 dex.
  Olivine is dropped: its near-infrared RMS was 0.20–0.21 at p = −1.5.
- **Inclination.** 70° lowered 12–20 µm by 0.06–0.08 dex at every carbon level
  with no effect on the 9.7-µm shape, at a cost of +0.018 near-infrared RMS. At
  60° every projected row would sit about 0.07 dex higher at 18–20 µm, beyond
  the tolerance. One inclination keeps the profile test clean; results are
  conditional on it.
- **Ice by column.** Setting the ice from the column v2 measured replaces v2's
  ice ladder, which tripled its model count. It also turns finding F5 into a
  prediction tested at new profiles.
- **Seed replicas.** Success margins here are a few hundredths of a dex, and no
  campaign has measured seed-to-seed scatter. The replicas sit where a model is
  most likely to land on a tolerance boundary.

### Declared comparisons

Every comparison is fixed in advance and reported regardless of fit.

| Comparison | Holds fixed | Count |
|---|---|---:|
| Profile: each shallower configuration against the anchor | silicate, carbon | 18 |
| Column at p = −1.0: c = 0.66 against 0.83 | silicate, carbon | 6 |
| Carbon: 15 and 25% against 20% | silicate, configuration | 16 |
| Silicate: pyroxene against Draine | carbon, configuration | 12 |
| Seed: 43002 against 43001 | everything else | 2 |

**26 runs and 54 declared comparisons in total.**

### Declared checks

- **Ice-column prediction.** F5 and the v1 offset predict raw 2.95- and 3.11-µm
  residuals within 0.15 dex in all 24 models. The spread is reported against
  configuration.
- **Near-infrared re-match.** The column factors predict a near-infrared RMS
  below 0.20 dex at p = −1.25 and in at least one of the two p = −1.0 columns,
  for every composition. A miss is reported per configuration before any
  success is discussed.
- **Bridge to v2.** The Draine, 20%-carbon anchor (9.16% ice) lies in
  composition between v2 models 15 and 17 (8% and 12% ice, same inclination and
  profile). Its 9.7-µm shape, 18.45/19.8-µm residuals and near-infrared RMS
  should fall within the range those two models span, widened by the measured
  seed scatter. A failure means the campaigns are not directly comparable. It is
  reported before any other interpretation.

## How results will be judged

Diagnostics are those of v2: the 9.7-µm core over the 8.4–11.0-µm power-law
continuum; the raw 2.95- and 3.11-µm bands; the raw 18.45- and 19.8-µm bands;
and the RMS of the nine near-infrared continuum residuals. Normalised ice depths
and the 3.36-µm red wing are reported as diagnostics. The 45/30/25 weighted
score and its nine scenarios remain descriptive only.

### Success, declared in advance

**One model** meeting all four at once:

| Diagnostic | Tolerance |
|---|---|
| 9.7-µm core/continuum against the frozen observed value | within 0.10 dex |
| 18.45- and 19.8-µm bands, raw | each within 0.15 dex |
| 2.95- and 3.11-µm bands, **raw flux** | each within 0.15 dex |
| Near-infrared continuum | RMS below 0.20 dex |

This settles the ambiguity in v2's design: the ice tolerance applies to raw
flux. A success lying within the measured seed scatter of any tolerance is
reported as marginal.

### What each outcome would mean

The 18–20-µm profile response is the median, over the declared profile
comparisons for one configuration, of the mean change in the 18.45- and
19.8-µm residuals.

| Outcome | Reading |
|---|---|
| One model meets all four | The first tested configuration that fits the near-infrared continuum, ice band, 9.7-µm feature and 18–20 µm together with standard MCFOST and the self-supplied H₂O constants. Envelope structure is part of the answer. This is one configuration, not a unique solution or a posterior. |
| 18–20 µm falls by ≥ 0.10 dex in a shallower configuration, but the near-infrared or ice fails wherever the other criteria pass | The excess tracks the radial distribution of warm dust, but a global power law cannot redistribute it without disturbing the outer envelope. Next is the inner-flattened envelope below, not more global exponents. |
| 18–20 µm falls by 0.05–0.10 dex | A partial response, reported as such. Whether the contingency is worth running is decided with the collaborators. |
| 18–20 µm falls by < 0.05 dex in every shallower configuration | Radial redistribution within this template does not reach the excess. The silicate search stops. 12–20 µm is reported as an open tension beyond this template. |
| Raw ice leaves the 0.15-dex band at the target column | The ice-column result does not extend to these profiles. The ice statement for the collaborators is narrowed to p = −1.5 and −1.75. |
| The 9.7-µm shape stays too deep at p = −1.0 even with 25% carbon | The carbon needed exceeds the sampled range at that profile. The limit is reported. |
| Seed scatter exceeds 0.03 dex in any criterion band | Responses of a few hundredths of a dex, including v2's silicate-material effect, are not resolved. All tolerances are read with that margin. |

A null or partial result is reported as such. No abundance is fitted to force
agreement.

## Contingency: inner-flattened envelope (not part of this run)

Only if the second outcome occurs. The envelope would get two zones sharing the
same dust: an inner zone from 1 AU to a break radius with a shallow exponent,
and an outer zone at p = −1.5 out to 3000 AU. That removes warm inner dust
without adding the outer mass that extincts the scattered near-infrared light.
It needs a builder extension (a second envelope zone and its normalisation), a
choice of break radius, and its own declared design. None of that is specified
here.

## Before launch

The collaborators have said the standard pipeline left artifacts in the MIRI
data. The v2 excess is broad (+0.2–0.3 dex across 12.7–19.8 µm), not fine
structure, and it far exceeds the 1.74% sub-band offsets measured earlier. One
question should still go to them first: does their artifact-corrected reduction
change the 12–20-µm level? If it moves by more than about 0.1 dex, the v2
residuals and the projection above must be recomputed before this grid is
frozen.

## Fixed inputs

Everything outside the axes above is carried over unchanged from v2:

- Separate, non-porous DHS species with vmax = 0.1; amin 0.03 µm, amax 0.4 µm,
  q = 2.75, 50 size bins for silicate, carbon and ice alike. Carbon shares the
  silicate size distribution by assumption.
- Self-supplied H₂O 30 K ice at density 0.94 g cm⁻³, table unchanged, with the
  accepted native treatment outside its 0.25–19.94-µm range.
- Cavity half-opening 17.5°; star 4000 K, 2.5 R☉, 0.5 M☉; accretion rate
  7.3011756855 × 10⁻⁸ M☉ yr⁻¹; the disk and its 70/30 silicate/carbon Mie
  species unchanged; envelope radii 1–3000 AU.
- 512,000 photon packets for temperature, SED and image; seed 43001, except the
  two replicas at 43002; a fresh temperature solution for every run; density
  grid 100 × 70 with 20 inner radial cells.
- 6001 × 6001 pixels over 6000 AU; image Method 2; signed total-I photometry in
  a source-centred 1″ aperture after the existing Gaussian PSF convolution;
  model and comparison distances 140 and 147 pc.
- 98 image wavelengths integrated into the same 19 bands, plus the three
  unweighted checks.
- `aperture_v2` quality policy: boundary diagnostics recorded; finite positive
  signed image and aperture flux and aperture support mandatory.

**Implementation prerequisites:**

1. The catalogue key must include the photon seed. v2's builder merges cases
   with identical physical parameters, which would drop the replicas.
2. The column factor c multiplies the reference mass before composition and
   profile matching. Its values, their derivation and the hashes of the v1
   results they come from are recorded in the experiment file.
3. The ice fraction becomes a derived, frozen quantity, computed from Stage 0
   like the masses. The closed form, both targets and the hash of
   `validation/silicate_search_v2_512k/ice_column_crossings.csv` are recorded.
4. p = −1.0 needs no new numerical setting. Its dust masses are 6.6–8.4 times
   those of the matching p = −1.5 anchors, and their outer-envelope temperatures
   rest on the same 512,000 packets.

## Cost

| | v2 | This design |
|---|---:|---:|
| Runs | 34 | 26 (24 physical + 2 replicas) |
| Images (98 per run) | 3,332 | 2,548 |
| Fresh temperature solves | 34 | 26 |
| Photon packets | 512,000 | 512,000 |

About 76% of v2's nominal workload, plus four Stage 0 dust-property
calculations. Resources as before: 64 CPUs and 160 GB per task, at most 16
concurrent tasks, MCFOST memory setting 112 GB, a 24-hour task limit and a
4-hour subprocess timeout. These are limits, not runtime predictions; the
p = −1.0 models may run longer.

## Input tables

| File | Source | Density (g cm⁻³) | SHA-256 |
|---|---|---:|---|
| `Draine_Si_sUV.dat` | Draine astronomical silicate, MCFOST distribution | 3.5 | `3a6af751…1dea80` |
| `ac_opct.dat` | amorphous carbon, Rouleau & Martin (1991) | 1.8 | `ef8f0984…a44de2` |
| `H2O_30K_Leiden_mcfost.dat` | self-supplied H₂O ice, 30 K | 0.94 | `29b572f8…8fb00c` |
| `Pyroxene_Mg05Fe05SiO3_Dorschner1995_mcfost.dat` | pyroxene Mg50, Dorschner et al. (1995) | 3.20 | `05c1ebf1…8b44a2` |

All four were rehashed against the copies in `runs/silicate_search_v2_512k/inputs`
on 2026-10-04 and match v2's frozen values. No new optical constants are
required.

## Limits

- A single global power law moves mass outward at all radii. It is the crudest
  radial redistribution. A rotating-infall structure would redistribute mass
  differently.
- The column factors come from v1 models without carbon, interpolated between
  sampled columns; v1's p = −1.5 optimum lies at the edge of its sampled range.
  The p = −1.0 factors are an extrapolation. The near-infrared re-match check
  tests them.
- The composition and column matchings do not certify identical sightline
  optical depths, scattered-light paths or aperture-integrated extinction on the
  discrete grid.
- At p = −1.0 the matched masses are 0.59–1.11 × 10⁻³ M☉ of envelope dust.
  Before a success at p = −1.0 is interpreted, that mass should be compared with
  independent envelope-mass estimates for TMC1A.
- The ice targets come from log-linear interpolation of v2's ladders.
  Pyroxene's target rests on two crossings, both at 60°. Holding the ice column
  at lower column factors is supported at p = −1.25 only.
- Only 70° is sampled. The spectrum does not constrain inclination.
- The projection assumes additive, non-interacting responses measured in other
  models. It motivates the test only.
- Two seed replicas estimate scatter at one configuration, not a covariance.
  The tolerances remain descriptive, not statistical acceptance regions.
- Out of scope: the carbon size distribution or other carbon tables; porosity;
  disk changes; the 3.36-µm red wing; an abundance posterior; any claim beyond
  this source.

## Reproducing the numbers in this document

- v2 evidence (ladder, opacity ratios, lever table, ice crossings, the
  near-infrared-matched step and projection):
  `python -B scripts/review_silicate_search_v2.py`, which writes
  `validation/silicate_search_v2_512k/FINDINGS.md` and its CSVs.
- Planning values (recovered DHS species opacities, column factors, ice targets
  and fractions, profile factors, matched masses, run count):
  `python -B scripts/silicate_search_v3_plan.py`. It recovers the v2 species
  opacities from the frozen catalogue to a maximum relative residual of 9 × 10⁻¹⁶.
