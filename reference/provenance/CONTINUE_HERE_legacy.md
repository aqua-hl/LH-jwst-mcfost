# TMC1A MCFOST — continue here

Current navigation (18 September 2026): see [START_HERE.md](START_HERE.md) and the
[dated retrospective](reports/restart_review_2026-09-18/REPORT.md). The old cluster
is unavailable. A clean local design/reference workspace exists at
`../tmc1a_restart`; no new model run was launched. The remaining sections below
are historical execution notes, including superseded current-state claims.

Last updated: **2026-07-27 16:20 CEST**.

## COMPLETE: the best dust scored on the band — a three-way tension

```text
15174381  COMPLETED  band-run authorization
15174383  COMPLETED  85/85 band tasks, all 0:0, zero failures
15174385  COMPLETED  fail-closed completion index, 2m08s
completion  57d611f38161f284e8915a0c3c235f5e60cfef1266450e45222a1cf4a4d1260d
tau scores  output/h2o_dust_axis_band_v1_tau.json
1-arcsec    output/h2o_band_1arcsec_dust_axis_v02_vmax0p1_amax0p4_v1.ecsv
figure      output/h2o_best_model_vs_jwst_full_spectrum.png
```

The near-IR winner `v02` (vmax 0.1, amax 0.4, i 70, x 0.04) was scored on the
frozen 0.35-arcsec observable.  The field question was settled before launching:
the earlier band smoke overflowed the n2401 field with *grown grains*, and the
sealed near-IR measurements attribute that to `amax` (x11.33 at 2.5464 microns)
not `vmax` (x1.34).  Projected worst-node edge 2.706e-04; measured 2.104e-04,
i.e. 4.8x headroom.  The projection was accurate and conservative.

```text
                          near-IR RMS (dex)   band RMS (tau)
accepted original dust          0.086            -- (no ice)
previous ice dust vmax 0.8      0.379            0.192
best dust vmax 0.1              0.133            0.287
```

**No variant is good at both.**  And the band degradation is NOT a uniform depth
offset -- an earlier reading of the mean signed residual (-0.2119) said it was,
and that was wrong.  Decomposed:

```text
            blue side (<3.08 um)        red side (>3.08 um)
best dust   mean -0.053  RMS 0.095      mean -0.410  RMS 0.417
prev dust   mean -0.066                 mean -0.250
```

The blue side and core are **slightly better** with the best dust.  The entire
RMS increase is the red side.  Consequently **re-optimising the H2O abundance
will not help**: tau scales the whole profile, the red side needs about 3.9x to
reach JWST, and applying that would grossly over-deepen a blue side that already
fits.  Do not run an abundance grid for this dust.

The measured three-way tension, all of it now RT-verified:

* the near-infrared continuum wants low `vmax` -- `v02` delivers it;
* the band blue side and core are fine either way, marginally better at low `vmax`;
* the red wing is the whole remaining problem and low `vmax` makes it *worse*,
  which the local opacity screen predicted (red/core 0.0845 at vmax 0.8 ->
  0.0778 at vmax 0.2, against an observed 0.2531).

### Next production run: ice optical constants, and only that

The one lever that moves the red wing without paying in the continuum is the ice
itself.  At unchanged grain size the local screen gives:

```text
                     red/core   nearIR penalty
H2O_30K (current)     0.0845        0.0975
ice_opct              0.0789        0.0916
waterice210K (warm)   0.1790        0.0961      <- doubles the red wing
JWST observed         0.2536
```

Hold the `v02` dust, vary the H2O optical constants over `waterice210K.dat` and
`ice_opct.dat`, score on the band.  Both files exist in `MCFOST_UTILS` on the
cluster; `porousice_james.dat` is there too but it is an 85-byte stub and
misbehaved locally, so exclude it.  Roughly `2 x (85 band + 5 near-IR checks) =
180 tasks`.

One code change is required first: swapping the optical-constants *filename* is
not a width-preserving numeric edit, which is what the renderer and its
line-by-line diff check are built on.  It needs a contained extension to admit
that one line, validated by filename.

Capacity note: the 384-core reservation is effectively full (we hold 192, four
other users hold 176, and the 16 idle cores are split 8+8 so no whole 16-core
slot fits).  The same `Cascade-GPU` partition has ~1036 idle cores outside the
reservation, including three fully idle nodes; a `--test-only` probe places a
non-reservation job on `s92node202`.  Its start-time estimates are not
trustworthy here, so the only real test is to submit and watch.

## COMPLETE: the near-IR slope is fixed by DHS vmax, not by grain growth

The chain is terminal and green.  The queue is empty.

```text
15172895  COMPLETED  dust-axis grid authorization
15172897  COMPLETED  104/104 near-IR tasks, all 0:0, zero failures
15172899  COMPLETED  fail-closed completion index, 8m25s
```

```text
manifest        6ddb1c8aefb2a9c20112c912f40f003c58a2b53898949e470d1d99861024925c
task table      00d57a97fca368bf5d265cd32114f61fa1ec2aaa2bfdc55606b04799fc43c3ee
authorization   3bd5021b13d1c755c4a3781252a03bb8dada61fe2ef27e0de16f637f91840a72
completion      55cae076bcb99aa1fd911d1576e77539022970ba4220044d455323225af38bda
scores          output/h2o_dust_axis_grid_v1_scores.json
```

Census exactly 8 anchors x 13 variants.  Model minus observed, dex, at the five
near-IR anchors; `slopeRMS` is their RMS and `blue-red` is the 1.2028 minus
2.5464 tilt, which is insensitive to overall normalization:

```text
variant                       vmax  amax  incl   1.203    1.499    1.759    2.192    2.546  slopeRMS  blue-red
ACCEPTED original dust           -     -  70.0  -0.085   +0.022   +0.031   -0.034   -0.165    0.0861    +0.080
sealed selected dust x=0.05    0.8   0.4  70.0  +0.522   +0.545   +0.415   +0.116   -0.143    0.3937    +0.665

v02_vmax0p1_amax0p4            0.1   0.4  70.0  +0.036   +0.154   +0.153   +0.136   -0.149    0.1334    +0.184
v01_vmax0p4_amax0p4            0.4   0.4  70.0  +0.180   +0.294   +0.263   +0.055   -0.151    0.2070    +0.331
v00_base                       0.8   0.4  70.0  +0.497   +0.524   +0.402   +0.106   -0.150    0.3786    +0.646
v09_base_incl60                0.8   0.4  60.0  +0.584   +0.593   +0.461   +0.147   -0.120    0.4341    +0.705
v10_base_incl50                0.8   0.4  50.0  +0.673   +0.656   +0.508   +0.172   -0.110    0.4866    +0.784
v03_vmax0p8_amax1p0            0.8   1.0  70.0  +0.167   -0.176   -0.466   -0.812   -1.055    0.6401    +1.222
v06_vmax0p8_amax2p0            0.8   2.0  70.0  +0.900   +0.297   -0.183   -0.807   -1.134    0.7575    +2.034
v07_vmax0p4_amax2p0            0.4   2.0  70.0  +0.618   +0.014   -0.476   -1.049   -1.371    0.8472    +1.989
v04_vmax0p4_amax1p0            0.4   1.0  70.0  -0.192   -0.573   -0.848   -1.138   -1.292    0.9000    +1.101
v12_vmax0p1_amax1p0_incl50     0.1   1.0  50.0  -0.154   -0.570   -0.872   -1.152   -1.298    0.9077    +1.144
v08_vmax0p1_amax2p0            0.1   2.0  70.0  +0.426   -0.160   -0.631   -1.192   -1.510    0.9281    +1.937
v11_vmax0p1_amax1p0_incl60     0.1   1.0  60.0  -0.253   -0.654   -0.938   -1.207   -1.347    0.9638    +1.094
v05_vmax0p1_amax1p0            0.1   1.0  70.0  -0.332   -0.721   -0.995   -1.259   -1.397    1.0156    +1.065
```

**`vmax` is the fix.**  Dropping the DHS vacuum fraction 0.8 -> 0.4 -> 0.1 at
unchanged grain size takes the five-anchor RMS from `0.379` to `0.207` to
`0.133` dex and the 1.2028-micron residual from `+0.497` to `+0.036`.  It fixes
the *slope* as well as the level: `blue-red` falls from `+0.646` to `+0.184`
against the accepted dust's `+0.080`.  That closes about two thirds of the gap
between the screening dust (`0.394`) and the accepted original dust (`0.086`),
by the mechanism the opacity ladder predicted — extinction, not albedo.

**Grain growth is actively harmful.**  Every `amax` cell is worse than the
control.  Growth fixes the level at 1.2 microns and destroys the slope
(`v03` runs `-0.18` to `-1.06` dex redward, `blue-red +1.222`).  This is a
measured conflict with the published TMC1A polarization constraint of `>~10`
micron envelope grains, not a suspicion.

**The published inclination does not help.**  70 -> 60 -> 50 degrees brightens
the near infrared monotonically (about `+0.09` dex per 10 degrees at 1.2
microns) and the model was already too bright, so RMS goes `0.379 -> 0.434 ->
0.487`.  It is a level knob, not a slope knob: `blue-red` moves only `+0.646`
to `+0.784` across 20 degrees.  The same shift *helps* the over-extinguished
`amax=1.0` variants (`1.016 -> 0.908`), which is the identical physics with the
sign flipped.


Figure: `output/h2o_dust_axis_grid_v1_scores.png` (script `plot_h2o_dust_axis_grid_v1_scores.py`, provenance sidecar
alongside).  Panel A is the residual spectrum, B the complete 3x3 dust
plane, C the inclination response.

### What to do next

1. **Pin `vmax`.**  The winner sits on the grid edge: `0.1` is the lowest value
   tested and the trend has not turned over, so the optimum may be at
   `vmax -> 0`, essentially Mie.  The opacity screen says little room remains
   (`kappa_ext/accepted` 0.940 at `vmax=0.05` versus 0.962 for Mie), so one
   8-task run at `vmax=0.05` or a Mie disjoint variant would settle it.
2. **Unblock the 3-micron band** (see below).  `v02` has to be scored on the
   band before it can be adopted; the near-IR is only half the objective.
3. **Note the geometry tension.**  The residual wants *more* line-of-sight
   extinction, i.e. higher inclination, but 70 degrees is already above the
   published 50--65 range.  Do not "fix" this by raising `i` further without a
   reason independent of the fit.
4. **The `-0.15` dex at 2.546 microns is not a dust problem.**  The accepted
   original dust shows the same `-0.165` there.  It is common to both and needs
   a separate explanation.

### What it varies and why

`x = 0.04` (the abundance-grid optimum) everywhere.  Block A is the complete
3 x 3 plane of DHS `vmax` (0.8 / 0.4 / 0.1) x envelope `amax` (0.4 / 1.0 / 2.0
um) at i = 70.  Block B probes the published TMC1A inclination (60 and 50
degrees) at both ends of that plane, so geometry and dust cannot be confounded.
All 104 tasks are 1-arcsec broadband anchors w001--w008 on the
`full_domain_n4801` 6000-au field; `w009` is the original-dust control and is
never re-run with this dust.

### The finding that motivated it, now confirmed in RT

The near-IR excess is an **extinction** effect from DHS `vmax`, not an albedo
effect.  A local `-dust_prop` decomposition (`local_opacity_screen/decompose.py`,
seconds per variant on the MacBook) attributes the three bundled dust changes at
1.2028 um as `kappa_ext` relative to the accepted dust:

```text
coated -> disjoint            0.987
ice_opct -> H2O_30K           0.962
Mie -> DHS(vmax=0.8)          0.736     <- the whole effect
```

The screening dust has *lower* albedo and *lower* kappa_sca than the accepted
dust yet is 3.3x brighter, because DHS removes a quarter of the near-IR
extinction and scattered starlight leaks out.  **Correction to earlier guidance
in this file: lowering `vmax` does NOT make the near-IR worse.**  Albedo rises
slightly (0.810 -> 0.851) but `kappa_ext` rises 31 percent, and in an optically
thick envelope extinction wins exponentially.  Score `kappa_ext` relative to the
accepted dust, never albedo.

`vmax` moves the near-IR continuum and leaves the 3-um band core alone
(1.977 -> 1.963 across the ladder); `amax` moves the red wing (excess-opacity
red/core 0.084 at 0.4 um, 0.214 at 1.0, 0.583 at 3.0, against an observed
0.254).  MCFOST **refuses DHS with coating**, so the accepted (Mie/coated) and
screening (DHS/disjoint) representations cannot be interpolated.

A single smoke task already confirms the direction in RT: at 1.2028 um, growing
`amax` 0.4 -> 1.0 um at unchanged `vmax` gives `7.1283e-16` W m-2 against the
sealed `1.6123e-15`, i.e. **+0.167 dex above JWST instead of +0.522**, with every
strict-v2 gate green.  One task moved three quarters of the discrepancy.  It is
one point, not a result.

### The 3-micron band is blocked, deliberately

The grid originally carried 255 band tasks at 0.35 arcsec.  Its smoke ray-traced
cleanly and then failed on finite-field truncation (plane-0 edge positive
fraction `1.6439e-03` against the `1e-3` limit, plus red component
convolution-loss gates) — the grown grains overflow the 3003.75-au field.  The
fix is a bigger field, and **that is not available**: the frozen v3 H2O kernel
admits only `legacy_n1201` / `expanded_n2401` and *requires* `expanded_n2401`
for phase `full_h2o`.

So scoring the band for grown grains needs that sealed kernel extended to a
larger field, with its own field-pair validation exactly as the guard-field
recovery got.  That was not done silently, and the 255 tasks were not spent
failing closed.  This is the top follow-up.

### Two defects the smokes caught, both fixed pre-launch

1. The band measurement call was transcribed without the `spatial_config`
   keyword the frozen v3 kernel requires.  A test now compares the runner's
   literal keyword arguments against both kernel signatures by introspection.
2. Workflow-specific phase names (`nearir` / `band`) were rejected by the frozen
   band kernel, which reads `task_identity["phase"]` and admits only
   `field_pilot` or `full_h2o`.  Phases are now the v3 execution phases,
   preserved from the source row rather than invented.

Both superseded input freezes are preserved byte-identically under
`runs/h2o_3um_dust_axis_grid_v1/_superseded_prelaunch_*/` because `freeze_bytes`
is write-once, and the smoke products under `_smoke_15172826_15172828/`.

## COMPLETE: the abundance optimum is drawn against JWST

`output/h2o_grid_x0p04_vs_jwst_full_spectrum.png` (script
`plot_h2o_grid_x0p04_vs_jwst_full_spectrum.py`) shows the grid optimum
`x = 0.04` beside the original fixed 5 percent, on the same three-panel
construction as before: A and B are 1-arcsec observables, C is the scored
0.35-arcsec optical depth, never pooled.

To put the band on the continuum axis, the 85 frozen `x = 0.04` products are
re-measured at 1 arcsec by `remeasure_h2o_abundance_grid_1arcsec_v1.py`
(read-only, no ray tracing; 85/85 gates green, max edge fraction `2.016e-04`,
table `output/h2o_band_1arcsec_grid_x0p04_v1.ecsv`, sha
`e0560f19dd555ebf01bca882a2c2a1bec6defe90d15c262c8c07af87f06b0a51`).

What it shows: at 1 arcsec the band bottom rises from `0.41x` the JWST median to
`0.71x`; at 0.35 arcsec the nine-anchor RMS falls from `0.311` to `0.192` and
the profile now tracks JWST from 2.86 to 3.22 um within `0.01--0.26` in tau.
The residual is concentrated at one anchor: JWST `0.59` versus model `0.17` at
3.2962 um.  The near-IR continuum is untouched by `x` and is the remaining
problem — `+0.52 / +0.54 / +0.41` dex above JWST at 1.20 / 1.50 / 1.76 um where
the accepted original dust sat within `-0.08 / +0.02 / +0.03`.  The nine
broadband anchors were ray-traced at `x = 0.05` only; a local `-dust_prop` check
gives the same near-IR albedo to three decimals and `kappa_sca` within 0.9
percent at `x = 0.04`, so they stand for both.

## COMPLETE: H2O abundance screening grid — the abundance axis is exhausted

The queue is empty.  The 255-task grid ran green end to end and is sealed.

```text
15172412  COMPLETED  grid authorization
15172414  COMPLETED  255/255 abundance-grid tasks, all 0:0, zero failures
15172416  COMPLETED  fail-closed completion index
```

```text
grid manifest        695d8214cbc87a4979971351ae78f931d11eab6ecc26796e5c4aac27e3c1c993
grid task table      0b00d97cc1b5e3e8e4627d9e83fad026d5eb64bfbb0f7b249b4f166a6c03dce4
grid authorization   673f159b1173eebe51b1ca4934d4c103ed983e35b8cbeb81b93632e4d6c162a2
completion index     c8435f0f56a6b14c9507f68d2d92af137cefdb3fd6373b13db5d47f42e2b3049
```

```text
   x     core tau    9-anchor RMS    mean signed
  0.02    1.1557        0.8734         -0.8258
  0.03    1.7144        0.5243         -0.5030
  0.04    2.2721        0.1916         -0.1475
  0.05    2.8296        0.3109         +0.1750   (reused post-pilot point)
```

Optical depth is linear in abundance: `tau = 0.0402 + 55.80 x`, maximum residual
`0.0004`.  Observed core tau is `2.349529`, giving a **screening bracket of
`x ~ 0.041--0.045`** — core-tau match `0.04139`, RMS minimum `0.04236`,
mean-signed-zero `0.04457`.  Those three disagree by 7.7 percent precisely
because the profile is wrong; a correct model would have them coincide.

**The key result is the residual that survives.**  Nine-anchor RMS falls only
from `0.3109` to `0.1916` (computed at `x=0.04`; about `0.179` interpolated).
The depth-scale diagnostic predicted an ideal pure-rescaling floor of `0.1671`,
so the grid reaches essentially that floor and no further.  About 58 percent of
the original RMS survives and all of it is band shape.  The sharpest symptom is
the red anchor `w010` at 3.2962 microns: observed tau about `0.59`, while every
modelled abundance collapses to `0.1--0.2`.

Do not treat `x ~ 0.041--0.045` as an ice abundance.  It is degenerate with
envelope column, the continuum temperature is not self-consistent for the
changed dust, and no likelihood is published.

Details in `H2O_3UM_ABUNDANCE_GRID_V1.md`; curve in
`output/h2o_abundance_grid_v1_curve.json`, figure
`output/h2o_abundance_grid_v1_curve.png`.  The curve builder
`analyze_h2o_3um_abundance_grid_v1.py` reuses the frozen neutral operator and
re-derives the sealed 5-percent point on every run, aborting if it drifts by
more than `1e-9` (observed drift `1.665e-16`).

### Next axis: composition and geometry, not abundance

A local envelope-only opacity screen (`-dust_prop`, seven variants in **four
seconds** on the MacBook; `local_opacity_screen/`) tested the near-IR continuum
levers and overturned two expectations:

```text
baseline envelope    kappa_sca 4879 cm2/g at 1.2 um, albedo 0.810
+20% carbon          albedo 0.649  -> reflectance ~0.65x   (helps, but grey)
amax 0.4 -> 2 um     albedo 0.868  -> reflectance ~1.19x   (WORSE)
vmax 0.8 -> 0.3      albedo 0.839  -> reflectance ~1.09x   (WORSE)
```

**SUPERSEDED 2026-07-27 — the conclusion drawn from those albedos was wrong.**
Growing grains and reducing porosity do raise the near-IR *albedo*, but albedo
is not the observable.  Both also raise `kappa_ext` far more, and in an
optically thick envelope extinction wins exponentially, so both make the near-IR
*better*, not worse.  See the dust-axis section at the top of this file for the
attribution and for the RT smoke that confirms it.  The one line that survives
unchanged is that this points at extinction or column rather than albedo.

Published TMC1A values worth adopting as priors (see the same document):
distance **141.8 pc** (model adopts 147), inclination **50--65 degrees** (model
uses **70**, outside the published range), measured **A(1.644 um) = 3.57--4.23
mag, A_V ~ 17--20** with a negative 1--2 micron slope, disk `a_max ~ 0.12 mm` or
`~4 mm` (degenerate), envelope grains `>~10 microns` from polarization.

Cheapest next steps, in order: check the model's own line-of-sight A_V against
the measured 17--20; move inclination into the published 50--65 degree range;
then a joint opacity screen over (carbon fraction x amax x vmax) scored on *both*
near-IR albedo and the 3-micron band shape — the Stage-1 selection scored band
shape only and never tested the continuum, which is how the near-IR regression
got in.

Two scheduler facts will bite again: `MinJobAge = 300`, so `afterok:<id>` fails
once the dependency job leaves the live controller five minutes after
completing; and `freeze_bytes` is write-once, so re-authorizing must preserve the
previous authority rather than overwrite it.  Authorize and submit in one
uninterrupted cluster-side pass.

A presentation figure covering the published screen is at
`output/h2o_screen_vs_jwst_full_spectrum.png` (script
`plot_h2o_screen_vs_jwst_full_spectrum.py`, provenance sidecar alongside).
Panels A and B are a 1-arcsec observable and panel C is the scored 0.35-arcsec
observable; they are never pooled.

To draw the ice band on the continuum axis, the 85 frozen full_h2o products are
re-measured in the 1-arcsec aperture with the unchanged strict-v2 kernel by
`remeasure_h2o_band_1arcsec_v1.py` (read-only, no ray tracing, all 85 gates
green, max edge fraction 2.0e-4).  Its table is
`output/h2o_band_1arcsec_remeasurement_v1.ecsv`.

Two cautions this figure exposed, both worth remembering:

- `output/continuum_sed_R100.ecsv` contains two **overlapping grating segments
  with duplicated wavelengths** (2.871--3.173 microns appears twice).  Plotting
  it in table order draws the line backwards through the ice band.  Use
  `output/tmc1a_sed_unstitched.ecsv` and plot one polyline per `segment`.
- The accepted original-dust model was only ray-traced at the nine broadband
  anchors, so it has **no prediction inside the ice band** and its curve stops
  there.  Running `parameters/original_dust_guardrail_rt.para` at the same 85
  band wavelengths would fill that in, at a cost of 85 further RT tasks.

At 1 arcsec the modelled band bottoms at `1.0533e-14` W m-2 versus a JWST median
of `2.5569e-14` over 3.00--3.15 microns, i.e. `0.41x` — too deep and too narrow,
the same sign as the 0.35-arcsec optical-depth result, so the two apertures
agree on the direction even though they are never pooled.

## COMPLETED: the H2O 3-micron screen is published

Everything below this section is history.  The production chain is terminal and
the Slurm queue is empty.

```text
15171420  COMPLETED  85/85 H2O-profile Method-2 tasks, all 0:0
15171421  TERMINAL   9 broadband guardrails: 7 x 0:0, rows 0 and 1 red
15171434  FAILED     obsolete all-green finalizer; published nothing (correct)
15171917  FAILED     guard-recovery inventory, first attempt (bug, see below)
15171922  COMPLETED  guard-recovery inventory after two pre-freeze corrections
15171926  COMPLETED  guard-recovery authorization
15171928  COMPLETED  2/2 n4801 guard recoveries, all 0:0
15171931  COMPLETED  guard-recovery finalizer; green seal
15172346  COMPLETED  composite mixed-lineage publication
```

The published science bundle is
`output/h2o_3um_method2_postpilot_recovered_v1/` (ten products plus SHA-256
sidecars), copied back locally and verified.  The guard-recovery bundle is
`output/h2o_3um_method2_postpilot_guard_field_recovery_v1/`.

```text
composite seal   a437cd4a0a3de8c2ea1ee1062827ffb2a9335dd9582ac69e90ed08a21d28fe56
recovery seal    6744c9b58b3bb9a11d02720b90aedc87685aa88186220ae1f1a68525051f1d86
recovery manifest ac544fa8aaf9c0b182c844495d26977cda41a2ef794880cd342dacd38cccdef7
recovery policy  001ef252001ce2bb552816a9ce8d20b06df154789593ba09c521ce5eef0ac6e1
terminal inventory f5448e015628ac3d19568ea3fca11c6cdb23feb98a5f9ce6ba0cadcafd673268
recovery authorization 4a3915e49e5a4f37afab836b0af7b74587a8dba2a659b9fa0799856c7ce5af7b
```

### Scientific result

The fixed 95:5 silicate:H2O screening model reproduces the 3-micron band
*shape* well but is systematically **too deep**.  It is not an abundance fit.

```text
unweighted RMS tau residual, nine anchors   0.3109497762878781
mean signed tau residual                    0.17499898371648642 (model too deep)
peak-normalized nine-anchor shape RMS       0.07534004547633222
model peak tau   / observed peak tau        2.82964617520377 / 2.3495291661829056
both peak at                                3.0764795490679795 um
effective depth scale alpha (diagnostic)    0.861177815684671
RMS after that scaling                      0.16714717638478654
largest residual                            +0.5279984627600212 at w005, 3.0088 um
only negative residual                      -0.3905623687094077 at w010, 3.2962 um
blue side w002--w005 RMS                    0.31255060502054505
red side  w007--w010 RMS                    0.2494714438645019
unscored wing closure w001/w011 RMS         0.0024833154558481663
```

The residual is positive at every scored anchor except the red 3.2962-micron
one.  Peak position is exact and the peak-normalized shape RMS is only 0.0753,
so the band *profile* is close; the amplitude is not.  The diagnostic depth
scale `alpha=0.8612` is the least-squares factor that would bring the model
onto the observation and cuts the RMS from 0.3109 to 0.1671; it is explicitly
**not** an abundance or a fitted parameter.  Read together with the red-side
sign flip, this says the screening dust is somewhat too absorbing *and*
slightly too broad on the red side, so the next step is a real
abundance/temperature grid rather than a rescaling of this run.

All five predeclared observation-only gain scenarios keep the same picture:
the nine-anchor RMS moves only over `0.3056--0.3179` and the conclusion is
unchanged.  No formal likelihood, chi-square, or information criterion is
published, because the supplied flux-only errors omit shared continuum
covariance.

The nine broadband guardrails are report-only and were never pooled with the
H2O objective.  The `w009` 27.513-micron original-dust control reproduces the
accepted continuum prediction to `1.872e-4 dex`, which is the intended
reproducibility check.  The selected H2O dust is brighter than the accepted
original-dust prediction at short wavelengths (`+0.607 dex` at 1.2028,
`+0.523 dex` at 1.4988, `+0.383 dex` at 1.7589) and slightly fainter in the
mid-infrared (`-0.218 dex` at 13.7999).  These are screening-only and use a
reused, non-self-consistent temperature.

### Two pre-freeze corrections were required

The prepared recovery workflow could not run as written.  Both bugs were fixed
before any freeze and before any recovery ray tracing; neither relaxes any
acceptance gate.  Full evidence is in
`H2O_3UM_METHOD2_POSTPILOT_GUARD_FIELD_RECOVERY_V1.md`.

1. The inventory's diagnostic relaxation monkeypatched
   `strict_v1.COEVAL_TOLERANCE`.  Strict-v1 also compares that constant for
   **exact equality** against the closure report's recorded `tolerance`
   provenance field, so relaxing it made a clean `1e-05` closure look forged
   and job `15171917` aborted.  `COEVAL_TOLERANCE` is now excluded from the
   relaxation.
2. The admission signature named only the `scattered_star` component edge
   sub-gate and required all six top-level strict-v2 checks green, so it could
   never admit the rows it existed for.  The component roster was in fact
   correct: `scattered_thermal` exceeds 1e-3 numerically but is **inactive**,
   so its sub-gate is vacuously green.  The real omission is the **top-level**
   plane-0 edge check, which has no active-component qualifier and goes red
   (0.0026655 and 0.0015189) because plane 0 is dominated by the truncated
   scattered-star map.  Admission now names exactly that one top-level check.

The recovery itself confirms the diagnosis: on the 6000-au 4801-square field
the plane-0 edge fractions drop to `6.039e-05` and `2.219e-05`, every strict-v2
gate is green, and the signed 1-arcsec aperture flux changes by only
`0.026067` and `0.044497` percent, inside the preferred 1-percent limit.

### Suggested next step

Define a small H2O abundance (and, if affordable, dust-temperature) grid around
the fixed 5-percent point rather than rescaling this screen.  The amplitude
mismatch is coherent and the shape is already good, so a one-dimensional
abundance scan at the same geometry should bracket it quickly.  Keep the
guardrails report-only and keep the 0.35-arcsec/1.0-arcsec aperture distinction.

## Superseded: H2O Method-2 field-size pilot v3 is GREEN

The replacement H2O 3-micron Method-2 workflow is implemented at
`runs/h2o_3um_method2_screen_v3/`.  The first cluster affinity smoke
(`15171274`) was green, but the first lambda smoke (`15171276`) exposed an
incorrect validator assumption and blocked before any science array.  MCFOST
completed normally: it read the exact request into binary32 for the SED, then
serialized the image header from that value with CFITSIO `%.7G`.  The old
validator incorrectly required SED and image `WAVE` values to be equal.

The complete rejected attempt, old manifest, all 18 old workflow files,
runtime logs, and FITS products are preserved under
`runs/h2o_3um_method2_screen_v3/_failed_attempts/prescience_manifest05ac_lambda_serialization_job15171276/`
and are evidence only.  The corrected contract binds exact request -> SED
binary32 -> seven-significant-digit image serialization; it never widens the
exact-request gate and never uses a rounded value for flux conversion.  All
84 unique wavelengths remain uniquely serialized, while six would have been
falsely rejected by the old direct-image gate.  The corrected freeze and
50-test suite are green locally and on the cluster.  The matching cluster
attempt is archived.  Fresh corrected-manifest smokes are green: affinity job
`15171278`, seal `684fce86ed7d1f420ba5070fb37d7d0bb860ac6e9874745beaf70c6edf85e877`;
lambda-contract job `15171280`, seal
`75845357853a7fb6615ad0f094492a04383b9fab087573e835d864826be7cf1b`.

The 44-task field-size pilot completed as Slurm array `15171282`: all 44 tasks
are `COMPLETED 0:0` with exactly 16 CPUs and immutable completion/measurement
records.  Fail-closed finalizer `15171295` completed `0:0` and published the
automatic green seal
`f031983db2157b5211d568ae1c5b5d060d5747283de77d21f2ee88f47af30628`.
The user queue is empty.

All 22 expanded-field tasks hard-pass.  All 22 field pairs pass; the maximum
1201-to-2401 symmetric aperture difference is 0.043555 percent versus the
exclusive 2-percent limit.  The maximum seed node-flux difference is 0.4103
percent versus 20 percent, the core-continuum difference is 0.09378 percent
versus 5 percent, and the core-tau absolute difference is 0.0009586 versus
0.10.  Legacy records contain 18 hard passes and four exact predeclared
scattered-star map-support diagnostics only.  Adopt n2401 plus 128k for later
H2O screening; never mix or use n1201 as the science geometry.

The fixed five-percent H2O pilot gives core tau 2.8306--2.8316 at 3.07648
microns, while the frozen JWST value is 2.34953.  Its local continuum is close
(about 0.1794 versus 0.1842 Jy) but its core is too deep (about 0.01057 versus
0.01757 Jy).  This is a preliminary direction toward lower effective H2O
abundance/column, not a fitted abundance result.  Full details are in
`H2O_3UM_METHOD2_FIELD_PILOT_V3_RESULT.md`; compact seals, 44 JSON record
pairs, a summary and the inspected diagnostic plot are in
`output/h2o_3um_method2_screen_v3_field_pilot/`.

The green seal cryptographically authorizes the fresh 85-task n2401 H2O
screen, but **do not launch it or the nine guardrails yet**.  The user asked to
report convergence and wait for a new direction at this decision point.

The pilot is two image fields (`legacy_n1201`, `expanded_n2401`) times two
128k-packet seeds (52001 and 52002) times 11 spectral cells.  The cells are
all five R=2400 Gauss--Hermite nodes for core target `h2o_3um_v3_w006`, plus
the central node of each of the six continuum supports.  Both fields preserve
the same 1.251039932-au pixel scale.  All 22 expanded-field tasks must pass
unchanged physical gates.  The legacy field may reproduce only the exact
known scattered-star convolution-support diagnostic; any total-plane failure
or any additional failure is fatal.  Every matched signed 0.35-arcsec
aperture pair must agree to less than two percent.  The expanded-field seed
gates remain 20 percent for target/support node flux, five percent for the
derived core continuum, and 0.10 absolute for core optical depth.

The full-screen estimator converts each exact Method-2 `lambda F_lambda`
node to `F_nu` in Jy before the frozen five-node aggregation and uses the
neutral observation-v3 operator directly.  The nine one-arcsec broadband
guardrails are report-only; eight use the selected H2O dust and one retains
the original dust as a control.  The selected dust overlay is the frozen
disjoint-population DHS (`vmax=0.8`) mixture: 95-percent Draine silicate and
five-percent 30-K H2O optical constants, with `q=2.75` and
`a=0.03--0.40 micron`.  It reuses the accepted continuum temperature and is
therefore an explicitly screening-only, non-self-consistent dust experiment.

Local verification is green: comparison-only freeze, 50 tests, Python
compilation, and all four Slurm shell syntax checks pass.  Current authorities:

```text
input manifest                     6faf36dd78b9cab8776d1634f6e2859782c2e4c2a8ba9b5e78aabfbc19999c79
44-task field pilot table          bcec528a3e252602b76ac3923ed56c9a2ea819718af8976e90e524c535e2317a
85-task full H2O table             fadcef787ccbaf05920ed2a4bc0800e24886fbb966e2343f8d45b6e74582cf98
9-task guardrail table             b9755f852f622536d62803514e3892b92f71a6d71b4c939790fea92e1fab034a
accepted continuum parameter       cee9777 (full SHA is pinned in the manifest)
accepted Temperature               1bd212f2f54dd8c1701a6064e129ba7fdb5721cb3fe296472433bbdd3b79ec1c
custom MCFOST binary                d96bbc7cb71c09fce93edcef5908e6af1eb0e90e40193a3afdd259c2a7786401
```

## Current handoff: corrected stellar H2O observation v3 frozen

The v2 observation contract failed independent preflight and is quarantined
byte-identically.  Do not score a model with any
`output/h2o_3um_observation_v2` table.  The four failures were an observed
`F_nu` versus model `lambda F_lambda` representation mismatch, a non-scale-
invariant legacy Huber fit, an artificial seam veto that moved B90 by 1.9635
R=2400 FWHM, and an undersampled spectral-degradation kernel.  Full evidence
and exact preserved hashes are in
`H2O_3UM_OBSERVATION_V2_FAILED_PREFLIGHT.md`.

The corrected local-only v3 contract is frozen at
`output/h2o_3um_observation_v3/`.  It uses separate native valid runs, cubic
interpolation to a uniform log grid and a 16-times-finer grid, and a full added
Gaussian.  Twelve sub-pixel phases recover R=2389.25--2400.01, inside the
declared one-percent R=2400 screening gate.  The actual-run plus/minus-five-
sigma gate allows the corrected B90 at 3.0087755423 microns to cross the 3.01-
micron stitch seam; every selected row is within 0.2126 FWHM of its ideal
target.

The observation and model now share one executable neutral operator in
`h2o_observation_measurement_kernel_v3.py`: MCFOST Method-2 plane-0
`lambda F_lambda` nodes are converted to `F_nu`/Jy before five-point
Gauss--Hermite aggregation, then the same scale-invariant six-support Huber
continuum and `-ln(F/continuum)` are applied to both sides.  Model ingestion
requires the exact hash-bound requested wavelength, the product FITS `WAVE`,
and the flux.  It validates the known seven-significant-digit FITS
serialization while using the exact request in the unit conversion; wrong,
adjacent, reordered, truncated, or fake-parent products abort.

V3 contains six continuum supports, 11 logical targets, nine scored feature
anchors, 30 support LSF nodes, and 55 target LSF nodes.  Its canonical result
is `N(H2O)=3.94415309e18 cm-2`, peak `tau=2.34952917` at 3.07647955 microns.
The v3 overlap-normalization sensitivity is plus/minus 1.9744 percent; the v2
1.9345-percent value is retained only as a predecessor comparison.  An
optional independent-segment calibration screen uses plus/minus 4.2426
percent relative gain.  All are fixed observation-only refit scenarios under
a zero-sum log-gain gauge; none is fitted or added per anchor as variance.

Validation is green: 17/17 v3 tests, 10/10 hardened continuum-kernel tests,
byte-identical `--check`, all seven SHA-256 sidecars, and visual inspection of
the diagnostic PNG pass.  Primary hashes are:

```text
v3 manifest                       519203acdb1271c02f10c927af7225e404d43856ae4bda1dd4a7dab9fba4a9e7
v3 diagnostic PNG                 c331d83a42db9bb15982586c2c22544095080445da59359ef326116595d668de
v3 dense profile                  abeec0f827f875a38d2511807efe0057529f42851b5e98e41acd4bd805fada24
v3 11 logical targets             87e06fad6ba261bdc56bd2100be5606e5754d4ff472f0aa3c4cf37275a5580e7
v3 six supports                   85ed565e146e08a2e4c101f00506c860c8f23f85e30a67327a33f81560a08b02
measurement kernel                a90f104f252ba97c6a1f0f78d703807c8b3aea3ca3503d12c0146eb37747ea72
scale-invariant continuum kernel  b14f05cfd4b46262ac4fe65f12807b1b14405191c240f81557293921dd30a521
```

The offline Stage-1 opacity rescore is now also frozen and independently
audited at `runs/h2o_3um_opacity_screen_v2/`.  It reuses the six existing
hash-closed `-dust_prop` products and scores only their R=2400 opacity shapes;
it does not fit abundance or amplitude and did not run MCFOST.  The complete
five-candidate order is unchanged in all 15 predeclared continuum/relative-
gain cases.  The selected representation is disjoint refractory and H2O
populations with DHS (`vmax=0.8`):

```text
canonical score                    0.10061788702175113
runner-up (disjoint / Mie)         0.12845956772582914
winner-to-runner-up gap            0.027841680704078012
opacity manifest SHA               00da86fb1fdc4c253d4ea74faa81c47fe0cb867c6c9a1d0d80837b98cdc779ab
selection SHA                      500532e7d927f87c519762e7c3f69ff54a99acd6653c3f322a152994588da73d
handoff SHA                        f517a2e6a7f65a2bd387b4ae97273b1058e10ed728ad2fa933099b2d9f4155da
```

The rescore passes 22/22 tests, byte-identical `--check`, compilation,
sidecar verification, visual inspection, and a separate independent score
reconstruction.  The generic-ice curve remains an unranked control.  Details
are in `H2O_3UM_OPACITY_RESCORE_V2.md`.

The next gate is an isolated v3 Method-2 field-size pilot, not the full H2O
production screen.  Its planned matrix is 44 executions: two image fields
(1201 and 2401 square), two 128k-packet seeds (52001/52002), and 11 spectral
cells (all five GH nodes of core target `h2o_3um_v3_w006` plus the central
node of each of six continuum supports).  The 1201-square products are
diagnostic comparators and may reproduce only the already-isolated
`<1e-4` PSF-convolution-support boundary failure; all 22 2401-square products
must pass every unchanged physical gate, and each matched signed 0.35-arcsec
aperture pair must agree within two percent.  Only a green cryptographic pilot
may authorize a fresh, uniform 85-execution n2401 H2O screen.  No cluster job
is currently running or authorized by the opacity handoff.

## Historical handoff: v1 archived; v2 frozen before failed preflight

There are currently **no jobs in the cluster queue**.  The focused H2O
convergence array `15170919` and its impossible-to-green finalizer `15170932`
were cancelled deliberately after four paired-seed cells exposed a uniform
finite-map-support problem.  This was not a wavelength-file, MCFOST-physics,
or photon-count failure.  Tasks 10, 13, 14, and 16 exceeded only the frozen
scattered-star PSF-convolution loss gate, by less than one percent of that
gate:

```text
hard exclusive limit             1.00000000e-4
task 10                          1.00638819e-4
task 13                          1.00277225e-4
task 14                          1.00892551e-4
task 16                          1.00685027e-4
paired passing seeds             9.82959e-5 -- 9.84759e-5
```

All four failed cells otherwise have exact coeval SED/map closure, zero
negative flux, passing edge-flux fractions, and valid apertures.  Seed-dependent
pass/fail switching at the same wavelength proves that the present 1201-square
map lies too close to the support threshold.  Do not reroll seeds until they
pass and do not relax the threshold.  Slurm recorded 14 completed cells and
four failures; 15 immutable completion seals exist because task 15 sealed just
before the cancellation race.  No convergence or wider-run seal was produced.
The complete attempt is now preserved in the independently copied and
checksum-verified cluster archive
`runs/h2o_3um_method2_screen_v1/_failed_attempts/convergence_array15170919_psf_convolution_support_boundary_sensitivity_v1`.
It contains 568.63 MiB plus checksum metadata; both the retained live tree and
the archive contain all 15 seals.

The planned rank-blind recovery is an expanded-support pilot, but it must be
rebuilt as v2 against the new observation contract rather than launched from
the old v1 wavelength table:

```text
image grid                         2401 x 2401
zoom                               1.9975051980918284
pixel scale                        1.251039932 au/pixel (unchanged)
field of view                      3003.74687672 au (doubled)
photons / seeds                    128k; 52001 and 52002 (unchanged)
matched pilot cells                22 at 1201 square + 22 at 2401 square
```

The logical subset remains all five LSF nodes of one core H2O target for both
seeds, plus the central LSF node of each of the other six convergence supports
for both seeds.  Because the observation, R=2400 LSF, anchor wavelengths and
aperture origin are changing, the formal field-size comparison needs the same
22 new-v2 cells at both 1201 and 2401 square.  Keep every physical gate and
require 1201-versus-2401 signed 0.35-arcsec aperture agreement below two
percent (one percent preferred).  If green, refreeze and rerun **all** H2O
cells uniformly at 2401 square; never mix per-cell map sizes.

A replacement observational package is now present at
`../DATA/reproducible_ice_extraction`.  Its 15 cube hashes and script/config/lab
provenance match, but it cannot be adopted wholesale:

- the useful `stellar_center` H2O extraction is the like-for-like target; its
  reproducible G395H 3.70--3.90-micron centroid is only 0.04660 arcsec from
  the old ice-package position and its G235H/G395H overlap passes at
  2.09-percent scatter;
- its delivered stellar flux was empirically corrected from a 0.35-arcsec
  aperture to 1.20 arcsec, whereas the MCFOST H2O estimator integrates an
  uncorrected 0.35-arcsec aperture;
- reconstructing the operator-matched, gap-safe R=2400 raw 0.35-arcsec profile
  gives `N(H2O)=3.91206e18 cm-2` and peak `tau=2.33903` with the dense
  unweighted reference continuum; applying the package algorithm to raw flux
  with random-only weights and using the newly reselected six-support Huber
  continuum give `3.89331e18` and `3.90519e18 cm-2`, a total scientifically
  usable estimator spread of only 0.48 percent; background subtraction changes
  the raw-column result by only 0.01884 percent;
- all three continuum estimators retain the same Stage-1 opacity ordering,
  headed by `disjoint_populations__dhs_vmax0p8`; the v2 sparse score will use
  the newly reselected six-support Huber operator identically on JWST and
  MCFOST, while the dense unweighted and random-only-weighted continua remain
  correlated-systematic checks.  The exact delivered package product is the
  aperture-corrected/total-error result (`3.70031e18 cm-2`), not the raw
  `3.89331e18` diagnostic;
- the new stellar centroid and old notebook coordinate differ by 0.121823
  arcsec.  At 0.35 arcsec this changes the observed column by 2.58 percent,
  so v2 identifies the new centroid with MCFOST `(0,0)` and keeps a separately
  re-extracted old-coordinate profile as report-only sensitivity.  The older
  1-arcsec continuum result remains valid because the corresponding model
  aperture change is at most 0.035 percent;
- the newly defined off-source `ice_peak`, 1.35363 arcsec from the star, is
  boundary-sensitive, has a 20.7-percent grating-overlap warning, 103/1260
  integration samples below S/N=3, and a nominal tau=5.193 spike measured at
  only S/N=0.16; do not fit it quantitatively;
- the CO2 4.27-micron homogenized products interpolate across a real G395H
  gap (0/120 native ice-peak blue-continuum samples become 116/116 nominally
  valid), so the 4.27 columns and 15.2/4.27 ratios are invalid pending a
  gap-safe redo.

The full send-back audit is in `ICE_ANALYSIS_REANALYSIS_REQUEST.md`.  The new
stellar-centred observation contract is now frozen at
`output/h2o_3um_observation_v2/`.  It is generated from the native central CSVs
using `flux_raw_Jy`, no local-annulus subtraction, no 1.20-arcsec correction,
and `error_spectral_covariance_Jy` as the random term.  The correlated
3-percent NIRSpec calibration floor remains a segment nuisance and is not
averaged down by the spectral kernel.  Its exact cardinalities are six
continuum supports, 11 logical targets (nine scored feature anchors), 55
target and 30 support constant-R=2400 Gauss--Hermite nodes.  All full
plus/minus-5-sigma gap/seam gates pass; 14/14 tests and byte-identical
`--check` pass.  Key hashes are:

```text
observation manifest              45d0397d551b10bc774f499bc5cb565cfaae2bd7f58f36fa5cf24d8b4c347f93
dense R=2400 profile              994fa6a0cfeff818e151998e082e7abac0d4524fc64808b06beaad0b72aac89d
11 logical targets                28cd6182160c1a1cca6eb3e599aabc7e7257b8692f14fd68bbbf2e2f47bb6306
six continuum supports            403dd774774b7ab9ed6c8b9e7a14159ccf1882629a28f17d6352411ff2693b31
```

The old support wavelengths produce a 6.85-percent column shift on the new
spectrum and are forbidden in v2.  The next local gate is a separately sealed
rescore of the six existing Stage-1 opacity curves against this observation;
no MCFOST dust-property calculation or Method-2 ray tracing is being launched
while that rescore and an independent observation-contract audit run.

## Completed handoff: source-128k Method-2 continuum screen

The uniform source-only remeasurement, independent seal, and canonical
ranking are complete.  Independent grid-contract and ranking-result audits
both pass with zero assertion failures.  There are no active jobs in the
user's Slurm queue.  The finished chain is:

```text
inventory                        15168723  COMPLETED/0:0  00:26:17
405-cell v2 remeasurement       15168724  405/405 COMPLETED/0:0
four-process independent seal   15168725  COMPLETED/0:0  01:06:16
four-process canonical analysis 15169142  COMPLETED/0:0  00:01:21
rank-free measurement index SHA 58703693e83d38717b800edb0e0947e399c25356b1952695ae53ddeef8103d8e
```

The canonical screen uses all 405 original seed-41001, 128k-packet,
Method-2 values and no high-N replacement.  It scores exactly `w001..w009`;
`w000` is absent.  The comparison figure contains 312 unscored grey JWST
R~100 bins and nine scored black anchors.

The nominal primary rank 1 is
`m2p25_c17p5_a0p40_q2p75_i0p05`: envelope-mass factor 2.25, cavity
half-opening 17.5 degrees, envelope size exponent q=2.75, amax=0.40 microns,
ice fraction 0.05, and inclination 70 degrees.  Its primary score is
0.1093125353 dex.  This is a screening minimum, not a unique physical fit:
rank 2 (`m2p25_c20p0_a0p40_q2p75_i0p05`) is only 0.0004482660 dex worse,
eight primary models fall inside the operational 0.0108579276-dex resolution
tolerance, and the two boundary-handling sensitivity scores choose rank-2's
20-degree cavity instead.  The defensible continuum region is therefore
mass factor 2.25, q=2.75, and cavity approximately 17.5--20 degrees, with a
secondary mass-2.50/q=3.25 family still degenerate.  Because mass and q sit
on lower grid boundaries, this screen does not bracket their continuum
minimum.

All 64 post-canonical Cartesian substitutions of the three available high-N
diagnostics preserve the winner within each scoring scenario.  Their largest
model-score change is only 0.00353 dex.  Task 371's red high-N pair remains
diagnostic and was never promoted into the canonical partition.  This removes
the three convergence diagnostics as the cause of the rank ordering, but it
does not turn the tiny rank-1/rank-2 gap into a resolved physical distinction.

Canonical outputs:

```text
output/strict_continuum_cavity_mass_q_source128k_uniform_v2_v1.png
output/strict_continuum_cavity_mass_q_source128k_uniform_v2_v1_predictions.ecsv
output/strict_continuum_cavity_mass_q_source128k_uniform_v2_v1_scores.ecsv
output/strict_continuum_cavity_mass_q_source128k_uniform_v2_v1_regions.ecsv
output/strict_continuum_cavity_mass_q_source128k_uniform_v2_v1_summary.json
output/strict_continuum_cavity_mass_q_source128k_uniform_v2_v1_validation.json
output/strict_continuum_cavity_mass_q_source128k_uniform_v2_v1_highn_sensitivity.json
```

The next scientific handoff should not call the nominal row a unique winner.
For the ice pilot, the top two rows are particularly convenient: they share
the same mass factor, q, amax, ice fraction, and inclination and differ only
in cavity angle.  Use the 17.5-degree row as an explicitly labelled nominal
screen baseline and propagate the 20-degree row as a geometry sensitivity
unless a lower-mass/lower-q continuum refinement is requested first.  The
ice pilot remains unlaunched pending the Stage-1 opacity choice and its 1x16
cluster affinity smoke test.

## Completed execution record: canonical source-128k Method-2 remeasurement

The original 45-model cavity--mass--q grid has all 405 raw Method-2
SED/image/closure triples, but its old v1 compact measurements were not a
uniform basis for ranking.  The source-only workflow therefore remeasured all
45 models x `w001..w009` from the original seed-41001/128k products with the
same v2 signed-aperture estimator.  No high-photon value is canonical; task
371 remains red and high-N products are optional post-ranking diagnostics
only.

The complete 19-file workflow passed an independent audit and 66/66 tests
locally and on the cluster.  Its secure v3 input adapter now performs two
complete reads through the same open descriptor and requires byte-for-byte
agreement before hashing/parsing, which closes the NFS timestamp-cache issue
found by the first cluster test.  A second preflight caught an over-literal
task-46 failure-log check: the original v1 execution stopped first at
`Intrinsic separated-component gate failed: scattered_star`, so its stderr
cannot also contain the later aperture diagnostic.  The corrected contract
separates the exact terminal exception from the later rank-blind diagnostic
labels and was independently re-audited.  The first frozen manifests were
preserved, not deleted, under:

```text
runs/_failed_preflight_contracts/source128k_uniform_v2_v1_task46_log_semantics_policy98de8631
```

Current audited/frozen identities:

```text
19-file scaffold bundle          6569f79198e8923af362f268de4abd46d9a95ac51d29c9e9f7e7fbea62446739
secure v3 kernel                 791ad275a72736245de157fac5a58dfa7f6b08d1316d00f021b8df02a7a15754
remeasurement policy             98de86310a95057fcbd2c56be2014d825e6335312bb78a652ac90461c50f4393
remeasurement workflow           c5042ca7f77e9394249e4d8cc441e990f358c4fb34aaab0075767974fa4dfcbb
canonical analysis policy        17cea8e952081153960ce9783705ae49a067b16847e44a325dea6e9fcc52a201
optional high-N policy           cf1763dcce442d99973c3b1e17d7f75fdf93fd11a99225fb4b812906c1805bd6
launch receipt                   1212bac7864a9b16165b659168a947cf8b055eb1e89d9d65aa27dceb7edb6e1c
inventory job                    15168723
405-cell one-CPU array           15168724 (maximum 16 simultaneous)
four-process independent seal    15168725
```

At launch there were 168 idle CPUs across reserved nodes 196--199.  This
chain post-processed existing FITS; it did not perform new MCFOST ray tracing.
The inventory passed before the array was released by dependency, and the seal
independently remeasured all 405 cells before ranking was allowed.  The
analysis was submitted with the externally captured index SHA and the
canonical/optional policy SHAs above.  Its final continuum figure shows all
312 unscored JWST bins in grey plus exactly nine black scored anchors.

For the following ice stage, the anchor design contains 44 evaluated
optical-depth targets per aperture: 11 each around H2O 3.0, CO2 4.27, CO 4.67,
and CO2 15.2 microns.  Of these, 36 feature-shape points enter the ice score
and eight BC/RC wings constrain local continua; 16 additional targets support
the local continua, and the nine broadband continuum anchors remain a
separate observable.  The missing full-44 scorer, exact 239-task wider
finalizer/index, fixed `sigma_tau,sys=0.10` screening floor, tau/residual plots,
312-grey+9-black continuum plot, run-local anchor-PNG provenance, task-row
semantic rederivation, and core/wider Huber fail-closed gates are now
implemented and independently audited: 57/57 tests pass.  Nothing was frozen
or launched.  The intentional launch blockers are the canonical continuum
winner, Stage-1 opacity selection, and a 1x16 cluster affinity smoke test.

## Terminal handoff: task 371 recovery failed closed by 0.01006 percentage points

The prospectively frozen task-371 recovery is terminal.  Both 2.048M-packet
Method-2 realizations completed on `s92node198`, passed every physical and
component-closure gate, and have exact launch/Slurm provenance.  Their signed
1-arcsec aperture fluxes at 140 pc are:

```text
seed 42001             1.1981590348472051e-14 W m-2
seed 42002             1.1684572414401619e-14 W m-2
symmetric difference   0.025100641540113086 = 2.5100641540113086 percent
frozen maximum         0.025 = 2.5 percent
```

The pair therefore misses the frozen gate by `0.000100641540113086` in
fractional units, or **0.0100641540113086 percentage points**.  Adjudicator job
`15168721` exited `2:0` deliberately and published a red decision: no pair mean
was adopted, no ranking was performed, and its required next action is
`stop_without_retry_and_publish_no_task371_replacement_flux`.  Do not relabel
this as green or silently relax the threshold after seeing the result.

```text
array task 0             COMPLETED/0:0, 16 CPUs, 01:25:10
array task 1             COMPLETED/0:0, 16 CPUs, 01:36:14
both physical records    green
terminal accounting      green
launch transaction       green
adjudication SHA-256      2dd718dd4a947a554a704770423ba28070558bcfc57e003300487ccb5fea6e08
```

The small adjudication, completion, measurement, closure, immutable-launch,
and Slurm-log records have been copied back locally with matching SHA-256
digests; the large FITS products remain on the cluster.  The new
`uniform_v2_v2` workflow must remain unfrozen and unsubmitted because its
prospective contract requires an independently green task-371 replacement.
Wait for a new scientific direction before defining any new recovery policy,
excluding a model, or changing the convergence threshold.

The locally repaired `uniform_v2_v2` scaffold has now passed its independent
code/provenance re-audit: 28/28 adversarial tests pass, including real inode
substitution and mid-read mutation rejection, and analysis policies are hashed
and parsed from one stable byte snapshot.  The exact readiness check with the
red adjudication SHA exits 2 as intended with no freeze.  A representative
serial seal benchmark projects about 0.26 hours for 405 cells locally, so keep
the current deterministic implementation and benchmark the real cluster seal
before adding process-level parallelism.

## Live handoff: task 371 recovery before uniform-v2-v2

The 405-cell cavity--mass--q source array is now terminal with exactly 402
valid v1 completions and three raw-only v1 quality failures: logical tasks
`[2,46,371]`.  Task 371 is model
`m2p75_c17p5_a0p40_q3p25_i0p05`, anchor w003 at
1.7588530996181813 microns.  It used 128k packets and 16 CPUs for 01:21:42.
MCFOST and the coeval SED/image closure completed successfully; the sole v1
failure was `direct_aperture_clip_symmetric_fraction_le_1e-3`.

A read-only v2 diagnostic was run as one 1-CPU Slurm job and wrote no
measurement, score, or rank product:

```text
diagnostic job                   15168717  COMPLETED/0:0 in 00:00:48
signed flux at 140 pc            1.0220288193235858e-14 W m-2
signed-vs-clipped discrepancy    0.0568376750786695 (diagnostic only)
full-plane negative fraction     0.19414671467963057 (diagnostic only)
physical v2 gates                all PASS
raw SED SHA-256                  c3006f025cc8881590946a8412afb3b15459db20ed7b26b4d707a7911a747dd5
raw image SHA-256                e3e28cf5f9d0bfcfc99eafd47f0c6abd107331d65158cfd5b25365b99d6c8c11
closure SHA-256                  eaf83a2c6fd0d6a94fe2e49c72dd538b887f56d4c4ca5e22668310a1223885b7
```

The single 128k v2 value is audit-only, not canonical.  The same-wavelength
task-2 precedent shows why: passing the relaxed diagnostic does not prove the
signed estimator has converged.  A separately named, prospectively frozen and
independently audited task-371 pair is now running with seeds 42001/42002,
2.048M packets each, one 16-core Method-2 process per seed, and the pre-existing
2.5-percent repeatability rule:

```text
workflow                         strict_continuum_cavity_mass_q_task371_highn_convergence_v2
source audit SHA-256             fac8c188121624b7e2a1880c20bfaaabf665022931cf8a757eeddd97933473ab
input manifest SHA-256           3b4156aef2e5a31aad17a219c8155941ed33befc4dc04a3286731439deefadb1
policy manifest SHA-256          b4a8ca7cc1972cb8dc3a9329d4ed8e175f825b3b9fb2af450fee201431480c50
workflow manifest SHA-256        040e78fd2a9c119ce26bf83a1816d8408ff771d8cb3b16e0a7c0019d94306605
Slurm array                      15168720 (tasks 0/1, each 1 x 16 CPUs)
Slurm adjudicator                15168721 afterany:15168720
launch intent SHA-256            b29501fdb9a917449e82c1990372eaf5418233bb2a47c6320d86084dca665877
array-submission SHA-256         a46d9ae328a8b1e33beb7d263e5507d727d2e623cbd90789ba7e0c4af12474a2
final receipt SHA-256            84d5f2cbad850bb5d95d4d54b7d1cb3d491f6fa7b09bc817da68e47d1143105c
```

At launch, nine whole 16-core slots were idle across the reserved nodes; the
two scientific tasks both started on s92node198, consuming 32 CPUs.  The
transaction is no-requeue/no-retry and binds the array job/task IDs, source
hashes, immutable receipt, and later `sacct` rows.  Adopt the arithmetic mean
only if both physical records and the pair decision are green; otherwise stop
without retry.

At the 03:21 CEST monitor, both recovery realizations were still healthy at
38 minutes on `s92node198`, each using about 15.7 of its 16 allocated CPUs and
about 4.0 GB RAM.  No completion or failure marker existed; the adjudicator
was correctly dependency-held.  Repeated Slurm IPMI messages concern only the
node energy meter and do not affect the calculation.

The old `strict_continuum_cavity_mass_q_uniform_v2_v1` chain correctly failed
closed because its frozen source partition required 403 completes plus only
failures `[2,46]`:

```text
source RT array                  15164641  terminal: 402 complete, failed [2,46,371]
old terminal inventory           15168712  FAILED CLOSED; no inventory published
old remeasurement array          15168713  CANCELLED; no element ran
old rank-free seal               15168714  CANCELLED; never ran
old uniform-v2 policy SHA-256    931bca7c6ee4d131fa3a7b86535e8e59fa29429b8fa447cfa0641eae8110e105
old uniform-v2 workflow SHA-256  c400a143e1f8f6f69e963a48fdd511c8d71b20534ed99bc87b21eb909991e963
```

Preserve that failed-closed audit trail; do not patch or retry v1.  If the new
task-371 pair is green, freeze a new `uniform_v2_v2` contract with exact
terminal partition 402 + failures `[2,46,371]`, uniformly remeasure all 405 raw
triples under v2, and canonically substitute only the three independently
green high-N pair means.  Ranking remains forbidden until the new 405-entry
index is sealed.  The later final figure must contain all 312 collapsed JWST
R~100 bins as unscored grey context behind exactly nine black scored anchors.

The corrected `uniform_v2_v2` scaffold remains deliberately unfrozen, unsynced,
and unsubmitted, but the provenance repairs are now implemented locally.  Its
source inventory is externally producer-anchored, every raw path is required
to equal the deterministic source path, and the held remeasurement transaction
binds one native array ID plus exact 405-row successful one-CPU accounting.
All 405 v2 measurements are independently recomputed from their bound FITS
triples at sealing; existing measurements cannot be reused.  The canonical
partition is described correctly as 402 source logical values plus three
high-N pair means derived from six realizations.  Twenty-eight uniform-v2-v2
contract/analysis tests and the independent re-audit pass.  Execution remains
correctly blocked because the now-terminal task-371 adjudication is red, so
this prospective contract cannot be frozen without a new scientific policy.

For the following ice pilot, keep the nine continuum anchors as a separate
continuum check and add a frozen ice-anchor table inside H2O 3.0, CO2 4.27,
CO 4.67, and CO2 15.2 microns.  Each band needs continuum wings, shoulders,
and trough/core samples; ice scoring is through matched 0.35-arcsec optical
depth rather than anonymous pooling with the continuum score.  Launch its
Method-2 RT array with one 16-core process per cell and all whole 16-core
reservation slots measured immediately before submission.

The rank-blind ice-anchor observation contract is now frozen and independently
audited.  It contains exactly 44 native-profile samples: 11 per band, ordered
as blue local-continuum wing; blue 25/50/75/90-percent shoulders/core; exact
peak; red 90/75/50/25-percent core/shoulders; red local-continuum wing.  The
eight local wings are nuisance constraints for each ice profile and are not
the nine broadband continuum anchors.  The latter remain a separate 1-arcsec
flux observable; the ice rows remain a 0.35-arcsec optical-depth observable.

```text
builder       build_ice_feature_anchor_contract.py
table         output/ice_feature_anchor_contract_v1.ecsv
manifest      output/ice_feature_anchor_contract_v1.json
diagnostic    output/ice_feature_anchor_contract_v1.png
table SHA     7625211b864188bf251e724a7ce73ef16c37f3b726bfce187896d7eb611b3de0
manifest SHA  1e12b459f1748d42aca555926ae64aa23af3b8b8e9dca47c3666a38b79ccadc1
tests         6/6 PASS; independent audit PASS
```

Do not launch the preliminary 42-row hand-selected pilot table: an independent
audit found that it diverged from the authoritative contract and crossed some
source masks.  The pilot scaffold is being corrected to consume the exact
44-row ECSV.  Any dense spectral samples required for an actual JWST LSF
convolution must be carried as explicitly unscored support wavelengths; 44
isolated monochromatic samples alone cannot implement a literal LSF operator.

The corrected local `ice_feature_method2_pilot_v2` scaffold now hash-binds the
exact 44 ice targets as 11 per H2O 3.0, CO2 4.27, CO 4.67, and CO2 15.2 micron
band.  Five LSF cells are generated per target, with 16 separately labelled
local-continuum targets (80 more LSF cells).  These ice-band anchors are the
explicit scored information requested for the next run: their five-node fluxes
are converted to matched 0.35-arcsec optical depths through a six-point local
continuum.  The nine broadband `w001`--`w009` anchors remain a separate 1-arcsec
continuum check; discarded `w000` is still absent.

The five audit blocker classes have now been repaired locally.  The continuum
adapter independently checks the exact 45-model score/rank permutation and
unique argmin, validates the authoritative temperature completion, and binds
the same temperature hash through all nine winner RT tasks.  Stage 1 is the
exact six-cell representation-family x Mie/DHS space at `amin=0.03`, the
continuum winner's `amax` and q, and exact band-specific 95/5 material rosters;
each score binds four terms and distinct kappa/lambda products.  Launches use a
single-use pre-submit lock, held intent/job receipt, exact command, and one
array identity.  The core `afterany` finalizer preserves terminal accounting
even on failure but only seals an index when all 140 elements are
`COMPLETED/0:0/16 CPU` on the reserved partition and nodes.  Sealing then
independently recomputes all 140 dual-aperture measurements from the bound FITS
in parallel.  Fast validation cannot launder the original raw-input-rehash
attestation, and wider tasks now check the immutable green digest without
repeating the 140-FITS analysis.  Launch aborts never infer a job ID from
malformed `sbatch` text.  Terminal-accounting capture is resumable across all
audited timeout, query-error, and interrupted TSV/state publication cases;
the seal independently re-queries Slurm and requires the same canonical 140
terminal rows.  Thirty-four adversarial tests and the final independent audit
pass.  The anchor manifest SHA-256 is
`73b078b4830aa26a39281ea6bd6580ea6c147215bae0b3aaf8be26b5e28431e7`.

This is still deliberately **not launchable yet** because the actual sealed
uniform-v2-v2 continuum winner and rank-blind Stage-1 opacity ledger/report do
not exist, so no ice input/task manifest has been frozen.  Revalidate the
actual handoffs after those upstream products exist, and run a 1x16 PSMN
affinity smoke test before submission.  Gaussian LSF and reused temperature
remain screening approximations, not a final production likelihood.

## New ice-input package audit

`DATA/mcfost_ice_inputs` has been inventoried and a staged integration plan is
recorded in `ICE_INPUTS_INTEGRATION_PLAN.md`. The package is usable, but is not
a complete TMC1A RT model: its `Dust/*.dat` files are laboratory H2O, CO2, and
CO optical constants; its parameter files are opacity-only placeholders; and
its `constraints/` products are JWST-derived fit targets.

The ice constraints use a 0.35-arcsec aperture at
`04:39:35.2258 +25:41:44.0392`, whereas the continuum fit uses a 1.0-arcsec
aperture at `04:39:35.221679 +25:41:44.181522`. The centers differ by 0.15283
arcsec. Keep these as two explicit observation operators: retain the 1-arcsec
continuum likelihood and initially score the supplied ice optical-depth
profiles through a matched 0.35-arcsec Method-2/PSF/LSF aperture calculation.

Do not use the package's absolute flux yet. A direct 0.35-arcsec re-extraction
is within about 3.6--5.5 percent in the NIRSpec ranges but differs by roughly a
factor of three near the MIRI 15-micron band. The named package generator and
its aperture-correction/background rules are missing. Normalized optical-depth
profiles are safer but remain aperture sensitive at roughly the 1--6 percent
integrated-depth level in a preliminary comparison.

Do not paste the ice-only mass fractions into the current envelope dust block.
They are normalized among the three measured ices, while the production model
uses a 95/5 silicate-core/generic-ice-mantle volume mixture. CO constants cover
only 4.547--4.996 microns and cannot support a broadband equilibrium run.

After the current continuum array is terminal and the continuum winner is
frozen, run a local opacity-only material/shape screen at the winning `amax`
and q, then a one-model sparse-wavelength Method-2 ice pilot with a new
versioned 0.35-arcsec kernel and two-seed convergence checks at the deep CO2
4.27-micron trough and narrow CO core. Only then define a small ice-abundance
and grain-mixing grid; do not launch an ice MCMC directly.

The reproducible package audit and shape-only plot are complete:

```text
audit_ice_input_package.py
plot_ice_opacity_vs_observed.py
output/mcfost_ice_inputs_audit.json                  PASS
output/ice_opacity_vs_observed_diagnostic.png
output/ice_opacity_vs_observed_diagnostic.json       PASS
package tree SHA-256  f270d65cae26ee4f124745594739cdc4bd43096af68214106fc10330996a29a3
```

The audit found 71 files, validated all 33 FITS files, reproduced the exact
negative-k clipping counts, and found no errors. The normalized diagnostic
confirms that the package's pure-ice Mie baseline is only a starting point:
CO is much too narrow and the 15.2-micron CO2 opacity is too blue and narrow.

## Critical live update: high-N v2 repair is green; source array continues

The original Method-2 source array `15164641` is still running at four
concurrent `1 x 16`-CPU cells.  At the 20:39 CEST snapshot it had 223 atomic
`complete.json` records and exactly 223 `measurement.json` records, the same
two known failed cells (tasks 2 and 46), and four active cells.  Query the live
state rather than treating this count as final.

The separately versioned, rank-blind high-photon convergence experiment for
those two failures has completed successfully:

```text
run root                 runs/strict_continuum_cavity_mass_q_highn_convergence_v2
array / adjudication     15165367 / 15165368, all exit 0
design                   tasks 2 and 46 x seeds 42001 and 42002
per realization          2,048,000 packets, Method 2, 1 x 16 OpenMP
input manifest           28486f9127c0e20c298cfe2fac5c7501899ec5f23169c6dde2d95730c15cacd8
policy manifest          422fb0650859271c2bd5af3d8c3b8f5751c08fc47648e80018cddd23f568f2e0
workflow manifest        af2553119555580276e5c84aa62afe4e0aed4421ee70468cd98099d943368671
adjudication             5daa1532f63dbb75a57fbea8f6f2d302cbfe47143bf408e77bd9ece547c05e67
```

Both pair decisions are green under the prospectively frozen 2.5%
repeatability limit:

```text
task 2  signed140 = 1.59470044898e-14, 1.59056890811e-14 W m-2
        repeatability = 0.00259415478 (0.2594%)
        adopted mean140 / mean147 = 1.59263467855e-14 / 1.44456660186e-14 W m-2

task 46 signed140 = 5.52094908575e-15, 5.51967190063e-15 W m-2
        repeatability = 0.000231361102 (0.02314%)
        adopted mean140 / mean147 = 5.52031049319e-15 / 5.00708434756e-15 W m-2
```

The v2 estimator remains the signed full-Mueller Method-2 total-I plane-0
aperture flux.  Zero-clipped and component-wise fluxes were not used.
Native-pixel clipping/negativity are diagnostics only; finite-flux, coeval
closure, geometry, support, edge, convolution, and physical component gates
all remained mandatory.  No extra seeds or retry-until-pass were used.

The next eligible action, after the source array is terminal, is to freeze and
run a separately versioned *uniform v2 reprocessing policy*: remeasure all
existing source products under v2 and substitute only the prospectively
adjudicated high-N means for logical tasks 2 and 46.  No grid ranking has yet
been produced from this repaired path.

## Historical v1 state: the original grid is not rankable by itself

The original Method-2 array `15164641` is still running and should be allowed
to finish for its valid scientific evidence.  Its throttle is back at four
concurrent `1 x 16`-CPU cells.  Query current atomic marker counts; do not use
the older progress snapshot below.

Two source cells have failed closed so far:

```text
task 2   model m2p25_c15p0_a0p40_q2p75_i0p05   w003 1.7588531 um
task 46  model m2p25_c15p0_a0p40_q3p00_i0p05   w002 1.4987957 um
```

Both completed MCFOST and coeval Method-2 closure, but the unchanged
measurement kernel rejected their faint signed maps.  No source
`measurement.json` or `complete.json` was published for either cell.

Task 2 then used its one prospectively frozen recovery attempt: independent
seed 41002, 128k packets, unchanged physics/kernel/gates, job `15165201`.
That attempt also failed closed.  Its exact diagnostics are:

```text
recovery signed / clipped aperture difference    0.03364285038 (limit 0.001)
source--recovery signed-flux difference           0.08674766906 (limit 0.025)
all closure/component/geometry/support gates      pass
measurement / adjudication / completion markers  absent
```

The frozen policy permits no second seed, retry-until-pass, clipped-flux
substitution, dropped anchor, or 512k fallback.  Therefore a complete
405-cell canonical index and ranking cannot be published from this grid as
currently defined.  Keep the valid source array running, but do not let the
old `afterok` sealer `15164642` or handoff `15164744` be treated as a valid
analysis route.  A later scientific decision is required: redesign the
near-IR estimator/quality contract or run a newly justified production
experiment.  Do not silently relax the frozen run.

Task 46's original v1 recovery did not run.  Two source-verification wrappers
exposed provenance-only decimal/sign mismatches before writing recovery
output.  A corrected v3 source-evidence package was later used only to bind
the independent high-N v2 experiment described above; the obsolete v1
recovery route remains unused.

Prospective terminal-accounting policy (frozen before recovery results):

```text
runs/strict_continuum_cavity_mass_q_v1_terminal_recoveries_v1/
aggregation_policy_manifest.json
SHA-256 a2aceb2f19619758004c0f851928afbd537dcc36eb8fc72b27a93a0ccffa950b
```

## Live cavity--mass--q production run (current state)

The previous 3-D-q recovery below completed successfully and selected
`m2p50_a0p40_q3p25_i0p05` as the native-nine-anchor provisional best
(`0.1159514582 dex`).  The next fixed, rank-blind screen is now live and
varies the two supported directions plus cavity geometry:

```text
mass factor                  2.25, 2.50, 2.75
envelope size exponent q    2.75, 3.00, 3.25
cavity half-angle state     15, 17.5, ~20, 22, 24 deg
models                       3 x 3 x 5 = 45
fixed amax / ice fraction   0.40 um / 0.05
fixed inclination           70 deg
photons / seed              128000 / 41001
scientific anchors          w001--w009 only (nine per model)
RT cells                     405, Method 2, 2401 x 2401
```

Frozen roots and hashes:

```text
runs/strict_continuum_cavity_mass_q_v1
  design config       f49c95212f0039dda207797d046211f9406c8e925620e911ed3e0e4736411a8d
  input manifest      2b4424009fcf9ad5e8f25a672d0a76ff84f835c1770beb52b0cd53be20bb7ed2
measurement kernel    8cfb9aa37e1a4eb6410d84ee2b66b58a00d574f43d94720ab0e513d69a4d951a
analysis policy       39ff20e66423a5b054ab79b38653ec1b4524aef99ed8420fa8dd5e59c9c0d638
```

Pre-launch validation passed: all 114 frozen files and 945 table references,
all 450 task mappings, the exact four-process analyzer, and exact-product
legacy/kernel/FITS-adapter parity on both 2401-pixel reference products.
The final parity report is
`output/strict_continuum_measurement_kernel_v1_cluster_parity.json`, SHA-256
`db13d0322c976a5a94109c81cac157f8cddd28150add0d41c89a0ca464c6fd82`.

Original Slurm chain:

```text
temperature array        15164640   COMPLETED 45/45
Method-2 RT array         15164641   RUNNING, 405 cells
completion-index job      15164642   CANCELLED intentionally (stale v1 route)
post-index handoff        15164744   CANCELLED intentionally (stale v1 route)
```

At the real submission snapshot the reserved nodes exposed exactly four whole
16-CPU slots (`[1,2,1,0]` by node), so the uncapped array launched at
`4 x 16 = 64` CPUs.  Launch receipt SHA-256:
`39d98869baec71ad8124fcfaee2f9be3548e3bf2bef22d427d5b73ab1bdf166b`.
At 11:02 CEST the RT grid had 3/405 atomic completion records, exactly three
measurements, four active 16-core tasks, and no failures.  Query current state
rather than treating that progress count as final.

The original `launch.py --submit-analysis` hook omitted the required analysis
policy environment hash.  Do not use it for this run.  The tested
orchestration-only job `15164744` fixes the handoff without modifying any
frozen scientific input: it validates the 45x9 completion-index seal and all
pinned hashes, then submits the existing four-process `job_analysis.sh` with
the policy, manifest, and index hashes and writes an atomic receipt.

Expected final products (currently absent) are
`output/strict_continuum_cavity_mass_q_v1_{predictions,scores,regions}` plus
the validation JSON, summary JSON, and PNG.  This is still a provisional
fixed-normalization screen, not a formal likelihood or final MCMC result.

## Current 3-D-q production state (supersedes the old stop state below)

The frozen 3-D screen sampled 36 logical models in
`(mass factor, envelope amax, grain-size exponent q)`. Thirty temperatures and
all 300 fresh one-wavelength Method-2 RT tasks completed successfully; 60
additional cells are hash-pinned reuses, giving a planned 360-cell matrix.

Canonical analysis job `15161574` deliberately stopped before predictions or
ranks because exactly one cell failed an unchanged numerical gate:

```text
task / model     24 / m2p00_a0p40_q3p75_i0p05
anchor           w004, 2.1916659199557045 um
failed metric    signed-vs-zero-clipped aperture symmetric fraction
measured / limit 0.0133539125 / 0.001
flux impact      0.0057996167 dex
```

All component maps, coeval Method-2 closure, edge support, and PSF-convolution
checks passed. The original result and logs are preserved under
`output/provenance/strict_continuum_3d_q_v1_failed_job15161574/`. The stale
dependent primary-nine job `15162152` was cancelled.

The rank-blind recovery policy is fixed before seeing another result:

1. Run exactly one separate seed-41002, 128k, 16-core realization of task 24;
   reuse the exact temperature and change no physics or sampling setting.
2. Keep the signed total-I plane-0 estimator and every original quality gate.
3. Also require the recovered signed aperture flux to agree with seed 41001
   within symmetric fraction 0.025 (the established 0.01085793-dex operational
   margin).
4. If it passes, replace only logical cell `(grid model 2, w004)` in a new,
   versioned canonical provenance chain. If it fails, stop: no seed shopping,
   no clipped/component substitution, and no task-specific gate relaxation.
5. The scientific score is the **nine-anchor** score using only `w001--w009`.
   `w000` is diagnostic context only and must not be ray-traced or scored in
   the next production grid.

After the recovery-aware nine-anchor result is validated, choose the companion
parameter from its landscape and make the next production grid vary cavity
opening angle. Immediately before that submission, count the reservation's
whole 16-CPU slots and set array concurrency to use all of those slots (not
fragmented remainders). At the latest check the four reserved nodes had
24/40/8/0 idle CPUs, i.e. three whole 16-core slots (48 usable cores).

## Stop state

- There are **no active, pending, or held Slurm jobs** for `uhlo01`.
- No new production calculation has been launched.
- The current result is a validated **Method-2 strict-continuum screen**, not a
  final parameter measurement.
- Wait for a new parameter-grid decision before submitting more work.

## Scientific decisions that are now fixed

1. Use Cartesian **MCFOST ray-tracing method 2** for both the SED and image.
   Do not return to the stock log-polar Method-1 SED or the old mixed
   `SED_total * image aperture fraction` estimator.
2. Use the strict continuum mask first. Ice and silicate feature regions stay
   excluded until the continuum geometry and dust-size distribution are
   settled.
3. Use the signed plane-0 flux in the exact **1.0 arcsec radius aperture** at
   CRPIX, with the established component-wise/coeval Method-2 quality checks.
4. Keep absolute flux normalization fixed; there is no fitted multiplicative
   amplitude.
5. Model distance is 140 pc and predictions are rescaled to 147 pc by the
   fixed inverse-square factor.

The ten strict anchors are:

```text
0.995, 1.203, 1.499, 1.759, 2.192,
2.546, 3.914, 5.391, 13.800, 27.513 um
```

## Current validated result

The support-remediated refinement is complete and valid:

```text
canonical cells                    90/90 valid
original unaffected cells          72/72 valid
expanded-support cells              18/18 valid
old/new support comparisons         18/18 pass <2%
worst aperture difference           0.05344%
```

Provisional ranking:

```text
rank  model                    mass factor   envelope amax   score
1     m2p25_a0p50_i0p05       2.25          0.50 um         0.134462755 dex
2     m2p50_a0p50_i0p05       2.50          0.50 um         0.155557892 dex
gap                                                         0.021095137 dex
```

`i0p05` means **ice fraction 0.05**, not inclination.

All three boundary-weight sensitivity cases keep the same rank-1 model. The
runner-up remains in the promotion union. The result is still provisional
because the predeclared independent high-photon confirmation has not run.

Main result bundle:

```text
output/strict_continuum_refinement_v1_support_v1.png
output/strict_continuum_refinement_v1_support_v1_scores.ecsv
output/strict_continuum_refinement_v1_support_v1_predictions.ecsv
output/strict_continuum_refinement_v1_support_v1_regions.ecsv
output/strict_continuum_refinement_v1_support_v1_summary.json
output/strict_continuum_refinement_v1_support_v1_validation.json
```

Important hashes:

```text
support input manifest
  a281f4ef937067662f9a1dd76d583a09732dc11176c923b1e3eef164bc0c538a
analysis manifest
  a2d6ad7d1d2aeb2f3e309025a067cbd85e9140ed7184725bb89bc84da0e7454d
summary JSON
  5f7fa38a8101f27272fe3d5c6b67af63b68881177409627f327136549538d578
validation JSON
  84fa5a05c06c2ad0c71c5508c661810d7e7762d4326a61c660678340c146ef47
fit figure
  0da7e03cc52f408aca16bafde47f2e569e1bc600dfa74177b60b525c62319e53
```

## What has actually been fitted

The scientifically valid strict Method-2 work is currently **two-dimensional**:

| Parameter | Values sampled |
|---|---|
| Envelope dust-mass factor | 1.25, 1.50, 1.75, 2.00, 2.25, 2.50 |
| Envelope maximum grain size | 0.10, 0.20, 0.35, 0.50, 0.60, 0.70 um |

The initial strict grid sampled mass `1.25--2.00` and `amax=0.10--0.50 um`.
The refinement sampled mass `2.00, 2.25, 2.50` crossed with
`amax=0.50, 0.60, 0.70 um`.

Held fixed in the strict fit:

```text
envelope ice/core fraction     0.05 / 0.95
envelope amin                  0.03 um
grain-size exponent q          3.5
inclination                    70 deg
total source luminosity        1.9 Lsun
cavity                         tokens 274.7, 100, 1.0 (~20 deg convention)
disk/star/envelope geometry    fixed
foreground extinction         not fitted
flux amplitude                 not fitted
```

Older exploratory work varied mass, `amax`, ice fraction, inclination,
luminosity, and cavity opening. Those grids used superseded Method-1 or mixed
estimators and are **directional history only**, not constraints.

## Why the present two-parameter fit is not final

For the nominal best model, the continuum residual changes sign through the
near-IR:

```text
model too bright   around 1.2--1.8 um
model too faint    around 2.2--5.4 um
close              near 1.0 and 27.5 um
```

The forced diagnostic reduced chi-square is about `11.6` for ten anchors and
two fitted coordinates. It is not a formal likelihood because the scorer uses
region-balanced log residuals, systematic floors, and no calibrated MCFOST
Monte-Carlo covariance. Nevertheless, the structured residual shows that
mass or luminosity alone cannot repair the continuum shape.

## Required next calculation

Before expanding the physics, confirm both promoted models independently:

```text
models                 m2p25_a0p50_i0p05 and m2p50_a0p50_i0p05
anchors                all ten strict anchors for each model
seed                   41002
packets per anchor     512000
quality gates          unchanged
repeat acceptance      pointwise symmetric difference <=2.5%
```

This is 20 Method-2 RT cells plus any required pinned temperatures. It is a
numerical confirmation, not another parameter dimension.

## Direction after confirmation

Do not launch a broad MCMC yet. Use staged screens:

1. Locally bracket the coupled mass--`amax` valley. A sensible range is mass
   factor `2.10--2.55` and `amax=0.40--0.60 um`, concentrated near
   `2.25 / 0.50 um`. The mass direction is bracketed, but the joint optimum is
   not, because the winner is on the refinement's lower-`amax` boundary.
2. Add one dust-shape coordinate: grain-size exponent `q`, initially around
   `3.25, 3.50, 3.75, 4.00`, while keeping `amin=0.03 um`. This makes the fit
   three-dimensional: `(mass, amax, q)`. It directly tests the proposed
   small-grain explanation for the sharp short-wavelength rise.
3. If dust-only changes cannot fix the shape, run a small
   literature-constrained inclination/cavity sensitivity grid and constrain it
   jointly with JWST morphology. Do not let the SED alone fit all geometry.
4. Restore ice and silicate bands only after the continuum is satisfactory,
   then fit composition/feature parameters separately.

## Cleanup performed on 2026-07-24

### Local

`jwst_mcfost` decreased from approximately **13 GiB to 554 MiB**. Compact
plots, ECSV/JSON summaries, scripts, manifests, input spectra/masks, and all
current Method-2 provenance were retained.

Deleted obsolete raw run trees:

```text
runs/baseline_R100
runs/feature_grid_screen
runs/feature_grid_high_photon_sed
runs/full_R100_finalists
runs/full_R100_aperture_images
runs/dust_grid_screen
runs/dust_refinement_screen
runs/cavity_opening_test
runs/envelope_mass_test
runs/inclination_test
runs/literature_boundary_extension
runs/literature_geometry_luminosity_screen
```

Python caches and `.DS_Store` files under `jwst_mcfost` were also removed.

Retained locally:

```text
runs/strict_continuum_screen
runs/strict_continuum_refinement_v1
runs/strict_continuum_refinement_v1_support_v1
runs/matched_method2_closure_4p5um
runs/matched_method2_nearir_0p994676
runs/aperture_normalization_convergence
runs/aperture_spatial_convergence*
runs/component_rt2_pilot
runs/strict_continuum_task114_plane0_recovery
```

The small task-22 recovery manifests remain locally because the historical
strict-screen amendment pins their hashes, although the incomplete cluster
work directory was removed.

### Cluster

Cluster root:

```text
/scratch/Cral/uhlo01/jwst_mcfost
```

The project decreased from approximately **22 GiB to 17 GiB**. Removed:

```text
runs/feature_grid_screen
runs/feature_grid_high_photon_sed
runs/literature_boundary_extension
runs/literature_geometry_luminosity_screen
runs/_input_variants
runs/aperture_normalization_convergence/_failed_attempts
runs/component_rt2_pilot/_failed_attempts
runs/strict_continuum_refinement_v1/_failed_attempts
runs/strict_continuum_task22_recovery_512k
all remaining runs/**/*.tmp and runs/**/*.partial files
```

Verified after deletion:

```text
queued jobs for uhlo01          0
tmp/partial files in runs       0
strict screen                   retained
strict refinement               retained
support-remediated refinement   retained
matched Method-2 tests          retained
```

No `scancel` was necessary because the queue was already empty. The failed
analyzer jobs `15157823` and cancelled retry `15159551` are accounting history;
the exact retry `15159553` completed successfully and produced the valid bundle.

## Cluster environment

```text
SSH alias       PSMN_sr650node230
partition       Cascade-GPU
reservation     Houches2026
project root    /scratch/Cral/uhlo01/jwst_mcfost
MCFOST          /home/uhlo01/mcfost/bin/mcfost
MCFOST_UTILS    /home/uhlo01/mcfost/utils
Python          /home/uhlo01/.venvs/jwst_mcfost/bin/python
binary SHA-256  342bb88df2bc7232d7d28a4472047620f8357bfe5e5a561ca3750df7a62f1950
```

Normal execution is one MCFOST process with OpenMP threads. Parallelism should
come from independent Slurm array cells, while respecting currently idle
reservation capacity. Recheck `sinfo`, `squeue`, and visible CPUs immediately
before any new launch.

## Minimal restart checklist

1. Read this file.
2. Open `output/strict_continuum_refinement_v1_support_v1.png` and its summary
   JSON.
3. Confirm the Slurm queue is empty and inspect reservation capacity.
4. Freeze manifests and analysis policy before inspecting any new ranking.
5. Run the two-model seed-41002/512k confirmation.
6. Stop if it fails the unchanged quality or repeatability gates; otherwise
   proceed to the staged `(mass, amax, q)` exploration above.

## H2O 3-micron Method-2 screen (active 2026-07-26)

The accepted continuum geometry is now held fixed while the 3-micron H2O
feature is screened.  This is a versioned, screening-only calculation, not a
formal abundance posterior or a thermally self-consistent refit after changing
the dust mixture.

Fixed continuum baseline:

```text
envelope mass factor       2.25
cavity half-opening angle  17.5 deg
density exponent q         2.75
inclination                70 deg
amin / amax                0.03 / 0.40 um
H2O fraction               0.05
photon packets             128k
```

Both the stellar-centred and observed-ice-offset 0.35-arcsec apertures are
measured independently.  Treat the offset aperture as a positional
sensitivity unless the MCFOST image axes are explicitly registered to the
observed position angle.  The canonical estimator is signed total-I from the
full-Mueller Method-2 output.  Nine H2O anchors are scored in optical depth,
with two continuum wings and four local continuum support points used only for
the reconstruction.  Nine broadband points are report-only guardrails; the
27.51-um point retains the original dust because it lies outside the selected
H2O optical-constant coverage.

### Sealed inputs and opacity choice

```text
continuum handoff SHA-256
  4f7878eb995991811eddce0055deb8a0f09589edc54b302c20792fb04751a33a

selected opacity shape
  disjoint_populations__dhs_vmax0p8

opacity handoff SHA-256
  6b5a0753f35d2f63c7c3f1e393cd4aa57b6b2cb0960c1dd4af890e2099bc2e7e

screen manifest SHA-256
  4a29d3432cbfa5a0729c1b2ae8d5e0ce503dbcccd17d5ad715108de75ddd1b41

convergence task-table SHA-256
  996224cb6efaa31bf3027c37f87cfc0ffcd053816d4050d9c4e30f6a8bdf5029

wider task-table SHA-256
  c030cbfdae2182eaa4ed22db464c00aa247137fa1e74cda60c684699a44fa127

1x16 affinity seal SHA-256
  72fa1010462750048e275099433e4a033b1421c5d6bca349107e538956ebe38a

single-wavelength Method-2 smoke seal SHA-256
  9565b151762ec6d63f00247a88ee085d2388bd19750789d8e868c88067af93bb
```

The six-candidate opacity screen selected the disjoint-population DHS shape
with score `0.0946627844`; the Mie version was runner-up at `0.1230373358`.
The `0.0283745515` gap is a screening preference, not evidence that this is a
unique physical dust composition.

Two invalid/inadequate Stage-1 attempts were preserved under
`runs/h2o_3um_method2_screen_v1/_failed_attempts`: one combined DHS with a
coating, which MCFOST forbids, and one did not provide enough wavelength margin
for the endpoint LSF kernels.  The final valid opacity grid spans 2.54--3.66 um
while the scored band remains 2.55--3.65 um.

Three later execution-preflight failures are also preserved with complete
checksums under `_failed_attempts`; none is a failed science model.  The first
attempt encoded a numeric wavelength count as if it were a header, so MCFOST
correctly interpreted it as an unwanted 1-um wavelength.  The next smoke
attempt omitted the nested `paired/seed=52099/data_th` parent.  The final
failed smoke produced a valid one-wavelength SED and 33x33 RT image, but the
Python validator incorrectly treated the coeval Method-2 marker as a second
wavelength-table row.  The runner now writes a one-record MCFOST wavelength
list, creates and verifies the exact product parent before launch, and accepts
both real MCFOST completion formats (`24.84s` and `0h 9m 30s`) using a strict
ASCII Fortran-number grammar.  Both authoritative and independent parsers
reject extra rows, truncation, reordering, duplicate/prefixed diagnostics,
non-finite/negative times, and Python-only underscore numerals.  Local and
cluster suites pass 91/91; the corrected dynamic smoke job `15170917`
completed `0:0` and independently revalidated the resulting FITS products.

### Cluster executable and active jobs

This workflow uses the custom Method-2 SED executable, not the older binary
listed in the general cluster section above:

```text
MCFOST /home/uhlo01/mcfost-custom/v4.1.13-rt_sed2-92bf1b1e086e/bin/mcfost
SHA    d96bbc7cb71c09fce93edcef5908e6af1eb0e90e40193a3afdd259c2a7786401
Git    92bf1b1e086e433a2b18e16b70e428ff4b5a1bd2
```

Current Slurm chain:

```text
15170915  completed  1x16 affinity smoke
15170917  completed  custom-MCFOST one-wavelength Method-2 smoke
15170919  running    70-element convergence array, throttle 11x16 CPUs
15170932  pending    afterany convergence finalizer
```

The live logs explicitly report `SED spatial sampler = 2 Cartesian adaptive`
and `Parallelized code on 16 processors`.  At the first checked wave, CPU
accumulation corresponded to roughly 12--13 busy cores per model and memory was
about 1.9--2.0 GiB per model.  All first-wave tasks were advancing without an
error signature.  The IPMI energy-accounting warning on node 198 is benign.

### Exact continuation gate

1. Let all 70 convergence elements terminate; do not launch the wider phase
   from partial products.
2. Require finalizer `15170932` to create and validate the green convergence
   seal.  If it is red, inspect seed-to-seed and LSF reconstruction gates and
   repair/re-run only with preserved provenance.
3. If green, submit the frozen 49-element wider phase through
   `launch_h2o_3um_method2_screen_v1.py`, using all currently available whole
   16-core slots and binding the exact green-seal SHA.
4. Submit the analysis job with an `afterany` dependency on that wider array;
   it must independently require all 49 terminal-success records and revalidate
   the 70-cell convergence seal.
5. Copy back the compact ECSV/JSON/PNG bundle.  The canonical H2O plot already
   includes the full JWST H2O profile as grey unscored context plus the 11
   anchors (nine scored).

A presentation-only view of the sealed Stage-1 screen is now available at
`output/h2o_stage1_six_opacity_shapes_vs_observed.png`, with its provenance in
the adjacent `.provenance.json`.  It shows all six H2O opacity shapes, the
observed profile, and the selected disjoint-population DHS curve.  The required
sealed reports and twelve small opacity FITS arrays were copied back under the
local `runs/h2o_3um_method2_screen_v1` tree, so the plot can be reproduced by:

```bash
python plot_h2o_stage1_opacity_shape_presentation.py \
  --source-root runs/h2o_3um_method2_screen_v1
```

This figure is explicitly shape-only and independently peak-normalized; it is
not an abundance, radiative-transfer, aperture, or absolute-depth fit.  A
presentation-only grey-point variant of the final Method-2 result can likewise
be made after the sealed analysis without changing the frozen workflow.

The observational ice-package audit and a ready-to-send reproducibility
request are recorded in `ICE_ANALYSIS_REANALYSIS_REQUEST.md`.  The package is
for TMC1A, but its 0.35-arcsec ice-peak extraction is not interchangeable with
the accepted 1.0-arcsec stellar continuum extraction; all physical inferences
from this screen remain conditional on that distinction and the documented
segment/PSF/uncertainty issues.

### Prepared post-analysis validation

Two non-invasive local helpers are ready; neither is part of, nor changes, the
frozen v1 cluster manifest:

```text
sync_h2o_3um_method2_screen_v1_compact.py
supplement_h2o_3um_method2_screen_v1.py
```

After the sealed cluster analysis exists, the sync helper first performs a
read-only inventory and then requires explicit `--apply`.  Its exact 78-file
allowlist includes the final compact products plus the small parameter,
temperature, accepted-continuum, Stage-1 provenance, and the complete
six-product/staged-input lambda-contract smoke evidence needed for independent
validation.  It excludes every scientific task tree and all other RT/SED
products, follows no symlinks, never deletes, and never overwrites differing
local files.

The supplement then checks the exact 405-source continuum, 5+1 opacity screen
and selected `disjoint_populations__dhs_vmax0p8` representation, fixed
geometry, all 70+49 task semantics, two aperture coordinates, gate thresholds,
guardrail identities, the dynamic one-wavelength log/SED/RT proof, final table
arithmetic, and source/output hashes.  It
will make separately sealed grey-point H2O, convergence, and broadband
guardrail figures with the positional-sensitivity limitations stated on the
plots.  The complete local H2O handoff/opacity/workflow/supplement/sync suite
currently passes 91/91 tests.  Do not run the build step until the final analysis bundle has
been copied back.

## 2026-07-26 current state: v3 pilot green; postpilot production active

The older v1 chain described above was superseded by the immutable
`h2o_3um_method2_screen_v3` field pilot and the separate
`h2o_3um_method2_postpilot_v1` production amendment.  Do not resume the old
v1 wider-phase instructions.

The v3 field pilot completed 44/44 cells and passed both gates.  Its largest
field-size flux change was 0.0436 percent; its largest seed-to-seed flux
change was 0.4103 percent; and its largest core optical-depth seed difference
was 0.0009586.  The accepted execution is therefore Method 2, 128000 packets,
2401-by-2401 pixels, inclination 70 degrees, the fixed continuum geometry
`m2p25_c17p5_a0p40_q2p75_i0p05`, and the selected opacity representation
`disjoint_populations__dhs_vmax0p8`.

```text
v3 input manifest
  6faf36dd78b9cab8776d1634f6e2859782c2e4c2a8ba9b5e78aabfbc19999c79

v3 field-pilot green seal
  f031983db2157b5211d568ae1c5b5d060d5747283de77d21f2ee88f47af30628

postpilot input manifest
  25a96777ed3363a6581a58d5bf8b8fc3daf307865ca81f26284a87f0c40ba9ad

postpilot authorization
  0d35c1309ff48a225ccf0ec43870ea9d8687303e49cc74e79fdf2fae62e49dea
```

An earlier postpilot manifest (`a20d5f...`) was rejected before any science
submission because Slurm had purged historical finalizer job `15171295` from
its live controller.  The exact rejected package is preserved under
`runs/h2o_3um_method2_postpilot_v1/_failed_attempts/`.  The corrected
authorization has no stale scheduler dependency: it recomputed all 44 green
cells, required the archived `15171295 | COMPLETED | 0:0` accounting row, and
then completed itself as job `15171417` with exit `0:0`.  Production arrays
remain strictly dependent on `afterok:15171417`.

Active production chain:

```text
15171420  running   85 H2O-profile tasks, throttle 10 x 16 CPUs
15171421  running    9 broadband guardrails, throttle 1 x 16 CPUs
15171434  pending   afterany final analysis/fail-closed publication
```

This run ray-traces one Stage-1-selected dust representation; it is not a
six-candidate RT reranking or a free-abundance MCMC.  It measures the
aperture-matched 3-micron profile at 85 LSF nodes, scores targets w002--w010,
keeps w001/w011 as unscored wings, derives only an effective depth-scale
diagnostic, and reports nine broadband anchors as non-likelihood guardrails.
The finalizer must see exactly 85+9 successful immutable completions before it
atomically publishes the tables, score report, lineage, and grey-JWST plots.

## 2026-07-26 finite-field guard recovery prepared; production still active

The original all-green publication path above is no longer valid.  MCFOST
itself completed broadband guard rows 0 and 1 with exact SED/image closure,
zero negative-flux fraction, and acceptable convolution loss, but strict-v2
correctly rejected their active `scattered_star` maps because the positive
edge fractions were `0.00266566` and `0.00152273`, above the frozen strict
`<0.001` gate.  The n2401, zoom-about-1.9975 source field is only about
3003.75 au wide at approximately 1.25 au per pixel, so this is a finite-field
sampling failure at the two shortest broadband wavelengths, not a failed H2O
profile model or a reason to relax the measurement gate.

Last verified live state before cluster polling became unavailable to Codex
was 25/85 full-profile rows complete with ten active and no failures.  Original
guard rows 2, 3, and 4 were complete, rows 0 and 1 had only the edge failure
above, row 5 was active, and the rest were queued.  Job `15171434` is the now
obsolete all-green finalizer; cancel it before it wakes.  Its fail-closed code
cannot publish from this mixed terminal roster, but it should not consume a
node merely to rediscover that fact.

An additive workflow is ready as
`h2o_3um_method2_postpilot_guard_field_recovery_v1`.  It must not be frozen or
launched until all nine rows of source array `15171421` are terminal.  It will:

1. Freeze the exact 0--8 source accounting roster and admit only failures that
   independently reconstruct to the single known active-scattered-star edge
   signature.
2. Re-run only those failed wavelengths, with all science inputs unchanged,
   on the expanded 6000-au, 4801-by-4801 field (`zoom=1`) while nearly
   preserving the source pixel scale.  Each task remains Method 2, seed 52001,
   128000 packets, one process by 16 OpenMP threads, 48 GiB, and four hours.
3. Require every standard strict-v2 gate to be green.  The central n2401 crop
   is diagnostic only and can never rescue a red result.
4. Require the old-red versus new-standard signed one-arcsecond aperture flux
   to differ by strictly less than 2 percent (less than 1 percent preferred).
5. Reuse all 85 original H2O rows and every successful original guard row;
   substitute an accepted n4801 result only for its exact failed guard index.

The recovery and mixed-lineage analysis code currently passes 31/31 targeted
tests, Python compilation, shell syntax checks, two independent audits, and a
synthetic eight-product atomic publication/readback.  The older v3 and
postpilot authorities still validate byte-for-byte at:

```text
v3 manifest       6faf36dd78b9cab8776d1634f6e2859782c2e4c2a8ba9b5e78aabfbc19999c79
postpilot manifest 25a96777ed3363a6581a58d5bf8b8fc3daf307865ca81f26284a87f0c40ba9ad
postpilot auth     0d35c1309ff48a225ccf0ec43870ea9d8687303e49cc74e79fdf2fae62e49dea
```

The continuation order is exact: wait for full terminal source accounting;
cancel `15171434`; sync the recovery workflow and the three hash-pinned
Stage-4b authority files; create
`runs/h2o_3um_method2_postpilot_guard_field_recovery_v1`; submit the inventory
phase from that root; freeze/check the resulting policy and manifest; submit
authorization; then submit the capacity-gated recovery with `--with-finalizer`.
Capture every printed job ID and each `.sha256` sidecar.  Only after both the
original full85 array and the recovery finalizer are terminal green may the
composite phase be submitted.  The composite deliberately has no historical
Slurm dependency, because old jobs can be purged from the live controller;
both its launcher and runtime independently re-query exact full85 and
recovery-finalizer accounting and validate the immutable recovery seal.

Do not use the original output directory as a partial result.  Accepted
recovery products publish atomically under
`output/h2o_3um_method2_postpilot_guard_field_recovery_v1`, and the final
mixed-lineage science bundle publishes under
`output/h2o_3um_method2_postpilot_recovered_v1`.

## 2026-07-31 silicate demonstration anchors added to the H2O summary figure

Panel A of `output/h2o_best_model_vs_jwst_full_spectrum.png` now contains four
new, real MCFOST predictions at 8.3, 9.7, 11.0, and 12.5 microns, together with
the aperture-matched JWST values.  They use the unchanged current best H2O
model, reused continuum temperature, standalone image Method 2, the MIRI PSF,
and the signed plane-0 one-arcsecond aperture estimator at 147 pc.  They are
presentation-only: no refit, rescaling, likelihood contribution, or ranking
change was made.  The table and lineage are in
`output/h2o_silicate_demo_anchors_v1.ecsv` and its provenance JSON.

The model/JWST residuals are -0.566, -1.034, -0.556, and -0.160 dex in
wavelength order.  Thus the illustrative 9.7-micron point is only 0.092 times
the observed flux: this current dust prescription predicts a substantially
deeper silicate absorption minimum than JWST.  Straight green segments in
panel A only connect the four computed anchors; they are not a densely sampled
silicate profile.

The first attempted coeval SED+image execution stopped without products in the
high-opacity silicate complex.  It was superseded by the previously validated
standalone image-Method-2 path.  The latter deliberately has no coeval SED
closure and is labelled that way in the figure.  Array `15205086` produced all
four images; its 11.0-micron task stopped at 98 percent and was preserved under
the remote `_failed_attempts` directory, then task-only recovery `15205092`
completed cleanly.  All four accepted images pass the pixel-scale, aperture
support, edge-flux, convolution-loss, finite-flux, and positive-flux gates.
