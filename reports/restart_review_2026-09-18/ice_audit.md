# Ice and dust audit for the restart review

Audited 2026-09-18 from local frozen tables, completion indexes, parameter files, and scripts. No MCFOST jobs or remote queries were run for this audit. Paths below are relative to `jwst_mcfost/`. Numerical data for new summary plots are in `ice_audit.json` beside this report.

## Current endpoint

The latest locally recorded ice/dust science product is the **July 31 silicate demonstration**, added after the July 27 dust and H2O grids. The `CONTINUE_HERE.md` header still says July 27 but its final section records July 31. No locally published August/September ice results were found. There is **no completed joint continuum + H2O + silicate fit and no ice/dust posterior**.

The most recent model, `v02_vmax0p1_amax0p4`, uses the nominal continuum geometry `m2p25_c17p5_a0p40_q2p75_i0p05`, but changes the envelope dust to separate silicate and H2O grain populations with DHS shape parameter `vmax=0.1`, `amax=0.4 µm`, and H2O species mass fraction `0.04`. Its inclination remains 70°. It improves the near-IR over the first ice-screen dust, but worsens the H2O red wing and produces far too little flux near the silicate minimum.

| Dust model | Five near-IR anchor RMS, dex, 1″ radius | Nine H2O anchor RMS, optical depth, 0.35″ radius | Meaning |
|---|---:|---:|---|
| Accepted original Mie coated dust | 0.086054 | Not measured with this H2O contract | Best of these three for the original near-IR observable |
| Disjoint DHS, vmax 0.8, H2O mass fraction 0.04 | 0.378639 | 0.191586 | Best *computed abundance-screen* H2O value |
| Latest `v02`, DHS vmax 0.1, H2O mass fraction 0.04 | 0.133425 | 0.286721 | Near-IR improvement; H2O red-wing regression |

These scores have different units and apertures and must not be added or treated as one ranking. They are unweighted screening metrics, not reduced chi-square values.

## Important correction to the existing figure

`plot_h2o_best_model_vs_jwst_full_spectrum.py:211` labels the accepted original dust **“no ice”**. That is incorrect. The frozen original parameter `runs/h2o_3um_method2_screen_v3/parameters/original_dust_guardrail_rt.para` explicitly contains a single Mie envelope species with mixing rule 2 (coating), a Draine silicate core at **0.95 volume fraction**, and an `ice_opct.dat` mantle at **0.05 volume fraction**. The original dust was not ice-free; it had not been independently fitted to the new H2O observable. Use “original Mie coated dust (5% mantle volume)” in the restart report. The original 5% volume and the subsequent 4–5% species mass fractions are different physical quantities.

## What was actually explored

| Stage | Explored space | Execution/result |
|---|---|---|
| Initial supplied ice package audit | H2O 3 µm, CO2 4.27 µm, CO 4.67 µm, CO2 15.2 µm opacity shapes | Package opacities and observed profiles compared; **not full source RT fits** |
| Corrected observation v3 | Stellar-centered raw 0.35″ extraction; R=2400; six continuum supports; 11 logical targets, nine scored | Frozen common observation/model operator. Earlier v2 is quarantined |
| Stage-1 opacity representation screen | Mie coating; Mie/DHS Bruggeman EMT; Mie/DHS disjoint populations; generic-ice Mie coating control | Five eligible candidates plus one unranked control; all at amin 0.03 µm, amax 0.4 µm, q 2.75 and nominal 95/5 composition |
| Method-2 field/seed pilot v3 | Two fields × two seeds × 11 spectral nodes | 44 executions; expanded n2401 field passes; 22 matched field pairs |
| Fixed-abundance RT H2O screen | Selected DHS vmax 0.8, 95/5 silicate/H2O by species mass | 85 H2O spectral-node tasks + nine continuum guardrails; two guardrails recovered on larger fields |
| Abundance grid | H2O species mass fraction 0.02, 0.03, 0.04; reference 0.05 reused | 255 new tasks = three physical models × 85 nodes; all completed |
| Dust-axis grid | Full vmax 0.8/0.4/0.1 × amax 0.4/1/2 µm at i=70°; i=60/50° at (vmax,amax)=(0.8,0.4) and (0.1,1.0) | 13 physical variants × eight broadband wavelengths = 104 tasks; all completed. H2O fraction fixed 0.04 |
| Best-dust H2O band test | v02 only | 85 nodes; all completed; H2O score worsens to 0.286721 |
| Silicate demonstration | Same v02 at 8.3, 9.7, 11.0, 12.5 µm | Four actual image Method-2 predictions; unscored presentation anchors |

The initially planned multi-band `ice_feature_method2_pilot_v1/v2` folders contain anchor scaffolds/templates; this audit found no corresponding completed all-band source fit. CO and CO2 **were not fitted in the completed H2O production**. The gas molecule `co.dat` in the parameter file is not a CO-ice optical-constants component.

Task counts are wavelength/LSF computations, not independent physical models. Repeated failed attempts, pilots, and field recoveries must not inflate the explored physical parameter-space count.

## Opacity screen and abundance result

The corrected Stage-1 scores (`runs/h2o_3um_opacity_screen_v2/selected_representation_v2.json`) are:

| Representation | Normalized opacity-shape score |
|---|---:|
| Disjoint populations, DHS vmax 0.8 | 0.1006178870 |
| Disjoint populations, Mie | 0.1284595677 |
| Bruggeman EMT, DHS vmax 0.8 | 0.1509104632 |
| Bruggeman EMT, Mie | 0.2051217856 |
| Two-component coating, Mie | 0.2391977338 |

All 15 predeclared continuum/relative-gain scenarios keep that candidate order. This is an **opacity-shape selection**, not a continuum or radiative-transfer ranking. MCFOST rejects coating+DHS for the same grain, hence that Cartesian combination is absent. The nominal 95/5 is volume fraction for composite/coated grains but species mass fraction for disjoint populations; normalized shape screening deliberately does not establish a common absolute ice abundance across representations.

At fixed selected dust geometry and reused temperature, the abundance grid gives:

| H2O mass fraction | Model core τ | Nine-anchor RMS Δτ | Mean Δτ |
|---:|---:|---:|---:|
| 0.02 | 1.155725 | 0.873408 | −0.825769 |
| 0.03 | 1.714400 | 0.524259 | −0.502963 |
| 0.04 | 2.272147 | 0.191586 | −0.147480 |
| 0.05, reference | 2.829646 | 0.310950 | +0.174999 |

Observed core τ=2.349529. Interpolation gives a core match at x=0.041388 and a minimum RMS at x=0.042359 (RMS 0.179004). These are conditional screening values, not a physical abundance posterior. A single abundance adjustment does not repair the red-wing/core mismatch in this tested family. It is too strong to say abundance can never help under a different geometry, thermal solution, or ice material.

## Dust-axis result and physical interpretation

All 13 exact rows are supplied in `ice_audit.json`; the source is `output/h2o_dust_axis_grid_v1_scores.json`. Reducing vmax 0.8 → 0.4 → 0.1 at fixed amax 0.4 µm improves the near-IR RMS 0.378639 → 0.207007 → 0.133425 dex. The local opacity decomposition attributes the original regression primarily to extinction: at 1.2028 µm Mie → DHS(vmax=0.8) lowers κext to approximately 0.736 of the accepted-dust value, whereas coating → disjoint and changing generic ice to H2O_30K have smaller effects (0.987 and 0.962 respectively). Albedo alone predicted the wrong flux direction in earlier notes; full radiative transfer is the authority.

Every amax=1 or 2 µm point in this particular fixed-geometry/reused-temperature grid has a worse near-IR RMS than the amax=0.4 µm control. This is **conditional evidence**, not proof that the real TMC1A envelope lacks larger grains. Grain-size constraints from other wavelengths/locations/dust representations must be compared with that limitation. Inclination 70→60→50 brightens the already overbright vmax0.8/amax0.4 branch, worsening RMS 0.378639→0.434080→0.486601. The same lower inclinations partly help the overly faint vmax0.1/amax1 branch. Geometry and dust remain degenerate; only two dust branches were tested at 50/60°.

The latest v02 H2O core is τ=2.238778, only 0.110751 below JWST, but at 3.296169 µm it gives τ=0.152755 versus observed 0.594608. Its red/core ratio is 0.068232 versus observed 0.253075. Thus a simple increase in depth would overcorrect the core before fixing the red wing. The increase in total H2O RMS from 0.191586 to 0.286721 arises primarily on the red side, not a uniform offset.

## Silicate anchors and figure interpretation

| Wavelength, µm | JWST λFλ, W m⁻² | Model λFλ, W m⁻² | log10(model/JWST), dex |
|---:|---:|---:|---:|
| 8.3 | 4.666080e−13 | 1.268036e−13 | −0.565820 |
| 9.7 | 1.457749e−13 | 1.346878e−14 | −1.034355 |
| 11.0 | 2.685093e−13 | 7.457231e−14 | −0.556382 |
| 12.5 | 4.220591e−13 | 2.918943e−13 | −0.160148 |

At 9.7 µm the model supplies 0.092394 times the JWST flux: an approximately 10.8-fold deficit. These predictions use the same latest v02 dust and reused temperature. They were generated through standalone `-img -rt2`, with no coeval SED closure after the custom coeval path stopped without products in this band. They pass the recorded finite-flux, pixel-scale, aperture-support, edge-flux and convolution-loss tests. They illustrate the discrepancy; four connected points are not a dense silicate spectrum or silicate fit. Source: `output/h2o_silicate_demo_anchors_v1.ecsv` and its provenance file.

`output/h2o_best_model_vs_jwst_full_spectrum.png` is the most useful existing summary. Panels A/B use 1″ radius, including H2O images remeasured at 1″; C uses the separately scored 0.35″ H2O observable. The continuous-looking green line outside the band is connected sparse model anchors. The 27.5 µm control used in earlier guardrail plots remains an **original-dust** calculation; do not relabel it as the new H2O dust model.

## Exact dust and geometry definitions

The common disk is two separate Mie populations: Draine silicate 70% and Rouleau–Martin amorphous carbon (`ac_opct.dat`) 30% by species mass; amin=0.03 µm, amax=1000 µm, q=3.5, 50 size bins. Disk dust mass is 1e−4 M⊙, radius 0.3–100 au, gas/dust=100.

The continuum reference envelope has one Mie-coated silicate/generic-ice population (95/5 **volume**), amin=0.03 µm, amax=0.4 µm, q=2.75. The H2O replacement has separate pure-silicate and pure-H2O DHS populations distributed through the same envelope, initially 95/5 and then 96/4 **mass**. Both have the same amin/amax/q. The ice constants are `H2O_30K_Leiden_mcfost.dat`; “30K” identifies the laboratory constants, not a solved uniform envelope temperature. DHS vmax is the maximum hollow-sphere shape parameter; the separate material-porosity field remains zero. Calling vmax the envelope's directly measured bulk porosity would be misleading.

The frozen continuum geometry has envelope dust mass 2.25e−4 M⊙ (factor 2.25), inner/outer radii 1/3000 au, radial density exponent −1.5, cavity half-opening 17.5°, and viewing inclination 70°. **q=2.75 is the grain-size power-law exponent, not the radial envelope-density exponent.** The filenames' `i0p05` token denotes the original ice fraction, not inclination. The parameter file runs at 140 pc and the aperture/flux operator converts to 147 pc. The star is fixed at 4000 K, 2.5 R⊙, 0.5 M⊙. Dust sublimation is switched off in these parameter files. No separate ice survival/chemistry boundary was fitted.

## Reliability limits and reusable work

All completed ice/abundance/dust-axis experiments reuse the accepted original-dust continuum temperature, hash `1bd212f2f54dd8c1701a6064e129ba7fdb5721cb3fe296472433bbdd3b79ec1c`. Changing dust opacities without recomputing temperature is the central physical approximation. The published products explicitly set `temperature_self_consistent_for_changed_dust=false` and `formal_likelihood=false`.

The observation contract is nevertheless valuable: both observation and model are measured in Fν after the model λFλ conversion, use the same R=2400 Gaussian/5-node Gauss–Hermite operator, six-support robust continuum, and τ=−ln(F/continuum). V3 fixed v2's units, scale-invariance, seam-veto and spectral-sampling defects. Keep v3; quarantine v2.

The new H2O stellar centroid is 0.121823″ from the accepted continuum extraction center, and 1.353637″ from the new extended ice peak. The H2O observable is the stellar 0.35″ raw aperture, not the ice-peak aperture or the provider's aperture-corrected/background-subtracted spectrum. The 1″ versus 0.35″ data were deliberately never pooled. A restart should explicitly adopt/propagate a consistent center/offset and extraction operator across instruments. The current H2O flux-error columns omit the full shared continuum covariance, so a formal likelihood requires an additional uncertainty model.

The 44-task field pilot found max matched-field difference 0.043555%, max two-seed node difference 0.410306%, and core Δτ=0.000959. These validate that frozen pilot's numerical observable; they do not establish convergence for every new dust/geometry variant. Grown grains exposed field truncation in the proposed wider dust-band grid, which was not executed; the accepted v02 band passes on the established n2401 field.

Keep the raw corrected data package, v3 observation contract, sparse anchor/LSF tables, measurement kernels, immutable parameter/opacity files, canonical source128k continuum ranking, final ice/dust score tables, completion indexes, and a compact set of representative FITS/temperature products and provenance. The historical diagnostics constrain useful response directions and identify dust/geometry degeneracies, even when the preferred old parameter values are not adopted as priors.

## Prepared but not completed

`h2o_3um_selfconsistent_grid_v1_common.py` and its launcher/runner/test scripts exist (July 27). They define **three ice files × three inclinations**: H2O_30K, waterice210K, ice_opct at 70°,60°,50°, keeping vmax0.1/amax0.4/x0.04; each variant would recompute temperature before eight continuum and 85 H2O-node evaluations. No corresponding `runs/h2o_3um_selfconsistent_grid_v1` manifest or published output is present locally. Treat this as a prepared design, not a completed or current production run. Warm-ice opacity shapes looked promising locally, but there is no full RT result here proving that they fix the red wing or are chemically appropriate throughout the envelope.

## Primary source paths

- `output/h2o_dust_axis_grid_v1_scores.json`
- `output/h2o_dust_axis_band_v1_tau.json`
- `output/h2o_abundance_grid_v1_curve.json`
- `output/h2o_silicate_demo_anchors_v1.ecsv`
- `output/h2o_best_model_vs_jwst_full_spectrum.provenance.json`
- `output/h2o_3um_method2_postpilot_recovered_v1/h2o_profile_scores_v3.json`
- `output/h2o_3um_method2_screen_v3_field_pilot/h2o_method2_field_pilot_v3_summary.json`
- `output/h2o_3um_observation_v3/h2o_3um_observation_contract_v3.json`
- `output/reproducible_ice_h2o_audit_v1/reproducible_ice_h2o_audit_v1.json`
- `runs/h2o_3um_opacity_screen_v2/selected_representation_v2.json`
- `runs/h2o_3um_method2_screen_v3/input_manifest.json`
- `runs/h2o_3um_method2_screen_v3/parameters/original_dust_guardrail_rt.para`
- `runs/h2o_3um_dust_axis_grid_v1/parameters/dust_v02_vmax0p1_amax0p4_rt.para`
- `runs/h2o_3um_abundance_grid_v1/completion_index.json`
- `runs/h2o_3um_dust_axis_grid_v1/completion_index.json`
- `runs/h2o_3um_dust_axis_band_v1/completion_index.json`
- `runs/h2o_silicate_demo_image_anchors_v1/input_manifest.json`
