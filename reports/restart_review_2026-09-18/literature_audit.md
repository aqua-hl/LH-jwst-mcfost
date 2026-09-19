# TMC1A literature audit for the local restart review

Audit date: 2026-09-18. This compares published constraints with the preserved local model; it does not launch a fit or use the retired cluster. Values from different tracers are not a joint posterior.

## What should change before restarting

The strongest reason to reopen the geometry is inclination. The literature does **not** establish 70° as a uniquely confirmed inclination. The major results below favor roughly 53–65°, with differing model assumptions. A useful restart should compare this interval with the old 70° solution while allowing luminosity and envelope/cavity dust to respond. A poor score at a lower inclination with all other parameters frozen does not rule out that inclination.

| Quantity | Preserved local assumption/result | Published comparison | Restart implication |
|---|---|---|---|
| Inclination, 0° face-on | Usually fixed at 70° in the late grids | 55° ± 10° [H14]; about 65° [A15]; 53° with formal −0.3/+0.2° intervals [A21] | Explore about 50–65° first; retain 70° as a comparison. Do not impose A21's very small formal error on a different geometry. |
| Stellar mass | 0.5 M⊙ fixed | 0.53 (+0.20/−0.10) M⊙ [H14], 0.68 [A15], 0.72 [A21]; 0.45 adopted in the JWST study [H23] | 0.5 M⊙ is not an obvious contradiction. Match the mass and inclination assumptions when using dynamical measurements. |
| Distance | Model at 140 pc, aperture operator rescaled to 147 pc | Recent ALMA/VLA study adopts 141.8 pc [A24]; older papers use 140 pc | Adopt one distance throughout. Recompute both angular aperture mapping and flux, not just one. |
| Heating luminosity | 4000 K, 2.5 R⊙ imply about 1.44 L⊙ photospheric luminosity | Published bolometric values include 2.3 L⊙ [A24] and 2.7 L⊙ [A15, J25] | Check stellar plus accretion heating and far-IR energy output. Lbol is not automatically equal to photospheric Lstar or the luminosity inside the JWST aperture. |
| Envelope dust mass | 2.25 × 10⁻⁴ M⊙ inside 3000 au | H14 gives about 0.12 M⊙ envelope mass from a submillimetre continuum estimate with a different aperture | Direct comparison needs matched radius, opacity, gas/dust ratio and spatial filtering. With gas/dust=100, the local mass would correspond to about 0.0225 M⊙ gas; this conversion is an assumption. |
| Disk radius | Inherited disk model | Keplerian radius about 80–100 au [H14, A15] | Match a kinematic radius to the appropriate disk parameter; Gaussian widths and continuum extents are different quantities. |
| Cavity half-opening angle | 17.5° nominal late solution | Historical outflow full opening 30–40° [C96]; a recent H₂ S(1) study lists an opening angle of 39° [F26] | The old half-angle is not clearly inconsistent with the historical full angle. Forward-model the tracer and projection before using a tight cavity prior. |

The local entries above are documented in [ice_audit.md](ice_audit.md) and [continuum_audit.json](continuum_audit.json). They describe conditional grid solutions, not a simultaneous measurement of all listed quantities.

## Dust and extinction: comparisons that must remain separate

The ALMA/VLA analysis favors a disk grain-size branch near amax ≈ 0.12 mm over a branch near 4 mm when polarization is included [A24]. Those are disk constraints. Our envelope amax=0.4 µm, constrained mainly by scattered near-IR continuum under a fixed dust representation, is not measuring that same population. Use separate disk and envelope dust rather than forcing either measurement onto both zones.

Assani et al. obtain AV approximately 17–20 toward four blueshifted jet positions from [Fe II] ratios. Their inference of enhanced MIR extinction depends on whether the NIR and MIR lines have the same excitation conditions. The introductory ≥10 µm grain argument concerns generic RAT alignment models of protostellar envelopes; it is **not** a TMC1A-specific measured lower bound on envelope amax [AS25]. It motivates a competing dust hypothesis, not a hard prior.

The JOYS overview lists AV=33 for TMC1A using the silicate optical depth and AV=18.5 τ9.7 [J25]. That conversion is uncertain and refers to a different diagnostic/sightline from the jet ratios. Neither number should be added as a foreground screen without checking which envelope attenuation MCFOST already computes. An aperture-integrated scattered continuum also need not behave as a single extinguished background source.

## Calibration and extraction provenance

I independently read the primary FITS headers of all 15 local files in `../../../cubes/` on 2026-09-18:

| Local products | Count | CAL_VER | CRDS_CTX | PROGRAM header |
|---|---:|---|---|---|
| NIRSpec `jw04537-…_g140h/g235h/g395h…s3d.fits` | 3 | 1.20.2 | jwst_1464.pmap | 02104 |
| MIRI `TMC1A_ch1…ch4-…_s3d.fits` | 12 | 1.16.1 | jwst_1303.pmap | 01290 |

The NIRSpec filename prefix and PROGRAM header differ; retain both as provenance and investigate before relabeling. The local cubes are **not** identified by their headers as pipeline 1.9.6. A filename such as `JWST_ver196_center.mrt` cannot establish equivalence to these cubes by itself. For restart, regenerate/check the extraction against these actual products and record aperture, center, background, masks, PSF, calibration and wavelength stitching explicitly. The header check does not prove which file versions an earlier cached extraction used.

## Suggested restart objective

Build a reproducible, self-consistent disk/envelope model that matches the same-aperture JWST continuum and selected images within calibration plus numerical uncertainties, while remaining compatible with external geometry and energy constraints. Only then infer ice properties, checking H₂O and silicate simultaneously. Preserve the old successful numerical operator, spectra and conditional parameter trends as benchmarks; treat the old best parameters as starting hypotheses.

Before production sampling, compare a small set of literature-oriented inclinations and heating luminosities, verify dust composition and temperature consistency, and examine residuals at held-out wavelengths. Use the published ranges as documented hypotheses, not intersected exact bounds. In particular, the listed luminosity values do not form a measured confidence interval, and neither cavity measurements nor jet AV give an exact prior for our adopted analytic envelope.

## Sources and locations checked

- **H14:** Harsono et al. 2014, *Rotationally-supported disks around Class I sources in Taurus*. Section 4.3.1, Tables 2 and 5. [Primary manuscript](https://arxiv.org/pdf/1312.5716), [journal DOI](https://doi.org/10.1051/0004-6361/201322646). The alternative line-transfer mass in §4.3.1 is 0.8 ± 0.3 M⊙, illustrating method dependence.
- **A15:** Aso et al. 2015, *ALMA Observations of the Transition from Infall Motion to Keplerian Rotation around the Late-phase Protostar TMC-1A*. [Primary manuscript](https://arxiv.org/pdf/1508.07013), [DOI](https://doi.org/10.1088/0004-637X/812/1/27).
- **A21:** Aso et al. 2021, *Multi-scale Dust Polarization and Spiral-like Stokes-I Residual in the Class I Protostellar System TMC-1A*. §4.2 for inclination; §4.3 for mass. [Primary manuscript](https://arxiv.org/pdf/2107.10646). Inclination intervals are 10th–90th percentiles conditional on their axisymmetric continuum model.
- **H23:** Harsono et al. 2023, *JWST Peers into the Class I Protostar TMC1A*. Introduction and §3. [Primary manuscript](https://arxiv.org/pdf/2306.08380).
- **A24:** Aso et al., *Grain Size in the Class I Protostellar System TMC-1A Constrained with ALMA and VLA Observations*, manuscript posted 2024. Introduction, §§V.1 and V.4. [Primary manuscript](https://arxiv.org/html/2411.13044v1). Distance is adopted from Krolikowski et al. 2021, not a new parallax measurement by this paper.
- **AS25:** Assani et al. 2025, *Mid-infrared extinction curve for protostellar envelopes from JWST-detected embedded jet emission: the case of TMC1A*. Introduction, §4 and Table 3. [Primary manuscript](https://arxiv.org/html/2504.02136v1), [published DOI](https://doi.org/10.1051/0004-6361/202555016).
- **J25:** van Dishoeck et al. 2025, *JWST Observations of Young protoStars (JOYS): Overview of program and early results*, Table 2. [Journal paper](https://www.aanda.org/articles/aa/pdf/2025/07/aa54444-25.pdf), [institutional full text](https://research.chalmers.se/publication/547708/file/547708_Fulltext.pdf).
- **C96:** Chandler et al. 1996, *Compact Outflows Associated with TMC-1 and TMC-1A*, §4.3. [DOI](https://doi.org/10.1086/177971), [ADS paper scan](https://ui.adsabs.harvard.edu/scan/manifest/1996ApJ...471..308C).
- **F26:** Francis et al. 2026, *JOYS+: A JWST/MIRI survey of the evolution of H₂ winds and jets from low-mass protostars*. §3.2 and Appendix D. [Primary manuscript](https://arxiv.org/html/2604.13773v1). Its Table 1 adopts rounded Lbol=3 L⊙ and distance=142 pc; these support an approximate luminosity scale, not exact bounds.
