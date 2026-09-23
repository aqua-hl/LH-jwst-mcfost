# Laboratory silicate inputs for MCFOST

These are alternative **material optical constants**, not optical constants recovered by fitting TMC1A. The published n,k rows are preserved without smoothing, clipping, interpolation, scaling or fitted adjustment. No preferred composition has been selected.

## Files to transfer

Copy this directory, especially its `Dust/` folder:

| MCFOST material file | Material | Bulk density [g cm⁻³] | Tabulated wavelength [µm] |
|---|---|---:|---:|
| `Dust/Olivine_MgFeSiO4_Dorschner1995_mcfost.dat` | Amorphous olivine-stoichiometry MgFeSiO₄ glass | 3.71 | 0.20–500 |
| `Dust/Pyroxene_Mg05Fe05SiO3_Dorschner1995_mcfost.dat` | Amorphous pyroxene-stoichiometry Mg₀.₅Fe₀.₅SiO₃ glass | 3.20 | 0.20–500 |

Use **one** as a replacement for `Draine_Si_sUV.dat` in the intended dust population. The files are isotropic amorphous glasses, not crystalline olivine/enstatite tensors and not a preselected mineral mixture.

Primary source: [Jena laboratory database](https://www2.astro.uni-jena.de/Laboratory/OCDB/amsilicates.html), [olivine table](https://www2.astro.uni-jena.de/Laboratory/OCDB/data/silicate/amorph/olmg50.txt), [pyroxene table](https://www2.astro.uni-jena.de/Laboratory/OCDB/data/silicate/amorph/pyrmg50.txt). Cite **Dorschner et al. (1995), Astronomy & Astrophysics 300, 503**, *Steps toward interstellar silicate mineralogy. II. Study of Mg-Fe-silicate glasses of variable composition*. The [author-associated methodology record](https://www2.mpia-hd.mpg.de/HJPDOC/PAPERS/4-dor95.html) explains the laboratory reflectance/transmittance analysis and sampled optical constants. These are not 109 independent direct measurements. The source requests citation; no explicit redistribution licence was verified, so this package does not assign a new licence to those data.

The source pages establish the material densities. A laboratory measurement temperature was not verified for these particular tables: **do not call them “30 K silicates.”** The appended MCFOST header uses an adopted **sublimation temperature of 1500 K**, following the historical silicate setup; that is a modelling choice, not the measurement temperature. Density changes opacity per gram of dust and must not be silently retained at the Draine value.

`manifest.json` records source URLs, hashes, conversion assumptions and companion files. `raw/` preserves the downloaded tables. `laboratory_vs_Draine_nk.png/.pdf` compares the input indices; `laboratory_opacity_controls.png/.pdf` compares computed grain extinction. Neither plot is a source-spectrum fit.

## Current numerical checks

Three small isolated-dust calculations were performed with local MCFOST 4.1.13: 517 wavelengths from 0.1 to 3000 µm, including all exact zero-k knots in the existing water table and neighbouring samples. Each used at most two threads and a 120-second timeout. No source temperature solution, image, SED or cluster job was run.

| Prescription | Result |
|---|---|
| Laboratory olivine core + generic `ice_opct.dat` mantle | Passed: finite nonnegative extinction, absorption and scattering |
| Laboratory pyroxene core + generic `ice_opct.dat` mantle | Passed: finite nonnegative extinction, absorption and scattering |
| Laboratory olivine core + unchanged `H2O_30K_Leiden_mcfost.dat` mantle | **Failed during dust initialization: SIGSEGV, return code −11** |

All used coated Mie grains, amin=0.03 µm, amax=0.4 µm, q=2.75 and a 5%-volume ice mantle. The disk was excluded from these numerical checks so the opacity is unambiguously that of the envelope population. Successful initialization is not physical validation of the extrapolated tails or convergence of a TMC1A model.

Compact tables, logs, parameters and the failed-case evidence are in `checks/`. The original complete attempts remain in `validation/laboratory_silicate_checks_v1/` in the project. No failed attempt was retried with secretly altered optical constants.

## Important limitation of the existing laboratory H₂O file

The unchanged `H2O_30K_Leiden_mcfost.dat` is included for provenance and future work, **not endorsed for joint production use yet**. Its coverage is only 0.250040–19.940540 µm, and it contains 206 zero-k rows, predominantly at 9.008670–10.236450 µm. MCFOST's unguarded logarithmic interpolation makes these zeros numerically hazardous. The crash does not uniquely prove that zeros caused it; other interpolation/coating effects require diagnosis.

The stand-alone examples with `H2O30K_UNVALIDATED` in their name are configuration illustrations, not approved executable production models. The pyroxene/H₂O combination has not been tested; success cannot be inferred from the silicate-only controls. A zero-k policy and the unmeasured water-ice tails need an explicit decision and validation before combining both laboratory replacements. The original H₂O file has **not** been changed or floored.

## Use with an independent MCFOST installation

Keep the ordinary complete `MCFOST_UTILS` installation and point the optional custom override at this package:

```bash
export MY_MCFOST_UTILS=/absolute/path/to/laboratory_silicates_v1
```

The package contains `Dust/`; it is not a replacement for the complete MCFOST utilities tree. Alternatively, install the selected material files in that installation's `Dust/` directory under their unique names. Do not overwrite `Draine_Si_sUV.dat` globally.

In the intended `.para` population, change the silicate filename. For example, our coated **envelope** optical-index rows become:

```text
Olivine_MgFeSiO4_Dorschner1995_mcfost.dat  0.95
ice_opct.dat                            0.05
```

These are only the two component rows inside a valid `Mie 2 2` species definition; the fractions are **volume fractions**, not independent species mass fractions. Full standalone parameter examples are under `examples/`. They retain the disk's Draine silicate/carbon mixture and the historical geometry; replacing disk silicate too is a separate physical change. They do not reproduce the old two-population DHS H₂O setup.

For a manual source calculation, a new temperature solution is required. From the examples directory, the historical control also requires the following physical CLI options (not a new recommended fit):

```bash
mcfost olivine_mg50_generic_ice.para \
  -cavity 317.159480236321 100 1 \
  -Mdot 1 7.3011756855e-08 \
  -root_dir olivine_control_output
```

This command **would start a source calculation**; it has not been run for this package. Set OpenMP/resource limits appropriate to your machine first. Do not reuse an old temperature file after changing dust.

## Use with this project's grid runner

The runner now supports explicit envelope-material selection. Initial controls that change the silicate only are:

- `config/grid.lab_olivine.example.json`
- `config/grid.lab_pyroxene.example.json`

Each contains one physical model, using the same 73-wavelength observation contract. They retain generic ice to isolate the silicate change. Prepare one only when ready:

```bash
python -B workflow.py prepare config/grid.lab_olivine.example.json
python -B workflow.py slurm runs/lab_olivine_control_v1 \
  --machine config/machine.cluster.json
```

Preparation and script generation do not submit jobs. The normal run-local/Slurm/analyze commands work with the frozen selected material. The new `dust` configuration accepts `envelope_core_file` and `envelope_mantle_file` paths and records their exact bytes, densities and wavelength coverage. The existing disk is unchanged. Material choice is fixed per prepared run; the two examples are separate controlled experiments, not one automatically combined composition posterior.

`allow_optical_extrapolation: true` explicitly acknowledges the finite wavelength coverage. In MCFOST 4.1.13, n and k are constant below the first table wavelength and independently power-law extrapolated from their last two samples above the last wavelength. Those extended values are not laboratory measurements. This flag is not evidence that the adopted tails yield correct heating/cooling.

Zero-k custom inputs are rejected unless `allow_zero_k: true` is explicitly set. That override is an experimental acknowledgement, **not a numerical fix or a recommendation to bypass the failed H₂O check**. No joint laboratory-H₂O production config is supplied.

## Reproduce the conversion offline

From `tmc1a_restart/`, with the vendored raw tables present:

```bash
PYTHONPATH=src python -B -m mcfost_grid.laboratory_silicates --examples
python -B -m unittest discover -s tests -v
```

Conversion checks pinned download hashes and preserves the numeric rows. Rebuilding the package does not execute the dust-property checks. All original `reference/` files remain unchanged. Include this data directory, active `src/`, configs and tests when transferring through Git; no Git commit or push was performed here.
