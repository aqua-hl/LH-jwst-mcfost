# Portable MCFOST grid workflow

This is a new, small workflow, separate from the immutable reference code and the legacy archive. One JSON file describes the science grid; a second describes the computer. The same prepared model catalogue runs locally or as a Slurm job array. No Python command submits jobs, queries a cluster, or changes the historical report.

## Quick local start

Run from `tmc1a_restart` with NumPy, SciPy, Astropy and Matplotlib installed. The laptop example uses Python 3.12 in the existing `/opt/anaconda3/bin/python` environment. For the new cluster, use a Python 3.12 environment with `pip install -r requirements-workflow.txt`. The executable and complete MCFOST utilities are separate external dependencies.

```bash
python -B workflow.py run-local runs/laptop_smoke_v4 --machine config/machine.local.json --workers 1
python -B workflow.py status runs/laptop_smoke_v4
python -B workflow.py analyze runs/laptop_smoke_v4 --workers 2
```

The supplied smoke run already completed: these commands reuse its results. To test different values, edit `config/grid.smoke.json`, choose a **new** `run_name`, then run `python -B workflow.py prepare config/grid.smoke.json` and use the newly printed run path. Preparation deliberately refuses to overwrite an existing experiment. The larger `continuum_grid_example` has also been prepared for script inspection only; it has not been simulated or submitted.

The smoke example contains two inclinations and two wavelengths, 2,000 photon packets, a reduced density/grain grid, and 151 × 151 images. It uses one model at a time, two OpenMP threads, and a five-minute limit per simulator command. It is a software demonstration, **not a converged physical fit**, and plots carry that label. The 128,000-packet full-resolution example is not automatically run.

## Choose parameter values

Edit `config/grid.example.json`. Lists under `grid` are combined as a Cartesian product:

```json
"grid": {
  "inclination_deg": [50, 60, 70],
  "envelope_dust_mass_msun": [0.00015, 0.000225],
  "cavity_half_opening_deg": [17.5, 25]
}
```

This creates 3 × 2 × 2 = 12 independent physical models. Values under `fixed` remain fixed unless overridden by a grid axis. Every model's complete parameter dictionary is saved in `manifest.json`; filenames are short IDs, not the definition of the physics.

For a selected list rather than every combination, replace `grid` with:

```json
"models": [
  {"inclination_deg": 60, "envelope_dust_mass_msun": 0.00015},
  {"inclination_deg": 70, "cavity_half_opening_deg": 20}
]
```

Do not include both `grid` and `models`. Unknown names, invalid values and duplicate physical models are rejected before simulation.

| Parameter | Definition |
|---|---|
| `envelope_dust_mass_msun` | Envelope **dust** mass, in solar masses, supplied to the density model; not total gas mass or the old dimensionless multiplier |
| `envelope_amax_um` | Maximum envelope grain radius, micrometres; disk grains unchanged |
| `envelope_size_exponent` | Positive q in dn/da ∝ a⁻q; not the envelope density exponent |
| `cavity_half_opening_deg` | Conical half-opening from the polar axis; full opening is twice this value |
| `inclination_deg` | Viewing angle from the pole, 0–90° |
| `envelope_ice_volume_fraction` | Mantle **volume** fraction of the coated Mie envelope grains, strictly between zero and one |
| `distance_pc` | Simulation distance; photometry converts both flux and angular scale to `observations.target_distance_pc` |
| `stellar_temperature_k`, `stellar_radius_rsun`, `stellar_mass_msun` | Stellar photosphere; MCFOST selects an atmosphere from its utilities |
| `accretion_rate_msun_per_year` | Stellar `-Mdot 1` prescription, not an independently fitted bolometric luminosity |

The first implementation supports the preserved two-zone **coated Mie** template: 70/30 silicate/carbon disk species and silicate-core/generic-ice-mantle envelope. It does not silently reinterpret the later disjoint DHS H₂O species mass fractions as mantle volume fractions. The latter dust family and its 0.35″ H₂O continuum/LSF operator need a dedicated extension before that fit is restarted. Arbitrary geometry templates, dust families and MCMC are outside this first implementation.

For heating experiments, varying stellar radius changes the photospheric luminosity at fixed temperature. Changing `-Mdot` changes another source component; neither is simply the observed Lbol. Inspect the MCFOST log for the resulting luminosities. A changed atmosphere can require additional files from the full utility installation.

## Choose observations and numerical resolution

`observations.preset = "strict9"` preserves the nine continuum anchors, excludes the old ≈0.995 μm anchor, and uses the 1″-radius observation. Detector overlaps are combined into one anchor using the original geometric mean and common calibration floor, rather than treated as independent points. The full grey JWST R100 spectrum remains in every spectral plot, including masked bands as context.

Omit `anchor_ids` for all nine. An explicit subset is useful for a smoke run. Six-region-balanced log RMS is the default screening metric; no flux normalization is fitted. With a subset, only its represented regions contribute, so its score is not comparable to a full-nine score. The diagnostic chi-square uses diagonal supplied errors and is neither a reduced chi-square nor a posterior likelihood.

For new anchors, provide `observations.anchors_file` as CSV or ECSV with columns:

```text
id,wavelength_um,instrument,region,flux_jy,uncertainty_jy,score
```

`instrument` is NIRSpec or MIRI; `score` is true/false. Optional fields are `image_npix`, `image_size_au`, `psf_fwhm_arcsec`. Diagnostic anchors with `score=false` are plotted but not used in the objective. Supply observations extracted in the configured aperture, and explicitly set `observations.spectrum_file` to the matching ECSV spectrum (`wavelength_um`, `flux_jy` columns) for the grey background. Custom anchors without this field are rejected. A different aperture is not obtained by simply changing the radius against the old data. No interpolation across detector boundaries or band-continuum extraction is performed automatically.

The production example uses 128,000 packets and the previous wavelength-dependent image fields/pixel counts. Explicit `numerics.image_npix` and `image_size_au` override those presets. Photon counts, density-grid dimensions and grain bins can also be set. Changing numerical resolution is not evidence of convergence by itself.

## What one task does

1. Check frozen inputs and the simulator's capabilities; acquire an exclusive per-model lock.
2. Compute a fresh radiative-equilibrium temperature for that physical model. A changed mass, grain population or source never borrows another model's temperature.
3. Ray-trace each requested wavelength with full polarized component output, using that temperature.
4. Immediately measure signed total-I plane 0 in the distance-scaled aperture after the inherited Gaussian JWST PSF approximation.
5. Save a compact measurement JSON, timing and provenance. Then proceed to the next wavelength.

Each attempt has its own log/output folder. A failed attempt is retained; a rerun skips successful temperatures and already measured wavelengths. No automatic retry changes photons, geometry or physical parameters. Never edit a frozen run's inputs; prepare a new run instead. The active entry point dispatches run/analysis commands to the code snapshot stored inside the prepared run.

The default `image_method2` backend uses standard MCFOST `-img ... -rt2`. Images use the Cartesian adaptive spatial sampler; `-rt2` separately chooses the source-function method. **No method-1 SED is generated or used to normalize aperture flux.** Optional `coeval_method2` requires the custom `-rt-sed-method 2` build and checks its same-realization SED/image plane-0 closure. Missing support fails explicitly; it never falls back to method 1. The laptop binary does not provide that custom extension.

## Why analysis should be quicker

The expensive operation is radiative transfer, not ranking a few hundred numerical fluxes. Image measurement now happens once inside the model task. For a symmetric Gaussian convolution G, aperture sum ⟨G I,A⟩ equals ⟨I,G A⟩: the code caches the smoothed aperture weights and evaluates a dot product against each signed image. Only boundary pixels need fractional-circle subsampling.

The final analysis reads small per-model JSON files, not hundreds of multi-component compressed FITS images. `--workers` parallelizes compact-file reads; the small score and plotting calculations are serial. Python BLAS pools are limited to one thread to prevent oversubscription. The image reader still reads the eight-plane FITS product once to validate it; this is not a promise that all I/O or simulation cost disappears. See the recorded synthetic benchmark for measured kernel timings.

Results are written under each run's `results/`:

- `ranking.ecsv`: parameters, score, diagnostic chi-square and exclusion reason;
- `predictions.ecsv`: observed/model fluxes and residuals at actual evaluated wavelengths;
- `summary.json`: best model, completion counts and limitations;
- `spectrum_comparison.png` / `.pdf`: grey JWST spectrum, anchors, best-model points and residuals;
- `parameter_space.png` / `.pdf`: tested coordinates and profile minima, not a posterior.

Only complete, valid models with the same scored-anchor contract enter the ranking. Missing models remain visible as excluded. Sparse model points are not presented as a simulated dense spectrum.

## Future cluster transfer

Nothing has been submitted. Edit `config/machine.cluster.example.json` for the new cluster: executable, complete utilities, Python environment, partition, optional account/reservation and module setup. Paths in the science JSON are relative to that JSON. Runtime machine paths can be absolute or use environment variables. No ENS-Lyon node names or reservation are embedded.

After transferring the workspace (or at least the prepared run plus a machine JSON), generate scripts:

```bash
python -B workflow.py prepare config/grid.example.json
python -B workflow.py slurm runs/continuum_grid_example --machine config/machine.cluster.example.json
```

Generation does not require `sbatch` or a local copy of the cluster binary. Review `job_array.sh`, `job_analysis.sh`, and `submit.sh` in the run directory. The **future** submission command is:

```bash
bash runs/continuum_grid_example/submit.sh /absolute/cluster/path/to/machine.cluster.json
```

If you transferred only a prepared run, call its saved entry point instead: `python RUN/code/workflow.py slurm RUN --machine /path/to/machine.cluster.json`. The prepared run already contains the code and extracted observations it needs; it still requires MCFOST and complete utilities on the target computer.

One array element is one physical model: one process with `cpus-per-task` OpenMP threads. The example is 16 threads × at most four simultaneous elements = at most 64 CPUs. It does not use `mpirun`, demand one node per model, or reserve exclusive nodes. Slurm can place several elements on a sufficiently large node. `max_parallel` controls the array throttle; edit/regenerate it deliberately rather than automatically consuming idle resources. Never request threads exceeding the actual allocation.

Absolute run and machine paths are exported at submission. The batch job does not derive the project location from Slurm's spool copy of `$0`, avoiding the earlier `/var/lib/slurm` failure. `submit.sh` runs `sbatch` from the run directory, where `logs/` already exists. The analysis job depends on `afterany`, so it reports failures as well as successes; incomplete models cannot become winners. Slurm's array identity and throttle follow its [official job-array documentation](https://slurm.schedmd.com/job_array.html).

## Checks retained and checks deferred

Retained inexpensive checks: parameter names/ranges, frozen input hashes, executable/backend identity, process exit and completion markers, expected output files, valid temperature, FITS units/shape/wavelength, distance and aperture geometry, finite positive signed aperture flux, field edges and PSF loss, and optional coeval closure. Locks prevent duplicate workers; incomplete outputs are never ranked.

Not run automatically: the old multi-seed campaigns, 512k-photon tests, repeated resolution grids, extensive per-component audit trees, complete data re-extraction, or repeated full-FITS analysis. Unit tests exercise code in seconds, not hours. The lightweight checks do **not** certify Monte Carlo convergence, new dust physics, the Gaussian PSF approximation, or the unresolved observation-centre offset. Before interpreting a changed model family or new cluster build, a small targeted control remains appropriate; passing a laptop smoke test is not permission to skip all numerical validation.

MCFOST's [overview](https://mcfost.readthedocs.io/en/latest/overview.html) describes the Monte Carlo temperature/source-function calculation followed by ray tracing. Raw cubes are not read by this workflow; it uses the preserved extracted observations. The historical report, archived simulations and reference sources are unchanged.

## Developer checks

```bash
python -B -m unittest discover -s tests -v
```

The new tests cover parameter rendering, aperture parity/units/distances, grid expansion, incomplete-result exclusion and script generation. They are separate from the preserved historical tests under `reference/source`.
