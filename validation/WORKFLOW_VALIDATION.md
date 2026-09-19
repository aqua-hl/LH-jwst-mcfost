# New workflow validation — 2026-09-18

Scope: local software integration only. No cluster connections, job submissions, full-resolution parameter grids, or renewed convergence campaigns were performed.

## Real MCFOST execution

- Run: `runs/laptop_smoke_v4`, with its complete frozen manifest/code/parameters.
- MCFOST 4.1.13, reported source SHA `fc14f46491616224d993503b2a3293c1dab7b47a`; the binary and full utility-tree fingerprints are in the run-level runtime binding.
- Two independent models, inclinations 60° and 70°; two anchors at 2.1916659199557045 and 13.799908320389276 μm.
- Fresh temperatures for both, followed by four stock Method-2 polarized images.
- 2,000 packets, 40 × 24 density grid, eight grain bins, 151 × 151 images over 6,000 au. One worker, two OpenMP threads.
- Both models completed and all four signed-aperture measurements passed the lightweight unit/geometry/flux checks.
- Total local-run elapsed time: **29.578 s**. Image calls took 5.672–8.519 s each; aperture measurement including small FITS reads took 0.0043–0.0052 s each.
- Analysis and PNG/PDF plotting via the frozen command-line entry point took **1.357 s** wall time. Both figures were visually inspected; grey observed spectrum, black anchors, green actual model samples, residuals and smoke labels are present.
- Repeating `run-local` returned both models as cached in **0.763 s**; attempt counts remained two temperature attempts and four image attempts. A direct array-style `task --index 0` also returned cached.

The coarse two-point score is not a new best fit and is not comparable to the previous nine-anchor scientific ranking. This test does not establish physical or Monte Carlo convergence.

## Automated checks

Final active source: **60 tests passed in 2.147 s** (warm plotting cache):

```bash
MPLCONFIGDIR="$PWD/.mplconfig" python -B -m unittest discover -s tests -v
```

Coverage includes parameter rendering and validation, exact original strict-nine observed flux/error agreement, fractional-aperture parity, Gaussian adjoint/direct parity, flux and angular distance scaling, malformed FITS rejection, optional coeval closure, missing/invalid model exclusion, stale/mixed provenance, changed frozen inputs/products/utilities, failed-task resume, allocation limits, and Slurm generation. The patched coeval backend was exercised with fixtures/mocks, **not with a real patched executable**, which is absent locally.

The successful smoke snapshot predates additional preparation checks: equivalent integer/float parameter values are now deduplicated, the model-count limit must be a positive integer, and custom anchors now require an explicit matching grey-spectrum file. Those refinements are covered by the final tests and the prepared cluster example. Simulation/photometry/analysis code in the successful smoke is otherwise the final implementation.

## Cluster preparation, not execution

`runs/continuum_grid_example` contains 12 parameter rows × nine strict anchors. Its generated scripts request one process × 16 CPUs and `--array=0-11%4` (maximum 64 concurrent CPUs), followed by an `afterany` analysis job. All three shell scripts passed `bash -n`. There are no simulation status records for these 12 models. Cluster-specific paths/partition remain explicit placeholders.

The full-resolution 128k-packet workflow has **not** been performance-benchmarked on the new cluster. Runtime and memory there remain to be measured. The real laptop test confirms stock-image execution, not real Slurm scheduling or the custom coeval extension.

## Preservation

All 109 original reference files still pass their existing SHA-256 verification. Historical report and archive bytes were not edited. The three earlier failed local development attempts and their fixes are documented in `workflow_development/README.md`, separate from current runs. Raw cubes and external DATA inputs remain untouched.

The isolated production-size aperture benchmark is in `docs/photometry_benchmark.json`; its kernel-only speedups must not be interpreted as whole-MCFOST speedups.
