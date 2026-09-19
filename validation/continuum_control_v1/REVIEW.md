# Continuum control review

## Downloaded results

The downloaded summary matches the local frozen manifest. All nine prediction rows report passing quality checks; their identities, observed fluxes/errors, residuals and score reproduce from the saved tables. One model is ranked, none excluded, and the analysis reports no warnings.

Regional log RMS: **0.08687054 dex**. Archived nominal control: **0.10931254 dex**. Diagnostic diagonal chi-square: **17.896088** (not reduced chi-square or a posterior likelihood).

| Wavelength (µm) | New − observed (%) | New − archived (%) |
|---:|---:|---:|
| 1.2028 | -16.01 | +2.10 |
| 1.4988 | +4.87 | -0.20 |
| 1.7589 | -1.36 | -8.19 |
| 2.1917 | -5.39 | +2.34 |
| 2.5464 | -23.95 | +11.20 |
| 3.9144 | -10.92 | +6.42 |
| 5.3906 | -28.02 | +14.46 |
| 13.7999 | +30.38 | +2.73 |
| 27.5131 | +7.91 | +1.02 |

## Scientific interpretation

The new nominal model differs from the archived nominal flux by up to 14.46%; 7/9 anchors exceed the provisional 2% restart target. The lower score alone does not establish improved physics: the physical setting is unchanged. This is a functioning continuum calculation, with archive reproduction and convergence still unresolved.

Archived images used seed 41001; the active workflow does not explicitly pass a seed. The reported simulator build has also changed, and this run recomputes temperature. These are possible contributors, not an attribution of the flux differences. Do not merge old and new scores into one ranking.

The four strongest fractional continuum residuals are at 13.80 µm (+30.38%), 5.39 µm (−28.02%), 2.55 µm (−23.95%) and 1.20 µm (−16.01%). A nine-point continuum run does not assess the H2O band or the silicate trough.

## Verification limits

Only the results directory was downloaded. Raw images, temperatures, per-model measurements, runtime identities, command records and timings are absent here. The table review cannot independently recheck image photometry, establish the actual binary/utilities hashes or benchmark runtime/memory. Do not rerun the ordinary analyzer locally over this partial download: it would replace the valid imported results with missing-model output.

Preserve and download `runtime_binding.json`, `models/*/{runtime_binding.json,status.json,measurements.json,temperature_complete.json}`, simulator `command.json`/`execution.json`/`mcfost.log` files and Slurm logs before making a reproducibility claim. FITS products are needed for independent photometry checks.

See `review.json` for input hashes and `anchor_comparison.csv` for the numerical comparison. Reproduce this review with `python -B scripts/review_continuum_control.py`.
