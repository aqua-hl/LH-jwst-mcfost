# Final resolution comparison at approximately 1 AU/pixel

This is the last numerical investigation requested for the 2.546353-micron,
50-degree reference model. It reuses the saved seed-42004 temperature byte for
byte and 2,048,000 image photons. Five matched image seeds (42001–42005) each
run two new images. Only pixel count changes within each pair:

| Geometry | Field width | Pixels per side | AU per pixel |
| --- | ---: | ---: | ---: |
| coarse | 6000 AU | 4801 | 1.249740 |
| fine | 6000 AU | 6001 | 0.999833 |

The odd 6001-pixel grid preserves a source-centred image. Both geometries use
the same 112-GB MCFOST memory setting. The prior 4801-pixel products are not
reused because their 64-GB setting would add a second changed variable.
There are no new temperature solves and no photon-budget arm.

Each of five array tasks requests 64 CPUs and 160 GB, runs its two images
sequentially, and has a 24-hour allocation (four-hour limit per MCFOST call).
The peak is 320 CPUs. The 6001-pixel, eight-plane, 64-thread image array alone
is approximately 74 GB in the audited default-real build; the allocation
includes additional working space. Analysis requests two CPUs and 16 GB.

The measurement remains a source-centred 1-arcsecond aperture, the original
Gaussian PSF, model distance 140 pc and measurement distance 147 pc. The
`aperture_v2` policy records whole-image boundary losses as warnings while
retaining finite-flux and central aperture/PSF-support requirements.

## Build and run on the cluster

After committing, pushing, and pulling the code, use the shared environment:

```bash
cd ~/LH-jwst-mcfost
git pull --ff-only

~/.venvs/mcfost-v11/bin/python -B scripts/build_extinction_ice_crossover_v2.py \
  --source runs/extinction_ice_numerics_v2_git \
  --final-resolution --output runs/final_resolution_1au_v1

final_run=runs/final_resolution_1au_v1
~/.venvs/mcfost-v11/bin/python -B "$final_run/code/configure_cluster.py" \
  "$final_run" --machine "$final_run/machine.template.json"
bash "$final_run/submit.sh" "$PWD/$final_run/machine.template.cluster.json"
```

Use `runs/extinction_ice_numerics_v2` as `--source` if that is the original
cluster run's name. The source must contain the original numerical experiment
and raw saved temperature plus its provenance, not only compact summaries.
The builder validates and copies those inputs and refuses to replace an
existing output directory. Preparation itself launches neither MCFOST nor
Slurm. Set any required partition, account or reservation in the new
`machine.template.json` before configuration. The two production/diagnostic
submissions are independent; this comparison is not a production dependency.

## Record the outcome and stop

The declared 1% reference rule uses total-I aperture flux. Both geometries'
95% Gaussian scatter upper bounds must be below 1%, and the complete 95%
paired-change interval must lie inside -1% to +1%. Intervals must be
nondegenerate and all five pairs must run on matching hosts. Five seeds do
not establish distributional assumptions or exclude rare outliers; the rule
is a local numerical reference, not a certification of all production models.

The dependent analysis remeasures all images and writes
`results/summary.json`, `seed_fluxes.csv`, `scatter.csv`, `spatial_changes.csv`,
`crossover_diagnostics.png/.pdf`, and `REVIEW.md`.
`final_reference_result.outcome` records `met`, `not_met`, or
`not_evaluated` for an incomplete/invalid run. A complete valid comparison
sets `investigation_closed: true` regardless of whether the 1% reference is
met. No further experiment is scheduled and the outcome does not veto
production. `production_convergence_certified` remains false. Execution or
integrity errors remain explicit instead of becoming a scientific result.

Download `results/`, `logs/`, `experiment.json`, `experiment.sha256`, and
`runtime_binding.json`; retain the complete run on the cluster for raw audits.
Read file status with:

```bash
~/.venvs/mcfost-v11/bin/python -B "$final_run/code/crossover_task.py" \
  "$final_run" --status
```

Use `squeue -r -j JOBID` and
`sacct -j JOBID -X --format=JobID%22,State%20,NodeList%20,Elapsed,ExitCode`
for live scheduler/accounting information. The diagnostic identifier is
`fixed_temperature_final_resolution_1au_v1`; the tools continue supporting
the separate, earlier 2401/4801-pixel diagnostic without changing its design.
