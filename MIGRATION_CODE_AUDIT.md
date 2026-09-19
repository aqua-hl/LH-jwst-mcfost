# Local migration and dependency audit

Date: 18 September 2026. Scope: the local `jwst_mcfost` tree, the new `tmc1a_restart` reference bundle, and directly referenced local resources. No cluster access or scientific simulation was performed.

## Inventory and preservation requirements

After relocation of the report and the earlier cache cleanup, the audited legacy tree contained **3,451 regular files and 1,094 directories**, with no symbolic links or special files. Its largest components were `runs/` (983,202,414 logical bytes), `output/` (87,409,158 bytes), and root-level files (8,981,432 bytes). These counts are the audit snapshot, not a promise that subsequent migration metadata has identical counts.

Hidden files include scientific `.sed_th.fits.gz` products and a partial image log, not just caches. The preservation archive must include hidden files and empty directories, retain relative paths and executable modes, and verify every regular file against a SHA-256 manifest before the legacy directory is removed. A complete archive is preferable to guessing which historical scripts or failed-run records can be discarded. Frozen manifests, selected convergence images, model inputs, temperature products, compact predictions, and recovery evidence are all retained this way.

The archive preserves the **available local evidence**, not all products formerly present on the retired cluster. Most production image arrays were not copied locally. The archive also does not supply a portable MCFOST executable or the full utility distribution.

## Code review

All **420 Python files** in the legacy tree, including **387 root-level modules**, parsed successfully with Python's abstract-syntax-tree parser. This is a syntax check, not a claim that every script executes successfully, that scientific algorithms are correct, or that the legacy workflow is portable.

The reusable measurement code has these required local dependencies:

| Reference module | Required sibling modules |
|---|---|
| `strict_continuum_measurement_kernel_v3.py` | `strict_continuum_measurement_kernel_v1.py`, `strict_continuum_measurement_kernel_v2.py` |
| `strict_continuum_measurement_kernel_v2.py` | `strict_continuum_measurement_kernel_v1.py` |
| `h2o_observation_measurement_kernel_v3.py` | `h2o_sparse_continuum_kernel_v2.py` |

The initial 37-file reference bundle omitted the second continuum kernel and the sparse H2O continuum kernel. These are required additions for an importable measurement-code reference set. Relevant standalone regression tests are `test_strict_continuum_measurement_kernel_v1.py`, `test_strict_continuum_measurement_kernel_v3.py`, and `test_h2o_sparse_continuum_kernel_v2.py`. Tests can be preserved independently of old cluster submission scripts.

Scientific Python dependencies include NumPy, SciPy, Astropy, Matplotlib, pandas, and photutils. Rebuilding the report also requires its PDF tooling. There is no verified, complete environment lock for every historical script; a new-cluster environment should be recorded and tested explicitly.

Many legacy modules derive `runs/` and `output/` from `Path(__file__)`; cluster runners additionally pin `/home/uhlo01/mcfost-custom/.../bin/mcfost`, `/home/uhlo01/mcfost/utils`, executable hashes, and workflow-source hashes. Their literal paths are provenance, not valid defaults for a new cluster. They should remain unchanged inside the archive. A new production implementation requires explicit configuration rather than editing historical manifests in place.

`extract_tmc1a_sed.py` expects the original project layout, with `cubes/` one directory above its script directory. A byte-identical copy under `reference/source/` is therefore **reference code, not a runnable relocated extractor**. Similarly, the H2O observation freezer depends on `DATA/reproducible_ice_extraction/`, the NIRSpec cubes, and sibling audit/measurement modules. The corrected v3 observation tables are already preserved as reference data.

Moving `build_review.py` and `inventory_workspace.py` changes the meaning of their original `HERE.parents[1]` root calculation. Their relocated versions must use explicitly preserved inputs or configurable roots. The historical inventory should not be overwritten by accidentally inventorying the new workspace.

## External resources

The following remain outside the directory being deleted:

| Resource | Role and migration requirement |
|---|---|
| `cubes/` (approximately 716 MiB allocated) | Original NIRSpec and MIRI inputs; needed for independent re-extraction, not for reading the saved report or predictions. |
| `DATA/reproducible_ice_extraction/` (approximately 33 MiB) | Supplied corrected ice-extraction code, native spectra, coordinates, and provenance; preserve with the larger project. |
| `DATA/mcfost_ice_inputs/` (approximately 11 MiB) | Laboratory optical constants, initial constraints, and opacity-only experiments; preserve without treating the older constraints as the adopted v3 observation contract. |
| `DATA/JWST_ver196_center.mrt` | Separate central-spectrum comparison; not interchangeable with the adopted aperture extraction. |
| `mcfost-tuto/` | Historical tutorial parameter file and baseline products referenced by early comparison scripts. |
| `/Users/harrisonlo/mcfost/utils/` (approximately 239 MiB) | Full local MCFOST utility distribution; install/transfer a verified compatible distribution before a new run. |

The following locally available files have hashes exactly matching the frozen production runtime contract in `h2o_3um_method2_screen_v3_common.py`:

| Path relative to local `mcfost/utils/` | SHA-256 |
|---|---|
| `Dust/Draine_Si_sUV.dat` | `3a6af751429f6f8dc948e60675b13be01ffcbed0662eec62934ec6d9c61dea80` |
| `Dust/ac_opct.dat` | `ef8f09849580a0859f29f61216486e7e249103759f851bcfa7abe11c09a44de2` |
| `Dust/ice_opct.dat` | `44625b4918b7015e97d3c92b3f19e9e06eb57d209f3b1a946ee8f141215c8697` |
| `Stellar_Spectra/lte4000-3.5.NextGen.fits.gz` | `a626d82c52e07dec12b7fe141aeb2c9af974516ffe39178a4c9a2094c2a6d87d` |

Draine silicate is already in the selected bundle. The other three small files should be copied and hashed with it so that the historical carbon, generic ice, and stellar spectrum can be transferred as explicit references. These selected assets are not a replacement for the full MCFOST utility distribution. The H2O optical constants and two custom MCFOST source patches are already present in the reference bundle.

## Cross-project references and deletion criterion

A recursive text scan of other project Python, shell, Markdown, JSON, notebook, TOML, YAML, text, and workspace-configuration files (including hidden configuration directories, excluding the two audited workspaces, bytecode, notebook checkpoints, and Git internals) found **no literal `jwst_mcfost` references**. No changes to unrelated project scripts are indicated by that scan. This does not cover external applications, shell history, arbitrary binary files, or dynamically constructed paths.

Removal of the old directory is justified only after: (1) its full preservation archive and per-file manifest have been independently verified; (2) report generation succeeds without consulting the old directory; (3) the selected reference bundle and required kernel dependencies pass integrity checks; and (4) documentation explains external data/runtime requirements and archive restoration. This preserves the old code and data without claiming the new workspace is already an executable production fit.
