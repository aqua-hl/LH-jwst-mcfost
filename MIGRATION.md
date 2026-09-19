# Local consolidation and preservation

The requested scope is the former `jwst_mcfost` directory only. Parent-project cubes, DATA, tutorials, notebooks and unrelated modelling projects are not deleted or modified.

## Preserved material

The complete report was moved to `reports/restart_review_2026-09-18/`. Every remaining regular file and directory in the legacy tree was inventoried, archived with its original relative path, and checked against a SHA-256 ledger. No distinction between apparently obsolete and useful scientific files was used for deletion: all locally available legacy material is recoverable from the archive.

- [Legacy archive](archive/jwst_mcfost_legacy_2026-09-18.tar.gz)
- [Per-file and directory ledger](archive/legacy_manifest.json)
- [Operation receipt and archive checksum](archive/migration_receipt.json)
- [Reusable reference manifest](reference/manifest.json)
- [Code and dependency review](MIGRATION_CODE_AUDIT.md)
- [Portability and test receipt](validation/pre_deletion_validation.json)

The compact active references were expanded to include missing kernel imports and test dependencies, the exact nominal continuum parameter file, anchor wavelengths, task/model manifests, and locally available production-matching dust and stellar-spectrum files. Copies retain the original scientific bytes. Full legacy code and selected convergence products are retained in the archive.

The original cache-cleanup and inventory helper scripts are saved under `archive/report_helpers_pre_migration/` for provenance only. Their relative-path assumptions describe the previous layout; they are not active maintenance commands. The frozen report audits describe the pre-migration tree. Rerunning an inventory against the new tree would answer a different question.

## Deletion safeguards

`migrate_legacy.py --delete-verified-legacy` requires an exact match between the archive and its ledger, rechecks every file in the source tree, verifies the active reference bundle, and requires a successful portability-validation receipt. It then removes only the adjacent directory named `jwst_mcfost`. A changed source or checksum blocks deletion. The receipt records whether deletion actually occurred.

## Restoration

From this workspace, restore into a new, empty directory rather than overwriting active work:

```bash
python -B migrate_legacy.py --verify
mkdir legacy_restore
tar -xzf archive/jwst_mcfost_legacy_2026-09-18.tar.gz -C legacy_restore
```

The historical tree appears as `legacy_restore/jwst_mcfost/`. Its files retain original bytes and within-tree paths. Absolute old workstation/cluster paths in scripts and manifests remain historical and require deliberate porting; restoring files does not make old jobs runnable. The moved, revised report remains in the active `reports/` directory rather than in the legacy archive.

## Completeness boundary

“Complete” here means complete preservation of the locally available old folder and successful reproduction of the current report. It does not mean that all former cluster outputs were downloaded, or that a new simulation stack has been installed. Missing production images/temperatures are recorded in the scientific report. Raw cubes, external ice inputs and the complete MCFOST utilities remain separate dependencies. The supplied kernel tests check measurement behaviour; they do not establish convergence for a new physical model.
