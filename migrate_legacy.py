#!/usr/bin/env python3
"""Preserve, verify, and (only explicitly) remove the retired local workspace.

No models, network calls, cluster submissions, or external deletions are performed.
The archive retains original paths and bytes, including hidden files and old logs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import stat
import tarfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
OLD = PROJECT / 'jwst_mcfost'
ARCHIVE = ROOT / 'archive/jwst_mcfost_legacy_2026-09-18.tar.gz'
LEDGER = ROOT / 'archive/legacy_manifest.json'
RECEIPT = ROOT / 'archive/migration_receipt.json'
REF = ROOT / 'reference'


def sha(path):
    with path.open('rb') as stream:
        return stream_sha(stream)


def stream_sha(stream):
    h = hashlib.sha256()
    for block in iter(lambda: stream.read(1024 * 1024), b''):
        h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def census(root):
    records = []
    for path in sorted(root.rglob('*')):
        st = path.lstat()
        if stat.S_ISLNK(st.st_mode) or not (stat.S_ISREG(st.st_mode) or stat.S_ISDIR(st.st_mode)):
            raise RuntimeError(f'Unsupported entry; manual preservation required: {path}')
        records.append({'path': path.relative_to(root).as_posix(),
                        'kind': 'file' if path.is_file() else 'directory',
                        'bytes': st.st_size if path.is_file() else 0,
                        'mode': stat.S_IMODE(st.st_mode),
                        'sha256': sha(path) if path.is_file() else None})
    return records


def copy_references():
    manifest_path = REF / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    entries = {r['path']: r for r in manifest['files']}
    additional = {
        'source/strict_continuum_measurement_kernel_v2.py': 'strict_continuum_measurement_kernel_v2.py',
        'source/h2o_sparse_continuum_kernel_v2.py': 'h2o_sparse_continuum_kernel_v2.py',
        'source/analyze_strict_continuum_refinement_v1.py': 'analyze_strict_continuum_refinement_v1.py',
        'source/parity_strict_continuum_measurement_kernel_v1.py': 'parity_strict_continuum_measurement_kernel_v1.py',
        'parameters/continuum_nominal.para': 'runs/strict_continuum_cavity_mass_q_v1/models/m2p25_c17p5_a0p40_q2p75_i0p05/disk_envelope_m2p25_c17p5_a0p40_q2p75_i0p05_rt.para',
        'provenance/CONTINUE_HERE_legacy.md': 'CONTINUE_HERE.md',
        'provenance/README_legacy.md': 'README.md',
    }
    for name in ['test_strict_continuum_measurement_kernel_v1.py',
                 'test_strict_continuum_measurement_kernel_v3.py',
                 'test_h2o_sparse_continuum_kernel_v2.py']:
        additional[f'source/{name}'] = name
    for run in ['strict_continuum_screen', 'strict_continuum_refinement_v1',
                'strict_continuum_3d_q_v1', 'strict_continuum_cavity_mass_q_v1',
                'h2o_3um_method2_screen_v3', 'h2o_3um_abundance_grid_v1',
                'h2o_3um_dust_axis_grid_v1', 'h2o_3um_dust_axis_band_v1',
                'h2o_silicate_demo_image_anchors_v1']:
        for path in (OLD / 'runs' / run).glob('*'):
            if path.is_file() and path.suffix in {'.json', '.tsv'}:
                additional[f'provenance/{run}/{path.name}'] = path.relative_to(OLD).as_posix()
    for path in (OLD / 'runs/strict_continuum_cavity_mass_q_v1/anchors').rglob('*.lambda'):
        additional['parameters/continuum_anchors/' + path.parent.name + '/' + path.name] = path.relative_to(OLD).as_posix()
    for path in (OLD / 'runs/h2o_silicate_demo_image_anchors_v1/inputs').glob('*.lambda'):
        additional['parameters/silicate_anchors/' + path.name] = path.relative_to(OLD).as_posix()
    pairs = [(REF / dest, OLD / source, 'jwst_mcfost/' + source) for dest, source in additional.items()]
    utils = Path.home() / 'mcfost/utils'
    for name in ['Dust/ac_opct.dat', 'Dust/ice_opct.dat', 'Stellar_Spectra/lte4000-3.5.NextGen.fits.gz']:
        pairs.append((REF / 'runtime_assets' / name, utils / name, str(utils / name)))
    for dest, source, origin in pairs:
        if not source.is_file() or source.is_symlink():
            raise RuntimeError(f'Missing reference: {source}')
        digest = sha(source)
        if dest.exists() and sha(dest) != digest:
            raise RuntimeError(f'Refusing to overwrite changed reference: {dest}')
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copy2(source, dest)
        if sha(dest) != digest:
            raise RuntimeError(f'Copy verification failed: {dest}')
        key = dest.relative_to(REF).as_posix()
        entries[key] = {'path': key, 'source': origin, 'sha256': digest, 'bytes': dest.stat().st_size}
    manifest.update(files=list(entries.values()), total_bytes=sum(r['bytes'] for r in entries.values()),
                    source_modules_are_reference_copies_not_self_contained_package=False,
                    measurement_kernel_import_closure_included=True,
                    extraction_and_simulation_workflows_require_porting=True,
                    complete_legacy_source_archive='archive/jwst_mcfost_legacy_2026-09-18.tar.gz')
    write_json(manifest_path, manifest)
    verify_references()


def verify_references():
    manifest = json.loads((REF / 'manifest.json').read_text())
    for row in manifest['files']:
        path = REF / row['path']
        if not path.is_file() or path.is_symlink() or sha(path) != row['sha256']:
            raise RuntimeError(f'Reference verification failed: {path}')
    return len(manifest['files'])


def verify_archive():
    ledger = json.loads(LEDGER.read_text())
    expected = {r['path']: r for r in ledger['entries']}
    found = set()
    with tarfile.open(ARCHIVE, 'r:gz') as archive:
        for member in archive:
            prefix = 'jwst_mcfost/'
            if member.name == 'jwst_mcfost':
                continue
            if not member.name.startswith(prefix):
                raise RuntimeError(f'Unexpected archive member: {member.name}')
            rel = member.name[len(prefix):].rstrip('/')
            if rel not in expected or rel in found:
                raise RuntimeError(f'Unexpected/duplicate archive entry: {rel}')
            row = expected[rel]
            if member.mode & 0o7777 != row['mode']:
                raise RuntimeError(f'Archive file-mode mismatch: {rel}')
            if row['kind'] == 'file':
                if not member.isfile() or member.size != row['bytes']:
                    raise RuntimeError(f'Wrong member type/size: {rel}')
                with archive.extractfile(member) as stream:
                    if stream_sha(stream) != row['sha256']:
                        raise RuntimeError(f'Archive checksum mismatch: {rel}')
            elif not member.isdir():
                raise RuntimeError(f'Wrong directory type: {rel}')
            found.add(rel)
    if found != set(expected):
        raise RuntimeError('Archive is missing entries')
    return ledger


def prepare():
    if not OLD.is_dir() or OLD.is_symlink():
        raise RuntimeError('Expected regular legacy directory not found')
    copy_references()
    if ARCHIVE.exists() or LEDGER.exists():
        raise RuntimeError('Archive already exists; use --verify, never overwrite silently')
    records = census(OLD)
    write_json(LEDGER, {'original_directory': str(OLD), 'archive_prefix': 'jwst_mcfost',
                       'entries': records, 'moved_reports': 'reports/restart_review_2026-09-18',
                       'files': sum(r['kind'] == 'file' for r in records),
                       'bytes': sum(r['bytes'] for r in records)})
    with tarfile.open(ARCHIVE, 'w:gz', compresslevel=6) as archive:
        archive.add(OLD, arcname='jwst_mcfost', recursive=True)
    verify_archive()
    if census(OLD) != records:
        raise RuntimeError('Legacy source changed during preservation; deletion is blocked')
    write_json(RECEIPT, {'prepared_utc': datetime.now(timezone.utc).isoformat(),
                        'archive': ARCHIVE.relative_to(ROOT).as_posix(), 'archive_sha256': sha(ARCHIVE),
                        'archive_bytes': ARCHIVE.stat().st_size, 'legacy_files': len([r for r in records if r['kind'] == 'file']),
                        'all_archive_member_hashes_verified': True, 'legacy_deleted': False,
                        'report_moved_to': 'reports/restart_review_2026-09-18'})
    print('Archive and reference copies verified. Legacy directory has NOT been deleted.')


def delete():
    receipt = json.loads(RECEIPT.read_text())
    if sha(ARCHIVE) != receipt['archive_sha256']:
        raise RuntimeError('Archive changed; deletion is blocked')
    ledger = verify_archive()
    verify_references()
    if OLD.resolve() != PROJECT.resolve() / 'jwst_mcfost' or OLD.is_symlink():
        raise RuntimeError('Unexpected deletion target')
    if census(OLD) != ledger['entries']:
        raise RuntimeError('Legacy tree changed since archival; deletion is blocked')
    gate = ROOT / 'validation/pre_deletion_validation.json'
    validation = json.loads(gate.read_text())
    if validation.get('passed') is not True:
        raise RuntimeError('Report/reference portability validation has not passed')
    for name, digest in validation['active_file_sha256'].items():
        if sha(ROOT / name) != digest:
            raise RuntimeError(f'Active workspace changed since validation: {name}')
    shutil.rmtree(OLD)
    receipt.update(legacy_deleted=True, deleted_utc=datetime.now(timezone.utc).isoformat(),
                   validation_receipt=gate.relative_to(ROOT).as_posix(),
                   validation_sha256=sha(gate), reference_files_verified=verify_references())
    write_json(RECEIPT, receipt)
    print('Removed only jwst_mcfost after archive, reference and portability verification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--verify', action='store_true')
    mode.add_argument('--delete-verified-legacy', action='store_true')
    args = parser.parse_args()
    if args.prepare:
        prepare()
    elif args.verify:
        ledger = verify_archive()
        print(f"Verified {ledger['files']} archived files and {verify_references()} active reference files")
    else:
        delete()
