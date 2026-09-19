#!/usr/bin/env python3
"""Verify reference copies without accessing the retired folder or a cluster."""
import json
from pathlib import Path
from migrate_legacy import verify_references

if __name__ == '__main__':
    count = verify_references()
    manifest = json.loads((Path(__file__).resolve().parent/'reference/manifest.json').read_text())
    print(f'Verified {count} reference files ({manifest["total_bytes"]/2**20:.2f} MiB).')
    print('No legacy directory, cluster access or model execution is required.')
