#!/usr/bin/env python3
"""Validate and rank the 3x3 strict-continuum method-2 refinement v1.

The expensive refinement contains eight fresh models (80 one-wavelength RT
tasks) and reuses the already validated ``m2p00_a0p50_i0p05`` corner.  Fresh
predictions are measured directly from polarized total-I plane 0 after the
same Gaussian PSF and deterministic nominal 1-arcsec-radius fractional-pixel
aperture used by the original strict screen.  Zero-clipped flux is a uniform
quality diagnostic only; the signed aperture flux remains canonical.

This program is deliberately a provisional screen, not a formal likelihood
or a production parameter estimate.  Real analysis requires the externally
supplied SHA-256 of a separately frozen, rank-blind analysis manifest.  That
manifest pins this analyzer, the refinement input manifest, and every reused
source artifact before any refinement output or rank is read.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.ndimage import gaussian_filter


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
RUN_ROOT = HERE / "runs/strict_continuum_refinement_v1"
OUTPUT_ROOT = HERE / "output"
MANIFEST = RUN_ROOT / "input_manifest.json"
ANALYSIS_MANIFEST = RUN_ROOT / "analysis_manifest.json"
MODELS_TABLE = RUN_ROOT / "models.tsv"
RT_TASKS_TABLE = RUN_ROOT / "rt_tasks.tsv"
ANALYZER = Path(__file__).resolve()

INPUT_MANIFEST_SHA256 = (
    "580b43d6d3e154ec5e6db508172fb753bb33c2b78ee5218ea3b6c997a4bdde10"
)

SOURCE_PREDICTIONS = OUTPUT_ROOT / "strict_continuum_method2_fit_v1_predictions.ecsv"
SOURCE_PREDICTIONS_SHA256 = (
    "f6a7ce537a9ab4adbfdae55dc868dc5fb36017266897aff6bdc94f9bb89d0568"
)
SOURCE_FIT_VALIDATION = OUTPUT_ROOT / "strict_continuum_method2_fit_v1_validation.json"
SOURCE_FIT_VALIDATION_SHA256 = (
    "90aabddb0342783941fb5c34e2d9a32b0bf1deaabd89709c80eb36cc9edf0814"
)
SOURCE_V3_VALIDATION = OUTPUT_ROOT / "strict_continuum_screen_validation_polarized_plane0_v3.json"
SOURCE_V3_VALIDATION_SHA256 = (
    "efa02d1c64ff462a9d35f7f39e389d51ad2ef557ff790c80203b2040841375ff"
)
SOURCE_MASK = OUTPUT_ROOT / "continuum_fit_mask_v1_major_ice_both_silicates_retained.ecsv"
SOURCE_MASK_SHA256 = (
    "a8d5c28ddeeec6727e5c68e800ddba920d7b475170066fa6a5493d4a01c3375b"
)
SOURCE_FIT_POLICY = (
    HERE / "runs/strict_continuum_screen/amendments/fit_policy_v1/policy.json"
)
SOURCE_FIT_POLICY_SHA256 = (
    "efb248506438f316756a29ead33b1d963968a32702d08f6eff5105cb45784f8e"
)
SOURCE_FIT_WORKFLOW = (
    HERE
    / "runs/strict_continuum_screen/amendments/fit_policy_v1/workflow_manifest.json"
)
SOURCE_FIT_WORKFLOW_SHA256 = (
    "06c4afb4d7ea1c4557aaf3ee9caef400a0f78de5971c925a0d2319eca6bf6946"
)

OUTPUT_PREFIX = "strict_continuum_refinement_v1"
VALIDATION_OUTPUT = OUTPUT_ROOT / f"{OUTPUT_PREFIX}_validation.json"
PREDICTIONS_OUTPUT = OUTPUT_ROOT / f"{OUTPUT_PREFIX}_predictions.ecsv"
SCORES_OUTPUT = OUTPUT_ROOT / f"{OUTPUT_PREFIX}_scores.ecsv"
REGIONS_OUTPUT = OUTPUT_ROOT / f"{OUTPUT_PREFIX}_regions.ecsv"
SUMMARY_OUTPUT = OUTPUT_ROOT / f"{OUTPUT_PREFIX}_summary.json"
PLOT_OUTPUT = OUTPUT_ROOT / f"{OUTPUT_PREFIX}.png"
DECLARED_OUTPUTS = (
    VALIDATION_OUTPUT,
    PREDICTIONS_OUTPUT,
    SCORES_OUTPUT,
    REGIONS_OUTPUT,
    SUMMARY_OUTPUT,
    PLOT_OUTPUT,
)

GRID = (
    (0, "m2p00_a0p50_i0p05", 2.00, 0.50, "upstream_reuse"),
    (1, "m2p00_a0p60_i0p05", 2.00, 0.60, "fresh"),
    (2, "m2p00_a0p70_i0p05", 2.00, 0.70, "fresh"),
    (3, "m2p25_a0p50_i0p05", 2.25, 0.50, "fresh"),
    (4, "m2p25_a0p60_i0p05", 2.25, 0.60, "fresh"),
    (5, "m2p25_a0p70_i0p05", 2.25, 0.70, "fresh"),
    (6, "m2p50_a0p50_i0p05", 2.50, 0.50, "fresh"),
    (7, "m2p50_a0p60_i0p05", 2.50, 0.60, "fresh"),
    (8, "m2p50_a0p70_i0p05", 2.50, 0.70, "fresh"),
)

# (bin, exact wavelength, region, spatial configuration, boundary compromise)
ANCHORS = (
    (2, "0.9946761707696127", "near_ir_short", "full_n2401", False),
    (21, "1.2028117593024095", "near_ir_short", "full_n2401", False),
    (43, "1.4987957445436022", "near_ir_short", "full_n2401", False),
    (59, "1.7588530996181813", "near_ir_rise", "full_n2401", False),
    (81, "2.1916659199557045", "near_ir_rise", "full_n2401", False),
    (96, "2.546352514425122", "near_ir_rise", "full_n2401", True),
    (139, "3.9143995604061015", "window_3p9", "zoom_n1201", True),
    (171, "5.3906283153397085", "window_5p4", "zoom_n1201", False),
    (265, "13.799908320389276", "window_13p8", "zoom_n1201", False),
    (334, "27.513091575690556", "window_27p5", "zoom_n1201", False),
)
REGIONS = (
    "near_ir_short",
    "near_ir_rise",
    "window_3p9",
    "window_5p4",
    "window_13p8",
    "window_27p5",
)
REGION_COUNTS = {
    "near_ir_short": 3,
    "near_ir_rise": 3,
    "window_3p9": 1,
    "window_5p4": 1,
    "window_13p8": 1,
    "window_27p5": 1,
}
SCENARIOS = (
    ("primary", 1.0, "all ten anchors at primary weight"),
    (
        "boundary_half",
        0.5,
        "bins 96 and 139 jointly at half their primary contribution",
    ),
    (
        "boundary_drop",
        0.0,
        "bins 96 and 139 jointly excluded from the score numerator",
    ),
)
BOUNDARY_BINS = frozenset((96, 139))

MODEL_DISTANCE_PC = 140.0
TARGET_DISTANCE_PC = 147.0
DISTANCE_FACTOR = (MODEL_DISTANCE_PC / TARGET_DISTANCE_PC) ** 2
APERTURE_RADIUS_ARCSEC = 1.0
APERTURE_SUBPIXELS = 64
PSF_TRUNCATION_SIGMA = 6.0
COEVAL_TOLERANCE = 1.0e-5
NEGATIVE_FRACTION_MAX = 1.0e-3
EDGE_FRACTION_MAX = 1.0e-3
CONVOLUTION_LOSS_MAX = 1.0e-4
APERTURE_CLIP_MAX = 1.0e-3
COMPONENT_SIGNIFICANT_FRACTION = 1.0e-8
COMPONENT_ACTIVE_FRACTION = 0.01
COMPONENT_APERTURE_FRACTION_MIN = -1.0e-6
COMPONENT_APERTURE_FRACTION_MAX = 1.0001
OPERATIONAL_SCORE_MARGIN_DEX = 0.010857927588208284
EXPECTED_SEED = 41001
EXPECTED_PACKETS = 128000
EXPECTED_THREADS = 16

COMPONENTS = {
    4: "direct_star",
    5: "scattered_star",
    6: "direct_thermal",
    7: "scattered_thermal",
}
CLOSURE_LABELS = {
    0: "total_I",
    4: "direct_star",
    5: "scattered_star",
    6: "thermal_dust",
    7: "scattered_thermal",
}

REQUIRED_RT_FIELDS = {
    "task_index",
    "fresh_model_index",
    "grid_model_index",
    "model_label",
    "mass_factor",
    "envelope_amax_um",
    "anchor_ordinal",
    "anchor_tag",
    "region",
    "bin_index",
    "wavelength_um",
    "wavelength_token",
    "instrument",
    "boundary_compromise",
    "spatial_config",
    "zoom_token",
    "nx",
    "ny",
    "expected_fov_au",
    "expected_pixel_au",
    "expected_cdelt_deg",
    "psf_fwhm_arcsec",
    "seed",
    "sed_photons",
    "parameter",
    "parameter_sha256",
    "wavelength_file",
    "wavelength_file_sha256",
    "temperature_task_index",
    "temperature_workdir",
    "task_workdir",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_regular(path: Path, label: str) -> str:
    if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
        raise ValueError(f"Missing regular non-symlink {label}: {path}")
    return sha256(path)


def require_hash(path: Path, expected: str, label: str) -> str:
    actual = require_regular(path, label)
    if actual != expected:
        raise ValueError(f"{label} SHA-256 changed: {actual}; expected {expected}")
    return actual


def resolve(value: object) -> Path:
    path = Path(str(value))
    return path.resolve() if path.is_absolute() else (PROJECT_ROOT / path).resolve()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))


def read_json(path: Path, label: str) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{label} root is not an object")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [dict(row) for row in csv.DictReader(stream, delimiter="\t")]
    if not rows:
        raise ValueError(f"Empty TSV: {path}")
    return rows


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def atomic_table(table: Table, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".ecsv", dir=path.parent
    )
    os.close(descriptor)
    try:
        table.write(temporary, format="ascii.ecsv", overwrite=True)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def atomic_figure(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".png", dir=path.parent
    )
    os.close(descriptor)
    try:
        figure.savefig(temporary, dpi=220, bbox_inches="tight")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def psf_fwhm_arcsec(instrument: str, wavelength_um: float) -> float:
    if instrument == "NIRSpec":
        return 0.033 * wavelength_um
    if instrument == "MIRI":
        return 0.033 * wavelength_um + 0.106
    raise ValueError(f"Unsupported instrument: {instrument}")


def symmetric_fraction(first: float, second: float) -> float:
    denominator = 0.5 * (abs(first) + abs(second))
    if denominator == 0.0:
        return 0.0 if first == second else math.inf
    return abs(first - second) / denominator


def relative_difference(value: float, reference: float) -> float:
    if reference == 0.0:
        return 0.0 if value == 0.0 else math.inf
    return (value - reference) / reference


def fractional_circle_weights(
    shape: tuple[int, int], centre_x: float, centre_y: float, radius: float
) -> np.ndarray:
    """Exact parity copy of the frozen 64x64 midpoint boundary algorithm."""
    ny, nx = shape
    weights = np.zeros(shape, dtype=np.float64)
    x_min = max(0, int(np.floor(centre_x - radius - 1.0)))
    x_max = min(nx - 1, int(np.ceil(centre_x + radius + 1.0)))
    y_min = max(0, int(np.floor(centre_y - radius - 1.0)))
    y_max = min(ny - 1, int(np.ceil(centre_y + radius + 1.0)))
    x_index = np.arange(x_min, x_max + 1)
    y_index = np.arange(y_min, y_max + 1)
    offsets = (
        (np.arange(APERTURE_SUBPIXELS, dtype=np.float64) + 0.5)
        / APERTURE_SUBPIXELS
        - 0.5
    )
    dx = x_index[:, None] + offsets[None, :] - centre_x
    for start in range(0, y_index.size, 8):
        selected = y_index[start : start + 8]
        dy = selected[:, None] + offsets[None, :] - centre_y
        inside = dy[:, None, :, None] ** 2 + dx[None, :, None, :] ** 2 <= radius**2
        weights[np.ix_(selected, x_index)] = inside.mean(axis=(2, 3))
    return weights


def _exact_grid_mapping() -> dict[int, tuple[str, float, float, str]]:
    return {
        index: (label, mass, amax, execution)
        for index, label, mass, amax, execution in GRID
    }


def _analysis_pinned_inputs() -> dict[str, str]:
    return {
        relative(SOURCE_PREDICTIONS): SOURCE_PREDICTIONS_SHA256,
        relative(SOURCE_FIT_VALIDATION): SOURCE_FIT_VALIDATION_SHA256,
        relative(SOURCE_FIT_POLICY): SOURCE_FIT_POLICY_SHA256,
        relative(SOURCE_FIT_WORKFLOW): SOURCE_FIT_WORKFLOW_SHA256,
        relative(SOURCE_V3_VALIDATION): SOURCE_V3_VALIDATION_SHA256,
        relative(SOURCE_MASK): SOURCE_MASK_SHA256,
    }


def require_declared_outputs_absent() -> None:
    """Refuse a rerun that could overwrite or mix an earlier analysis bundle."""
    existing = [
        relative(path)
        for path in DECLARED_OUTPUTS
        if path.exists() or path.is_symlink()
    ]
    if existing:
        raise ValueError(
            "Refusing analysis with pre-existing declared outputs: "
            + ", ".join(existing)
        )


def load_analysis_manifest(expected_sha256: str) -> dict[str, Any]:
    """Load the prospective, rank-blind analysis contract fail-closed."""
    require_hash(ANALYSIS_MANIFEST, expected_sha256, "refinement analysis manifest")
    analysis = read_json(ANALYSIS_MANIFEST, "refinement analysis manifest")
    expected_identity = {
        "schema_version": 1,
        "analysis_id": "strict_continuum_refinement_v1_analysis",
        "classification": "strict_continuum_refinement_v1_analysis_frozen_rank_blind",
        "screen_only": True,
        "provisional": True,
        "formal_likelihood": False,
        "production_final": False,
        "final_winner_claim_allowed": False,
        "refinement_outputs_or_ranks_read_before_freeze": False,
        "source_screen_used_to_define_refinement_grid": True,
        "source_scores_or_ranks_consumed_by_refinement_analyzer": False,
        "entry_point": relative(ANALYZER),
        "input_manifest": relative(MANIFEST),
        "input_manifest_sha256": INPUT_MANIFEST_SHA256,
        "outputs_absent_at_freeze_verified": True,
        "real_run_requires_declared_outputs_absent": True,
        "real_run_requires_explicit_analysis_manifest_sha256": True,
    }
    changed = [
        key for key, expected in expected_identity.items() if analysis.get(key) != expected
    ]
    if changed:
        raise ValueError(f"Analysis manifest identity/policy changed: {changed}")
    if analysis.get("exact_matrix") != {
        "models": 9,
        "anchors_per_model": 10,
        "expected_cells": 90,
        "fresh_rt_tasks": 80,
        "reused_cells": 10,
    }:
        raise ValueError("Analysis manifest exact matrix changed")
    if analysis.get("canonical_observable") != (
        "seed41001_128k_signed_direct_coeval_method2_polarized_total_i_plane0_"
        "after_same_gaussian_psf_and_nominal_1arcsec_radius_64x64_fractional_"
        "pixel_aperture_scaled_140pc_to_147pc"
    ):
        raise ValueError("Analysis manifest canonical observable changed")
    if analysis.get("uniform_quality_policy") != {
        "significant_component_coeval_relative_difference_max": COEVAL_TOLERANCE,
        "component_significant_fraction": COMPONENT_SIGNIFICANT_FRACTION,
        "component_active_fraction": COMPONENT_ACTIVE_FRACTION,
        "active_component_negative_fraction_max_exclusive": NEGATIVE_FRACTION_MAX,
        "positive_edge_fraction_max_exclusive": EDGE_FRACTION_MAX,
        "psf_convolution_flux_loss_fraction_max_exclusive": CONVOLUTION_LOSS_MAX,
        "component_aperture_fraction_min": COMPONENT_APERTURE_FRACTION_MIN,
        "component_aperture_fraction_max": COMPONENT_APERTURE_FRACTION_MAX,
        "aperture_clip_symmetric_fraction_max": APERTURE_CLIP_MAX,
        "retain_intrinsic_separated_component_product_integrity_gates": True,
        "plane0_vs_sum_planes4to7_gate_applied": False,
        "aggregate_componentwise_flux_gate_applied": False,
        "signed_flux_used_for_scoring": True,
        "zero_clipped_flux_used_for_scoring": False,
        "full_plane_negative_fraction_role": "diagnostic_only",
    }:
        raise ValueError("Analysis manifest uniform quality policy changed")
    if analysis.get("ranking_policy") != {
        "metric": "unchanged_six_region_balanced_fixed_normalization_log10_rms",
        "free_normalization_fitted": False,
        "boundary_scenarios": {
            "primary": 1.0,
            "boundary_half": 0.5,
            "boundary_drop": 0.0,
        },
        "boundary_bins": [96, 139],
        "promotion_margin_dex": OPERATIONAL_SCORE_MARGIN_DEX,
        "promotion_set": "union_of_rank1_and_near_ties_in_all_three_scenarios",
        "independent_seed41002_512k_confirmation_required_before_production_claim": True,
    }:
        raise ValueError("Analysis manifest ranking policy changed")
    analyzer_digest = str(analysis.get("analyzer_sha256", ""))
    if len(analyzer_digest) != 64:
        raise ValueError("Analysis manifest lacks an analyzer SHA-256")
    require_hash(ANALYZER, analyzer_digest, "analysis-manifest-pinned analyzer")
    if analysis.get("pinned_inputs") != _analysis_pinned_inputs():
        raise ValueError("Analysis manifest pinned-input table changed")
    for name, digest in _analysis_pinned_inputs().items():
        require_hash(resolve(name), digest, f"analysis-manifest input {name}")
    declared = [relative(path) for path in DECLARED_OUTPUTS]
    if analysis.get("outputs_absent_at_freeze") != declared:
        raise ValueError("Analysis manifest prospective output declaration changed")
    require_hash(
        ANALYSIS_MANIFEST,
        expected_sha256,
        "refinement analysis manifest post-read",
    )
    return analysis


def load_manifest(expected_sha256: str) -> dict[str, Any]:
    if expected_sha256 != INPUT_MANIFEST_SHA256:
        raise ValueError("Analysis manifest does not pin the expected refinement input")
    require_hash(MANIFEST, expected_sha256, "refinement input manifest")
    manifest = read_json(MANIFEST, "refinement manifest")
    execution = manifest.get("execution", {})
    reuse = manifest.get("reuse", {})
    expected = {
        "schema_version": 1,
        "workflow_id": "strict_continuum_refinement_v1",
        "classification": "frozen_parallel_strict_method2_refinement_screen",
        "screen_only": True,
        "production_final": False,
        "post_analysis_included": False,
    }
    if any(manifest.get(key) != value for key, value in expected.items()):
        raise ValueError("Refinement manifest identity/classification changed")
    if not (
        manifest.get("grid")
        == {
            "mass_factors": [2.0, 2.25, 2.5],
            "envelope_amax_um": [0.5, 0.6, 0.7],
            "model_count": 9,
        }
        and execution.get("reused_model_count") == 1
        and execution.get("fresh_model_count") == 8
        and execution.get("temperature_task_count") == 8
        and execution.get("rt_task_count") == 80
        and execution.get("anchors_per_model") == 10
        and execution.get("seed") == EXPECTED_SEED
        and execution.get("accepted_packets") == EXPECTED_PACKETS
        and execution.get("openmp_threads") == EXPECTED_THREADS
        and execution.get("one_wave_coeval_method2") is True
        and reuse.get("model_label") == GRID[0][1]
        and reuse.get("upstream_model_index") == 15
        and reuse.get("upstream_rt_task_indices") == list(range(150, 160))
        and reuse.get("copy_or_symlink_performed") is False
        and reuse.get("fit_validation_sha256") == SOURCE_FIT_VALIDATION_SHA256
    ):
        raise ValueError("Refinement manifest grid/execution/reuse semantics changed")
    for group in ("pinned_upstream", "frozen_files"):
        values = manifest.get(group)
        if not isinstance(values, dict) or not values:
            raise ValueError(f"Refinement manifest lacks {group}")
        for name, digest in values.items():
            require_hash(resolve(name), str(digest), f"manifest {group} file {name}")
    require_hash(SOURCE_PREDICTIONS, SOURCE_PREDICTIONS_SHA256, "reused predictions")
    require_hash(SOURCE_FIT_VALIDATION, SOURCE_FIT_VALIDATION_SHA256, "reused fit validation")
    require_hash(SOURCE_V3_VALIDATION, SOURCE_V3_VALIDATION_SHA256, "reused v3 validation")
    require_hash(SOURCE_MASK, SOURCE_MASK_SHA256, "strict Scenario-B mask")
    return manifest


def validate_tables(manifest: dict[str, Any]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    models = read_tsv(MODELS_TABLE)
    rows = read_tsv(RT_TASKS_TABLE)
    if len(models) != 9 or [int(row["grid_model_index"]) for row in models] != list(range(9)):
        raise ValueError("models.tsv is not the exact nine-model ordered grid")
    grid = _exact_grid_mapping()
    for row in models:
        index = int(row["grid_model_index"])
        label, mass, amax, execution = grid[index]
        fresh = "" if index == 0 else str(index - 1)
        if not (
            row["model_label"] == label
            and math.isclose(float(row["mass_factor"]), mass, abs_tol=1e-14)
            and math.isclose(float(row["envelope_amax_um"]), amax, abs_tol=1e-14)
            and row["execution"] == execution
            and row["fresh_model_index"] == fresh
        ):
            raise ValueError(f"Grid model row {index} changed")
        require_hash(resolve(row["temperature_parameter"]), row["temperature_parameter_sha256"], f"model {label} T parameter")
        require_hash(resolve(row["rt_parameter"]), row["rt_parameter_sha256"], f"model {label} RT parameter")

    if len(rows) != 80 or set(rows[0]) != REQUIRED_RT_FIELDS:
        missing = sorted(REQUIRED_RT_FIELDS.difference(rows[0] if rows else set()))
        extra = sorted(set(rows[0] if rows else set()).difference(REQUIRED_RT_FIELDS))
        raise ValueError(f"rt_tasks.tsv schema/count changed; missing={missing}, extra={extra}")
    if [int(row["task_index"]) for row in rows] != list(range(80)):
        raise ValueError("Fresh RT task indices are not exactly 0..79")
    seen: set[tuple[int, int]] = set()
    for row in rows:
        task = int(row["task_index"])
        fresh = int(row["fresh_model_index"])
        grid_index = int(row["grid_model_index"])
        ordinal = int(row["anchor_ordinal"])
        if not (fresh == grid_index - 1 and task == fresh * 10 + ordinal):
            raise ValueError(f"Malformed fresh/grid/task identity at task {task}")
        label, mass, amax, execution = grid[grid_index]
        anchor = ANCHORS[ordinal]
        bin_index, token, region, spatial, compromise = anchor
        expected_fwhm = psf_fwhm_arcsec(row["instrument"], float(token))
        if not (
            execution == "fresh"
            and row["model_label"] == label
            and math.isclose(float(row["mass_factor"]), mass, abs_tol=1e-14)
            and math.isclose(float(row["envelope_amax_um"]), amax, abs_tol=1e-14)
            and int(row["bin_index"]) == bin_index
            and row["anchor_tag"] == f"w{ordinal:03d}"
            and row["region"] == region
            and row["spatial_config"] == spatial
            and (row["boundary_compromise"] == "true") is compromise
            and row["wavelength_um"] == token
            and row["wavelength_token"] == token
            and int(row["seed"]) == EXPECTED_SEED
            and int(row["sed_photons"]) == EXPECTED_PACKETS
            and math.isclose(float(row["psf_fwhm_arcsec"]), expected_fwhm, abs_tol=1e-14)
            and int(row["temperature_task_index"]) == fresh
        ):
            raise ValueError(f"Fresh task {task} physical/anchor identity changed")
        require_hash(resolve(row["parameter"]), row["parameter_sha256"], f"task {task} RT parameter")
        require_hash(resolve(row["wavelength_file"]), row["wavelength_file_sha256"], f"task {task} wavelength")
        pair = (grid_index, ordinal)
        if pair in seen:
            raise ValueError(f"Duplicate refinement cell {pair}")
        seen.add(pair)
    if seen != {(grid_index, ordinal) for grid_index in range(1, 9) for ordinal in range(10)}:
        raise ValueError("Fresh tasks are not the exact eight-by-ten matrix")
    # Bind again after parsing to catch concurrent replacement.
    frozen = manifest["frozen_files"]
    require_hash(MODELS_TABLE, frozen[relative(MODELS_TABLE)], "models.tsv post-read")
    require_hash(RT_TASKS_TABLE, frozen[relative(RT_TASKS_TABLE)], "rt_tasks.tsv post-read")
    return models, rows


def combine_observation_rows(rows: Iterable[Any]) -> dict[str, Any]:
    selected = list(rows)
    if not selected:
        raise ValueError("Cannot combine an empty observation bin")
    flux = np.asarray([float(row["lambda_f_lambda_w_m2"]) for row in selected])
    error = np.asarray([float(row["lambda_f_lambda_error_w_m2"]) for row in selected])
    floor = np.asarray([float(row["systematic_fraction"]) for row in selected])
    wavelength = np.asarray([float(row["wavelength_um"]) for row in selected])
    if (
        np.any(~np.isfinite(flux))
        or np.any(flux <= 0.0)
        or np.any(~np.isfinite(error))
        or np.any(error <= 0.0)
        or np.any(~np.isfinite(floor))
        or np.any(floor <= 0.0)
        or not np.allclose(wavelength, wavelength[0], rtol=0.0, atol=1e-13)
        or not np.allclose(floor, floor[0], rtol=0.0, atol=1e-12)
    ):
        raise ValueError("Observation rows have invalid flux/error/wavelength/floor")
    fractional = error / flux
    if np.any(fractional + 1e-14 < floor):
        raise ValueError("An observation error omits its systematic floor")
    common_fraction = float(floor[0])
    independent_variance = np.maximum(fractional**2 - floor**2, 0.0)
    combined_fraction = math.sqrt(
        common_fraction**2
        + float(np.sum(independent_variance, dtype=np.float64)) / len(selected) ** 2
    )
    combined_flux = float(10.0 ** np.mean(np.log10(flux)))
    return {
        "wavelength_um": float(wavelength[0]),
        "n_detector_rows": len(selected),
        "segments": ";".join(sorted({str(row["segment"]) for row in selected})),
        "observed_flux_geometric_w_m2": combined_flux,
        "observed_error_w_m2": combined_flux * combined_fraction,
        "observed_fractional_error": combined_fraction,
        "observed_log10_sigma_dex": math.log10(1.0 + combined_fraction),
    }


def load_observations() -> dict[int, dict[str, Any]]:
    table = Table.read(SOURCE_MASK, format="ascii.ecsv")
    required = {
        "instrument",
        "segment",
        "bin_index",
        "wavelength_um",
        "lambda_f_lambda_w_m2",
        "lambda_f_lambda_error_w_m2",
        "systematic_fraction",
    }
    if not required.issubset(table.colnames):
        raise ValueError("Strict mask observation schema changed")
    if not (
        table.meta.get("fit_mask_scenario_id") == "defer_major_ice_and_both_silicates"
        and table.meta.get("fit_mask_weighting_applied") is False
        and table.meta.get("segment_rescaling_applied") is False
        and table.meta.get("background_subtraction_applied") is False
        and math.isclose(float(table.meta.get("target_resolution", math.nan)), 100.0)
    ):
        raise ValueError("Strict Scenario-B metadata changed")
    bins = np.asarray(table["bin_index"], dtype=int)
    result: dict[int, dict[str, Any]] = {}
    for bin_index, token, region, _, compromise in ANCHORS:
        combined = combine_observation_rows(table[bins == bin_index])
        if not math.isclose(combined["wavelength_um"], float(token), abs_tol=1e-13):
            raise ValueError(f"Observed anchor wavelength changed at bin {bin_index}")
        combined.update(
            bin_index=bin_index,
            region=region,
            boundary_compromise=compromise,
        )
        result[bin_index] = combined
    if set(result) != {row[0] for row in ANCHORS}:
        raise ValueError("Strict mask did not yield exactly ten observations")
    require_hash(SOURCE_MASK, SOURCE_MASK_SHA256, "strict mask post-read")
    return result


def _completion_paths(row: dict[str, str]) -> dict[str, Path]:
    root = resolve(row["task_workdir"])
    data = root / "paired/seed=41001/data_th"
    return {
        "root": root,
        "completion": root / "complete.json",
        "sed_rt": data / "sed_rt.fits.gz",
        "image_rt": data / "RT.fits.gz",
        "closure": root / "closure.json",
        "log": root / "logs/paired_method2.log",
    }


def validate_completion(
    row: dict[str, str], manifest: dict[str, Any], manifest_sha: str
) -> tuple[dict[str, Path], dict[str, str]]:
    task = int(row["task_index"])
    paths = _completion_paths(row)
    completion_sha = require_regular(paths["completion"], f"task {task} completion")
    completion = read_json(paths["completion"], f"task {task} completion")
    expected = {
        "schema_version": 1,
        "workflow_id": "strict_continuum_refinement_v1",
        "kind": "rt",
        "task_index": task,
        "model_label": row["model_label"],
        "seed": EXPECTED_SEED,
        "accepted_packets": EXPECTED_PACKETS,
        "threads": EXPECTED_THREADS,
        "manifest_sha256": manifest_sha,
    }
    changed = [key for key, value in expected.items() if completion.get(key) != value]
    if changed:
        raise ValueError(f"Task {task} completion fields changed: {changed}")
    command = completion.get("command")
    expected_command = [
        str(manifest["custom_binary"]["path"]),
        Path(row["parameter"]).name,
        "-cavity",
        "274.7",
        "100",
        "1.0",
        "-Mdot",
        "1",
        "7.3011756855e-08",
        "-no_T",
        "-Tfile",
        "Temperature.fits.gz",
        "-zoom",
        row["zoom_token"],
        "-resol",
        row["nx"],
        row["ny"],
        "-rt2",
        "-rt-sed-method",
        "2",
        "-seed",
        row["seed"],
        "-root_dir",
        "../paired",
        "-no_backup",
    ]
    if not isinstance(command, list) or list(map(str, command)) != expected_command:
        raise ValueError(f"Task {task} completion command contract changed")
    products = completion.get("products")
    if not isinstance(products, dict) or not {"sed_rt", "image_rt", "closure", "log"}.issubset(products):
        raise ValueError(f"Task {task} completion product table changed")
    expected_hashes = {"completion": completion_sha}
    for name in ("sed_rt", "image_rt", "closure", "log"):
        item = products[name]
        if not isinstance(item, dict):
            raise ValueError(f"Task {task} malformed completion product {name}")
        recorded = resolve(item.get("path"))
        if recorded != paths[name].resolve():
            raise ValueError(f"Task {task} completion path changed for {name}")
        digest = str(item.get("sha256", ""))
        require_hash(paths[name], digest, f"task {task} {name}")
        expected_hashes[name] = digest
    marker = str(manifest["custom_binary"]["runtime_marker"])
    if marker not in paths["log"].read_text(encoding="utf-8", errors="replace"):
        raise ValueError(f"Task {task} coeval method-2 runtime marker is absent")
    return paths, expected_hashes


def assert_products_unchanged(paths: dict[str, Path], hashes: dict[str, str], task: int) -> None:
    for name in ("completion", "sed_rt", "image_rt", "closure", "log"):
        require_hash(paths[name], hashes[name], f"task {task} post-read {name}")


def validate_closure(
    closure_path: Path,
    sed_path: Path,
    image_path: Path,
    wavelength: float,
) -> dict[str, Any]:
    closure = read_json(closure_path, "coeval closure report")
    inputs = closure.get("inputs", {})
    if not (
        closure.get("schema_version") == 1
        and closure.get("passes_requested_closure_tolerance") is True
        and math.isclose(float(closure.get("tolerance", math.nan)), COEVAL_TOLERANCE, abs_tol=0.0)
        and math.isclose(
            float(closure.get("significant_fraction", math.nan)),
            COMPONENT_SIGNIFICANT_FRACTION,
            abs_tol=0.0,
        )
        and math.isfinite(
            float(closure.get("max_significant_custom_image_relative_difference", math.inf))
        )
        and float(closure.get("max_significant_custom_image_relative_difference", math.inf))
        <= COEVAL_TOLERANCE
        and resolve(inputs.get("custom_sed")) == sed_path.resolve()
        and resolve(inputs.get("image")) == image_path.resolve()
        and inputs.get("custom_sed_sha256") == sha256(sed_path)
        and inputs.get("image_sha256") == sha256(image_path)
        and inputs.get("stock_sed") is None
        and inputs.get("stock_sed_sha256") is None
        and closure.get("stock_component_closure_relative") is None
        and math.isclose(float(closure.get("wavelength_um", math.nan)), wavelength, abs_tol=1e-6)
    ):
        raise ValueError("Coeval SED/map closure report failed")
    planes = closure.get("planes")
    if not isinstance(planes, list) or len(planes) != len(CLOSURE_LABELS):
        raise ValueError("Closure report plane table changed")
    selected = {
        int(row["plane"]): row
        for row in planes
        if isinstance(row, dict) and row.get("plane") in CLOSURE_LABELS
    }
    if set(selected) != set(CLOSURE_LABELS):
        raise ValueError("Closure report is not the exact plane 0,4,5,6,7 set")
    for plane, label in CLOSURE_LABELS.items():
        if selected[plane].get("label") != label:
            raise ValueError(f"Closure report plane-{plane} semantics changed")
    row = selected[0]
    if not (
        row.get("label") == "total_I"
        and row.get("significant") is True
        and abs(float(row.get("custom_relative_difference", math.inf))) <= COEVAL_TOLERANCE
    ):
        raise ValueError("Coeval plane-0 closure failed")
    return closure


def _closure_rows(closure: dict[str, Any]) -> dict[int, dict[str, Any]]:
    return {int(row["plane"]): row for row in closure["planes"]}


def measure_fresh_component_integrity(
    *,
    sed: np.ndarray,
    image: np.ndarray,
    closure: dict[str, Any],
    sigma_pixels: float,
    weights: np.ndarray,
    border_mask: np.ndarray,
) -> dict[str, Any]:
    """Independently retain the v3 intrinsic plane-4:7 product gates."""
    row_by_plane = _closure_rows(closure)
    np_totals = {
        plane: float(np.sum(image[plane, 0, 0], dtype=np.float64))
        for plane in CLOSURE_LABELS
    }
    total_scale = max(abs(float(sed[0])), abs(np_totals[0]))
    significance_cut = COMPONENT_SIGNIFICANT_FRACTION * total_scale
    significant_errors: dict[int, float] = {}
    for plane, label in CLOSURE_LABELS.items():
        sed_value = float(sed[plane])
        map_value = np_totals[plane]
        significant = plane == 0 or max(abs(sed_value), abs(map_value)) > significance_cut
        coeval = relative_difference(sed_value, map_value)
        reported = row_by_plane[plane]
        if not (
            reported.get("label") == label
            and reported.get("significant") is significant
            and math.isclose(
                float(reported.get("custom_sed_flux", math.nan)),
                sed_value,
                rel_tol=1.0e-13,
                abs_tol=0.0,
            )
            and math.isclose(
                float(reported.get("image_flux", math.nan)),
                map_value,
                rel_tol=1.0e-13,
                abs_tol=0.0,
            )
            and math.isclose(
                float(reported.get("custom_relative_difference", math.nan)),
                coeval,
                rel_tol=1.0e-10,
                abs_tol=1.0e-14,
            )
        ):
            raise ValueError(f"Independent closure report binding failed for plane {plane}")
        if significant:
            error = abs(coeval)
            if not math.isfinite(error) or error > COEVAL_TOLERANCE:
                raise ValueError(
                    f"Independent significant-plane closure failed for plane {plane}: {error}"
                )
            significant_errors[plane] = error
    recomputed_max = max(significant_errors.values(), default=math.inf)
    if not math.isclose(
        recomputed_max,
        float(closure["max_significant_custom_image_relative_difference"]),
        rel_tol=1.0e-10,
        abs_tol=1.0e-14,
    ):
        raise ValueError("Closure report maximum does not match independently reopened planes")

    component_l1 = float(np.sum(np.abs(sed[4:8]), dtype=np.float64))
    component_rows: dict[str, dict[str, Any]] = {}
    for plane, name in COMPONENTS.items():
        values = np.asarray(image[plane, 0, 0], dtype=np.float64)
        map_total = math.fsum(float(value) for value in values.ravel(order="C"))
        map_total_numpy = np_totals[plane]
        sed_component = float(sed[plane])
        positive = float(np.sum(values[values > 0.0], dtype=np.float64))
        negative = float(np.sum(np.abs(values[values < 0.0]), dtype=np.float64))
        negative_fraction = (
            negative / positive
            if positive > 0.0
            else (0.0 if negative == 0.0 else float(np.finfo(float).max))
        )
        edge_fraction = (
            float(
                np.sum(np.clip(values, 0.0, None)[border_mask], dtype=np.float64)
            )
            / positive
            if positive > 0.0
            else 0.0
        )
        active = (
            abs(sed_component) >= COMPONENT_ACTIVE_FRACTION * component_l1
            if component_l1 > 0.0
            else False
        )
        significant = plane in significant_errors
        coeval = relative_difference(sed_component, map_total_numpy)
        if positive > 0.0:
            convolved = gaussian_filter(
                values,
                sigma_pixels,
                mode="constant",
                cval=0.0,
                truncate=PSF_TRUNCATION_SIGMA,
            )
            convolved_total = float(np.sum(convolved, dtype=np.float64))
            map_aperture = float(np.sum(convolved * weights, dtype=np.float64))
            convolution_loss = abs(convolved_total - map_total) / max(
                abs(map_total), np.finfo(float).tiny
            )
            del convolved
        else:
            map_aperture = 0.0
            convolution_loss = (
                0.0 if map_total == 0.0 else float(np.finfo(float).max)
            )
        zero_map_contract = map_total != 0.0 or sed_component == 0.0
        aperture_fraction = 0.0 if map_total == 0.0 else map_aperture / map_total
        finite_map_contract = all(
            math.isfinite(value)
            for value in (
                sed_component,
                map_total,
                map_total_numpy,
                positive,
                negative,
                negative_fraction,
                edge_fraction,
                map_aperture,
                aperture_fraction,
                convolution_loss,
                coeval,
            )
        )
        aperture_fraction_contract = finite_map_contract and (
            COMPONENT_APERTURE_FRACTION_MIN
            <= aperture_fraction
            <= COMPONENT_APERTURE_FRACTION_MAX
        )
        active_morphology_contract = (not active) or (
            0.0 <= negative_fraction < NEGATIVE_FRACTION_MAX
            and 0.0 <= edge_fraction < EDGE_FRACTION_MAX
            and 0.0 <= convolution_loss < CONVOLUTION_LOSS_MAX
        )
        quality = (
            zero_map_contract
            and finite_map_contract
            and aperture_fraction_contract
            and active_morphology_contract
            and ((not significant) or abs(coeval) <= COEVAL_TOLERANCE)
        )
        if not quality:
            raise ValueError(f"Intrinsic separated-component gate failed: {name}")
        component_rows[name] = {
            "plane_zero_based": plane,
            "significant_ge_1e-8_total_scale": significant,
            "active_ge_1pct_component_l1": active,
            "sed_flux_140pc_w_m2": sed_component,
            "map_total_fixed_fsum_140pc_w_m2": map_total,
            "map_total_numpy_140pc_w_m2": map_total_numpy,
            "map_min_pixel_flux_w_m2": float(np.min(values)),
            "map_max_pixel_flux_w_m2": float(np.max(values)),
            "map_positive_flux_140pc_w_m2": positive,
            "map_negative_magnitude_140pc_w_m2": negative,
            "coeval_sed_minus_map_over_map": coeval,
            "negative_flux_fraction": negative_fraction,
            "edge_positive_flux_fraction": edge_fraction,
            "convolution_flux_loss_fraction": convolution_loss,
            "psf_map_aperture_flux_140pc_w_m2": map_aperture,
            "psf_map_aperture_fraction": aperture_fraction,
            "map_pixels_finite": True,
            "zero_map_requires_zero_sed_component_pass": zero_map_contract,
            "finite_map_contract_pass": finite_map_contract,
            "aperture_fraction_contract_pass": aperture_fraction_contract,
            "active_morphology_contract_pass": active_morphology_contract,
            "quality_pass": quality,
        }
        del values
    significant_components = [
        item
        for item in component_rows.values()
        if item["significant_ge_1e-8_total_scale"]
    ]
    active_components = [
        item for item in component_rows.values() if item["active_ge_1pct_component_l1"]
    ]
    return {
        "evidence_kind": "direct_recomputation_from_coeval_sed_and_image_products",
        "intrinsic_component_product_integrity_gates_applied": True,
        "plane0_vs_sum_planes4to7_gate_applied": False,
        "aggregate_componentwise_flux_gate_applied": False,
        "significant_fraction_of_total_scale": COMPONENT_SIGNIFICANT_FRACTION,
        "active_fraction_of_component_l1": COMPONENT_ACTIVE_FRACTION,
        "component_l1_sed_flux_140pc_w_m2": component_l1,
        "significant_component_count": len(significant_components),
        "active_component_count": len(active_components),
        "max_abs_significant_all_plane_coeval_relative_difference": recomputed_max,
        "max_abs_significant_component_coeval_relative_difference": max(
            (abs(float(item["coeval_sed_minus_map_over_map"])) for item in significant_components),
            default=0.0,
        ),
        "max_active_component_negative_flux_fraction": max(
            (float(item["negative_flux_fraction"]) for item in active_components),
            default=0.0,
        ),
        "max_active_component_edge_positive_flux_fraction": max(
            (float(item["edge_positive_flux_fraction"]) for item in active_components),
            default=0.0,
        ),
        "max_active_component_convolution_flux_loss_fraction": max(
            (float(item["convolution_flux_loss_fraction"]) for item in active_components),
            default=0.0,
        ),
        "components": component_rows,
        "quality_pass": True,
    }


def validate_reused_component_integrity(v3_task: dict[str, Any]) -> dict[str, Any]:
    """Normalize and revalidate pinned v3 component evidence for the reused corner."""
    measurement = v3_task.get("measurement", {})
    source_components = measurement.get("components", {})
    if not (
        v3_task.get("retained_intrinsic_component_product_integrity_gates") is True
        and v3_task.get("aggregate_componentwise_estimator_affects_validity") is False
        and measurement.get("intrinsic_component_product_integrity_gates_applied") is True
        and measurement.get("total_component_closure_gate_applied") is False
        and measurement.get("aggregate_componentwise_finite_positive_gate_applied") is False
        and set(source_components) == set(COMPONENTS.values())
    ):
        raise ValueError("Pinned upstream intrinsic component contract changed")
    component_l1 = math.fsum(
        abs(float(source_components[name]["sed_flux_140pc_w_m2"]))
        for name in COMPONENTS.values()
    )
    plane0 = v3_task["plane0_total_i"]
    total_scale = max(
        abs(float(plane0["coeval_plane0_custom_sed_flux_140pc_w_m2"])),
        abs(float(plane0["coeval_plane0_image_flux_140pc_w_m2"])),
    )
    significance_cut = COMPONENT_SIGNIFICANT_FRACTION * total_scale
    component_rows: dict[str, dict[str, Any]] = {}
    for plane, name in COMPONENTS.items():
        source = source_components[name]
        sed_component = float(source["sed_flux_140pc_w_m2"])
        map_total = float(source["map_total_140pc_w_m2"])
        map_aperture = float(source["psf_map_aperture_flux_140pc_w_m2"])
        aperture_fraction = float(source["psf_map_aperture_fraction"])
        negative_fraction = float(source["negative_flux_fraction"])
        edge_fraction = float(source["edge_positive_flux_fraction"])
        convolution_loss = float(source["convolution_flux_loss_fraction"])
        coeval = relative_difference(sed_component, map_total)
        significant = max(abs(sed_component), abs(map_total)) > significance_cut
        active = (
            abs(sed_component) >= COMPONENT_ACTIVE_FRACTION * component_l1
            if component_l1 > 0.0
            else False
        )
        finite_map_contract = all(
            math.isfinite(value)
            for value in (
                sed_component,
                map_total,
                map_aperture,
                aperture_fraction,
                negative_fraction,
                edge_fraction,
                convolution_loss,
                coeval,
            )
        )
        zero_map_contract = map_total != 0.0 or sed_component == 0.0
        aperture_fraction_contract = finite_map_contract and (
            COMPONENT_APERTURE_FRACTION_MIN
            <= aperture_fraction
            <= COMPONENT_APERTURE_FRACTION_MAX
        )
        active_morphology_contract = (not active) or (
            0.0 <= negative_fraction < NEGATIVE_FRACTION_MAX
            and 0.0 <= edge_fraction < EDGE_FRACTION_MAX
            and 0.0 <= convolution_loss < CONVOLUTION_LOSS_MAX
        )
        quality = (
            int(source.get("plane_zero_based", -1)) == plane
            and source.get("active_ge_1pct") is active
            and source.get("quality_pass") is True
            and finite_map_contract
            and zero_map_contract
            and aperture_fraction_contract
            and active_morphology_contract
            and ((not significant) or abs(coeval) <= COEVAL_TOLERANCE)
        )
        if not quality:
            raise ValueError(f"Pinned upstream component integrity changed: {name}")
        component_rows[name] = {
            "plane_zero_based": plane,
            "significant_ge_1e-8_total_scale": significant,
            "active_ge_1pct_component_l1": active,
            "sed_flux_140pc_w_m2": sed_component,
            "map_total_fixed_fsum_140pc_w_m2": map_total,
            "coeval_sed_minus_map_over_map": coeval,
            "negative_flux_fraction": negative_fraction,
            "edge_positive_flux_fraction": edge_fraction,
            "convolution_flux_loss_fraction": convolution_loss,
            "psf_map_aperture_flux_140pc_w_m2": map_aperture,
            "psf_map_aperture_fraction": aperture_fraction,
            "map_pixels_finite": True,
            "map_pixels_finite_evidence": "pinned_upstream_v3_validation",
            "zero_map_requires_zero_sed_component_pass": zero_map_contract,
            "finite_map_contract_pass": finite_map_contract,
            "aperture_fraction_contract_pass": aperture_fraction_contract,
            "active_morphology_contract_pass": active_morphology_contract,
            "quality_pass": quality,
        }
    significant_components = [
        item
        for item in component_rows.values()
        if item["significant_ge_1e-8_total_scale"]
    ]
    active_components = [
        item for item in component_rows.values() if item["active_ge_1pct_component_l1"]
    ]
    source_closure_max = float(v3_task.get("coeval_max_relative_difference", math.inf))
    if not math.isfinite(source_closure_max) or source_closure_max > COEVAL_TOLERANCE:
        raise ValueError("Pinned upstream significant-plane closure changed")
    return {
        "evidence_kind": "pinned_upstream_v3_intrinsic_component_validation",
        "intrinsic_component_product_integrity_gates_applied": True,
        "plane0_vs_sum_planes4to7_gate_applied": False,
        "aggregate_componentwise_flux_gate_applied": False,
        "significant_fraction_of_total_scale": COMPONENT_SIGNIFICANT_FRACTION,
        "active_fraction_of_component_l1": COMPONENT_ACTIVE_FRACTION,
        "component_l1_sed_flux_140pc_w_m2": component_l1,
        "significant_component_count": len(significant_components),
        "active_component_count": len(active_components),
        "max_abs_significant_all_plane_coeval_relative_difference": source_closure_max,
        "max_abs_significant_component_coeval_relative_difference": max(
            (abs(float(item["coeval_sed_minus_map_over_map"])) for item in significant_components),
            default=0.0,
        ),
        "max_active_component_negative_flux_fraction": max(
            (float(item["negative_flux_fraction"]) for item in active_components),
            default=0.0,
        ),
        "max_active_component_edge_positive_flux_fraction": max(
            (float(item["edge_positive_flux_fraction"]) for item in active_components),
            default=0.0,
        ),
        "max_active_component_convolution_flux_loss_fraction": max(
            (float(item["convolution_flux_loss_fraction"]) for item in active_components),
            default=0.0,
        ),
        "components": component_rows,
        "quality_pass": True,
    }


def measure_fresh_task(
    row: dict[str, str], manifest: dict[str, Any], manifest_sha: str
) -> dict[str, Any]:
    task = int(row["task_index"])
    paths, hashes = validate_completion(row, manifest, manifest_sha)
    wavelength = float(row["wavelength_token"])
    closure = validate_closure(
        paths["closure"], paths["sed_rt"], paths["image_rt"], wavelength
    )
    closure_plane0 = _closure_rows(closure)[0]

    with fits.open(paths["sed_rt"], memmap=False) as hdus:
        sed_raw = np.array(hdus[0].data, dtype=np.float64, copy=True)
        sed_shape = tuple(sed_raw.shape)
        sed_header = hdus[0].header.copy()
        sed_bitpix = int(sed_header.get("BITPIX", 0))
        sed_unit = str(sed_header.get("BUNIT", "")).lower().replace(" ", "")
        wave = np.asarray(hdus[1].data, dtype=np.float64).ravel()
    sed = sed_raw[:, 0, 0, 0]
    sed_plane0 = float(sed[0])
    nx, ny = int(row["nx"]), int(row["ny"])
    with fits.open(paths["image_rt"], memmap=False) as hdus:
        image_raw = np.array(hdus[0].data, copy=True)
        image_shape = tuple(image_raw.shape)
        image_bitpix = int(hdus[0].header.get("BITPIX", 0))
        header = hdus[0].header.copy()
        image_unit = str(header.get("BUNIT", "")).lower().replace(" ", "")
    plane0 = np.asarray(image_raw[0, 0, 0], dtype=np.float64)
    expected_component_tokens = {
        5: ("direct", "star"),
        6: ("scattered", "star"),
        7: ("direct", "thermal"),
        8: ("scattered", "thermal"),
    }
    component_semantics_pass = all(
        all(token in str(header.get(f"FLUX_{header_plane}", "")).lower() for token in tokens)
        for header_plane, tokens in expected_component_tokens.items()
    )
    if not (
        sed_shape == (8, 1, 1, 1)
        and sed_bitpix == -64
        and wave.shape == (1,)
        and np.all(np.isfinite(sed_raw))
        and image_shape == (8, 1, 1, ny, nx)
        and image_bitpix == -32
        and np.all(np.isfinite(image_raw))
        and sed_unit in {"w.m-2", "w/m2"}
        and image_unit == "w.m-2.pixel-1"
        and "total flux" in str(header.get("FLUX_1", "")).lower()
        and component_semantics_pass
        and math.isclose(float(wave[0]), wavelength, rel_tol=1e-6, abs_tol=1e-8)
        and math.isclose(float(header.get("WAVE", math.nan)), float(np.float32(wavelength)), rel_tol=1e-6, abs_tol=1e-8)
    ):
        raise ValueError(f"Task {task} FITS shape/precision/unit/wavelength gate failed")
    for product_header, label in ((sed_header, "SED"), (header, "map")):
        if any(key in product_header for key in ("BSCALE", "BZERO")):
            raise ValueError(f"Task {task} {label} uses forbidden FITS scaling")
    if any(key in header for key in ("PC1_1", "PC1_2", "PC2_1", "PC2_2", "CD1_1", "CD1_2", "CD2_1", "CD2_2", "CROTA1", "CROTA2")):
        raise ValueError(f"Task {task} map WCS matrix/rotation overrides are forbidden")

    cdelt1 = float(header.get("CDELT1", math.nan))
    cdelt2 = float(header.get("CDELT2", math.nan))
    expected_cdelt = float(row["expected_cdelt_deg"])
    centre_x = float(header.get("CRPIX1", math.nan)) - 1.0
    centre_y = float(header.get("CRPIX2", math.nan)) - 1.0
    if not (
        math.isclose(abs(cdelt1), expected_cdelt, rel_tol=5e-6)
        and math.isclose(abs(cdelt2), expected_cdelt, rel_tol=5e-6)
        and math.isclose(centre_x, (nx - 1) / 2.0, abs_tol=1e-6)
        and math.isclose(centre_y, (ny - 1) / 2.0, abs_tol=1e-6)
    ):
        raise ValueError(f"Task {task} map scale/centre changed")
    target_pixel_arcsec = abs(cdelt1) * 3600.0 * MODEL_DISTANCE_PC / TARGET_DISTANCE_PC
    expected_target_pixel = expected_cdelt * 3600.0 * MODEL_DISTANCE_PC / TARGET_DISTANCE_PC
    radius_pixels = APERTURE_RADIUS_ARCSEC / target_pixel_arcsec
    fwhm = psf_fwhm_arcsec(row["instrument"], wavelength)
    sigma_pixels = fwhm / 2.354820045 / target_pixel_arcsec
    half_field = min(
        (centre_x + 0.5) * target_pixel_arcsec,
        (nx - 0.5 - centre_x) * target_pixel_arcsec,
        (centre_y + 0.5) * target_pixel_arcsec,
        (ny - 0.5 - centre_y) * target_pixel_arcsec,
    )
    support_minimum = APERTURE_RADIUS_ARCSEC + PSF_TRUNCATION_SIGMA * fwhm / 2.354820045
    if not (
        math.isclose(target_pixel_arcsec, expected_target_pixel, rel_tol=5e-6)
        and math.isclose(float(row["psf_fwhm_arcsec"]), fwhm, abs_tol=1e-14)
        and half_field >= support_minimum
    ):
        raise ValueError(f"Task {task} PSF/aperture support geometry failed")

    map_total = math.fsum(float(value) for value in plane0.ravel(order="C"))
    map_total_numpy = float(np.sum(plane0, dtype=np.float64))
    coeval_relative = relative_difference(sed_plane0, map_total)
    reported_relative = float(closure_plane0["custom_relative_difference"])
    if not (
        abs(coeval_relative) <= COEVAL_TOLERANCE
        and math.isclose(float(closure_plane0["custom_sed_flux"]), sed_plane0, rel_tol=1e-13, abs_tol=0.0)
        and math.isclose(float(closure_plane0["image_flux"]), map_total_numpy, rel_tol=1e-13, abs_tol=0.0)
        and math.isclose(reported_relative, relative_difference(sed_plane0, map_total_numpy), rel_tol=1e-10, abs_tol=1e-14)
    ):
        raise ValueError(f"Task {task} independently recomputed plane-0 closure failed")

    positive = float(np.sum(plane0[plane0 > 0.0], dtype=np.float64))
    negative = float(np.sum(np.abs(plane0[plane0 < 0.0]), dtype=np.float64))
    negative_fraction = negative / positive if positive > 0.0 else math.inf
    border = max(1, min(nx, ny) // 100)
    border_mask = np.ones((ny, nx), dtype=bool)
    border_mask[border:-border, border:-border] = False
    weights = fractional_circle_weights((ny, nx), centre_x, centre_y, radius_pixels)
    component_integrity = measure_fresh_component_integrity(
        sed=sed,
        image=image_raw,
        closure=closure,
        sigma_pixels=sigma_pixels,
        weights=weights,
        border_mask=border_mask,
    )
    del image_raw, sed_raw
    edge_fraction = (
        float(np.sum(np.clip(plane0, 0.0, None)[border_mask], dtype=np.float64)) / positive
        if positive > 0.0
        else math.inf
    )
    signed_convolved = gaussian_filter(
        plane0, sigma_pixels, mode="constant", cval=0.0, truncate=PSF_TRUNCATION_SIGMA
    )
    signed_convolved_total = float(np.sum(signed_convolved, dtype=np.float64))
    signed_aperture_140 = float(np.sum(signed_convolved * weights, dtype=np.float64))
    convolution_loss = abs(signed_convolved_total - map_total) / max(abs(map_total), np.finfo(float).tiny)
    if negative == 0.0:
        clipped_aperture_140 = signed_aperture_140
    else:
        clipped_convolved = gaussian_filter(
            np.clip(plane0, 0.0, None),
            sigma_pixels,
            mode="constant",
            cval=0.0,
            truncate=PSF_TRUNCATION_SIGMA,
        )
        clipped_aperture_140 = float(np.sum(clipped_convolved * weights, dtype=np.float64))
        del clipped_convolved
    dclip = symmetric_fraction(signed_aperture_140, clipped_aperture_140)
    checks = {
        "finite_positive_plane0_flux": math.isfinite(map_total) and map_total > 0.0 and positive > 0.0,
        "coeval_plane0_abs_relative_difference_le_1e-5": math.isfinite(coeval_relative) and abs(coeval_relative) <= COEVAL_TOLERANCE,
        "edge_positive_fraction_lt_1e-3": math.isfinite(edge_fraction) and 0.0 <= edge_fraction < EDGE_FRACTION_MAX,
        "convolution_flux_loss_lt_1e-4": math.isfinite(convolution_loss) and 0.0 <= convolution_loss < CONVOLUTION_LOSS_MAX,
        "finite_positive_signed_aperture": math.isfinite(signed_aperture_140) and signed_aperture_140 > 0.0,
        "direct_aperture_clip_symmetric_fraction_le_1e-3": math.isfinite(dclip) and 0.0 <= dclip <= APERTURE_CLIP_MAX,
        "aperture_plus_6sigma_support": half_field >= support_minimum,
    }
    if not all(checks.values()):
        failed = sorted(key for key, value in checks.items() if not value)
        raise ValueError(f"Task {task} canonical quality gate failed: {failed}")
    prediction_147 = signed_aperture_140 * DISTANCE_FACTOR
    assert_products_unchanged(paths, hashes, task)
    return {
        "fresh_task_index": task,
        "grid_model_index": int(row["grid_model_index"]),
        "model_label": row["model_label"],
        "anchor_ordinal": int(row["anchor_ordinal"]),
        "anchor_tag": row["anchor_tag"],
        "region": row["region"],
        "bin_index": int(row["bin_index"]),
        "wavelength_um": wavelength,
        "instrument": row["instrument"],
        "boundary_compromise": row["boundary_compromise"] == "true",
        "signed_aperture_flux_140pc_w_m2": signed_aperture_140,
        "signed_aperture_flux_147pc_w_m2": prediction_147,
        "zero_clipped_aperture_flux_140pc_w_m2": clipped_aperture_140,
        "aperture_clip_symmetric_fraction": dclip,
        "aperture_clip_evidence_kind": "direct_post_psf_nominal_1arcsec_symmetric_fraction",
        "full_plane_negative_fraction_diagnostic": negative_fraction,
        "plane0_edge_positive_flux_fraction": edge_fraction,
        "plane0_convolution_flux_loss_fraction": convolution_loss,
        "coeval_plane0_relative_difference": coeval_relative,
        "coeval_max_significant_relative_difference": float(
            closure["max_significant_custom_image_relative_difference"]
        ),
        "intrinsic_component_integrity": component_integrity,
        "target_pixel_arcsec": target_pixel_arcsec,
        "aperture_radius_pixels": radius_pixels,
        "psf_fwhm_arcsec": fwhm,
        "available_half_field_arcsec": half_field,
        "quality_checks": checks,
        "completion_sha256": hashes["completion"],
        "sed_rt_sha256": hashes["sed_rt"],
        "image_rt_sha256": hashes["image_rt"],
        "closure_sha256": hashes["closure"],
    }


def load_reused_measurements() -> list[dict[str, Any]]:
    predictions = Table.read(SOURCE_PREDICTIONS, format="ascii.ecsv")
    selected = predictions[np.asarray(predictions["model_index"], dtype=int) == 15]
    validation = read_json(SOURCE_FIT_VALIDATION, "source fit validation")
    v3 = read_json(SOURCE_V3_VALIDATION, "source v3 validation")
    qa = {
        int(row["task_index"]): row
        for row in validation.get("qa", [])
        if row.get("model_label") == GRID[0][1]
    }
    v3_tasks = {
        int(row["task_index"]): row
        for row in v3.get("tasks", [])
        if row.get("model_label") == GRID[0][1]
    }
    if not (
        len(selected) == 10
        and set(np.asarray(selected["task_index"], dtype=int)) == set(range(150, 160))
        and set(qa) == set(range(150, 160))
        and set(v3_tasks) == set(range(150, 160))
        and validation.get("valid") is True
    ):
        raise ValueError("Pinned upstream corner is not an exact validated ten-cell model")
    records: list[dict[str, Any]] = []
    selected.sort("anchor_ordinal")
    for ordinal, source in enumerate(selected):
        source_task = int(source["task_index"])
        evidence = qa[source_task]
        v3_task = v3_tasks[source_task]
        plane0 = v3_task["plane0_total_i"]
        component_integrity = validate_reused_component_integrity(v3_task)
        negative = float(evidence["full_plane_negative_fraction_diagnostic"])
        metric = float(evidence["aperture_clip_metric_or_upper_bound"])
        if not (
            int(source["anchor_ordinal"]) == ordinal
            and str(source["anchor_tag"]) == f"w{ordinal:03d}"
            and str(source["canonical_estimator"]) == "direct_polarized_total_i_plane0"
            and bool(source["recovery_flux_substituted"]) is False
            and evidence.get("uniform_quality_pass") is True
            and math.isclose(negative, 0.0, abs_tol=0.0)
            and math.isclose(metric, 0.0, abs_tol=0.0)
            and float(v3_task.get("coeval_max_relative_difference", math.inf)) <= COEVAL_TOLERANCE
            and math.isclose(
                float(source["model_aperture_flux_147pc_w_m2"]),
                float(plane0["prediction_total_i_plane0_aperture_flux_147pc_w_m2"]),
                rel_tol=1e-14,
                abs_tol=0.0,
            )
        ):
            raise ValueError(f"Reused source task {source_task} provenance/quality changed")
        records.append(
            {
                "fresh_task_index": -1,
                "upstream_task_index": source_task,
                "grid_model_index": 0,
                "model_label": GRID[0][1],
                "anchor_ordinal": ordinal,
                "anchor_tag": str(source["anchor_tag"]),
                "region": str(source["region"]),
                "bin_index": int(source["bin_index"]),
                "wavelength_um": float(source["wavelength_um"]),
                "instrument": str(source["instrument"]),
                "boundary_compromise": bool(source["boundary_compromise"]),
                "signed_aperture_flux_147pc_w_m2": float(source["model_aperture_flux_147pc_w_m2"]),
                "aperture_clip_symmetric_fraction": 0.0,
                "aperture_clip_evidence_kind": "exact_zero_from_upstream_no_negative_plane0_flux",
                "full_plane_negative_fraction_diagnostic": 0.0,
                "plane0_edge_positive_flux_fraction": float(plane0["plane0_edge_positive_flux_fraction"]),
                "plane0_convolution_flux_loss_fraction": float(plane0["plane0_convolution_flux_loss_fraction"]),
                "coeval_plane0_relative_difference": float(plane0["coeval_plane0_custom_relative_difference"]),
                "coeval_max_significant_relative_difference": float(v3_task["coeval_max_relative_difference"]),
                "intrinsic_component_integrity": component_integrity,
                "quality_checks": {"pinned_upstream_uniform_quality_pass": True},
                "source_predictions_sha256": SOURCE_PREDICTIONS_SHA256,
                "source_fit_validation_sha256": SOURCE_FIT_VALIDATION_SHA256,
                "source_v3_validation_sha256": SOURCE_V3_VALIDATION_SHA256,
            }
        )
    require_hash(SOURCE_PREDICTIONS, SOURCE_PREDICTIONS_SHA256, "reused predictions post-read")
    require_hash(SOURCE_FIT_VALIDATION, SOURCE_FIT_VALIDATION_SHA256, "reused fit validation post-read")
    require_hash(SOURCE_V3_VALIDATION, SOURCE_V3_VALIDATION_SHA256, "reused v3 validation post-read")
    return records


def component_integrity_extrema(
    measurements: list[dict[str, Any]],
) -> dict[str, Any]:
    """Summarize retained component gates with deterministic provenance witnesses."""
    component_entries: list[dict[str, Any]] = []
    task_entries: list[dict[str, Any]] = []
    for measurement in measurements:
        integrity = measurement["intrinsic_component_integrity"]
        grid_index = int(measurement["grid_model_index"])
        common = {
            "execution": GRID[grid_index][4],
            "grid_model_index": grid_index,
            "model_label": str(measurement["model_label"]),
            "fresh_task_index": int(measurement.get("fresh_task_index", -1)),
            "upstream_task_index": int(measurement.get("upstream_task_index", -1)),
            "anchor_ordinal": int(measurement["anchor_ordinal"]),
            "anchor_tag": str(measurement["anchor_tag"]),
            "wavelength_um": float(measurement["wavelength_um"]),
        }
        task_entries.append(
            {
                **common,
                "max_abs_significant_all_plane_coeval_relative_difference": float(
                    integrity[
                        "max_abs_significant_all_plane_coeval_relative_difference"
                    ]
                ),
            }
        )
        for name, component in integrity["components"].items():
            component_entries.append(
                {
                    **common,
                    "component": str(name),
                    "plane_zero_based": int(component["plane_zero_based"]),
                    "significant": bool(
                        component["significant_ge_1e-8_total_scale"]
                    ),
                    "active": bool(component["active_ge_1pct_component_l1"]),
                    "abs_coeval_relative_difference": abs(
                        float(component["coeval_sed_minus_map_over_map"])
                    ),
                    "negative_flux_fraction": float(
                        component["negative_flux_fraction"]
                    ),
                    "edge_positive_flux_fraction": float(
                        component["edge_positive_flux_fraction"]
                    ),
                    "convolution_flux_loss_fraction": float(
                        component["convolution_flux_loss_fraction"]
                    ),
                    "psf_map_aperture_fraction": float(
                        component["psf_map_aperture_fraction"]
                    ),
                }
            )

    identity_keys = (
        "execution",
        "grid_model_index",
        "model_label",
        "fresh_task_index",
        "upstream_task_index",
        "anchor_ordinal",
        "anchor_tag",
        "wavelength_um",
    )

    def witness(
        entries: list[dict[str, Any]],
        metric: str,
        *,
        largest: bool = True,
        threshold: float | None = None,
        threshold_relation: str | None = None,
    ) -> dict[str, Any]:
        if not entries:
            return {
                "value": None,
                "threshold": threshold,
                "threshold_relation": threshold_relation,
                "witness": None,
            }
        ordered = sorted(
            entries,
            key=lambda item: (
                -float(item[metric]) if largest else float(item[metric]),
                int(item["grid_model_index"]),
                int(item["anchor_ordinal"]),
                int(item.get("plane_zero_based", -1)),
            ),
        )
        selected = ordered[0]
        witness_keys = identity_keys + (
            ("component", "plane_zero_based")
            if "component" in selected
            else ()
        )
        return {
            "value": float(selected[metric]),
            "threshold": threshold,
            "threshold_relation": threshold_relation,
            "witness": {key: selected[key] for key in witness_keys},
        }

    def summarize(
        components: list[dict[str, Any]], tasks: list[dict[str, Any]]
    ) -> dict[str, Any]:
        significant = [item for item in components if item["significant"]]
        active = [item for item in components if item["active"]]
        return {
            "validated_task_count": len(tasks),
            "component_record_count": len(components),
            "significant_component_count": len(significant),
            "active_component_count": len(active),
            "max_abs_reported_all_significant_plane_coeval_relative_difference": witness(
                tasks,
                "max_abs_significant_all_plane_coeval_relative_difference",
                threshold=COEVAL_TOLERANCE,
                threshold_relation="<=",
            ),
            "max_abs_significant_component_coeval_relative_difference": witness(
                significant,
                "abs_coeval_relative_difference",
                threshold=COEVAL_TOLERANCE,
                threshold_relation="<=",
            ),
            "max_active_component_negative_flux_fraction": witness(
                active,
                "negative_flux_fraction",
                threshold=NEGATIVE_FRACTION_MAX,
                threshold_relation="<",
            ),
            "max_active_component_edge_positive_flux_fraction": witness(
                active,
                "edge_positive_flux_fraction",
                threshold=EDGE_FRACTION_MAX,
                threshold_relation="<",
            ),
            "max_active_component_convolution_flux_loss_fraction": witness(
                active,
                "convolution_flux_loss_fraction",
                threshold=CONVOLUTION_LOSS_MAX,
                threshold_relation="<",
            ),
            "min_all_component_psf_map_aperture_fraction": witness(
                components,
                "psf_map_aperture_fraction",
                largest=False,
                threshold=COMPONENT_APERTURE_FRACTION_MIN,
                threshold_relation=">=",
            ),
            "max_all_component_psf_map_aperture_fraction": witness(
                components,
                "psf_map_aperture_fraction",
                threshold=COMPONENT_APERTURE_FRACTION_MAX,
                threshold_relation="<=",
            ),
        }

    fresh_components = [item for item in component_entries if item["execution"] == "fresh"]
    reused_components = [
        item for item in component_entries if item["execution"] == "upstream_reuse"
    ]
    fresh_tasks = [item for item in task_entries if item["execution"] == "fresh"]
    reused_tasks = [item for item in task_entries if item["execution"] == "upstream_reuse"]
    return {
        "combined_90_cells": summarize(component_entries, task_entries),
        "fresh_80_cells": summarize(fresh_components, fresh_tasks),
        "pinned_upstream_reuse_10_cells": summarize(reused_components, reused_tasks),
    }


def prediction_record(
    measurement: dict[str, Any],
    observations: dict[int, dict[str, Any]],
) -> dict[str, Any]:
    grid_index = int(measurement["grid_model_index"])
    _, label, mass, amax, execution = GRID[grid_index]
    ordinal = int(measurement["anchor_ordinal"])
    observed = observations[int(measurement["bin_index"])]
    predicted = float(measurement["signed_aperture_flux_147pc_w_m2"])
    if not (math.isfinite(predicted) and predicted > 0.0):
        raise ValueError(f"Non-positive canonical prediction for grid model {grid_index}, anchor {ordinal}")
    residual = math.log10(observed["observed_flux_geometric_w_m2"] / predicted)
    component_integrity = measurement["intrinsic_component_integrity"]
    return {
        "cell_index": grid_index * 10 + ordinal,
        "fresh_task_index": int(measurement.get("fresh_task_index", -1)),
        "upstream_task_index": int(measurement.get("upstream_task_index", -1)),
        "grid_model_index": grid_index,
        "model_index": grid_index,
        "model_label": label,
        "mass_factor": mass,
        "envelope_amax_um": amax,
        "execution": execution,
        "anchor_ordinal": ordinal,
        "anchor_tag": str(measurement["anchor_tag"]),
        "region": str(measurement["region"]),
        "bin_index": int(measurement["bin_index"]),
        "wavelength_um": float(measurement["wavelength_um"]),
        "instrument": str(measurement["instrument"]),
        "boundary_compromise": bool(measurement["boundary_compromise"]),
        "n_detector_rows": int(observed["n_detector_rows"]),
        "segments": str(observed["segments"]),
        "observed_flux_geometric_w_m2": float(observed["observed_flux_geometric_w_m2"]),
        "observed_error_w_m2": float(observed["observed_error_w_m2"]),
        "observed_fractional_error": float(observed["observed_fractional_error"]),
        "observed_log10_sigma_dex": float(observed["observed_log10_sigma_dex"]),
        "model_aperture_flux_147pc_w_m2": predicted,
        "canonical_estimator": "signed_direct_coeval_method2_polarized_total_i_plane0",
        "observed_over_model": float(observed["observed_flux_geometric_w_m2"] / predicted),
        "log10_observed_over_model_dex": residual,
        "diagnostic_standardized_log_residual": residual / observed["observed_log10_sigma_dex"],
        "aperture_clip_symmetric_fraction": float(measurement["aperture_clip_symmetric_fraction"]),
        "aperture_clip_evidence_kind": str(measurement["aperture_clip_evidence_kind"]),
        "aperture_clip_threshold": APERTURE_CLIP_MAX,
        "full_plane_negative_fraction_diagnostic": float(measurement["full_plane_negative_fraction_diagnostic"]),
        "plane0_edge_positive_flux_fraction": float(measurement["plane0_edge_positive_flux_fraction"]),
        "plane0_convolution_flux_loss_fraction": float(measurement["plane0_convolution_flux_loss_fraction"]),
        "coeval_plane0_relative_difference": float(measurement["coeval_plane0_relative_difference"]),
        "coeval_max_significant_relative_difference": float(measurement["coeval_max_significant_relative_difference"]),
        "intrinsic_component_product_integrity_gates_applied": True,
        "significant_component_count": int(
            component_integrity["significant_component_count"]
        ),
        "active_component_count": int(component_integrity["active_component_count"]),
        "max_abs_significant_component_coeval_relative_difference": float(
            component_integrity[
                "max_abs_significant_component_coeval_relative_difference"
            ]
        ),
        "max_active_component_negative_flux_fraction": float(
            component_integrity["max_active_component_negative_flux_fraction"]
        ),
        "max_active_component_edge_positive_flux_fraction": float(
            component_integrity[
                "max_active_component_edge_positive_flux_fraction"
            ]
        ),
        "max_active_component_convolution_flux_loss_fraction": float(
            component_integrity[
                "max_active_component_convolution_flux_loss_fraction"
            ]
        ),
        "uniform_quality_pass": True,
        "signed_flux_used_for_scoring": True,
        "zero_clipped_flux_used_for_scoring": False,
        "empirical_mc_scatter_added": False,
        "formal_likelihood": False,
        "free_normalization_fitted": False,
        "screen_only": True,
        "operational_score_gap_tolerance_dex": OPERATIONAL_SCORE_MARGIN_DEX,
        "operational_tolerance_is_formal_uncertainty": False,
    }


def score_one_model(
    records: list[dict[str, Any]], scenario: str, boundary_weight: float
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if len(records) != 10:
        raise ValueError("A refinement model must contain exactly ten anchors")
    by_region: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        by_region[str(record["region"])].append(record)
    if {key: len(value) for key, value in by_region.items()} != REGION_COUNTS:
        raise ValueError("Six-region anchor allocation changed")
    region_rows: list[dict[str, Any]] = []
    score_squared = 0.0
    for region in REGIONS:
        selected = by_region[region]
        primary_count = REGION_COUNTS[region]
        terms = []
        effective_weights = []
        for record in selected:
            weight = boundary_weight if int(record["bin_index"]) in BOUNDARY_BINS else 1.0
            residual = float(record["log10_observed_over_model_dex"])
            terms.append(weight * residual**2)
            effective_weights.append(weight)
        contribution = math.fsum(terms) / primary_count
        effective_region_weight = math.fsum(effective_weights) / primary_count
        score_squared += contribution / len(REGIONS)
        region_rows.append(
            {
                "model_index": int(records[0]["model_index"]),
                "grid_model_index": int(records[0]["grid_model_index"]),
                "model_label": str(records[0]["model_label"]),
                "mass_factor": float(records[0]["mass_factor"]),
                "envelope_amax_um": float(records[0]["envelope_amax_um"]),
                "scenario": scenario,
                "region": region,
                "primary_anchor_denominator": primary_count,
                "effective_region_weight": effective_region_weight,
                "fixed_normalized_region_mean_square_dex2": contribution,
                "fixed_outer_region_weight": 1.0 / len(REGIONS),
                "renormalized_after_downweighting": False,
                "screen_only": True,
            }
        )
    score = math.sqrt(score_squared)
    residuals = np.asarray([float(row["log10_observed_over_model_dex"]) for row in records])
    summary = {
        "model_index": int(records[0]["model_index"]),
        "grid_model_index": int(records[0]["grid_model_index"]),
        "model_label": str(records[0]["model_label"]),
        "mass_factor": float(records[0]["mass_factor"]),
        "envelope_amax_um": float(records[0]["envelope_amax_um"]),
        "execution": str(records[0]["execution"]),
        "scenario": scenario,
        "boundary_anchor_weight": boundary_weight,
        "boundary_bins": "96;139",
        "n_physical_anchor_bins": 10,
        "n_equal_weight_regions": 6,
        "fixed_normalization_region_balanced_log_rms_dex": score,
        "fixed_normalization_mismatch_factor": 10.0**score,
        "primary_anchor_unweighted_log_rms_dex": float(np.sqrt(np.mean(residuals**2))),
        "median_observed_over_model": float(10.0 ** np.median(residuals)),
        "worst_absolute_log_residual_dex": float(np.max(np.abs(residuals))),
        "empirical_mc_scatter_added": False,
        "formal_likelihood": False,
        "free_normalization_fitted": False,
        "renormalized_after_downweighting": False,
        "diagnostic_rank": 0,
        "score_gap_to_scenario_rank1_dex": math.nan,
        "top_screen_diagnostic": False,
        "final_winner_claim_allowed": False,
        "screen_only": True,
    }
    return summary, region_rows


def score_predictions(
    predictions: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[int]]:
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in predictions:
        grouped[int(row["grid_model_index"])].append(row)
    if set(grouped) != set(range(9)) or any(len(rows) != 10 for rows in grouped.values()):
        raise ValueError("Predictions are not the exact nine-by-ten matrix")
    scores: list[dict[str, Any]] = []
    regions: list[dict[str, Any]] = []
    promoted: set[int] = set()
    for scenario, boundary_weight, _ in SCENARIOS:
        scenario_rows = []
        for model_index in range(9):
            summary, contribution = score_one_model(grouped[model_index], scenario, boundary_weight)
            scenario_rows.append(summary)
            regions.extend(contribution)
        scenario_rows.sort(
            key=lambda row: (
                float(row["fixed_normalization_region_balanced_log_rms_dex"]),
                int(row["grid_model_index"]),
            )
        )
        best = float(scenario_rows[0]["fixed_normalization_region_balanced_log_rms_dex"])
        for rank, row in enumerate(scenario_rows, 1):
            gap = float(row["fixed_normalization_region_balanced_log_rms_dex"]) - best
            unresolved = gap <= OPERATIONAL_SCORE_MARGIN_DEX
            row.update(
                diagnostic_rank=rank,
                score_gap_to_scenario_rank1_dex=gap,
                top_screen_diagnostic=rank == 1,
                operational_score_gap_tolerance_dex=OPERATIONAL_SCORE_MARGIN_DEX,
                unresolved_with_scenario_rank1=unresolved,
                promotion_candidate_in_this_scenario=unresolved,
                operational_tolerance_is_formal_uncertainty=False,
            )
            if unresolved:
                promoted.add(int(row["grid_model_index"]))
        scores.extend(scenario_rows)
    for row in scores:
        row["promotion_set_union_member"] = int(row["grid_model_index"]) in promoted
    return scores, regions, sorted(promoted)


def output_table(
    records: list[dict[str, Any]],
    product: str,
    manifest_sha: str,
    analysis_manifest_sha: str,
) -> Table:
    table = Table(rows=records)
    table.meta.update(
        {
            "product": product,
            "classification": "strict_continuum_refinement_v1_provisional_screen",
            "manifest_sha256": manifest_sha,
            "analysis_manifest_sha256": analysis_manifest_sha,
            "screen_only": True,
            "formal_likelihood": False,
            "production_final": False,
            "final_winner_claim_allowed": False,
            "canonical_estimator": "signed direct coeval method2 polarized total-I plane0",
            "mask": "strict Scenario B; ice and silicate features deferred",
            "aperture": "nominal 1arcsec radius; 64x64 fractional-pixel boundary quadrature",
            "zero_clipped_flux_used_for_scoring": False,
            "operational_tolerance_is_formal_uncertainty": False,
        }
    )
    return table


def make_plot(
    predictions: list[dict[str, Any]],
    scores: list[dict[str, Any]],
    promoted: list[int],
) -> None:
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in predictions:
        grouped[int(row["grid_model_index"])].append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: int(row["anchor_ordinal"]))
    primary = sorted(
        (row for row in scores if row["scenario"] == "primary"),
        key=lambda row: int(row["diagnostic_rank"]),
    )
    best_index = int(primary[0]["grid_model_index"])
    figure = plt.figure(figsize=(15.8, 8.8), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, width_ratios=(1.45, 1.0))
    sed_axis = figure.add_subplot(grid[0, 0])
    residual_axis = figure.add_subplot(grid[1, 0], sharex=sed_axis)
    heat_axis = figure.add_subplot(grid[0, 1])
    rank_axis = figure.add_subplot(grid[1, 1])

    for index, rows in grouped.items():
        wave = [row["wavelength_um"] for row in rows]
        flux = [row["model_aperture_flux_147pc_w_m2"] for row in rows]
        if index == best_index:
            continue
        colour = "#f4a261" if index in promoted else "0.82"
        sed_axis.plot(wave, flux, color=colour, lw=1.0, alpha=0.8)
    best = grouped[best_index]
    wave = np.asarray([row["wavelength_um"] for row in best])
    observed = np.asarray([row["observed_flux_geometric_w_m2"] for row in best])
    observed_error = np.asarray([row["observed_error_w_m2"] for row in best])
    best_flux = np.asarray([row["model_aperture_flux_147pc_w_m2"] for row in best])
    best_label = str(best[0]["model_label"])
    sed_axis.plot(wave, best_flux, "s-", color="#2166ac", lw=2.4, ms=5, label=f"nominal primary rank 1: {best_label}")
    sed_axis.errorbar(wave, observed, yerr=observed_error, fmt="o", color="black", capsize=2, label="strict continuum anchors")
    sed_axis.set(xscale="log", yscale="log", ylabel=r"$\lambda F_\lambda$ (W m$^{-2}$)")
    sed_axis.set_title("Absolute signed 1-arcsec-aperture SED; no fitted normalization")
    sed_axis.grid(alpha=0.2, which="both")
    sed_axis.legend(fontsize=8)

    for index in promoted:
        rows = grouped[index]
        residual_axis.plot(
            wave,
            [row["log10_observed_over_model_dex"] for row in rows],
            "o-",
            lw=2.2 if index == best_index else 1.2,
            ms=4,
            color="#2166ac" if index == best_index else "#f4a261",
            label=str(rows[0]["model_label"]),
        )
    residual_axis.axhline(0.0, color="black", lw=0.8)
    residual_axis.set(
        xscale="log",
        xlabel=r"Wavelength ($\mu$m)",
        ylabel=r"$\log_{10}(F_{obs}/F_{model})$ (dex)",
        title="Residuals for the operational promotion-set union",
    )
    residual_axis.grid(alpha=0.2, which="both")
    residual_axis.legend(fontsize=7, ncol=2)

    masses = [2.0, 2.25, 2.5]
    amax_values = [0.5, 0.6, 0.7]
    heat = np.full((3, 3), np.nan)
    ranks = np.zeros((3, 3), dtype=int)
    for row in primary:
        iy = masses.index(float(row["mass_factor"]))
        ix = amax_values.index(float(row["envelope_amax_um"]))
        heat[iy, ix] = float(row["fixed_normalization_region_balanced_log_rms_dex"])
        ranks[iy, ix] = int(row["diagnostic_rank"])
    image = heat_axis.imshow(heat, origin="lower", aspect="auto", cmap="viridis_r")
    figure.colorbar(image, ax=heat_axis, label="primary score (dex)")
    heat_axis.set_xticks(range(3), [f"{value:g}" for value in amax_values])
    heat_axis.set_yticks(range(3), [f"{value:g}" for value in masses])
    heat_axis.set(xlabel=r"Envelope $a_{max}$ ($\mu$m)", ylabel="Envelope mass factor", title="Primary score / diagnostic rank")
    midpoint = 0.5 * (float(np.nanmin(heat)) + float(np.nanmax(heat)))
    for iy in range(3):
        for ix in range(3):
            heat_axis.text(ix, iy, f"{heat[iy, ix]:.3f}\n#{ranks[iy, ix]}", ha="center", va="center", fontsize=8, color="white" if heat[iy, ix] > midpoint else "black")

    labels = [str(row["model_label"]) for row in reversed(primary)]
    y = np.arange(len(labels))
    lookup = {
        (str(row["scenario"]), str(row["model_label"])): float(row["fixed_normalization_region_balanced_log_rms_dex"])
        for row in scores
    }
    styles = {"primary": ("#2166ac", "o"), "boundary_half": ("#f4a261", "^"), "boundary_drop": ("0.35", "x")}
    for scenario, _, _ in SCENARIOS:
        colour, marker = styles[scenario]
        rank_axis.scatter([lookup[(scenario, label)] for label in labels], y, color=colour, marker=marker, s=32, label=scenario)
    rank_axis.set_yticks(y, labels, fontsize=7)
    rank_axis.set(xlabel="six-region fixed-normalization log RMS (dex)", title="Primary and mask-boundary sensitivity")
    rank_axis.grid(axis="x", alpha=0.2)
    rank_axis.legend(fontsize=7)
    figure.suptitle("TMC1A method-2 strict-continuum refinement v1 — provisional screen", fontsize=13)
    atomic_figure(figure, PLOT_OUTPUT)
    plt.close(figure)


def static_preflight(
    analysis_manifest_sha: str, manifest_sha_override: str | None = None
) -> None:
    analysis = load_analysis_manifest(analysis_manifest_sha)
    require_declared_outputs_absent()
    manifest_sha = str(analysis["input_manifest_sha256"])
    if manifest_sha_override is not None and manifest_sha_override != manifest_sha:
        raise ValueError("Explicit input-manifest SHA differs from the analysis manifest")
    manifest = load_manifest(manifest_sha)
    models, rows = validate_tables(manifest)
    reused = load_reused_measurements()
    observations = load_observations()
    if not (len(models) == 9 and len(rows) == 80 and len(reused) == 10 and len(observations) == 10):
        raise AssertionError("Static preflight count changed")
    print(
        "Refinement-v1 analyzer static preflight passed: exact 9-model grid, "
        "80 fresh tasks, ten pinned reused predictions, strict observations, "
        "and rank-blind analysis-manifest binding."
    )


def self_test() -> None:
    if not math.isclose(
        OPERATIONAL_SCORE_MARGIN_DEX,
        math.log10((1.0 + 0.025 / 2.0) / (1.0 - 0.025 / 2.0)),
        abs_tol=1e-15,
    ):
        raise AssertionError("Operational score margin changed")
    if not (
        math.isclose(symmetric_fraction(1.0, 1.001), 0.000999500249875035, rel_tol=1e-12)
        and symmetric_fraction(0.0, 0.0) == 0.0
    ):
        raise AssertionError("Dclip symmetric-fraction formula changed")
    weights = fractional_circle_weights((31, 31), 15.0, 15.0, 7.25)
    if not math.isclose(float(np.sum(weights)), math.pi * 7.25**2, rel_tol=5e-4):
        raise AssertionError("Fractional-circle aperture geometry changed")
    synthetic_sed = np.zeros(8, dtype=np.float64)
    synthetic_sed[0] = 1.0
    synthetic_sed[4:8] = (0.6, 0.3, 0.1, 0.0)
    synthetic_image = np.zeros((8, 1, 1, 31, 31), dtype=np.float64)
    synthetic_image[0, 0, 0, 15, 15] = 1.0
    for plane in COMPONENTS:
        synthetic_image[plane, 0, 0, 15, 15] = synthetic_sed[plane]
    synthetic_closure_rows = []
    for plane, label in CLOSURE_LABELS.items():
        value = float(synthetic_sed[plane])
        synthetic_closure_rows.append(
            {
                "plane": plane,
                "label": label,
                "significant": plane != 7,
                "image_flux": value,
                "custom_sed_flux": value,
                "custom_relative_difference": 0.0,
            }
        )
    synthetic_component_integrity = measure_fresh_component_integrity(
        sed=synthetic_sed,
        image=synthetic_image,
        closure={
            "planes": synthetic_closure_rows,
            "max_significant_custom_image_relative_difference": 0.0,
        },
        sigma_pixels=0.5,
        weights=weights,
        border_mask=np.pad(
            np.zeros((29, 29), dtype=bool), 1, constant_values=True
        ),
    )
    if not (
        synthetic_component_integrity["quality_pass"] is True
        and synthetic_component_integrity["significant_component_count"] == 3
        and synthetic_component_integrity["active_component_count"] == 3
        and synthetic_component_integrity["plane0_vs_sum_planes4to7_gate_applied"]
        is False
        and synthetic_component_integrity["aggregate_componentwise_flux_gate_applied"]
        is False
    ):
        raise AssertionError("Intrinsic separated-component integrity contract changed")
    synthetic: list[dict[str, Any]] = []
    for model_index in range(9):
        for ordinal, (bin_index, token, region, _, compromise) in enumerate(ANCHORS):
            synthetic.append(
                {
                    "grid_model_index": model_index,
                    "model_index": model_index,
                    "model_label": GRID[model_index][1],
                    "mass_factor": GRID[model_index][2],
                    "envelope_amax_um": GRID[model_index][3],
                    "execution": GRID[model_index][4],
                    "anchor_ordinal": ordinal,
                    "region": region,
                    "bin_index": bin_index,
                    "log10_observed_over_model_dex": 0.02 * model_index,
                }
            )
    scores, regions, promoted = score_predictions(synthetic)
    primary = sorted(
        (row for row in scores if row["scenario"] == "primary"),
        key=lambda row: int(row["diagnostic_rank"]),
    )
    if not (
        len(scores) == 27
        and len(regions) == 162
        and primary[0]["grid_model_index"] == 0
        and promoted == [0]
        and math.isclose(primary[1]["score_gap_to_scenario_rank1_dex"], 0.02, abs_tol=1e-15)
    ):
        raise AssertionError("Nine-model scoring/promotion contract changed")
    print("Refinement-v1 analyzer synthetic self-test passed.")


def run_analysis(
    analysis_manifest_sha: str, manifest_sha_override: str | None = None
) -> None:
    analysis = load_analysis_manifest(analysis_manifest_sha)
    require_declared_outputs_absent()
    manifest_sha = str(analysis["input_manifest_sha256"])
    if manifest_sha_override is not None and manifest_sha_override != manifest_sha:
        raise ValueError("Explicit input-manifest SHA differs from the analysis manifest")
    manifest = load_manifest(manifest_sha)
    _, task_rows = validate_tables(manifest)
    observations = load_observations()
    reused = load_reused_measurements()
    measurements = list(reused)
    fresh_qa: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    for row in task_rows:
        try:
            measurement = measure_fresh_task(row, manifest, manifest_sha)
            measurements.append(measurement)
            fresh_qa.append(measurement)
        except Exception as exc:
            errors.append(
                {
                    "task_index": int(row["task_index"]),
                    "grid_model_index": int(row["grid_model_index"]),
                    "model_label": row["model_label"],
                    "anchor_tag": row["anchor_tag"],
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
            )
    require_hash(
        ANALYSIS_MANIFEST,
        analysis_manifest_sha,
        "refinement analysis manifest post-measurement",
    )
    require_hash(
        ANALYZER,
        str(analysis["analyzer_sha256"]),
        "analysis-manifest-pinned analyzer post-measurement",
    )
    valid = not errors and len(fresh_qa) == 80 and len(measurements) == 90
    component_extrema = component_integrity_extrema(measurements)
    validation = {
        "schema_version": 1,
        "classification": (
            "strict_continuum_refinement_v1_uniform_quality_valid"
            if valid
            else "strict_continuum_refinement_v1_incomplete_or_invalid"
        ),
        "valid": valid,
        "screen_only": True,
        "formal_likelihood": False,
        "production_final": False,
        "final_winner_claim_allowed": False,
        "analysis_manifest_path": relative(ANALYSIS_MANIFEST),
        "analysis_manifest_sha256": analysis_manifest_sha,
        "analysis_manifest_classification": analysis["classification"],
        "analyzer_sha256": analysis["analyzer_sha256"],
        "manifest_path": relative(MANIFEST),
        "manifest_sha256": manifest_sha,
        "fresh_validated_task_count": len(fresh_qa),
        "fresh_expected_task_count": 80,
        "reused_validated_task_count": len(reused),
        "combined_validated_cell_count": len(measurements),
        "combined_expected_cell_count": 90,
        "canonical_estimator": "signed direct coeval method2 polarized total-I plane0 after Gaussian PSF and nominal 1arcsec aperture at 147pc",
        "coeval_relative_difference_max": COEVAL_TOLERANCE,
        "aperture_clip_symmetric_fraction_definition": "abs(F_signed-F_zero_clipped)/(0.5*(abs(F_signed)+abs(F_zero_clipped)))",
        "aperture_clip_symmetric_fraction_max": APERTURE_CLIP_MAX,
        "signed_flux_used_for_scoring": True,
        "zero_clipped_flux_used_for_scoring": False,
        "full_plane_negative_fraction_role": "diagnostic_only",
        "retained_intrinsic_separated_component_product_integrity_gates": True,
        "component_significant_fraction": COMPONENT_SIGNIFICANT_FRACTION,
        "component_active_fraction": COMPONENT_ACTIVE_FRACTION,
        "active_component_negative_fraction_max_exclusive": NEGATIVE_FRACTION_MAX,
        "active_component_edge_positive_fraction_max_exclusive": EDGE_FRACTION_MAX,
        "active_component_convolution_loss_fraction_max_exclusive": CONVOLUTION_LOSS_MAX,
        "component_aperture_fraction_range_inclusive": [
            COMPONENT_APERTURE_FRACTION_MIN,
            COMPONENT_APERTURE_FRACTION_MAX,
        ],
        "plane0_vs_sum_planes4to7_gate_applied": False,
        "aggregate_componentwise_flux_gate_applied": False,
        "component_integrity_extrema": component_extrema,
        "reused_source_predictions_path": relative(SOURCE_PREDICTIONS),
        "reused_source_predictions_sha256": SOURCE_PREDICTIONS_SHA256,
        "errors": errors,
        "fresh_qa": fresh_qa,
        "reused_qa": reused,
    }
    atomic_json(VALIDATION_OUTPUT, validation)
    if not valid:
        raise SystemExit(
            f"Refinement validation failed: {len(fresh_qa)}/80 fresh tasks valid; "
            f"{len(errors)} error(s). No scoring outputs were written."
        )

    predictions = [prediction_record(item, observations) for item in measurements]
    predictions.sort(key=lambda row: int(row["cell_index"]))
    if [int(row["cell_index"]) for row in predictions] != list(range(90)):
        raise ValueError("Combined prediction identities are not exactly 0..89")
    scores, regions, promoted = score_predictions(predictions)
    products = (
        (predictions, "strict_continuum_refinement_v1_predictions", PREDICTIONS_OUTPUT),
        (scores, "strict_continuum_refinement_v1_scores", SCORES_OUTPUT),
        (regions, "strict_continuum_refinement_v1_region_contributions", REGIONS_OUTPUT),
    )
    for records, product, path in products:
        atomic_table(
            output_table(records, product, manifest_sha, analysis_manifest_sha), path
        )
    make_plot(predictions, scores, promoted)

    top_by_scenario: dict[str, Any] = {}
    for scenario, weight, description in SCENARIOS:
        ordered = sorted(
            (row for row in scores if row["scenario"] == scenario),
            key=lambda row: int(row["diagnostic_rank"]),
        )
        top_by_scenario[scenario] = {
            "boundary_anchor_weight": weight,
            "description": description,
            "grid_model_index": int(ordered[0]["grid_model_index"]),
            "model_label": str(ordered[0]["model_label"]),
            "score_dex": float(ordered[0]["fixed_normalization_region_balanced_log_rms_dex"]),
            "nominal_rank_only": True,
        }
    primary = {
        int(row["grid_model_index"]): row
        for row in scores
        if row["scenario"] == "primary"
    }
    summary = {
        "schema_version": 1,
        "classification": "strict_continuum_refinement_v1_provisional_screen_complete",
        "valid": True,
        "screen_only": True,
        "formal_likelihood": False,
        "production_final": False,
        "final_winner_claim_allowed": False,
        "exact_grid": "9_models_x_10_strict_continuum_anchors",
        "execution": "1 pinned upstream model plus 8 fresh models / 80 fresh RT tasks",
        "analysis_manifest": {
            "path": relative(ANALYSIS_MANIFEST),
            "sha256": analysis_manifest_sha,
        },
        "manifest": {"path": relative(MANIFEST), "sha256": manifest_sha},
        "canonical_estimator": "signed direct coeval method2 polarized total-I plane0 aperture flux",
        "component_integrity_extrema": component_extrema,
        "top_by_scenario": top_by_scenario,
        "primary_nominal_rank1": top_by_scenario["primary"],
        "operational_score_gap_tolerance_dex": OPERATIONAL_SCORE_MARGIN_DEX,
        "operational_tolerance_is_formal_uncertainty": False,
        "promotion_model_indices": promoted,
        "promotion_models": [
            {
                "grid_model_index": index,
                "model_label": str(primary[index]["model_label"]),
                "mass_factor": float(primary[index]["mass_factor"]),
                "envelope_amax_um": float(primary[index]["envelope_amax_um"]),
                "primary_rank": int(primary[index]["diagnostic_rank"]),
                "primary_score_dex": float(primary[index]["fixed_normalization_region_balanced_log_rms_dex"]),
                "primary_gap_dex": float(primary[index]["score_gap_to_scenario_rank1_dex"]),
            }
            for index in promoted
        ],
        "required_confirmation": {
            "independent_seed_512k_confirmation_required_before_production_claim": True,
            "confirmation_seed": 41002,
            "confirmation_accepted_packets_per_anchor": 512000,
            "confirmation_anchor_scope": "all ten strict anchors for every model in the promotion-set union",
            "pointwise_symmetric_fraction_max": 0.025,
            "all_uniform_quality_gates_required": True,
        },
        "validation": {"path": relative(VALIDATION_OUTPUT), "sha256": sha256(VALIDATION_OUTPUT)},
        "outputs": {
            path.name: {"path": relative(path), "sha256": sha256(path)}
            for path in (PREDICTIONS_OUTPUT, SCORES_OUTPUT, REGIONS_OUTPUT, PLOT_OUTPUT)
        },
    }
    atomic_json(SUMMARY_OUTPUT, summary)
    print("Validated and ranked exact 9 x 10 strict-continuum refinement v1.")
    print(
        f"Nominal primary rank 1: {top_by_scenario['primary']['model_label']} "
        f"(score={top_by_scenario['primary']['score_dex']:.6f} dex)"
    )
    print(
        "Promotion set: "
        + ", ".join(item["model_label"] for item in summary["promotion_models"])
    )
    print("No final winner claim; independent seed-41002/512k confirmation remains required.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--self-test", action="store_true")
    mode.add_argument("--static-preflight", action="store_true")
    parser.add_argument(
        "--analysis-manifest-sha256",
        help="externally supplied SHA-256 of the frozen rank-blind analysis manifest",
    )
    parser.add_argument(
        "--manifest-sha256",
        help="optional redundant refinement input-manifest SHA consistency check",
    )
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if (
        not args.analysis_manifest_sha256
        or len(args.analysis_manifest_sha256) != 64
    ):
        raise ValueError(
            "Static/real analysis requires --analysis-manifest-sha256 HEX64"
        )
    if args.manifest_sha256 is not None and len(args.manifest_sha256) != 64:
        raise ValueError("--manifest-sha256 must be HEX64 when supplied")
    if args.static_preflight:
        static_preflight(args.analysis_manifest_sha256, args.manifest_sha256)
        return
    run_analysis(args.analysis_manifest_sha256, args.manifest_sha256)


if __name__ == "__main__":
    main()
