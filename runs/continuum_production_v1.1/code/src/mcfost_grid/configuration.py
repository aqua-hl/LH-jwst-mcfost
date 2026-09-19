"""Freeze a readable experiment configuration into an explicit task catalogue."""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .physics import PHYSICAL_KEYS, render_parameter, physical_args

WORKSPACE = Path(__file__).resolve().parents[2]
DEFAULT_PARAMETERS = {
    "envelope_dust_mass_msun": 0.000225,
    "envelope_amax_um": 0.4,
    "envelope_size_exponent": 2.75,
    "envelope_ice_volume_fraction": 0.05,
    "cavity_half_opening_deg": 17.5,
    "inclination_deg": 70.0,
    "distance_pc": 140.0,
    "stellar_temperature_k": 4000.0,
    "stellar_radius_rsun": 2.5,
    "stellar_mass_msun": 0.5,
    "accretion_rate_msun_per_year": 7.3011756855e-8,
}
DEFAULT_NUMERICS = {
    "photons_temperature": 128000, "photons_image": 128000,
    "photons_sed": 128000, "image_npix": 2401, "image_size_au": 6000.0,
}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def load_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON: {x}")))


def resolve(base, value):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (Path(base) / path).resolve()


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _strict_keys(obj, allowed, name):
    if not isinstance(obj, dict):
        raise ValueError(f"{name} must be an object")
    extra = set(obj) - set(allowed)
    if extra:
        raise ValueError(f"Unknown {name} keys: {sorted(extra)}")


def expand_models(config):
    limit = config.get("max_models", 10000)
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("max_models must be a positive integer")
    fixed = {**DEFAULT_PARAMETERS, **config.get("fixed", {})}
    _strict_keys(fixed, PHYSICAL_KEYS, "physical parameter")
    if "grid" in config and "models" in config:
        raise ValueError("Use either grid (Cartesian product) or models (explicit overrides), not both")
    if "models" in config:
        overrides = config["models"]
        if not isinstance(overrides, list) or not overrides:
            raise ValueError("models must be a nonempty list of parameter objects")
    else:
        axes = config.get("grid", {})
        _strict_keys(axes, PHYSICAL_KEYS, "grid")
        for key, values in axes.items():
            if not isinstance(values, list) or not values or len({canonical_hash(x) for x in values}) != len(values):
                raise ValueError(f"grid.{key} must be a nonempty list without duplicate values")
        count = math.prod(len(values) for values in axes.values())
        if count > config.get("max_models", 10000):
            raise ValueError(f"Grid has {count} models; exceeds max_models")
        overrides = [dict(zip(axes, values)) for values in itertools.product(*axes.values())]
    if len(overrides) > config.get("max_models", 10000):
        raise ValueError("Too many models")
    models, seen = [], set()
    for i, override in enumerate(overrides):
        _strict_keys(override, PHYSICAL_KEYS, "model")
        parameters = {**fixed, **override}
        physical_args(parameters)  # Validate before numeric canonicalization.
        parameters = {key: float(value) for key, value in parameters.items()}
        key = canonical_hash(parameters)
        if key in seen:
            raise ValueError("Duplicate physical model in catalogue")
        seen.add(key)
        models.append({"index": i, "id": f"m{i:04d}_{key[:8]}", "parameters": parameters})
    return models


def read_anchors(config, base):
    from astropy.table import Table
    selection = config.get("observations", {})
    _strict_keys(selection, {"preset", "anchor_ids", "anchors_file", "spectrum_file",
                            "aperture_radius_arcsec", "target_distance_pc", "aperture_subpixels"}, "observations")
    if selection.get("preset", "strict9") != "strict9" and "anchors_file" not in selection:
        raise ValueError("Unknown preset; use strict9 or an explicit anchors_file")
    spectrum = resolve(base, selection["spectrum_file"]) if "spectrum_file" in selection else WORKSPACE / "reference/observations/continuum_sed_R100.ecsv"
    if "anchors_file" in selection:
        if "spectrum_file" not in selection:
            raise ValueError("Custom anchors require an explicit matching spectrum_file; do not silently reuse the 1-arcsec grey spectrum")
        path = resolve(base, selection["anchors_file"])
        table = Table.read(path, format="ascii.ecsv" if path.suffix == ".ecsv" else "ascii.csv")
        required = {"id", "wavelength_um", "instrument", "region", "flux_jy", "uncertainty_jy"}
        if not required <= set(table.colnames):
            raise ValueError(f"Anchor table requires {sorted(required)}")
        anchors = []
        for row in table:
            item = {key: (row[key].item() if hasattr(row[key], "item") else row[key]) for key in table.colnames}
            if "score" in item and isinstance(item["score"], str):
                if item["score"].lower() not in {"true", "false"}:
                    raise ValueError("CSV score must be true or false")
                item["score"] = item["score"].lower() == "true"
            item.setdefault("score", True)
            anchors.append(item)
    else:
        table = Table.read(spectrum)
        source = WORKSPACE / "reference/provenance/strict_continuum_cavity_mass_q_v1/anchors.tsv"
        anchors = []
        with source.open() as stream:
            for old in csv.DictReader(stream, delimiter="\t"):
                candidates = table[(table["bin_index"] == int(old["bin_index"])) &
                                   (table["instrument"] == old["instrument"])]
                if not len(candidates):
                    raise ValueError(f"Missing observed anchor {old['anchor_tag']}")
                # Detector overlaps are the same physical anchor, not independent
                # repeated constraints. Preserve the original geometric mean and
                # common (not divided down) calibration floor.
                fluxes = [float(row["flux_jy"]) for row in candidates]
                floors = [float(row["systematic_fraction"]) for row in candidates]
                errors = [float(row["fit_error_jy"]) for row in candidates]
                for value in fluxes + floors + errors:
                    _positive(value, "observed flux/error/floor")
                if any(abs(f - floors[0]) > 1e-12 for f in floors):
                    raise ValueError("Different calibration floors in an overlap")
                if any(abs(float(row["wavelength_um"]) - float(old["wavelength_um"])) > 1e-13 for row in candidates):
                    raise ValueError("Observed anchor wavelength changed")
                combined_flux = 10 ** (sum(math.log10(f) for f in fluxes) / len(fluxes))
                fraction = math.sqrt(floors[0]**2 + sum(max((e/f)**2 - s**2, 0.0)
                                      for f, e, s in zip(fluxes, errors, floors)) / len(fluxes)**2)
                anchors.append({"id": old["anchor_tag"], "wavelength_um": float(old["wavelength_um"]),
                                "instrument": old["instrument"], "region": old["region"],
                                "flux_jy": combined_flux, "uncertainty_jy": combined_flux * fraction,
                                "n_detector_rows": len(candidates),
                                "segments": ";".join(sorted(str(row["segment"]) for row in candidates)),
                                "score": True, "image_npix": int(old["nx"]),
                                "image_size_au": float(old["expected_fov_au"]),
                                "psf_fwhm_arcsec": float(old["psf_fwhm_arcsec"])})
        # These data are a specific 1-arcsec-radius extraction, not interchangeable apertures.
        if selection.get("aperture_radius_arcsec", 1.0) != 1.0:
            raise ValueError("strict9 observations require aperture_radius_arcsec=1.0; supply new data for another aperture")
    if "anchor_ids" in selection:
        ids = selection["anchor_ids"]
        if not isinstance(ids, list) or not ids or len(ids) != len(set(ids)):
            raise ValueError("anchor_ids must be a nonempty unique list")
        if set(ids) - {a["id"] for a in anchors}:
            raise ValueError("Unknown anchor_ids")
        anchors = [a for a in anchors if a["id"] in ids]
    seen = set()
    for item in anchors:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", str(item["id"])) or item["id"] in seen:
            raise ValueError("Anchor IDs must be unique safe filename tokens")
        seen.add(item["id"])
        for key in ("wavelength_um", "flux_jy", "uncertainty_jy"):
            _positive(item[key], key)
        if item["instrument"] not in ("NIRSpec", "MIRI") or not isinstance(item["score"], bool):
            raise ValueError("Anchor instrument/score invalid")
    if not anchors or not any(a["score"] for a in anchors):
        raise ValueError("Need at least one scored anchor")
    measurement = {"aperture_radius_arcsec": selection.get("aperture_radius_arcsec", 1.0),
                   "target_distance_pc": selection.get("target_distance_pc", 147.0),
                   "aperture_subpixels": selection.get("aperture_subpixels", 64)}
    for key, value in measurement.items():
        _positive(value, key)
    if not isinstance(measurement["aperture_subpixels"], int) or measurement["aperture_subpixels"] > 128:
        raise ValueError("aperture_subpixels must be integer <=128")
    return anchors, spectrum, measurement


def numerics_for_anchor(config, anchor):
    result = {**DEFAULT_NUMERICS, **config.get("numerics", {})}
    # Production preset uses the old wavelength-dependent map geometry. Explicit
    # numerical fields override it, which is useful for the coarse laptop smoke.
    for key in ("image_npix", "image_size_au"):
        if key not in config.get("numerics", {}) and key in anchor:
            result[key] = anchor[key]
    return result


def prepare_run(config_path):
    config_path = Path(config_path).resolve()
    config = load_json(config_path)
    _strict_keys(config, {"schema_version", "run_name", "output_dir", "template", "fixed", "grid", "models",
                         "numerics", "observations", "smoke_test", "max_models", "notes"}, "configuration")
    if config.get("schema_version") != 1:
        raise ValueError("schema_version must be 1")
    name = config.get("run_name", "")
    if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_-]+(?:\.[A-Za-z0-9_-]+)*", name):
        raise ValueError("run_name must be a safe nonempty filename token; dots may only separate nonempty components")
    if not isinstance(config.get("smoke_test", False), bool):
        raise ValueError("smoke_test must be boolean")
    base = config_path.parent
    run = resolve(base, config.get("output_dir", "../runs")) / name
    if run.exists():
        raise FileExistsError(f"Run already exists: {run}. Resume it, or choose a new run_name.")
    template = resolve(base, config["template"]) if "template" in config else WORKSPACE / "reference/parameters/continuum_nominal.para"
    template_text = template.read_text()
    models = expand_models(config)
    anchors, spectrum, measurement = read_anchors(config, base)
    rendered = {}
    for model in models:
        p = model["parameters"]
        physical_args(p)
        rendered[model["id"]] = {"temperature": render_parameter(template_text, p, numerics_for_anchor(config, {}), "temperature")}
        for anchor in anchors:
            n = numerics_for_anchor(config, anchor)
            rendered[model["id"]][anchor["id"]] = {
                "image": render_parameter(template_text, p, n, "image", anchor["wavelength_um"]),
                "coeval": render_parameter(template_text, p, n, "coeval", anchor["wavelength_um"]),
            }
    # All configuration validation precedes directory creation.
    run.mkdir(parents=True)
    inputs = run / "inputs"
    inputs.mkdir()
    shutil.copy2(template, inputs / "template.para")
    shutil.copy2(spectrum, inputs / "observed_spectrum.ecsv")
    atomic_json(inputs / "configuration.json", config)
    # Freeze the small active code and needed dust inputs, not the legacy tree.
    shutil.copytree(WORKSPACE / "src/mcfost_grid", run / "code/src/mcfost_grid",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copy2(WORKSPACE / "workflow.py", run / "code/workflow.py")
    dust_dir = inputs / "utils/Dust"
    dust_dir.mkdir(parents=True)
    for src in (WORKSPACE / "reference/dust/Draine_Si_sUV.dat",
                WORKSPACE / "reference/runtime_assets/Dust/ac_opct.dat",
                WORKSPACE / "reference/runtime_assets/Dust/ice_opct.dat"):
        shutil.copy2(src, dust_dir / src.name)
    stars = inputs / "utils/Stellar_Spectra"
    stars.mkdir()
    shutil.copy2(WORKSPACE / "reference/runtime_assets/Stellar_Spectra/lte4000-3.5.NextGen.fits.gz", stars)
    for model in models:
        directory = run / "models" / model["id"]
        directory.mkdir(parents=True)
        (directory / "temperature.para").write_text(rendered[model["id"]]["temperature"])
        for anchor in anchors:
            sub = directory / "anchors" / anchor["id"]
            sub.mkdir(parents=True)
            for stage in ("image", "coeval"):
                (sub / f"{stage}.para").write_text(rendered[model["id"]][anchor["id"]][stage])
            (sub / "wavelength.lambda").write_text(f"1\n{anchor['wavelength_um']:.17g}\n")
    files = [p for p in run.rglob("*") if p.is_file()]
    manifest = {"schema_version": 1, "run_id": name, "created_utc": now(), "configuration": config,
                "models": models, "anchors": anchors, "measurement": measurement,
                "observation_spectrum": "inputs/observed_spectrum.ecsv", "model_directory": "models/{id}",
                "input_hashes": {str(p.relative_to(run)): sha256(p) for p in files},
                "estimator": "signed_direct_method2_image_total_I", "formal_likelihood": False,
                "temperature_policy": "fresh equilibrium for every physical model"}
    atomic_json(run / "manifest.json", manifest)
    return run
