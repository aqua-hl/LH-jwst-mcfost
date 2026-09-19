#!/usr/bin/env python3
"""Prepare the 1-AU extinction/ice production run; never simulate or submit."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcfost_grid.configuration import atomic_json, canonical_hash, prepare_run, sha256
from production_observables import build_contract

NAME = "extinction_ice_production_1au_v1"
INCLINATIONS = (45., 50., 55., 60., 65., 70.)
MASSES = (.00015, .000225)
ICE_MASS_FRACTIONS = (0., .02, .04, .08)


def ice_volume_fraction(x):
    return (x / 1.2) / (x / 1.2 + (1-x) / 3.5)


def catalogue(parent):
    """Retain all 96 old points; add controlled ice/mass contrasts without fits."""
    from mcfost_grid.configuration import expand_models
    previous = expand_models(parent)
    if len(previous) != 96:
        raise ValueError("Expected the frozen 96-model v1.1 catalogue")
    records = []
    for model in previous:
        parameters = model["parameters"]
        f = parameters["envelope_ice_volume_fraction"]
        records.append(dict(index=len(records), role="parent_control", parent_index=model["index"],
                            parameters=parameters, ice_mass_fraction=f*1.2/(f*1.2+(1-f)*3.5),
                            contrast_group_id=None, dust_family="generic_ice_coated_Mie"))
    base = previous[0]["parameters"]
    for inclination, mass in itertools.product(INCLINATIONS, MASSES):
        group = f"i{inclination:g}_m{mass:.6f}"
        for x in ICE_MASS_FRACTIONS:
            records.append(dict(index=len(records), role="ice_ladder", parent_index=None,
                                parameters={**base, "inclination_deg": inclination,
                                            "envelope_dust_mass_msun": mass,
                                            "envelope_ice_volume_fraction": ice_volume_fraction(x)},
                                ice_mass_fraction=x, contrast_group_id=group,
                                dust_family="bare_silicate_Mie" if x == 0 else "generic_ice_coated_Mie"))
    if len(records) != 144 or len({canonical_hash(r["parameters"]) for r in records}) != 144:
        raise ValueError("Duplicate or unexpected production models")
    return records


def build_package(destination):
    from astropy.table import Table
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(f"Refusing to replace an existing run: {destination}")
    parent_path = ROOT / "config/grid.continuum-production-v1.1.json"
    parent = json.loads(parent_path.read_text())
    records = catalogue(parent)
    contract = build_contract(ROOT)
    sources = [parent_path, ROOT / "reference/parameters/continuum_original_dust.para",
               ROOT / "reference/observations/continuum_sed_R100.ecsv",
               ROOT / "reference/observations/tmc1a_sed_unstitched.ecsv",
               ROOT / "reference/observations/continuum_fit_mask_v1_manifest.json",
               ROOT / "docs/EXTINCTION_ICE_PRODUCTION_1AU_V1.md",
               ROOT / "requirements-cluster-recovery.txt", Path(__file__).resolve()]
    copies = {"extinction_ice_production_task.py": "production_task.py",
              "analyze_extinction_ice_production.py": "analyze_extinction_ice_production.py",
              "configure_extinction_ice_numerics_v2.py": "configure_cluster.py",
              "production_observables.py": "production_observables.py"}
    sources += [ROOT / "scripts" / name for name in copies]
    sources += list((ROOT / "src/mcfost_grid").glob("*.py")) + [ROOT / "workflow.py"]
    before = {str(p.relative_to(ROOT)): sha256(p) for p in sources}
    # Prepare with temporary source tables; all executable inputs are copied
    # into the run and pinned by the final manifest before any launch script.
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="mcfost-production-inputs-", dir=destination.parent) as tmp:
        temporary = Path(tmp)
        anchor_path = temporary / "probes.ecsv"
        Table(rows=contract["probes"]).write(anchor_path, format="ascii.ecsv")
        config = dict(schema_version=1, run_name=destination.name, output_dir=str(destination.parent),
                      template=str(ROOT / "reference/parameters/continuum_original_dust.para"),
                      numerical_only=True, smoke_test=False, max_models=144,
                      models=[r["parameters"] for r in records],
                      numerics={**parent["numerics"], "photons_temperature": 2048000,
                                "photons_image": 2048000, "photons_sed": 2048000,
                                "random_seed": 43001, "image_npix": 6001, "image_size_au": 6000.},
                      observations=dict(anchors_file=str(anchor_path),
                                        spectrum_file=str(ROOT / "reference/observations/continuum_sed_R100.ecsv"),
                                        aperture_radius_arcsec=1., target_distance_pc=147.,
                                        aperture_subpixels=64, quality_policy="aperture_v2"),
                      notes="Fresh equilibrium for every model including each ice composition. Numerical-only probes are unscored by the generic analyzer; the dedicated analyzer integrates them into the frozen 19-band observation contract. Fixed approximately 1 AU/pixel; no dependency on the final diagnostic passing. Single seed per physical model, no per-model convergence claim.")
        config_path = temporary / "grid.json"
        atomic_json(config_path, config)
        run = prepare_run(config_path)
        shutil.copy2(anchor_path, run / "inputs/production_probes.ecsv")
    for source, target in copies.items():
        shutil.copy2(ROOT / "scripts" / source, run / "code" / target)
    shutil.copy2(ROOT / "docs/EXTINCTION_ICE_PRODUCTION_1AU_V1.md", run / "README.md")
    shutil.copy2(ROOT / "requirements-cluster-recovery.txt", run / "requirements-cluster.txt")
    shutil.copy2(parent_path, run / "inputs/parent_grid_reference.json")
    atomic_json(run / "inputs/production_observation_contract.json", contract)
    manifest = json.loads((run / "manifest.json").read_text())
    # Remove now-dead preparation paths from the provenance copy. Runtime uses
    # frozen .para files, but this also leaves the configuration understandable
    # after the bundle moves to a different host.
    frozen_config = manifest["configuration"]
    frozen_config.update(output_dir="../..", template="template.para")
    frozen_config["observations"].update(anchors_file="production_probes.ecsv",
                                       spectrum_file="observed_spectrum.ecsv")
    atomic_json(run / "inputs/configuration.json", frozen_config)
    tasks = [dict(index=m["index"], model_id=m["id"]) for m in manifest["models"]]
    experiment = dict(schema_version=1, experiment_id=NAME, catalogue=records, tasks=tasks,
                      model_count=144, temperature_solves=144,
                      image_requests=144*len(contract["probes"]), probes_per_model=len(contract["probes"]),
                      observable_count=len(contract["bands"]),
                      selection_rule="All 96 previous physical combinations, plus 12 fixed inclination/mass backgrounds crossed with four ice mass fractions; no fit/rank selection",
                      resources=dict(cpus_per_task=64, max_concurrent_tasks=16, peak_cpus=1024,
                                     slurm_memory_per_task="160G", mcfost_memory_gb=112),
                      final_diagnostic=dict(run="final_resolution_1au_v1", dependency=False,
                                            investigation_policy="One final comparison; record outcome and close investigation. Production is authorized independently at approximately 1 AU/pixel, without claiming convergence passed."),
                      numerics=frozen_config["numerics"],
                      block_weights=contract["block_weights"],
                      dust=dict(core="Draine_Si_sUV.dat", mantle="ice_opct.dat", core_density_g_cm3=3.5,
                                mantle_density_g_cm3=1.2, ice_mass_fractions=list(ICE_MASS_FRACTIONS),
                                mass_convention="Total envelope dust mass including ice; outer grain size distribution held fixed in each matched ladder",
                                bare_limit="Exactly zero ice uses one pure-silicate Mie component, not a zero-thickness coated species",
                                H2O30K_used=False, material_extension_validated=False,
                                temperature_dependent_ice_survival=False),
                      quadrature_reference_fraction=.01, production_convergence_certified=False,
                      formal_likelihood=False, posterior_inference=False,
                      limits=["Conditional dust-family sensitivity and degeneracy experiment, not proof that the spectrum uniquely constrains all dust/ice properties.",
                              "Broad 1-arcsec H2O flux windows are distinct from the archived 0.35-arcsec R2400 optical-depth operator.",
                              "MIRI covariance and all common-gain scales are sensitivity assumptions; within-window structure is not a measured artifact amplitude.",
                              "Five-node versus three-node integration is an embedded diagnostic, not a bound on every quadrature error.",
                              "One production seed per model; the last numerical check samples one fixed temperature and wavelength only.",
                              "True bare controls are included, but no infinitesimal-mantle opacity-equivalence or physical ice-survival certification is inferred."],
                      source_sha256=before)
    for filename in ("experiment.json", "production_experiment.json"):
        atomic_json(run / filename, experiment)
    (run / "experiment.sha256").write_text(sha256(run / "experiment.json") + "  experiment.json\n")
    for p in run.rglob("*"):
        if p.is_file() and p.name not in {"manifest.json", "manifest.sha256", "machine.template.json"}:
            manifest["input_hashes"][str(p.relative_to(run))] = sha256(p)
    atomic_json(run / "manifest.json", manifest)
    (run / "manifest.sha256").write_text(sha256(run / "manifest.json") + "  manifest.json\n")
    machine = dict(schema_version=1, mcfost_executable="mcfost", mcfost_utils="$MCFOST_UTILS",
                   backend="image_method2", threads=64, timeout_seconds=14400, max_memory_gb=112,
                   slurm=dict(time="24:00:00", memory="160G", max_parallel=16,
                              python="/SET/BY/CONFIGURE_CLUSTER", setup_lines=[],
                              analysis_cpus=2, analysis_memory="8G", analysis_time="02:00:00"),
                   notes="Approximately 1 AU/pixel: 6001 pixels across6000AU. Four-hour per-command and24-hour per-model limits; cached stages resume after interruption. 64CPUs/160GB per model;16concurrent. Automatic integrated-band/contrast analysis runs afterany, and records incomplete results honestly. Use the shared mcfost-v11 environment to configure before submission.")
    atomic_json(run / "machine.template.json", machine)
    if before != {p: sha256(ROOT / p) for p in before}:
        raise RuntimeError("Source changed during preparation; do not launch this run")
    return dict(bundle=str(run), models=144, images=experiment["image_requests"],
                bands=len(contract["bands"]), probes_per_model=len(contract["probes"]),
                temperature_solves=144, mcfost_invoked=False, jobs_submitted=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "runs" / NAME)
    args = parser.parse_args()
    print(json.dumps(build_package(args.output), indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
