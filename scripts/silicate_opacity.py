#!/usr/bin/env python3
"""Grain opacities from supplied optical constants, for the silicate workspace.

Computes size-distribution-averaged absorption and scattering opacities from an
optical-constant table, so a candidate silicate can be screened in opacity space
before any radiative-transfer run. Mie theory (Bohren & Huffman) for compact
spheres, with optional vacuum inclusion through the Bruggeman effective medium.

The quantity the screen turns on is

    S = kappa_abs(9.7 um) / kappa_ext(2.2 um)

the 9.7-um opacity per unit near-infrared extinction. S compares optical depths
at fixed NIR extinction, requiring a mass adjustment when kappa_ext changes.
It is not a fixed-mass or emergent-aperture-flux prediction. The tabulated
flux_gain_dex_if_tau_* values are pure-screen illustrations with that fixed-NIR
normalization, not forecasts for the fresh-temperature DHS pilot.

This module computes opacities. It does not invent optical constants: every
material must come from a table on disk, recorded in the workspace manifest.

    python -B scripts/silicate_opacity.py --help
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONSTANTS = ROOT / "reference/dust/Draine_Si_sUV.dat"

# Production envelope size distribution, from continuum_original_dust.para.
AMIN_UM, AMAX_UM, QEXP, NGRAINS = 0.03, 0.4, 2.75, 50
NEAR_IR_UM, SILICATE_UM = 2.2, 9.7


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_constants(path: Path) -> tuple[float, np.ndarray]:
    """Return (bulk density g/cm3, array of [lambda_um, n, k]) from an MCFOST table."""
    rows: list[list[float]] = []
    density: float | None = None
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        if len(f) == 2 and density is None:
            density = float(f[0])
            continue
        if len(f) >= 3:
            try:
                rows.append([float(x) for x in f[:3]])
            except ValueError:
                continue
    if density is None or not rows:
        raise SystemExit(f"Could not parse optical constants from {path}")
    table = np.asarray(rows)
    return density, table[np.argsort(table[:, 0])]


def mie(m: complex, x: float) -> tuple[float, float, float]:
    """Qext, Qsca and the asymmetry parameter g for a compact sphere (BHMIE)."""
    if x < 1e-8:
        return 0.0, 0.0, 0.0
    y = m * x
    nmax = int(x + 4.0 * x ** (1.0 / 3.0) + 2.0)
    nmx = max(nmax, int(abs(y))) + 15

    # Downward recurrence for the logarithmic derivative of the Riccati-Bessel function.
    d = np.zeros(nmx + 1, dtype=complex)
    for n in range(nmx, 0, -1):
        d[n - 1] = n / y - 1.0 / (d[n] + n / y)

    psi0, psi1 = math.cos(x), math.sin(x)
    chi0, chi1 = -math.sin(x), math.cos(x)
    xi1 = complex(psi1, -chi1)

    qext = qsca = gsca = 0.0
    an1 = bn1 = 0j
    for n in range(1, nmax + 1):
        fn = (2 * n + 1) / (n * (n + 1))
        psi = (2 * n - 1) * psi1 / x - psi0
        chi = (2 * n - 1) * chi1 / x - chi0
        xi = complex(psi, -chi)

        da = d[n] / m + n / x
        db = d[n] * m + n / x
        an = (da * psi - psi1) / (da * xi - xi1)
        bn = (db * psi - psi1) / (db * xi - xi1)

        qext += (2 * n + 1) * (an.real + bn.real)
        qsca += (2 * n + 1) * (abs(an) ** 2 + abs(bn) ** 2)
        gsca += fn * (an * bn.conjugate()).real
        if n > 1:
            gsca += ((n - 1) * (n + 1) / n) * (
                (an1 * an.conjugate()).real + (bn1 * bn.conjugate()).real)

        an1, bn1 = an, bn
        psi0, psi1 = psi1, psi
        chi0, chi1 = chi1, chi
        xi1 = complex(psi1, -chi1)

    qext *= 2.0 / x ** 2
    qsca *= 2.0 / x ** 2
    gsca *= 4.0 / (x ** 2 * qsca) if qsca > 0 else 0.0
    return qext, qsca, gsca


def bruggeman(m_solid: complex, porosity: float) -> complex:
    """Effective index of a solid/vacuum mixture (Bruggeman rule)."""
    if porosity <= 0:
        return m_solid
    e1, e2, f = m_solid ** 2, complex(1.0, 0.0), porosity
    # Solve (1-f)(e1-e)/(e1+2e) + f(e2-e)/(e2+2e) = 0 for e.
    a = -2.0
    b = (3 * (1 - f) - 1) * e1 + (3 * f - 1) * e2
    c = e1 * e2
    disc = (b ** 2 - 4 * a * c) ** 0.5
    for root in ((-b + disc) / (2 * a), (-b - disc) / (2 * a)):
        if root.real > 0 and root.imag >= 0:
            return root ** 0.5
    raise SystemExit("Bruggeman mixing produced no physical root")


def size_grid(amin: float, amax: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Log-spaced radii (um) and their dn/da * da weights for a power law."""
    edges = np.logspace(math.log10(amin), math.log10(amax), n + 1)
    radii = np.sqrt(edges[:-1] * edges[1:])
    weights = radii ** (-QEXP) * np.diff(edges)
    return radii, weights


def opacities(constants: Path, wavelengths: np.ndarray, *, amax_um: float = AMAX_UM,
              porosity: float = 0.0, amin_um: float = AMIN_UM) -> dict:
    """Mass opacities per gram of grain material, averaged over the size distribution."""
    density, table = load_constants(constants)
    radii, weights = size_grid(amin_um, amax_um, NGRAINS)
    # Mass per unit dn/da weight; porosity lowers the mean density of a grain.
    solid_density = density * (1.0 - porosity)
    mass = weights * (4.0 / 3.0) * math.pi * (radii * 1e-4) ** 3 * solid_density
    total_mass = mass.sum()

    kabs = np.zeros_like(wavelengths)
    ksca = np.zeros_like(wavelengths)
    gbar = np.zeros_like(wavelengths)
    for i, lam in enumerate(wavelengths):
        n = float(np.interp(lam, table[:, 0], table[:, 1]))
        k = float(np.interp(lam, table[:, 0], table[:, 2]))
        m = bruggeman(complex(n, k), porosity)
        area_abs = area_sca = area_g = 0.0
        for a, w in zip(radii, weights):
            x = 2.0 * math.pi * a / lam
            qext, qsca, g = mie(m, x)
            geo = w * math.pi * (a * 1e-4) ** 2
            area_abs += geo * (qext - qsca)
            area_sca += geo * qsca
            area_g += geo * qsca * g
        kabs[i] = area_abs / total_mass
        ksca[i] = area_sca / total_mass
        gbar[i] = area_g / area_sca if area_sca > 0 else 0.0
    return {"wavelength_um": wavelengths, "kappa_abs": kabs, "kappa_sca": ksca,
            "g": gbar, "density_g_cm3": density, "porosity": porosity,
            "amax_um": amax_um, "amin_um": amin_um}


def feature_metrics(constants: Path, **kw) -> dict:
    """S = kappa_abs(9.7) / kappa_ext(2.2), plus the feature contrast in kappa_abs."""
    lam = np.array([NEAR_IR_UM, 8.5, SILICATE_UM, 12.5])
    o = opacities(constants, lam, **kw)
    kext_nir = o["kappa_abs"][0] + o["kappa_sca"][0]
    kabs_97 = o["kappa_abs"][2]
    shoulder = 0.5 * (o["kappa_abs"][1] + o["kappa_abs"][3])
    return {"amax_um": o["amax_um"], "porosity": o["porosity"],
            "kappa_ext_2p2": kext_nir, "kappa_abs_9p7": kabs_97,
            "S_9p7_per_nir": kabs_97 / kext_nir,
            "feature_contrast": kabs_97 / shoulder}


def self_test() -> None:
    """Check the Mie routine against the analytic Rayleigh limit."""
    m, x = complex(1.7, 0.03), 1e-3
    qext, qsca, _ = mie(m, x)
    e = (m ** 2 - 1) / (m ** 2 + 2)
    assert abs(qsca / ((8.0 / 3.0) * x ** 4 * abs(e) ** 2) - 1) < 1e-3, "Rayleigh Qsca failed"
    assert abs(qext / (4 * x * e.imag + (8.0 / 3.0) * x ** 4 * abs(e) ** 2) - 1) < 1e-2, \
        "Rayleigh Qext failed"
    # A large, strongly absorbing sphere approaches the geometric limit Qext -> 2.
    qext_big, _, _ = mie(complex(1.7, 0.1), 300.0)
    assert abs(qext_big - 2.0) < 0.1, f"geometric limit failed: {qext_big}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constants", type=Path, default=DEFAULT_CONSTANTS)
    parser.add_argument("--amax", type=float, nargs="*", default=[0.4, 1.0, 3.0, 10.0, 30.0])
    parser.add_argument("--porosity", type=float, nargs="*", default=[0.0, 0.3, 0.6])
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    self_test()
    base = feature_metrics(args.constants, amax_um=AMAX_UM, porosity=0.0)

    rows = []
    for porosity in args.porosity:
        for amax in args.amax:
            m = feature_metrics(args.constants, amax_um=amax, porosity=porosity)
            m["S_relative_to_baseline"] = m["S_9p7_per_nir"] / base["S_9p7_per_nir"]
            # If tau_9.7 scales with S at fixed near-infrared column, the emergent
            # flux brightens by 0.4343 * tau * (1 - S_rel) dex. tau is not known
            # from opacities alone, so report the gain for representative values.
            for tau in (1.0, 2.0, 3.0):
                m[f"flux_gain_dex_if_tau_{int(tau)}"] = (
                    0.4342944819 * tau * (1.0 - m["S_relative_to_baseline"]))
            rows.append(m)

    print(f"material: {args.constants.name}  sha256 {sha256(args.constants)[:16]}")
    print(f"baseline (amax {AMAX_UM} um, compact): "
          f"S = {base['S_9p7_per_nir']:.4f}, contrast = {base['feature_contrast']:.3f}\n")
    print(f"{'porosity':>8} {'amax/um':>8} {'S':>9} {'S/S_base':>9} {'contrast':>9} "
          f"{'gain@tau=2':>11}")
    for m in rows:
        print(f"{m['porosity']:8.2f} {m['amax_um']:8.2f} {m['S_9p7_per_nir']:9.4f} "
              f"{m['S_relative_to_baseline']:9.3f} {m['feature_contrast']:9.3f} "
              f"{m['flux_gain_dex_if_tau_2']:+11.2f}")

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(
            {"material": str(args.constants.relative_to(ROOT)),
             "sha256": sha256(args.constants), "baseline": base, "variations": rows,
             "definition": "S = kappa_abs(9.7 um) / kappa_ext(2.2 um)",
             "note": "Opacity-space screen only. It does not predict emergent flux, "
                     "which also depends on geometry, temperature and radiative transfer."},
            indent=2, default=float) + "\n")
        try:
            shown = args.json.resolve().relative_to(ROOT)
        except ValueError:
            shown = args.json
        print(f"\nwrote {shown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
