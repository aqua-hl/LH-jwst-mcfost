"""Validated edits of the preserved two-zone MCFOST 4.1 parameter template.

This is deliberately not a general MCFOST parameter parser.  All edited rows
are identified by their section/field labels, and unrecognised dust structures
are rejected.  Cavity geometry and stellar accretion are command-line physics,
not fields silently invented in the parameter file.
"""
from __future__ import annotations

import math
from collections.abc import Mapping


PHYSICAL_KEYS = frozenset({
    "envelope_dust_mass_msun", "envelope_amax_um", "envelope_size_exponent",
    "cavity_half_opening_deg", "inclination_deg", "distance_pc",
    "stellar_radius_rsun", "stellar_temperature_k", "stellar_mass_msun",
    "accretion_rate_msun_per_year", "envelope_ice_volume_fraction",
})
NUMERIC_KEYS = frozenset({
    "photons_temperature", "photons_image", "photons_sed", "grid_nr",
    "grid_ntheta", "grid_n_inner", "image_npix", "image_size_au", "grains",
})
STAGES = frozenset({"temperature", "image", "coeval"})


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number, not {value!r}")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return float(value)


def _positive(value: object, name: str) -> float:
    number = _number(value, name)
    if number <= 0:
        raise ValueError(f"{name} must be greater than zero")
    return number


def _integer(value: object, name: str, minimum: int = 1) -> int:
    number = _number(value, name)
    if number != math.floor(number) or number < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(number)


def _token(value: float | int) -> str:
    return format(value, ".15g")


def validate_parameters(parameters: Mapping[str, object]) -> None:
    """Validate named physical overrides without supplying hidden defaults."""
    if not isinstance(parameters, Mapping):
        raise ValueError("parameters must be a mapping of named physical quantities")
    unknown = set(parameters) - PHYSICAL_KEYS
    if unknown:
        raise ValueError(f"Unsupported physical parameters: {sorted(unknown)}; "
                         "dust families/compositions other than the coated template need a separate implementation")
    for name, value in parameters.items():
        number = _number(value, name)
        if name in {"inclination_deg", "accretion_rate_msun_per_year", "envelope_size_exponent"}:
            if number < 0:
                raise ValueError(f"{name} must be non-negative")
            if name == "inclination_deg" and number > 90:
                raise ValueError("inclination_deg must be between 0 and 90 degrees from the pole")
        elif name == "cavity_half_opening_deg":
            if not 0 < number < 90:
                raise ValueError("cavity_half_opening_deg must lie strictly between 0 and 90")
        elif name == "envelope_ice_volume_fraction":
            if not 0 < number < 1:
                raise ValueError("envelope_ice_volume_fraction must lie strictly between 0 and 1; "
                                 "zero-ice and pure-ice populations require different dust templates")
        elif number <= 0:
            raise ValueError(f"{name} must be greater than zero")


def physical_args(parameters: Mapping[str, object]) -> list[str]:
    """Return an explicit conical cavity and star-1 accretion prescription.

    The cavity is z = h0 (r/100 au)^1, with h0 = 100/tan(theta).
    Its requested half-angle is measured from the polar axis; the density grid
    discretises that geometry.  ``-Mdot`` is stellar, not viscous disk heating.
    Neither ``-rt2`` nor ``-mol`` belongs to these physical arguments.
    """
    validate_parameters(parameters)
    required = {"cavity_half_opening_deg", "accretion_rate_msun_per_year"}
    if required - set(parameters):
        raise ValueError(f"Explicit CLI physics is required: {sorted(required - set(parameters))}")
    angle = float(parameters["cavity_half_opening_deg"])
    tangent = math.tan(math.radians(angle))
    if tangent == 0:
        raise ValueError("cavity_half_opening_deg is too small to represent numerically")
    h0 = 100.0 / tangent
    if not math.isfinite(h0):
        raise ValueError("cavity_half_opening_deg produces a non-finite cavity height")
    return ["-cavity", _token(h0), "100", "1", "-Mdot", "1",
            _token(float(parameters["accretion_rate_msun_per_year"]))]


class _Template:
    def __init__(self, text: str):
        if not isinstance(text, str) or not text.strip():
            raise ValueError("template_text must be a non-empty MCFOST parameter file")
        self.lines = text.splitlines()
        first = next(line for line in self.lines if line.strip())
        if first.split()[0] != "4.1" or "mcfost version" not in first.lower():
            raise ValueError("Only labelled MCFOST 4.1 parameter templates are supported")

    def section(self, start: str, end: str | None = None) -> range:
        first = self.one(start) + 1
        last = self.one(end) if end else len(self.lines)
        if first >= last:
            raise ValueError(f"Unsupported section ordering: {start!r} / {end!r}")
        return range(first, last)

    def one(self, tag: str, scope: range | None = None) -> int:
        matches = [i for i in (scope if scope is not None else range(len(self.lines)))
                   if tag.casefold() in self.lines[i].casefold()]
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one labelled field {tag!r}; found {len(matches)}")
        return matches[0]

    def values(self, tag: str, n: int, scope: range | None = None) -> list[str]:
        index = self.one(tag, scope)
        line = self.lines[index]
        tokens = line[:line.casefold().index(tag.casefold())].split()
        if len(tokens) != n:
            raise ValueError(f"Field {tag!r} must have {n} values before its label; got {tokens}")
        return tokens

    def set(self, tag: str, values: list[object], scope: range | None = None) -> None:
        index = self.one(tag, scope)
        line = self.lines[index]
        comment = line[line.casefold().index(tag.casefold()):]
        tokens = [_token(v) if isinstance(v, (float, int)) else str(v) for v in values]
        self.lines[index] = "  " + "  ".join(tokens) + "    " + comment


def render_parameter(template_text: str, parameters: Mapping[str, object],
                     numerics: Mapping[str, object], stage: str,
                     wavelength_um: float | None = None) -> str:
    """Render temperature, stock-image, or custom-coeval parameter content.

    Unspecified supported quantities retain the template values.  ``image``
    needs runner arguments ``-img <wavelength> -rt2`` and the freshly computed
    temperature file.  ``coeval`` additionally needs a capability-checked
    custom executable, ``-rt_sed_method 2``, and a one-node ``wavelength.lambda``
    file.  Here ``-rt2`` is the source-function method, not that spatial option.
    """
    validate_parameters(parameters)
    if not isinstance(numerics, Mapping):
        raise ValueError("numerics must be a mapping")
    unknown = set(numerics) - NUMERIC_KEYS
    if unknown:
        raise ValueError(f"Unsupported numeric parameters: {sorted(unknown)}")
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {sorted(STAGES)}, not {stage!r}")
    if stage in {"image", "coeval"}:
        _positive(wavelength_um, "wavelength_um")
    elif wavelength_um is not None:
        raise ValueError("temperature stage must not have a monochromatic wavelength_um")
    for name, value in numerics.items():
        if name == "image_size_au":
            _positive(value, name)
        else:
            _integer(value, name, 2 if name in {"grid_nr", "grid_ntheta", "grains"} else 1)

    t = _Template(template_text)
    zone_count_scope = t.section("#Number of zones", "#Density structure")
    if t.values("needs to be 1 if", 1, zone_count_scope) != ["2"]:
        raise ValueError("Only the two-zone disk + envelope template is supported")
    density = t.section("#Density structure", "#Grain properties")
    z1 = t.one("ZONE 1 - zone type", density)
    z2 = t.one("ZONE 2 - zone type", density)
    if z1 >= z2 or t.values("ZONE 1 - zone type", 1, density) != ["1"] or t.values("ZONE 2 - zone type", 1, density) != ["3"]:
        raise ValueError("Expected zone 1 disk (type 1), zone 2 envelope (type 3)")
    envelope_density = range(z2, density.stop)

    grains = t.section("#Grain properties", "#Molecular RT settings")
    species_headers = [i for i in grains if "number of species" in t.lines[i].casefold()]
    if len(species_headers) != 2 or "ZONE 1" not in t.lines[species_headers[0]] or "ZONE 2" not in t.lines[species_headers[1]]:
        raise ValueError("Expected two explicitly zone-labelled dust species blocks")
    disk_grains = range(species_headers[0], species_headers[1])
    envelope_grains = range(species_headers[1], grains.stop)
    if t.values("Number of species", 1, disk_grains) != ["2"] or t.values("Number of species", 1, envelope_grains) != ["1"]:
        raise ValueError("Only two disk species and one coated envelope species are supported; DHS/disjoint populations are unsupported")
    envelope_type = t.values("Grain type", 6, envelope_grains)
    if envelope_type[:3] != ["Mie", "2", "2"] or float(envelope_type[3]) != 0.0 or float(envelope_type[4]) != 1.0:
        raise ValueError("Envelope must be a non-porous Mie species with two coating components and unit species mass fraction")
    core = t.values("Optical indices file", 2, range(species_headers[1], t.one("ice_opct.dat", envelope_grains)))
    if core[0] != "Draine_Si_sUV.dat":
        raise ValueError("Only a Draine_Si_sUV.dat core with ice_opct.dat mantle is supported")
    ice_row = t.one("ice_opct.dat", envelope_grains)
    ice = t.values("Optical indices file", 2, range(ice_row, ice_row + 1))
    if ice[0] != "ice_opct.dat" or not math.isclose(float(core[1]) + float(ice[1]), 1.0, abs_tol=1e-10):
        raise ValueError("Envelope core and mantle volume fractions must sum to one")
    if not (0 < float(core[1]) < 1 and 0 < float(ice[1]) < 1):
        raise ValueError("The supported coated template requires nonzero core and mantle fractions")
    optical_rows = [i for i in disk_grains if "optical indices file" in t.lines[i].casefold()]
    type_rows = [i for i in disk_grains if "grain type" in t.lines[i].casefold()]
    if len(optical_rows) != 2 or len(type_rows) != 2:
        raise ValueError("Unsupported disk dust species structure")
    if [t.values("Optical indices file", 2, range(i, i+1))[0] for i in optical_rows] != ["Draine_Si_sUV.dat", "ac_opct.dat"]:
        raise ValueError("Disk must retain the Draine silicate and amorphous-carbon species")
    for i in type_rows:
        kind = t.values("Grain type", 6, range(i, i+1))
        if kind[:3] != ["Mie", "1", "1"] or float(kind[3]) != 0:
            raise ValueError("Unsupported disk grain type, mixing rule or porosity")

    geometry = t.section("#Grid geometry and size", "#Maps")
    if t.values("1 = cylindrical", 1, geometry) != ["2"]:
        raise ValueError("Envelope calculations require the spherical geometry template")
    scattering = t.section("#Scattering method", "#Symmetries")
    if t.values("1=exact phase function", 1, scattering) != ["1"]:
        raise ValueError("Full polarized output requires the exact phase-function template, not HG")
    star_scope = t.section("#Star properties")
    if t.values("Number of stars", 1, star_scope) != ["1"]:
        raise ValueError("Only one stellar source is supported")
    stellar_tag = "Temp, radius (solar radius),M (solar mass)"
    star = t.values(stellar_tag, 7, star_scope)
    if star[6] != "T" or any(float(x) != 0.0 for x in star[3:6]):
        raise ValueError("Expected one centred star with automatic stellar-spectrum selection")

    photon_tags = {"photons_temperature": "nbr_photons_eq_th", "photons_sed": "nbr_photons_lambda",
                   "photons_image": "nbr_photons_image"}
    for key, tag in photon_tags.items():
        old = t.values(tag, 1)
        t.set(tag, [int(numerics[key])] if key in numerics else old)
    flags = {"temperature": ["T", "F", "T"], "image": ["F", "F", "T"], "coeval": ["F", "T", "F"]}
    t.values("compute temperature?", 3)
    t.set("compute temperature?", flags[stage])
    t.values("wavelength file (", 1)
    t.set("wavelength file (", ["wavelength.lambda"])
    t.values("separation of different contributions?", 2)
    t.set("separation of different contributions?", ["T", "T"])
    t.values("image symmetry", 1)
    t.set("image symmetry", ["F"])

    grid_tag = "n_rad (log distribution)"
    grid = t.values(grid_tag, 4, geometry)
    if grid[2] != "1":
        raise ValueError("Only the axisymmetric n_az=1 density template is supported")
    for key, index in {"grid_nr": 0, "grid_ntheta": 1, "grid_n_inner": 3}.items():
        if key in numerics:
            grid[index] = int(numerics[key])
    t.set(grid_tag, grid, geometry)
    maps = t.section("#Maps", "#Scattering method")
    size = t.values("grid (nx,ny), size [AU]", 3, maps)
    if "image_npix" in numerics:
        size[:2] = [int(numerics["image_npix"])] * 2
    if "image_size_au" in numerics:
        size[2] = float(numerics["image_size_au"])
    t.set("grid (nx,ny), size [AU]", size, maps)
    incl = t.values("RT: imin, imax, n_incl, centered ?", 4, maps)
    if incl[2:] != ["1", "F"] or float(incl[0]) != float(incl[1]):
        raise ValueError("Expected one explicitly sampled inclination, not an inclination range")
    if "inclination_deg" in parameters:
        incl[:2] = [float(parameters["inclination_deg"])] * 2
    t.set("RT: imin, imax, n_incl, centered ?", incl, maps)
    if "distance_pc" in parameters:
        t.values("distance (pc)", 1, maps)
        t.set("distance (pc)", [float(parameters["distance_pc"])], maps)

    if "envelope_dust_mass_msun" in parameters:
        mass = t.values("dust mass,", 2, envelope_density)
        mass[0] = float(parameters["envelope_dust_mass_msun"])
        t.set("dust mass,", mass, envelope_density)
    size_tag = "amin, amax [mum], aexp, n_grains"
    envelope_size = t.values(size_tag, 4, envelope_grains)
    amin = float(envelope_size[0])
    amax = float(parameters.get("envelope_amax_um", envelope_size[1]))
    if amax <= amin:
        raise ValueError(f"envelope_amax_um must exceed the template minimum grain size ({amin:g} micron)")
    if "envelope_amax_um" in parameters:
        envelope_size[1] = amax
    if "envelope_size_exponent" in parameters:
        envelope_size[2] = float(parameters["envelope_size_exponent"])
    t.set(size_tag, envelope_size, envelope_grains)
    if "envelope_ice_volume_fraction" in parameters:
        fraction = float(parameters["envelope_ice_volume_fraction"])
        core_row = t.one("Draine_Si_sUV.dat", envelope_grains)
        t.set("Optical indices file", ["Draine_Si_sUV.dat", 1-fraction], range(core_row, core_row+1))
        t.set("Optical indices file", ["ice_opct.dat", fraction], range(ice_row, ice_row+1))
    if "grains" in numerics:
        for i in grains:
            if size_tag in t.lines[i]:
                values = t.values(size_tag, 4, range(i, i+1))
                values[3] = int(numerics["grains"])
                t.set(size_tag, values, range(i, i+1))
    for key, index in {"stellar_temperature_k": 0, "stellar_radius_rsun": 1, "stellar_mass_msun": 2}.items():
        if key in parameters:
            star[index] = float(parameters[key])
    t.set(stellar_tag, star, star_scope)
    return "\n".join(t.lines) + "\n"
