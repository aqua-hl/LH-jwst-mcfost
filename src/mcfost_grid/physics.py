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
    "envelope_ice_mass_fraction",
})
NUMERIC_KEYS = frozenset({
    "photons_temperature", "photons_image", "photons_sed", "grid_nr",
    "grid_ntheta", "grid_n_inner", "image_npix", "image_size_au", "grains", "random_seed",
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
                         "dust families/compositions outside the supported templates need a separate implementation")
    if {"envelope_ice_volume_fraction", "envelope_ice_mass_fraction"} <= set(parameters):
        raise ValueError("Specify an envelope ice volume fraction or mass fraction, not both")
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
            if not 0 <= number < 1:
                raise ValueError("envelope_ice_volume_fraction must be >=0 and <1; "
                                 "zero renders a genuine single-component bare-silicate species")
        elif name == "envelope_ice_mass_fraction":
            if not 0 <= number < 1:
                raise ValueError("envelope_ice_mass_fraction must be >=0 and <1")
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


def numerical_args(numerics: Mapping[str, object]) -> list[str]:
    """Return optional simulator CLI numerics without changing its defaults.

    MCFOST's ``-seed`` also nests products under ``root_dir/seed=VALUE``;
    ``-Tfile`` remains relative to ``root_dir``. The seed is frozen numerical
    configuration, separate from the physical model and parameter-file rows.
    """
    if not isinstance(numerics, Mapping):
        raise ValueError("numerics must be a mapping")
    if "random_seed" not in numerics:
        return []
    seed = numerics["random_seed"]
    if isinstance(seed, bool) or not isinstance(seed, int) or not 1 <= seed <= 2**31 - 1:
        raise ValueError("random_seed must be an integer between 1 and 2147483647")
    return ["-seed", str(seed)]


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


def _disjoint_dhs_scopes(t: _Template, envelope: range) -> tuple[range, range]:
    """Validate the archived v02 pure-silicate + pure-H2O population layout.

    The ice population is fixed, including its size distribution and DHS
    shape. Only the silicate size distribution and the two species' mass
    weights can be edited by this deliberately bounded renderer.
    """
    headers = [i for i in envelope if "grain type" in t.lines[i].casefold()]
    if len(headers) != 2:
        raise ValueError("DHS/disjoint envelope requires exactly two labelled grain species")
    scopes = (range(headers[0], headers[1]), range(headers[1], envelope.stop))
    fractions = []
    for scope, filename in zip(scopes, ("Draine_Si_sUV.dat", "H2O_30K_Leiden_mcfost.dat")):
        kind = t.values("Grain type", 6, scope)
        if kind[:3] != ["DHS", "1", "1"] or float(kind[3]) != 0 or float(kind[5]) != .1:
            raise ValueError("DHS/disjoint species must be pure, non-porous DHS with Vmax=0.1")
        fraction = float(kind[4])
        if not math.isfinite(fraction) or not 0 < fraction < 1:
            raise ValueError("DHS/disjoint template requires positive species mass fractions below one")
        fractions.append(fraction)
        optical = t.values("Optical indices file", 2, scope)
        if optical[0] != filename or float(optical[1]) != 1:
            raise ValueError(f"DHS/disjoint species requires pure {filename} with unit within-species volume fraction")
        if t.values("Heating method", 1, scope) != ["1"]:
            raise ValueError("DHS/disjoint species must retain equilibrium LTE heating")
        size = [float(value) for value in t.values("amin, amax [mum], aexp, n_grains", 4, scope)]
        if any(not math.isfinite(value) for value in size) or size[0] != .03 or size[3] != 50:
            raise ValueError("DHS/disjoint species must retain amin=0.03 micron and 50 grain bins")
        if size[1] <= size[0] or size[2] < 0:
            raise ValueError("DHS/disjoint species has an invalid grain size distribution")
        if filename == "H2O_30K_Leiden_mcfost.dat" and size[1:3] != [.4, 2.75]:
            raise ValueError("DHS/disjoint H2O population must retain amax=0.4 micron and exponent=2.75")
    if not math.isclose(sum(fractions), 1.0, rel_tol=0, abs_tol=1e-10):
        raise ValueError("DHS/disjoint species mass fractions must sum to one")
    return scopes


def render_parameter(template_text: str, parameters: Mapping[str, object],
                     numerics: Mapping[str, object], stage: str,
                     wavelength_um: float | None = None) -> str:
    """Render temperature, stock-image, or custom-coeval parameter content.

    Unspecified supported quantities retain the template values.  ``image``
    needs runner arguments ``-img <wavelength> -rt2`` and the freshly computed
    temperature file.  ``coeval`` additionally needs a capability-checked
    custom executable, ``-rt_sed_method 2``, and a one-node ``wavelength.lambda``
    file.  Here ``-rt2`` is the source-function method, not that spatial option.

    ``envelope_amax_um`` and ``envelope_size_exponent`` refer to the coated
    population for Mie templates, and only to the silicate population for the
    v02 disjoint DHS template. Its separate H2O grain distribution stays fixed.
    """
    validate_parameters(parameters)
    if not isinstance(numerics, Mapping):
        raise ValueError("numerics must be a mapping")
    unknown = set(numerics) - NUMERIC_KEYS
    if unknown:
        raise ValueError(f"Unsupported numeric parameters: {sorted(unknown)}")
    numerical_args(numerics)
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {sorted(STAGES)}, not {stage!r}")
    if stage in {"image", "coeval"}:
        _positive(wavelength_um, "wavelength_um")
    elif wavelength_um is not None:
        raise ValueError("temperature stage must not have a monochromatic wavelength_um")
    for name, value in numerics.items():
        if name == "random_seed":
            continue  # Validated above; supplied on the CLI, not in .para.
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
    if t.values("Number of species", 1, disk_grains) != ["2"]:
        raise ValueError("Only two disk species are supported")
    envelope_count = t.values("Number of species", 1, envelope_grains)
    disjoint_dhs = envelope_count == ["2"]
    bare_template = False
    ice_row = None
    if disjoint_dhs:
        if "envelope_ice_mass_fraction" not in parameters:
            raise ValueError("DHS/disjoint template requires explicit envelope_ice_mass_fraction")
        if float(parameters["envelope_ice_mass_fraction"]) == 0:
            raise ValueError("DHS/disjoint template requires a positive ice mass fraction; a bare DHS null needs a separate template")
        if "grains" in numerics and numerics["grains"] != 50:
            raise ValueError("DHS/disjoint template must retain 50 grain bins for the fixed ice prescription")
        silicate_grains, ice_grains = _disjoint_dhs_scopes(t, envelope_grains)
    else:
        if envelope_count != ["1"]:
            raise ValueError("Only one Mie envelope species or the two-species DHS/disjoint template is supported")
        if "envelope_ice_mass_fraction" in parameters:
            raise ValueError("Mie envelope templates require an ice volume fraction, not an ice mass fraction")
        silicate_grains = envelope_grains
        envelope_type = t.values("Grain type", 6, envelope_grains)
        bare_template = envelope_type[:3] == ["Mie", "1", "1"]
        if (not bare_template and envelope_type[:3] != ["Mie", "2", "2"]) or float(envelope_type[3]) != 0.0 or float(envelope_type[4]) != 1.0:
            raise ValueError("Envelope must be non-porous coated or bare Mie with unit species mass fraction")
        optical_envelope = [i for i in envelope_grains if "optical indices file" in t.lines[i].casefold()]
        if len(optical_envelope) != (1 if bare_template else 2):
            raise ValueError("Envelope optical-component count disagrees with grain type")
        core_row = optical_envelope[0]
        core = t.values("Optical indices file", 2, range(core_row, core_row+1))
        if core[0] != "Draine_Si_sUV.dat":
            raise ValueError("Only a Draine_Si_sUV.dat core with ice_opct.dat mantle is supported")
        ice_row = None if bare_template else optical_envelope[1]
        if bare_template:
            if float(core[1]) != 1 or parameters.get("envelope_ice_volume_fraction", 0) != 0:
                raise ValueError("Bare template requires unit core fraction and zero ice; use the coated template to add ice")
        else:
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
    envelope_size = t.values(size_tag, 4, silicate_grains)
    amin = float(envelope_size[0])
    amax = float(parameters.get("envelope_amax_um", envelope_size[1]))
    if amax <= amin:
        raise ValueError(f"envelope_amax_um must exceed the template minimum grain size ({amin:g} micron)")
    if "envelope_amax_um" in parameters:
        envelope_size[1] = amax
    if "envelope_size_exponent" in parameters:
        envelope_size[2] = float(parameters["envelope_size_exponent"])
    t.set(size_tag, envelope_size, silicate_grains)
    if disjoint_dhs:
        fraction = float(parameters["envelope_ice_mass_fraction"])
        for scope, weight in ((silicate_grains, 1-fraction), (ice_grains, fraction)):
            species = t.values("Grain type", 6, scope)
            species[4] = weight
            t.set("Grain type", species, scope)
    if "envelope_ice_volume_fraction" in parameters:
        fraction = float(parameters["envelope_ice_volume_fraction"])
        core_row = t.one("Draine_Si_sUV.dat", envelope_grains)
        t.set("Optical indices file", ["Draine_Si_sUV.dat", 1-fraction], range(core_row, core_row+1))
        if ice_row is not None:
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
    # Remove the mantle completely at the null hypothesis, rather than asking
    # a coated-sphere solver to evaluate a zero-thickness mantle. Do this last
    # so deleting its optical-constant row cannot shift the validated scopes.
    if parameters.get("envelope_ice_volume_fraction") == 0 and not bare_template and not disjoint_dhs:
        envelope_type[:3] = ["Mie", "1", "1"]
        t.set("Grain type", envelope_type, envelope_grains)
        del t.lines[ice_row]
    return "\n".join(t.lines) + "\n"
