#!/usr/bin/env python3
"""Extract an unstitched TMC1A spectrum from the JWST IFU cubes.

This is the spectrum-extraction part of ``JWST_notebook.ipynb`` rewritten as
a small standalone first step for the JWST--MCFOST workflow.  It extracts a
fixed circular aperture from every NIRSpec and MIRI cube and converts the
pipeline surface brightness (MJy/sr) into integrated flux density (Jy).

The result is deliberately *unstitched*: no background subtraction, segment
rescaling, spectral binning, line masking, or aperture correction is applied.
Those choices should be inspected before fitting MCFOST models.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
CUBE_DIR = PROJECT_ROOT / "cubes"
OUTPUT_DIR = HERE / "output"

# Keep Matplotlib's cache inside this workflow, which is writable on the
# school installation and avoids warnings about ~/.matplotlib permissions.
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy import units as u
from astropy.wcs import FITSFixedWarning, WCS
from astropy.wcs.utils import proj_plane_pixel_scales, skycoord_to_pixel
from photutils.aperture import CircularAperture


# Astropy reports harmless DATE/observer-coordinate header normalizations when
# it constructs each WCS.  Suppress those so the 15-line progress report stays
# readable for an interactive beginner run.
warnings.filterwarnings("ignore", category=FITSFixedWarning)


# The on-source coordinate currently used in JWST_notebook.ipynb.
SOURCE_COORD = SkyCoord(
    "4h39m35.2216793062s",
    "25d41m44.1815219788s",
    frame="icrs",
)

# Start with the same physical region at every wavelength.  We will compare
# other apertures and check the MIRI centroid/background in the next step.
APERTURE_RADIUS_ARCSEC = 1.0

# JWST's DO_NOT_USE flag is bit 0.  Other flagged pixels are retained for now;
# their effect can be inspected rather than silently applying a broad DQ cut.
DO_NOT_USE = 1


@dataclass(frozen=True)
class Segment:
    name: str
    instrument: str
    filename: str


SEGMENTS = [
    Segment(
        "G140H_F100LP",
        "NIRSpec",
        "jw04537-o001_t001_nirspec_clean_jun26_ns_clean_g140h-f100lp_s3d.fits",
    ),
    Segment(
        "G235H_F170LP",
        "NIRSpec",
        "jw04537-o001_t001_nirspec_clean_jun26_ns_clean_g235h-f170lp_s3d.fits",
    ),
    Segment(
        "G395H_F290LP",
        "NIRSpec",
        "jw04537-o001_t001_nirspec_clean_jun26_ns_clean_g395h-f290lp_s3d.fits",
    ),
    Segment("CH1_short", "MIRI", "TMC1A_ch1-short_s3d.fits"),
    Segment("CH1_medium", "MIRI", "TMC1A_ch1-medium_s3d.fits"),
    Segment("CH1_long", "MIRI", "TMC1A_ch1-long_s3d.fits"),
    Segment("CH2_short", "MIRI", "TMC1A_ch2-short_s3d.fits"),
    Segment("CH2_medium", "MIRI", "TMC1A_ch2-medium_s3d.fits"),
    Segment("CH2_long", "MIRI", "TMC1A_ch2-long_s3d.fits"),
    Segment("CH3_short", "MIRI", "TMC1A_ch3-short_s3d.fits"),
    Segment("CH3_medium", "MIRI", "TMC1A_ch3-medium_s3d.fits"),
    Segment("CH3_long", "MIRI", "TMC1A_ch3-long_s3d.fits"),
    Segment("CH4_short", "MIRI", "TMC1A_ch4-short_s3d.fits"),
    Segment("CH4_medium", "MIRI", "TMC1A_ch4-medium_s3d.fits"),
    Segment("CH4_long", "MIRI", "TMC1A_ch4-long_s3d.fits"),
]


def wavelength_axis(header: fits.Header) -> np.ndarray:
    """Return the linear spectral WCS axis in micrometres."""
    pixel_fits = np.arange(header["NAXIS3"], dtype=float) + 1.0
    wavelength = (
        header["CRVAL3"]
        + (pixel_fits - header["CRPIX3"]) * header["CDELT3"]
    )
    unit = u.Unit(header.get("CUNIT3", "um"))
    return (wavelength * unit).to_value(u.um)


def aperture_weight_map(
    header: fits.Header,
    image_shape: tuple[int, int],
    source_coord: SkyCoord = SOURCE_COORD,
    aperture_radius_arcsec: float = APERTURE_RADIUS_ARCSEC,
) -> tuple[np.ndarray, float, float, float]:
    """Make the notebook's exact circular aperture at the source position."""
    celestial_wcs = WCS(header).celestial
    x_pixel, y_pixel = skycoord_to_pixel(
        source_coord, celestial_wcs, origin=0, mode="all"
    )

    scales_deg = proj_plane_pixel_scales(celestial_wcs)
    pixel_scale_arcsec = float(np.sqrt(np.abs(np.prod(scales_deg))) * 3600.0)
    radius_pixel = aperture_radius_arcsec / pixel_scale_arcsec

    aperture = CircularAperture((x_pixel, y_pixel), r=radius_pixel)
    mask = aperture.to_mask(method="exact")
    weight = mask.to_image(image_shape)
    if weight is None:
        raise ValueError("The requested aperture does not overlap this cube")
    weight = np.nan_to_num(weight, nan=0.0)
    return weight, float(x_pixel), float(y_pixel), pixel_scale_arcsec


def extract_segment(
    segment: Segment,
    source_coord: SkyCoord = SOURCE_COORD,
    aperture_radius_arcsec: float = APERTURE_RADIUS_ARCSEC,
) -> Table:
    """Extract one segment as aperture-integrated F_nu in Jy."""
    cube_path = CUBE_DIR / segment.filename
    if not cube_path.exists():
        raise FileNotFoundError(f"Missing cube: {cube_path}")

    # do_not_scale_image_data permits memory mapping of the unsigned DQ array.
    # Bit 0 has the same value in its raw signed representation.
    with fits.open(
        cube_path,
        memmap=True,
        do_not_scale_image_data=True,
    ) as hdul:
        sci_hdu = hdul["SCI"]
        header = sci_hdu.header
        if header.get("BUNIT", "").lower().replace(" ", "") != "mjy/sr":
            raise ValueError(
                f"{segment.name}: expected SCI in MJy/sr, got "
                f"{header.get('BUNIT')!r}"
            )

        sci = sci_hdu.data
        err = hdul["ERR"].data
        dq = hdul["DQ"].data
        wmap = hdul["WMAP"].data

        wavelengths_um = wavelength_axis(header)
        n_wave, ny, nx = sci.shape
        if len(wavelengths_um) != n_wave:
            raise ValueError(f"{segment.name}: WCS/data wavelength mismatch")

        weights_2d, x_pixel, y_pixel, pixel_scale = aperture_weight_map(
            header,
            (ny, nx),
            source_coord=source_coord,
            aperture_radius_arcsec=aperture_radius_arcsec,
        )
        aperture_weight = weights_2d.sum()
        if aperture_weight <= 0:
            raise ValueError(f"{segment.name}: empty source aperture")

        pixel_area_sr = header.get("PIXAR_SR")
        if pixel_area_sr is None:
            raise KeyError(f"{segment.name}: SCI header has no PIXAR_SR")
        mjy_sr_to_jy_pixel = float(pixel_area_sr) * 1.0e6

        flux_jy = np.full(n_wave, np.nan, dtype=float)
        error_jy = np.full(n_wave, np.nan, dtype=float)
        coverage = np.zeros(n_wave, dtype=float)

        # Work in chunks so each cube remains comfortably within memory.
        chunk_size = 256
        aperture_weights = weights_2d[None, :, :]
        for start in range(0, n_wave, chunk_size):
            stop = min(start + chunk_size, n_wave)
            sci_chunk = np.asarray(sci[start:stop], dtype=float)
            err_chunk = np.asarray(err[start:stop], dtype=float)
            dq_chunk = np.asarray(dq[start:stop])
            wmap_chunk = np.asarray(wmap[start:stop], dtype=float)

            valid = (
                np.isfinite(sci_chunk)
                & np.isfinite(err_chunk)
                & (err_chunk > 0.0)
                & ((dq_chunk & DO_NOT_USE) == 0)
                & np.isfinite(wmap_chunk)
                & (wmap_chunk > 0.0)
                & (aperture_weights > 0.0)
            )
            effective_weights = aperture_weights * valid

            coverage[start:stop] = (
                effective_weights.sum(axis=(1, 2)) / aperture_weight
            )
            summed_flux = np.sum(
                np.where(valid, sci_chunk, 0.0) * aperture_weights,
                axis=(1, 2),
            )
            summed_variance = np.sum(
                np.where(valid, err_chunk, 0.0) ** 2 * aperture_weights**2,
                axis=(1, 2),
            )

            has_data = coverage[start:stop] > 0.0
            flux_jy[start:stop][has_data] = (
                summed_flux[has_data] * mjy_sr_to_jy_pixel
            )
            error_jy[start:stop][has_data] = (
                np.sqrt(summed_variance[has_data]) * mjy_sr_to_jy_pixel
            )

    # nu F_nu = lambda F_lambda, useful for direct comparison with an SED.
    c_m_s = 299_792_458.0
    lambda_f_lambda = (
        c_m_s / (wavelengths_um * 1.0e-6) * flux_jy * 1.0e-26
    )
    lambda_f_lambda_error = (
        c_m_s / (wavelengths_um * 1.0e-6) * error_jy * 1.0e-26
    )

    table = Table()
    table["wavelength_um"] = wavelengths_um
    table["flux_jy"] = flux_jy
    table["error_jy"] = error_jy
    table["lambda_f_lambda_w_m2"] = lambda_f_lambda
    table["lambda_f_lambda_error_w_m2"] = lambda_f_lambda_error
    table["coverage_fraction"] = coverage
    table["instrument"] = np.repeat(segment.instrument, n_wave)
    table["segment"] = np.repeat(segment.name, n_wave)
    table["source_file"] = np.repeat(segment.filename, n_wave)
    table.meta.update(
        {
            "source": "TMC1A / IRAS 04365+2535",
            "ra_icrs_deg": source_coord.ra.deg,
            "dec_icrs_deg": source_coord.dec.deg,
            "aperture_radius_arcsec": aperture_radius_arcsec,
            "flux_definition": "raw aperture-integrated F_nu",
            "background_subtracted_here": False,
            "segment_rescaling_applied": False,
            "x_source_pixel": x_pixel,
            "y_source_pixel": y_pixel,
            "pixel_scale_arcsec": pixel_scale,
        }
    )
    return table


def plot_spectrum(table: Table, output_path: Path) -> None:
    """Plot every segment independently so overlaps remain visible."""
    fig, ax = plt.subplots(figsize=(11, 6.5))
    colors = plt.get_cmap("turbo")(np.linspace(0.05, 0.95, len(SEGMENTS)))

    for color, segment in zip(colors, SEGMENTS):
        keep_segment = np.asarray(table["segment"] == segment.name)
        wavelength = np.asarray(table["wavelength_um"])[keep_segment]
        sed = np.asarray(table["lambda_f_lambda_w_m2"])[keep_segment]
        coverage = np.asarray(table["coverage_fraction"])[keep_segment]
        good = np.isfinite(sed) & (sed > 0.0) & (coverage >= 0.5)
        ax.plot(
            wavelength[good],
            sed[good],
            lw=0.8,
            color=color,
            label=segment.name,
        )

    ax.set(xscale="log", yscale="log")
    ax.set_xlabel(r"Wavelength $\lambda$ ($\mu$m)")
    ax.set_ylabel(r"$\lambda F_\lambda$ (W m$^{-2}$)")
    ax.set_title(
        f'TMC1A: raw JWST spectrum in a {APERTURE_RADIUS_ARCSEC:.1f}" radius aperture'
    )
    ax.grid(alpha=0.2, which="both")
    ax.legend(ncol=3, fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(
        "TMC1A source: "
        f"RA={SOURCE_COORD.ra.to_string(unit=u.hour)}, "
        f"Dec={SOURCE_COORD.dec.to_string(unit=u.deg, alwayssign=True)}"
    )
    print(f'Aperture radius: {APERTURE_RADIUS_ARCSEC:.2f}"')

    tables = []
    for index, segment in enumerate(SEGMENTS, start=1):
        print(f"[{index:02d}/{len(SEGMENTS)}] Extracting {segment.name}")
        segment_table = extract_segment(segment)
        usable = np.count_nonzero(
            np.isfinite(segment_table["flux_jy"])
            & (segment_table["coverage_fraction"] >= 0.5)
        )
        print(
            f"       {segment_table['wavelength_um'][0]:.3f}--"
            f"{segment_table['wavelength_um'][-1]:.3f} um; "
            f"{usable}/{len(segment_table)} wavelengths >=50% covered"
        )
        tables.append(segment_table)

    combined = vstack(tables, metadata_conflicts="silent")
    combined.meta = {
        "source": "TMC1A / IRAS 04365+2535",
        "ra_icrs_deg": SOURCE_COORD.ra.deg,
        "dec_icrs_deg": SOURCE_COORD.dec.deg,
        "aperture_radius_arcsec": APERTURE_RADIUS_ARCSEC,
        "input_unit": "MJy/sr",
        "output_flux_unit": "Jy",
        "background_subtracted_here": False,
        "segment_rescaling_applied": False,
        "note": "Formal ERR propagation does not include cube covariance.",
    }

    table_path = OUTPUT_DIR / "tmc1a_sed_unstitched.ecsv"
    plot_path = OUTPUT_DIR / "tmc1a_sed_unstitched.png"
    combined.write(table_path, format="ascii.ecsv", overwrite=True)
    plot_spectrum(combined, plot_path)

    print(f"Saved table: {table_path}")
    print(f"Saved plot:  {plot_path}")
    print("This is intentionally unstitched; inspect aperture/background next.")


if __name__ == "__main__":
    main()
