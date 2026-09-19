"""Rebuild the VEQDB manuscript figures from the frozen release inputs.

The script does not alter reported numbers. Its conversion branch needs the
VEQPy compact-record and G-EQDSK projection interfaces; all inputs and accepted
records required for the rendering branch are local to this repository.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from contourpy import contour_generator
from figure_style import (
    BLUE,
    CONFIG,
    DOUBLE_FIGURE_WIDTH_IN,
    GREEN,
    GREY,
    LINE_STYLES,
    ORANGE,
    RED,
    SINGLE_FIGURE_WIDTH_IN,
    TEAL,
    YELLOW,
    apply_style,
    physical_layout,
    quantitative_grid,
    save_figure,
    top_legend_y,
)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import RectBivariateSpline, griddata

from veqpy import Equilibrium, Geqdsk
from gfile_projection.implementation import settings as fitting_settings
from gfile_projection.implementation.flux import setup as flux_setup

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PAPER = ROOT / "paper"
FIGURE_DATA = ROOT / "data" / "figure-data"
TCV_BATCH = (
    ROOT
    / "data"
    / "figure-inputs"
    / "formal_TCV_GAQ_2026-09-15"
    / "TCV"
)
START_BATCH = (
    ROOT
    / "data"
    / "figure-inputs"
    / "START_GAQ_test_2026-09-15_corrected2"
    / "START"
)
DEVICE_MANIFESTS = ROOT / "data" / "release"
EAST_SOURCE = (
    ROOT
    / "data"
    / "gfile"
    / "east_standard"
    / "standard"
    / "g048059.03650.geqdsk"
)
FREEGSNKE = ROOT / "data" / "gfile" / "freegsnke_reference"

DEVICE_COVERAGE_HEIGHT_IN = 4.80
TCV_GEOMETRY_WIDTH_IN = 468.0 / 72.27  # \textwidth from the 3p final layout.
TCV_GEOMETRY_HEIGHT_IN = 3.15
START_PROFILES_HEIGHT_IN = 4.75
REFERENCE_COMPARISON_HEIGHT_IN = 7.85
TOPOLOGY_FAMILY_HEIGHT_IN = 4.70


def _positive_fit_input(gfile: Geqdsk) -> Geqdsk:
    """Prepare a g-file in the single convention used by the fit comparison.

    The fit study treats plasma current as a positive scalar.  The complete
    flux-dependent payload is reflected when needed so that the axis-to-LCFS
    interval is positive as well.  The normalized flux map is invariant under
    this reflection, while the resulting compact record has positive ``Ip`` and
    ``psi_lcfs`` diagnostics in its axis-zero gauge.
    """

    span = float(gfile.psi_bound - gfile.psi_axis)
    if span == 0.0:
        raise ValueError("cannot canonicalise a zero g-file flux span")
    sign = 1.0 if span > 0.0 else -1.0
    return replace(
        gfile,
        Ip=abs(float(gfile.Ip)),
        psi_axis=float(gfile.psi_axis * sign),
        psi_bound=float(gfile.psi_bound * sign),
        psi=np.asarray(gfile.psi) * sign,
        P_psi=np.asarray(gfile.P_psi) * sign,
        FF_psi=np.asarray(gfile.FF_psi) * sign,
    )


def _active_interior_fit_parameter_count(budget: dict, topology: dict) -> int:
    """Count the interior geometry and radial-flux degrees of freedom.

    The comparison holds LCFS parameterization fixed, so its free boundary values
    are deliberately excluded from the labelled precision budget.
    """

    config = fitting_settings(topology, K_max=None, **budget)
    _, _, flux, family, _ = flux_setup(np.array((0.0, 1.0)), config.counts, config.powers)
    return int(family.size + flux.shape[1])


def _matches_reference_budget(candidate: Equilibrium, budget: dict, topology: dict) -> bool:
    """Reject cached fits that use a different LCFS or interior configuration."""

    config = fitting_settings(topology, K_max=None, **budget)
    return (
        candidate.Ip > 0.0
        and candidate.psi[-1] > 0.0
        and candidate.Nr == budget["Nr"]
        and candidate.Nt == budget["Nt"]
        and candidate.c_lcfs.size == budget["c_order"] + 1
        and candidate.s_lcfs.size == budget["s_order"]
        and candidate.h_coeffs.size == config.counts[0]
        and candidate.v_coeffs.size == config.counts[1]
        and candidate.kappa_coeffs.size == config.counts[2]
        and candidate.c_coeffs.shape[1] == max(config.counts)
        and candidate.s_coeffs.shape[1] == max(config.counts)
    )


def _close_curve(r: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.r_[r, r[0]], np.r_[z, z[0]]


def _horizontal_boundary_envelope(
    curves: np.ndarray, z_grid: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return inboard and outboard LCFS envelopes on horizontal cuts."""

    left = np.full((curves.shape[0], z_grid.size), np.nan)
    right = np.full_like(left, np.nan)
    for curve_index, curve in enumerate(curves):
        r, z = _close_curve(curve[:, 0], curve[:, 1])
        z0, z1 = z[:-1], z[1:]
        r0, r1 = r[:-1], r[1:]
        nonhorizontal = z0 != z1
        for level_index, level in enumerate(z_grid):
            crosses = nonhorizontal & (
                ((z0 <= level) & (level < z1))
                | ((z1 <= level) & (level < z0))
            )
            if not np.any(crosses):
                continue
            fraction = (level - z0[crosses]) / (z1[crosses] - z0[crosses])
            intersections = r0[crosses] + fraction * (r1[crosses] - r0[crosses])
            left[curve_index, level_index] = float(np.min(intersections))
            right[curve_index, level_index] = float(np.max(intersections))

    def across_family(values: np.ndarray, reducer: object) -> np.ndarray:
        result = np.full(z_grid.size, np.nan)
        valid = np.any(np.isfinite(values), axis=0)
        result[valid] = reducer(values[:, valid], axis=0)
        return result

    return (
        across_family(left, np.nanmin),
        across_family(left, np.nanmax),
        across_family(right, np.nanmin),
        across_family(right, np.nanmax),
    )


def _load_device_coverage() -> list[dict[str, float | str]]:
    devices: list[dict[str, float | str]] = []
    case_count = 0
    for path in sorted(DEVICE_MANIFESTS.glob("*.manifest.json")):
        manifest = json.loads(path.read_text(encoding="utf-8"))
        cases = manifest["cases"]
        case_count += len(cases)
        boundaries = [case["reconstruction"]["boundary"] for case in cases]
        targets = [case["reconstruction"]["targets"] for case in cases]
        major_radius = float(np.median([item["R0"] for item in boundaries]))
        minor_radius = float(np.median([item["a"] for item in boundaries]))
        devices.append(
            {
                "name": manifest["device"],
                "R0": major_radius,
                "a": minor_radius,
                "B0": float(np.median([abs(item["B0"]) for item in boundaries])),
                "Ip": float(
                    np.median([abs(item["Ip"]) for item in targets]) / 1.0e6
                ),
                "aspect": major_radius / minor_radius,
            }
        )
    if len(devices) != 267 or case_count != 13_291:
        raise RuntimeError(
            f"Expected 267 devices and 13,291 cases, found {len(devices)} "
            f"devices and {case_count} cases"
        )
    return devices


def _marker_area(aspect_ratio: np.ndarray) -> np.ndarray:
    """Map aspect ratio directly to marker area in points squared."""

    return 14.0 * np.asarray(aspect_ratio, dtype=float)


def build_device_coverage() -> None:
    devices = _load_device_coverage()
    r0 = np.asarray([item["R0"] for item in devices], dtype=float)
    a = np.asarray([item["a"] for item in devices], dtype=float)
    b0 = np.asarray([item["B0"] for item in devices], dtype=float)
    ip = np.asarray([item["Ip"] for item in devices], dtype=float)
    spherical = np.asarray([item["aspect"] for item in devices], dtype=float) < 2.0

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(DOUBLE_FIGURE_WIDTH_IN, DEVICE_COVERAGE_HEIGHT_IN),
    )
    fig.subplots_adjust(
        **physical_layout(
            DOUBLE_FIGURE_WIDTH_IN,
            DEVICE_COVERAGE_HEIGHT_IN,
            left_in=0.92,
            right_in=0.18,
            top_extra_in=0.20,
        ),
        wspace=0.25,
    )

    conventional = ~spherical
    encodings = (
        (conventional, "conventional", BLUE, "o"),
        (spherical, "spherical", YELLOW, "o"),
    )
    for axis, values, ylabel, title in (
        (
            axes[0],
            b0,
            r"Median toroidal field $|B_0|$ [T]",
            "(a) Toroidal field",
        ),
        (
            axes[1],
            ip,
            r"Median plasma current $|I_p|$ [MA]",
            "(b) Plasma current",
        ),
    ):
        for mask, _name, color, marker in encodings:
            axis.scatter(
                r0[mask],
                values[mask],
                s=_marker_area(r0[mask] / a[mask]),
                marker=marker,
                facecolor=color,
                edgecolor="white",
                linewidth=0.35,
                alpha=0.80,
                zorder=2,
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlim(0.05, 12.0)
        axis.margins(y=0.10)
        axis.set_xlabel(r"Median major radius $R_0$ [m]")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        quantitative_grid(axis)

    class_handles = [
        Line2D(
            [],
            [],
            marker=marker,
            linestyle="none",
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.35,
            markersize=5.2,
            label=label,
        )
        for label, color, marker in (
            (rf"conventional, $R_0/a\geq2$ ($N={int(np.count_nonzero(conventional))}$)", BLUE, "o"),
            (rf"spherical, $R_0/a<2$ ($N={int(np.count_nonzero(spherical))}$)", YELLOW, "o"),
        )
    ]
    fig.legend(
        handles=class_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, top_legend_y(DEVICE_COVERAGE_HEIGHT_IN)),
        ncol=2,
        handletextpad=0.45,
        columnspacing=1.4,
        borderaxespad=0.0,
        fontsize=CONFIG.super_legend_font_size,
    )
    save_figure(fig, HERE / "device_coverage")


def build_tcv_geometry() -> None:
    manifest = json.loads((TCV_BATCH / "manifest.json").read_text(encoding="utf-8"))
    ordered = sorted(manifest["cases"], key=lambda case: case["sample"]["delta"])
    selected = [ordered[index] for index in (0, 16, 33, 49)]

    fig, axes = plt.subplots(
        1,
        4,
        sharex=False,
        sharey=True,
        figsize=(TCV_GEOMETRY_WIDTH_IN, TCV_GEOMETRY_HEIGHT_IN),
    )
    # Preserve the negative-to-positive triangularity scan in reading order.
    # A 6.25 in canvas gives the four equal-aspect TCV panels enough room for
    # their titles and ticks without stretching the physical geometry.
    fig.subplots_adjust(
        left=0.62 / TCV_GEOMETRY_WIDTH_IN,
        right=1.0 - 0.10 / TCV_GEOMETRY_WIDTH_IN,
        bottom=0.48 / TCV_GEOMETRY_HEIGHT_IN,
        top=1.0 - 0.44 / TCV_GEOMETRY_HEIGHT_IN,
        wspace=0.09,
    )

    for index, (axis, case) in enumerate(zip(axes, selected, strict=True)):
        state = Equilibrium.load_json(
            TCV_BATCH / case["compact"]
        ).replace(Nt=512)
        r, z = _close_curve(np.asarray(state.R_lcfs), np.asarray(state.Z_lcfs))
        delta = case["sample"]["delta"]
        color = BLUE if delta < 0.0 else RED
        axis.plot(r, z, color=color, linewidth=1.25)
        title = axis.set_title(
            rf"$\delta={delta:+.2f}$" + "\n"
            rf"$\kappa={case['sample']['kappa']:.2f}$",
            pad=2.0,
        )
        title.set_linespacing(1.0)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlim(0.55, 1.20)
        axis.set_xticks((0.6, 1.1))
        axis.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        axis.set_xlabel(r"$R$ [m]")
        axis.grid(False)
    axes[0].set_ylabel(r"$Z$ [m]")
    save_figure(fig, HERE / "tcv_geometry")


def build_start_profiles() -> None:
    statistics = json.loads(
        (FIGURE_DATA / "figure_statistics_ch3.json").read_text(encoding="utf-8")
    )
    cases = statistics["start_profiles"]["cases"]

    fig, axes = plt.subplots(
        2,
        2,
        sharex="col",
        figsize=(DOUBLE_FIGURE_WIDTH_IN, START_PROFILES_HEIGHT_IN),
    )
    fig.subplots_adjust(
        **physical_layout(
            DOUBLE_FIGURE_WIDTH_IN,
            START_PROFILES_HEIGHT_IN,
            left_in=1.08,
            right_in=0.22,
            top_extra_in=0.08,
        ),
        wspace=0.26,
        hspace=0.10,
    )

    fields = (
        ("P", 1.0e-3, r"$P$ [kPa]"),
        (
            "FF_psi",
            1.0,
            r"$FF_\psi$ [T]",
        ),
        ("P_psi", 1.0e-6, r"$P_\psi$ [MA/m³]"),
        ("q", 1.0, r"$q$"),
    )
    colors = (BLUE, YELLOW, RED, GREEN)
    handles: list[Line2D] = []

    for case, color, line_style in zip(cases, colors, LINE_STYLES, strict=True):
        state = Equilibrium.load_json(
            START_BATCH / "cases" / f"{case['case_id']}.compact.json",
        ).replace(Nt=case["display_Nt"])
        sample = case["sample"]
        handles.append(
            Line2D(
                [],
                [],
                color=color,
                linestyle=line_style,
                linewidth=1.25,
                label=(
                    rf"$\beta_t={100.0 * case['beta_actual']:.1f}\%$, "
                    rf"$B_t={sample['Bt_T']:.3f}$ T, "
                    rf"$I_p={sample['Ip_A'] / 1000.0:.0f}$ kA"
                ),
            )
        )
        for axis, (field, scale, ylabel) in zip(axes.flat, fields, strict=True):
            axis.plot(
                np.asarray(state.psin),
                scale * np.asarray(getattr(state, field)),
                color=color,
                linestyle=line_style,
                linewidth=1.25,
            )
            axis.set_xlim(0.0, 1.0)
            axis.set_ylabel(ylabel)
            quantitative_grid(axis)

    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.55, top_legend_y(START_PROFILES_HEIGHT_IN)),
        ncol=2,
        columnspacing=1.35,
        handlelength=2.4,
        handletextpad=0.55,
        labelspacing=0.40,
        fontsize=CONFIG.super_legend_font_size,
    )
    for column in range(2):
        axes[0, column].tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )
        axes[1, column].set_xticks((0.0, 0.5, 1.0))
        axes[1, column].tick_params(
            axis="x", which="both", bottom=True, top=False, labelbottom=True
        )
        axes[1, column].set_xlabel(r"$\hat{\psi}$")
    save_figure(fig, HERE / "start_profiles")


def _square_window(reference: object, pad: float = 1.06) -> tuple[float, float, float, float]:
    boundary = np.asarray(reference.boundary)
    center_r = 0.5 * (boundary[:, 0].min() + boundary[:, 0].max())
    center_z = 0.5 * (boundary[:, 1].min() + boundary[:, 1].max())
    half = 0.5 * pad * max(float(np.ptp(boundary[:, 0])), float(np.ptp(boundary[:, 1])))
    return (
        max(reference.Rmin, center_r - half),
        min(reference.Rmax, center_r + half),
        max(reference.Zmin, center_z - half),
        min(reference.Zmax, center_z + half),
    )


def _flux_on_rect(
    state: object, r_grid: np.ndarray, z_grid: np.ndarray
) -> np.ndarray:
    radius = np.asarray(state.R)
    height = np.asarray(state.Z)
    normalized_flux = np.asarray(state.psin)[:, None] * np.ones_like(radius)
    rr, zz = np.meshgrid(r_grid, z_grid, indexing="ij")
    points = np.c_[radius.ravel(), height.ravel()]
    boundary = np.c_[radius[-1], height[-1]]
    inside = (
        MplPath(boundary).contains_points(points)
        & (rr.min() <= points[:, 0])
        & (points[:, 0] <= rr.max())
        & (zz.min() <= points[:, 1])
        & (points[:, 1] <= zz.max())
    )
    return griddata(
        points[inside],
        normalized_flux.ravel()[inside],
        (rr, zz),
        method="linear",
    )


def _surface_coordinates(state: object, normalized_flux: float) -> np.ndarray:
    """Return one compact-equilibrium surface at a prescribed flux value."""

    psin = np.asarray(state.psin)
    radius = np.asarray(state.R)
    height = np.asarray(state.Z)
    return np.column_stack(
        (
            [np.interp(normalized_flux, psin, radius[:, index]) for index in range(radius.shape[1])],
            [np.interp(normalized_flux, psin, height[:, index]) for index in range(height.shape[1])],
        )
    )


def _curve_rms_distance(points: np.ndarray, curve: np.ndarray) -> float:
    """Return the RMS closest-segment distance from points to a closed curve."""

    if np.allclose(curve[0], curve[-1]):
        curve = curve[:-1]
    closed = np.vstack((curve, curve[0]))
    start = closed[:-1]
    segment = closed[1:] - start
    denominator = np.sum(segment * segment, axis=1)
    offset = points[:, None, :] - start[None, :, :]
    fraction = np.clip(
        np.sum(offset * segment[None, :, :], axis=2) / denominator[None, :],
        0.0,
        1.0,
    )
    closest = start[None, :, :] + fraction[:, :, None] * segment[None, :, :]
    squared_distance = np.sum((points[:, None, :] - closest) ** 2, axis=2)
    return float(np.sqrt(np.mean(np.min(squared_distance, axis=1))))


def build_reference_comparison() -> None:
    cases = (
        (
            "EAST",
            EAST_SOURCE,
        ),
        (
            "MAST-U",
            FREEGSNKE / "MASTU_freegsnke.geqdsk",
        ),
        (
            "ITER",
            FREEGSNKE / "ITER_freegsnke_repaired.geqdsk",
        ),
    )
    levels = (0.1, 0.3, 0.5, 0.7, 0.9)

    fig = plt.figure(figsize=(DOUBLE_FIGURE_WIDTH_IN, REFERENCE_COMPARISON_HEIGHT_IN))
    grid = fig.add_gridspec(
        4,
        3,
        **physical_layout(
            DOUBLE_FIGURE_WIDTH_IN,
            REFERENCE_COMPARISON_HEIGHT_IN,
        left_in=0.54,
        right_in=0.12,
        top_extra_in=0.15,
        bottom_extra_in=0.08,
        ),
        wspace=0.16,
        hspace=0.32,
    )
    geometry_axes = [fig.add_subplot(grid[0, column_index]) for column_index in range(3)]
    surface_error_axes = [fig.add_subplot(grid[1, column_index]) for column_index in range(3)]
    pressure_axes = [fig.add_subplot(grid[2, column_index]) for column_index in range(3)]
    q_axes = [fig.add_subplot(grid[3, column_index]) for column_index in range(3)]
    # Sparse source dashes identify the G-EQDSK curve without masking closely
    # overlapping fits.  Drawing the precision sweep from high to low leaves
    # lower-budget departures visible instead of covering every colour.
    reference_style = (0, (1.1, 5.0))
    fit_styles = (LINE_STYLES[1], LINE_STYLES[2], LINE_STYLES[3], "-")
    pressure_errors: list[np.ndarray] = []
    q_errors: list[np.ndarray] = []
    surface_errors: list[np.ndarray] = []

    legend_handles = (
        Line2D([], [], color="#1F2937", linestyle=reference_style, linewidth=1.00,
               label="reference G-EQDSK"),
        Line2D([], [], color=BLUE, linestyle=LINE_STYLES[1], linewidth=1.15, label="low (20)"),
        Line2D([], [], color=TEAL, linestyle=LINE_STYLES[2], linewidth=1.15, label="medium (40)"),
        Line2D([], [], color=ORANGE, linestyle=LINE_STYLES[3], linewidth=1.15, label="high (60)"),
        Line2D([], [], color=RED, linestyle="-", linewidth=1.55, marker="o",
               markerfacecolor="white", markeredgewidth=0.85, markersize=3.4,
               label="ultra (100)"),
    )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, top_legend_y(REFERENCE_COMPARISON_HEIGHT_IN)),
        ncol=5,
        handlelength=2.4,
        columnspacing=1.1,
        borderaxespad=0.0,
        fontsize=CONFIG.super_legend_font_size,
    )

    for column_index, (device, source_path) in enumerate(cases):
        reference = _positive_fit_input(Geqdsk.load_geqdsk(source_path))
        # Every column receives the same interior-precision sweep.  The LCFS
        # parameterization remains fixed, so profile and geometry differences
        # are attributable to the interior fit budget.
        budgets = (
            ("p20", 20, dict(Nr=32, Nt=32, c_order=15, s_order=15),
             dict(h_count=3, v_count=3, kappa_count=3, c_counts=(4, 1, 1), s_counts=(1, 1, 1)), BLUE),
            ("p40", 40, dict(Nr=32, Nt=32, c_order=15, s_order=15),
             dict(h_count=6, v_count=4, kappa_count=4, c_counts=(6, 4, 3), s_counts=(4, 3, 2)), TEAL),
            ("p60", 60, dict(Nr=32, Nt=32, c_order=15, s_order=15),
             dict(h_count=8, v_count=5, kappa_count=5, c_counts=(8, 7, 6), s_counts=(6, 5, 4)), ORANGE),
            ("p100", 100, dict(Nr=32, Nt=32, c_order=15, s_order=15),
             dict(h_count=11, v_count=10, kappa_count=10, c_counts=(11, 10, 9), s_counts=(11, 10, 9)), RED),
        )
        fitted = []
        for tag, parameter_count, budget, topology, colour in budgets:
            actual_count = _active_interior_fit_parameter_count(budget, topology)
            if actual_count != parameter_count:
                raise RuntimeError(
                    f"{tag} labels {parameter_count} interior fit parameters, but config activates {actual_count}"
                )
            budget_path = FIGURE_DATA / f"crossdev_{device.lower()}_{tag}.compact.json"
            if budget_path.exists():
                candidate = Equilibrium.load_json(budget_path)
                if not _matches_reference_budget(candidate, budget, topology):
                    candidate = Equilibrium.from_geqdsk(
                        reference,
                        topology=topology,
                        **budget,
                    )
                    candidate.write_json(budget_path, header=f"{device}-{tag}-{parameter_count}", precision=12)
            else:
                candidate = Equilibrium.from_geqdsk(
                    reference,
                    topology=topology,
                    **budget,
                )
                candidate.write_json(budget_path, header=f"{device}-{tag}-{parameter_count}", precision=12)
            fitted.append((tag, parameter_count, candidate, colour))
        r_grid = np.linspace(reference.Rmin, reference.Rmax, 220)
        z_grid = np.linspace(reference.Zmin, reference.Zmax, 220)
        rr, zz = np.meshgrid(r_grid, z_grid, indexing="ij")
        spline = RectBivariateSpline(
            np.linspace(reference.Rmin, reference.Rmax, reference.NR),
            np.linspace(reference.Zmin, reference.Zmax, reference.NZ),
            reference.psi,
        )
        normalized_flux = (spline(r_grid, z_grid) - reference.psi_axis) / (
            reference.psi_bound - reference.psi_axis
        )
        boundary = np.asarray(reference.boundary)
        reference_minor_radius = 0.5 * float(np.ptp(boundary[:, 0]))
        region = (
            MplPath(boundary)
            .contains_points(np.c_[rr.ravel(), zz.ravel()])
            .reshape(rr.shape)
            & np.isfinite(normalized_flux)
            & (normalized_flux >= 0.0)
            & (normalized_flux <= 0.95)
        )
        contour_source = contour_generator(
            x=r_grid,
            y=z_grid,
            z=np.ma.masked_where(~region.T, normalized_flux.T),
        )
        reference_surfaces = {
            level: max(contour_source.lines(level), key=len) for level in levels
        }

        geometry = geometry_axes[column_index]
        dense_candidates = []
        for (tag, parameter_count, _candidate, colour), line_style in zip(
            fitted, fit_styles, strict=True
        ):
            dense_candidate = Equilibrium.load_json(
                FIGURE_DATA / f"crossdev_{device.lower()}_{tag}.compact.json",
            ).replace(Nt=1024)
            dense_candidates.append((tag, dense_candidate, colour, line_style))

        for tag, dense_candidate, colour, line_style in reversed(dense_candidates):
            candidate_flux = _flux_on_rect(dense_candidate, r_grid, z_grid)
            geometry.contour(
                r_grid, z_grid, np.where(region, candidate_flux, np.nan).T,
                levels=levels, colors=colour, linestyles=[line_style], linewidths=1.10, zorder=3,
            )
            comp_r, comp_z = _close_curve(
                np.asarray(dense_candidate.R_lcfs), np.asarray(dense_candidate.Z_lcfs)
            )
            ultra = tag == "p100"
            geometry.plot(
                comp_r,
                comp_z,
                color=colour,
                linestyle=line_style,
                linewidth=1.55 if ultra else 1.05,
                marker="o" if ultra else None,
                markerfacecolor="white" if ultra else None,
                markeredgewidth=0.85 if ultra else None,
                markersize=2.9 if ultra else None,
                markevery=64 if ultra else None,
                zorder=3,
            )
        # Overlay the source only as sparse dashes.  This retains a visible
        # reference at print scale without fabricating separation where curves
        # agree within a linewidth.
        geometry.contour(
            r_grid, z_grid, np.where(region, normalized_flux, np.nan).T,
            levels=levels, colors="#1F2937", linestyles=[reference_style], linewidths=0.80,
            zorder=4,
        )
        ref_r, ref_z = _close_curve(boundary[:, 0], boundary[:, 1])
        geometry.plot(
            ref_r,
            ref_z,
            color="#1F2937", linestyle=reference_style, linewidth=0.85, zorder=4,
        )
        geometry.plot(reference.Raxis, reference.Zaxis, "kx", markersize=4.5)
        r0, r1, z0, z1 = _square_window(reference)
        geometry.set_xlim(r0, r1)
        geometry.set_ylim(z0, z1)
        geometry.set_aspect("equal", adjustable="box")
        geometry.set_title(f"({chr(97 + column_index)}) {device}", loc="left", pad=4.0)
        geometry.set_xlabel(r"$R$ [m]")
        if column_index == 0:
            geometry.set_ylabel(r"$Z$ [m]")
        geometry.grid(False)

        surface_error = surface_error_axes[column_index]
        pressure = pressure_axes[column_index]
        safety_factor = q_axes[column_index]
        profile_x = np.linspace(0.0, 1.0, reference.P.size)
        reference_pressure = np.asarray(reference.P)
        pressure_scale = float(np.max(np.abs(reference_pressure)))
        reference_q = np.abs(np.asarray(reference.q))
        q_mask = profile_x <= 0.95
        for tag, dense_candidate, colour, line_style in dense_candidates:
            ultra = tag == "p100"
            geometric_error = np.asarray(
                [
                    _curve_rms_distance(
                        _surface_coordinates(dense_candidate, level),
                        reference_surfaces[level],
                    )
                    / reference_minor_radius
                    for level in levels
                ]
            )
            surface_errors.append(geometric_error)
            surface_error.plot(
                levels,
                np.where(geometric_error > 0.0, geometric_error, np.nan),
                color=colour,
                linestyle=line_style,
                linewidth=1.55 if ultra else 1.05,
                marker="o" if ultra else None,
                markerfacecolor="white" if ultra else None,
                markeredgewidth=0.85 if ultra else None,
                markersize=2.9 if ultra else None,
            )
            candidate_x = np.asarray(dense_candidate.psin)
            candidate_pressure = np.interp(
                profile_x, candidate_x, np.asarray(dense_candidate.P)
            )
            pressure_error = np.abs(candidate_pressure - reference_pressure) / pressure_scale
            pressure_errors.append(pressure_error)
            pressure.plot(profile_x, np.where(pressure_error > 0.0, pressure_error, np.nan),
                          color=colour, linestyle=line_style, linewidth=1.55 if ultra else 1.05,
                          marker="o" if ultra else None, markerfacecolor="white" if ultra else None,
                          markeredgewidth=0.85 if ultra else None, markersize=2.9 if ultra else None,
                          markevery=24 if ultra else None)
            candidate_q = np.interp(
                profile_x[q_mask], candidate_x, np.abs(np.asarray(dense_candidate.q))
            )
            q_error = np.abs(candidate_q - reference_q[q_mask]) / reference_q[q_mask]
            q_errors.append(q_error)
            safety_factor.plot(profile_x[q_mask], np.where(q_error > 0.0, q_error, np.nan),
                               color=colour, linestyle=line_style, linewidth=1.55 if ultra else 1.05,
                               marker="o" if ultra else None, markerfacecolor="white" if ultra else None,
                               markeredgewidth=0.85 if ultra else None, markersize=2.9 if ultra else None,
                               markevery=24 if ultra else None)
        for profile_axis, ylabel in (
            (surface_error, r"$d_{\mathrm{surf,RMS}}/a_{\mathrm{ref}}$"),
            (pressure, r"$|\Delta p|/p_{\mathrm{G\! -\! EQDSK}}(0)$"),
            (safety_factor, r"$|\Delta |q||/|q|_{\mathrm{G\! -\! EQDSK}}$"),
        ):
            profile_axis.set_xlim(0.0, 1.0)
            profile_axis.set_xlabel(r"$\hat{\psi}$")
            profile_axis.set_ylabel(ylabel if column_index == 0 else "")
            quantitative_grid(profile_axis)
            profile_axis.tick_params(top=True, right=True)
    surface_positive = np.concatenate([error[error > 0.0] for error in surface_errors])
    pressure_positive = np.concatenate([error[error > 0.0] for error in pressure_errors])
    q_positive = np.concatenate([error[error > 0.0] for error in q_errors])
    # The exact source endpoint can agree to machine precision.  Its isolated
    # roundoff tail should not consume the readable range of a log-error plot.
    surface_floor = float(np.quantile(surface_positive, 0.02))
    surface_limit = 1.15 * max(float(np.max(error)) for error in surface_errors)
    pressure_floor = float(np.quantile(pressure_positive, 0.02))
    pressure_limit = 1.15 * max(float(np.max(error)) for error in pressure_errors)
    q_floor = float(np.quantile(q_positive, 0.02))
    q_limit = 1.15 * max(float(np.max(error)) for error in q_errors)
    for axis in surface_error_axes:
        axis.set_ylim(surface_floor / 1.4, surface_limit)
        axis.set_yscale("log")
    for axis in pressure_axes:
        axis.set_ylim(pressure_floor / 1.4, pressure_limit)
        axis.set_yscale("log")
    for axis in q_axes:
        axis.set_ylim(q_floor / 1.4, q_limit)
        axis.set_yscale("log")

    save_figure(fig, HERE / "east_reference_comparison")


def build_same_topology_family() -> None:
    data = json.loads(
        (FIGURE_DATA / "east_same_topology_family.json").read_text(encoding="utf-8")
    )
    accepted_records = [
        record for record in data["records"] if record.get("accepted")
    ]
    if len(accepted_records) != 45:
        raise RuntimeError(
            f"Expected 45 accepted local-variation records, found "
            f"{len(accepted_records)}"
        )

    grouped_records: dict[tuple[float, float], list[dict[str, object]]] = {}
    for record in accepted_records:
        key = (float(record["dkappa"]), float(record["dR0_m"]))
        grouped_records.setdefault(key, []).append(record)
    if len(grouped_records) != 15 or any(
        len(records) != 3 for records in grouped_records.values()
    ):
        raise RuntimeError("Expected 15 geometry states with three source scalings each")

    invariant_fields = (
        "beta_t",
        "P_axis_Pa",
        "q_abs_min",
        "Ip_A",
        "kappa_measured",
        "R_center_m",
        "a_measured_m",
    )
    for records in grouped_records.values():
        for field in invariant_fields:
            values = np.asarray([float(record[field]) for record in records])
            if not np.allclose(values, values[0], rtol=1.0e-11, atol=1.0e-12):
                raise RuntimeError(
                    f"Common source scaling changed {field}; the records are no "
                    "longer duplicate physical states"
                )

    states = [
        min(records, key=lambda record: abs(float(record["dsource"])))
        for _, records in sorted(grouped_records.items())
    ]

    base = Equilibrium.load_json(FIGURE_DATA / "east.compact.json").replace(Nt=512)
    base_r = np.asarray(base.R_lcfs)
    base_z = np.asarray(base.Z_lcfs)

    fig = plt.figure(figsize=(DOUBLE_FIGURE_WIDTH_IN, TOPOLOGY_FAMILY_HEIGHT_IN))
    grid = fig.add_gridspec(
        1,
        2,
        **physical_layout(
            DOUBLE_FIGURE_WIDTH_IN,
            TOPOLOGY_FAMILY_HEIGHT_IN,
            left_in=0.80,
            right_in=0.95,
        ),
        width_ratios=(0.78, 1.45),
        wspace=0.24,
    )
    axes = (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]))

    geometry = axes[0]
    family_curves: list[np.ndarray] = []
    for record in states:
        r = base_r + float(record["dR0_m"])
        requested_kappa = float(record["kappa_required"])
        z = base.Z0 + (base_z - base.Z0) * requested_kappa / base.kappa_lcfs
        family_curves.append(np.c_[r, z])
    curves = np.asarray(family_curves)
    z_span = float(np.ptp(curves[:, :, 1]))
    z_grid = np.linspace(
        float(np.min(curves[:, :, 1])) + 1.0e-6 * z_span,
        float(np.max(curves[:, :, 1])) - 1.0e-6 * z_span,
        480,
    )
    left_min, left_max, right_min, right_max = _horizontal_boundary_envelope(
        curves, z_grid
    )
    geometry.fill_betweenx(
        z_grid,
        left_min,
        left_max,
        color=BLUE,
        alpha=0.22,
        linewidth=0.0,
    )
    geometry.fill_betweenx(
        z_grid,
        right_min,
        right_max,
        color=BLUE,
        alpha=0.22,
        linewidth=0.0,
    )
    seed_r, seed_z = _close_curve(base_r, base_z)
    geometry.plot(
        seed_r,
        seed_z,
        color=RED,
        linewidth=1.65,
    )
    geometry.set_xlim(
        float(np.min(curves[:, :, 0])) - 0.025,
        float(np.max(curves[:, :, 0])) + 0.025,
    )
    geometry.set_ylim(
        float(np.min(curves[:, :, 1])) - 0.025,
        float(np.max(curves[:, :, 1])) + 0.025,
    )
    geometry.set_aspect("equal", adjustable="box")
    geometry.set_xlabel(r"$R$ [m]")
    geometry.set_ylabel(r"$Z$ [m]")
    geometry.set_title("(a) Local LCFS family")
    geometry.legend(
        handles=(
            Line2D([], [], color=RED, linewidth=1.65, label="EAST seed"),
            Patch(
                facecolor=BLUE,
                edgecolor="none",
                alpha=0.22,
                label="15-state range",
            ),
        ),
        loc="center",
        fontsize=CONFIG.legend_font_size,
        handlelength=1.6,
    )
    geometry.grid(False)

    elongation = axes[1]
    pairs = sorted(
        {
            (float(record["kappa_required"]), float(record["kappa_measured"]))
            for record in states
        }
    )
    elongation.plot(
        [1.30, 1.90],
        [1.30, 1.90],
        color=GREY,
        linewidth=0.9,
        zorder=0,
        label="requested = measured",
    )
    elongation.plot(
        [pair[0] for pair in pairs],
        [pair[1] for pair in pairs],
        color=BLUE,
        marker="o",
        linestyle="none",
        markersize=4.2,
        label=r"realised $\kappa$",
    )
    elongation.set_xlabel(r"requested boundary elongation $\kappa$")
    elongation.set_ylabel(r"realised elongation $\kappa_{\rm LCFS}$", color=BLUE)
    elongation.tick_params(axis="y", colors=BLUE)
    elongation.set_title("(b) Stored-state response")
    quantitative_grid(elongation)

    separation = elongation.twinx()
    grouped_separation: dict[int, list[tuple[float, float]]] = {}
    for record in states:
        shift_cm = round(float(record["dR0_m"]) * 100.0)
        grouped_separation.setdefault(shift_cm, []).append(
            (
                float(record["kappa_required"]),
                100.0
                * (
                    float(record["R_axis_m"])
                    - float(record["boundary_center_r_m"])
                ),
            )
        )
    for (shift_cm, values), line_style, marker in zip(
        sorted(grouped_separation.items()), ("-", "--", ":"), ("o", "s", "^")
    ):
        values.sort()
        separation.plot(
            [item[0] for item in values],
            [item[1] for item in values],
            color=RED,
            linestyle=line_style,
            marker=marker,
            markersize=3.5,
            linewidth=1.05,
            label=(
                rf"$R_{{\rm axis}}-R_{{\rm c}}$, "
                rf"$\Delta R_0={shift_cm:+d}$ cm"
            ),
        )
    separation.set_ylabel(r"$R_{\rm axis}-R_{\rm c}$ [cm]", color=RED)
    separation.tick_params(axis="y", colors=RED)
    handles_left, labels_left = elongation.get_legend_handles_labels()
    handles_right, labels_right = separation.get_legend_handles_labels()
    elongation.legend(
        handles_left + handles_right,
        labels_left + labels_right,
        loc="lower right",
        fontsize=CONFIG.legend_font_size,
        handlelength=2.0,
    )
    save_figure(fig, HERE / "east_same_topology_family")


def main() -> None:
    apply_style()
    build_device_coverage()
    build_tcv_geometry()
    build_start_profiles()
    build_reference_comparison()
    build_same_topology_family()


if __name__ == "__main__":
    main()
