"""Shared publication style for the VEQDB manuscript figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler


class _FigureConfig:
    """Typography and colour constants fixed at the manuscript figure scale."""

    font_family = "DejaVu Sans"
    font_size = 9.25
    legend_font_size = 8.25
    super_legend_font_size = 9.25
    line_width = 1.15
    marker_size = 4.2
    save_dpi = 450
    grid_color = "#CBD5E1"
    grid_alpha = 0.72


CONFIG = _FigureConfig()

SINGLE_FIGURE_WIDTH_IN = 3.5
DOUBLE_FIGURE_WIDTH_IN = 7.2

# Layout is specified in physical inches rather than figure fractions.  This
# keeps the lower label reserve and the outer top/bottom whitespace identical
# after the common LaTeX scaling, independent of the figure height.
OUTER_VERTICAL_PAD_IN = 0.18
AXES_BOTTOM_MARGIN_IN = 0.86
AXES_TOP_MARGIN_IN = 0.62

NAVY = "#102A43"
BLUE = "#4A7BD0"
TEAL = "#58B368"
GREEN = TEAL
ORANGE = "#E7B23C"
RED = "#E45C3A"
YELLOW = "#E7B23C"
GREY = "#9A9A9A"
# Keep every precision distinguishable where fits nearly overlap the reference.
# In particular, ultra is dash-dotted rather than another solid curve.
LINE_STYLES = ((0, (8, 1.8, 1.5, 1.8)), (0, (5, 2.5)), (0, (2, 1.6)), ":")


def physical_layout(
    width: float,
    height: float,
    *,
    left_in: float,
    right_in: float,
    top_extra_in: float = 0.0,
    bottom_extra_in: float = 0.0,
) -> dict[str, float]:
    """Convert fixed physical margins to Matplotlib figure fractions."""

    return {
        "left": left_in / width,
        "right": 1.0 - right_in / width,
        "bottom": (AXES_BOTTOM_MARGIN_IN + bottom_extra_in) / height,
        "top": 1.0 - (AXES_TOP_MARGIN_IN + top_extra_in) / height,
    }


def top_legend_y(height: float) -> float:
    """Anchor a top legend below the common outer edge padding."""

    return 1.0 - OUTER_VERTICAL_PAD_IN / height


def apply_style() -> None:
    """Install a high-contrast, colour-blind-safe journal figure style."""

    mpl.rcParams.update(
        {
            "font.family": CONFIG.font_family,
            "font.size": CONFIG.font_size,
            "mathtext.fontset": "dejavusans",
            "axes.labelsize": CONFIG.font_size,
            "axes.titlesize": CONFIG.font_size,
            "axes.linewidth": 0.8,
            "axes.edgecolor": "#334E68",
            "axes.labelcolor": "#102A43",
            "axes.titleweight": "bold",
            "axes.prop_cycle": cycler(color=(BLUE, TEAL, ORANGE, RED, YELLOW)),
            "axes.axisbelow": True,
            "axes.unicode_minus": False,
            "xtick.labelsize": CONFIG.font_size,
            "ytick.labelsize": CONFIG.font_size,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 3.6,
            "ytick.major.size": 3.6,
            "xtick.minor.size": 2.1,
            "ytick.minor.size": 2.1,
            "legend.fontsize": CONFIG.legend_font_size,
            "legend.frameon": False,
            "lines.linewidth": CONFIG.line_width,
            "lines.markersize": CONFIG.marker_size,
            "lines.solid_capstyle": "round",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": CONFIG.save_dpi,
            "savefig.transparent": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def quantitative_grid(axis: plt.Axes) -> None:
    """Use a light major grid that supports numerical reading."""

    axis.grid(
        True,
        which="major",
        color=CONFIG.grid_color,
        alpha=CONFIG.grid_alpha,
        linewidth=0.55,
    )


def save_figure(figure: plt.Figure, stem: Path) -> None:
    """Export publication raster and vector variants at the native canvas size."""

    stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(stem.with_suffix(".png"), dpi=CONFIG.save_dpi)
    figure.savefig(stem.with_suffix(".pdf"))
    plt.close(figure)
