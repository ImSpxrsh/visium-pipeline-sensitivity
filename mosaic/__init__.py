"""MOSAIC: Multiverse Of Spatial Analysis — Instability & Concordance."""

__version__ = "0.1.0"


def _paper_style() -> None:
    """Larger default text so figures stay legible when shrunk to column width."""
    try:
        import matplotlib
    except ImportError:  # plotting is optional
        return
    matplotlib.rcParams.update({
        "font.size": 14, "axes.titlesize": 16, "axes.labelsize": 15,
        "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 12,
        "figure.titlesize": 17, "axes.titleweight": "normal",
    })


_paper_style()
