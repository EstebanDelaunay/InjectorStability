import numpy as np
import matplotlib as mpl
from cycler import cycler
import matplotlib.pyplot as plt


def figSize(scaleW: float, scaleH: float = 1.0) -> tuple:

    figWidthPt = 450.0  # Get this from LaTeX using \the\textwidth
    inchesPerPt = 1.0 / 72.27  # Convert pt to inch

    figWidth = figWidthPt * inchesPerPt * scaleW

    figHeight = figWidthPt * inchesPerPt * scaleH

    return figWidth, figHeight


stylePDF = {
    "pgf.texsystem": "pdflatex",
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 10,  # LaTeX default is 10pt font.
    "font.size": 10,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.titlesize": 11,
    "figure.figsize": figSize(0.8, 0.6),  # default fig size of 0.9 textwidth
    "pgf.preamble": r"\usepackage[utf8x]{inputenc}",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "grid.linestyle": ":",
    "axes.grid": True,
    "figure.constrained_layout.use": True,
    "lines.linewidth": 1.,
    "legend.borderaxespad" : 0.2,
    "legend.frameon" : False,
    "legend.edgecolor" : "0.0",
    "legend.fancybox" : False,
}

cyclerColor = cycler(color=["darkblue", "orangered", "mediumseagreen", "royalblue", "gold", "forestgreen"]) 

cyclerColorLine = (cycler(color=["darkblue", "orangered", "mediumseagreen", "royalblue", "gold", "forestgreen"]) 
                    + cycler(linestyle=["-", "-.", ":","-", "-.", ":"]))


def ApplyStyle(custom_cycler = cyclerColor) -> None:
    mpl.rcParams.update(stylePDF)
    plt.rc("axes", prop_cycle=custom_cycler)


def saveFig(figure: plt.Figure, fileName: str) -> None:
    figure.savefig(f"python/images/{fileName}.pdf")
