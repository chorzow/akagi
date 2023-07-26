import os
from typing import Optional, Tuple

import cooler
import cooltools
import hicstraw
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_cooler(clr: cooler.Cooler,
                region: Optional[Tuple] = None,
                vmin: Optional[float] = None,
                vmax: Optional[float] = None,
                balance: bool = True,
                cmap: str = 'Blues',
                ylabel: Optional[str] = None,
                ax: Optional[plt.Axes] = None,
                **kwargs
                ) -> plt.Axes:
    """Plots log-normalized contact matrix from a cooler.Cooler object.
    :arg clr: object of type cooler.Cooler.
    :arg region: tuple of format (chrom, start, end) where chrom is a string, start and end are int.
    Defines a genomic region to plot. Default is whole genome (can take a longer time for large genomes).
    :arg vmin: min value for log-normalization.
    :arg vmax: max value for log-normalization.
    :arg balance: whether matrix should be balanced. Default is True.
    :arg cmap: name of the colormap to use. Only matplotlib default cmaps are supported for now.
    :arg ylabel: ylabel.
    :arg ax: plt.Axes object to plot matrix on. If not specified, such object is created.
    :returns: plt.Axes object of a contact matrix extracted from clr."""

    if ax is None:
        _, ax = plt.subplots(**kwargs)

    resolution = clr.binsize
    norm = LogNorm(vmax=vmax, vmin=vmin)

    if region is None:
        plot_data = clr.matrix(balance=balance)[:]
    else:
        plot_data = clr.matrix(balance=balance).fetch(region)

    im = ax.matshow(plot_data, norm=norm, cmap=cmap)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax, label='corrected frequencies')

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x * resolution // 1000)} Kb'))

    if ylabel:
        ax.set_ylabel(ylabel)

    ax.xaxis.set_visible(False)

    return ax


def plot_hic(hic: np.ndarray,
             vmin: float = 0,
             vmax: Optional[float] = None,
             region: Optional[Tuple] = None,
             cmap: Optional[str] = 'Blues',
             ylabel: Optional[str] = None,
             ax: Optional[plt.Axes] = None,
             **kwargs
             ) -> plt.Axes:
    """Plots normalized contact matrix obtained from .hic file. Downstream processing of read_hic from utils.
    :arg hic: contact matrix.
    :arg region: tuple of format (chrom, start, end) where chrom is a string, start and end are int.
    Defines a genomic region to plot. Default is whole genome (can take a longer time for large genomes).
    :arg vmin: min value for log-normalization.
    :arg vmax: max value for log-normalization.
    :arg cmap: name of the colormap to use. Only matplotlib default cmaps are supported for now.
    :arg ylabel: ylabel.
    :arg ax: plt.Axes object to plot matrix on. If not specified, such object is created.
    :returns: plt.Axes object of a contact matrix extracted from hic."""

    if ax is None:
        _, ax = plt.subplots(**kwargs)
    pass
