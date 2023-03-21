import os

import cooler
import cooltools
import hicstraw
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from typing import Optional, List


def plot_cooler(clr: cooler.Cooler,
                region: Optional[List[int]],
                vmin=Optional[float],
                vmax=Optional[float],
                balance: bool = True,
                cmap: str = 'Blues',
                xlabel=Optional[str],
                ylabel=Optional[str],
                **kwargs
                ) -> plt.Axes:
    """Plots log-normalized contact matrix from a cooler.Cooler object.
    :arg clr: object of type cooler.Cooler.
    :arg region: list of what region of the matrix to show IN BINS (not genomic coordinates). Format: [start, stop].
    Default is whole genome (can take a long time for large genomes).
    :arg vmin: min value for log-normalization.
    :arg vmax: max value for log-normalization.
    :arg balance: whether matrix should be balanced. Default is True.
    :arg cmap: name of the colormap to use. Only matplotlib default cmaps are supported for now.
    :arg xlabel: xlabel.
    :arg ylabel: ylabel.
    :returns: a figure of a contact matrix extracted from clr."""

    resolution = clr.binsize
    nbins = len(clr.bins()[:])

    f, ax = plt.subplots(**kwargs)

    norm = LogNorm(vmax=vmax, vmin=vmin)

    if not region:
        im = ax.matshow(
            clr.matrix(balance=balance)[:],
            norm=norm,
            cmap=cmap
        )
    else:
        try:
            im = ax.matshow(clr.matrix()[slice(region[0], region[1])])
        except IndexError:
            raise IndexError('Wrong regions format. Accepted format: [start, stop] where start and stop are integers.')
        except Exception as e:
            raise e

    plt.axis([0, nbins, nbins, 0])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax, label='corrected frequencies')

    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.xaxis.set_visible(False)

    return ax


def plot_hic(hic: hicstraw.HiCFile,
             vmin=Optional[float],
             vmax=Optional[float],
             ):
    pass
