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

def scaling_curve(clr: cooler.Cooler,
        smooth: Optional[bool] = False,
        aggregate_smoothed: Optional[bool] = False,
        nproc: Optional[int] = 4,
        ax: Optional[plt.Axes] = None,
        **kwargs) -> plt.Axes:
    """
    Plots a scaling curve in log-log scale with its derivative in semi-log scale.
    Data used for plotting can be accessed via attributes of the function:
    scaling_curve.cvd = contact-versus-distance dataframe;
    scaling_curve.der = derivative of a P(s) curve.
    :param clr: cooler.Cooler object to calculate curves from.
    :param smooth: whether to perform smoothing of the P(s) curve. Default is False.
    :param aggregate_smoothed: whether to aggregate curves. Default is False.
    :param nproc: number of processors to use. Default is 4.
    :returns: a plt.Axes object with plotted P(s) curve and derivative.
    """

    resolution = clr.binsize

    cvd = cooltools.expected_cis(clr=clr,
                                 smooth=smooth,
                                 aggregate_smoothed=aggregate_smoothed,
                                 nproc=nproc)
    
    cvd['s_bp'] = cvd['dist'] * resolution

    if ax is None:
        _, ax = plt.subplots(**kwargs)

    if smooth and aggregate_smoothed:
        cvd['balanced.avg.smoothed'].loc[cvd['dist'] < 2] = np.nan
        cvd = cvd.drop_duplicates(subset=['dist'])
        cvd = cvd[3:]  # removes peak values at the beginning for prokaryotes

        ys = cvd['balanced.avg.smoothed.agg']
    
    if smooth and not aggregate_smoothed:
        ys = cvd['balanced.avg.smoothed']
    
    else:
        ys = cvd['balanced.avg']

    xs = cvd['s_bp']

    ax.loglog(xs, ys)

    ax.set_ylabel('IC contact frequency')

    scaling_curve.data = cvd

    return ax


def scaling_curve_der(
        clr: cooler.Cooler,
        c: Optional[str] = None,
        nproc: Optional[int] = 4,
        label: Optional[str] = None,
        ax: Optional[plt.Axes] = None,
        **kwargs) -> plt.Axes:
    """
    Plots a scaling curve derivative in semi-log scale.
    Data used for plotting can be accessed via attributes of the function:
    scaling_curve.data = derivative of a P(s) curve.
    :param clr: cooler.Cooler object to calculate curves from.
    :param c: color of the line.
    :param nproc: number of processors to use. Default is 4.
    :returns: a plt.Axes object with plotted P(s) curve and derivative.
    """

    resolution = clr.binsize

    cvd = cooltools.expected_cis(clr=clr,
                                 smooth=True,
                                 aggregate_smoothed=True,
                                 nproc=nproc
                                 )
    
    cvd['s_bp'] = cvd['dist'] * resolution

    if ax is None:
        _, ax = plt.subplots(**kwargs)

    cvd['balanced.avg.smoothed'].loc[cvd['dist'] < 2] = np.nan
    cvd = cvd.drop_duplicates(subset=['dist'])


    cvd['der'] = np.gradient(np.log(cvd['balanced.avg.smoothed.agg']),
                        np.log(cvd['s_bp'])
                        )
    cvd = cvd[3:]  # removes peak values at the beginning for prokaryotes

    xs = cvd['s_bp']
    ys = cvd['der']


    ax.semilogx(xs, ys, color=c, label=label)

    ax.set(xlabel='separation, bp',
           ylabel='slope'
           )

    scaling_curve_der.data = cvd

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
