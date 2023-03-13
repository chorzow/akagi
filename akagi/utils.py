import cooler
import pandas as pd
import os
from typing import Optional, Union
import numpy as np


def read_cooler(cool_path: Union[str, os.PathLike], resolution: Optional[int]) -> (cooler.Cooler, list):
    """Reads a .cool or .mcool file with given resolution.
    :param cool_path: input path.
    :param resolution: desired resolution as integer.
    :returns: a tuple with parsed .cool file and a list of chromosomes
              and their sizes retrieved from the same file."""

    try:
        if cool_path.endswith('.mcool'):
            cool_file = cooler.Cooler(f'{cool_path}::/resolutions/{resolution}')

        elif cool_path.endswith('.cool'):
            cool_file = cooler.Cooler(cool_path)

        else:
            raise ValueError('Provide a .cool or .mcool file')

        chromsizes = cool_file.chromsizes

    except FileNotFoundError:
        raise FileNotFoundError(f"No file found at {cool_path}")

    except Exception as e:
        raise e

    return cool_file, chromsizes


def read_inter(file: Union[str, os.PathLike]) -> pd.DataFrame:
    """Parse inter.txt file with Juicer Hi-C statistics.
    :param file: path to file.
    :returns: pandas DataFrame with processed file."""

    with open(file, 'r') as inf:
        lines = inf.readlines()

    lines = [x.strip() for x in lines if not x.startswith('WARN')]
    experiment_description = lines[0].split()
    genome = experiment_description[experiment_description.index('-g') + 1]
    lines = lines[1:]

    stats = dict()
    for line in lines:
        l = line.split(': ')
        stats[l[0]] = l[1].split(' ')
    for key, val in stats.items():
        stats[key] = list(filter(None, val))
    for key, val in stats.items():
        if key != "3' Bias (Long Range)" and key != 'Pair Type %(L-I-O-R)':
            stats[key] = [val[0]]
        else:
            stats[key] = [''.join(val)]

    stats = pd.DataFrame.from_dict(stats, orient='index', columns=[genome])
    stats[genome] = stats[genome].str.replace(',', '')
    stats[genome] = [int(i) if '%' not in i else i for i in stats[genome]]

    return stats


def poke_cool(input: Union[str, os.PathLike],
              n_diags: int,
              output: Union[str, os.PathLike],
              resolution: int=1000
              ) -> cooler.Cooler:
    """Poke out main diagonals of a contact matrix.
    :param input: path to a .cool file to operate with.
    :param n_diags: number of diagonals to poke out.
    :param resolution: specified resolution (if you are working with .mcool file).
    :param output: output file path.
    :returns: a cool file of specified resolution with n_diags diagonals nullified."""

    clr, csz = read_cooler(cool_path=input,
                           resolution=resolution)

    pixels = pd.DataFrame(clr.pixels()[:])
    bins = pd.DataFrame(clr.bins()[:])
    bins = bins[['chrom', 'start', 'end']]  # to remove previously calculated weights

    for i in range(n_diags):
        pixels['count'] = np.where(abs(pixels['bin1_id'] - pixels['bin2_id']) == i, 0, pixels['count'])
    # remove the last row and column (they are always masked in coolers we get)
    pixels['count'] = np.where(pixels['bin1_id'] == max(pixels['bin1_id']), 0, pixels['count'])
    pixels['count'] = np.where(pixels['bin2_id'] == max(pixels['bin2_id']), 0, pixels['count'])

    try:
        cooler.create_cooler(output, bins, pixels)
    except Exception as e:
        raise e
    c_poked = cooler.Cooler(output)

    try:  # re-balance
        with c_poked.open('r+') as f:
            f['bins'].create_dataset(
                'weight',
                data=cooler.balance_cooler(c_poked)[0],
                compression='gzip',
                compression_opts=6
            )
    except ValueError as e:
        print(e)

    return c_poked
