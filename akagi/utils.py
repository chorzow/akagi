import os
import subprocess
from typing import Optional, Union, List
import itertools

import pandas as pd
import numpy as np
from scipy import sparse

import cooler
import hicstraw


def cool_or_mcool(path: Union[str, os.PathLike]) -> str:
    """Checks if file is a single-resolution, multi-resolution cooler or none of them.
    :param path: path to file.
    :returns: 'sr' if single-resolution, 'mr' if multi-resolution."""

    if cooler.fileops.is_cooler(path):
        return 'sr'
    elif cooler.fileops.is_multires_file(path):
        return 'mr'
    else:
        raise ValueError('Provide a .cool or .mcool file')


def read_cooler(path: Union[str, os.PathLike], resolution: Optional[int] = None) -> cooler.Cooler:
    """Reads a .cool or .mcool file with given resolution.
    :param cool_path: input path.
    :param resolution: desired resolution as integer. if not provided and file is .mcool, reads the highest resolution.
    :returns: class cooler.Cooler."""

    if cool_or_mcool(path) == 'sr':
        cool_file = cooler.Cooler(path)
    elif cool_or_mcool(path) == 'mr':

        if resolution is None:
            resolutions_list = [int(i.split('/')[-1]) for i in cooler.fileops.list_coolers(path)]
            resolutions_list.sort()
            resolution = resolutions_list[0]

        cool_file = cooler.Cooler(f'{path}::/resolutions/{resolution}')

    return cool_file


def read_hic(path: Union[str, os.PathLike],
             resolution: int,
             norm: str = 'KR',
             datatype: str = 'observed'
             ) -> Union[np.ndarray, dict]:
    """Returns a contact matrices for the specified chromosomes.
    :arg path: path to a .hic file.
    :arg resolution: desired resolution in bp.
    :arg norm: normalization to get. Default is 'KR' (Knight-Ruiz). Supported: NONE, VC, VC_SQRT, KR, SCALE.
    :arg datatype: type of data to get. Default is 'observed'. Supported: 'observed', 'oe' (observed/expected).
    """

    if norm not in {'NONE', 'VC', 'VC_SQRT', 'KR', 'SCALE'}:
        raise ValueError(f"Wrong normalization {norm} provided. Supported: NONE, VC, VC_SQRT, KR, SCALE.")

    hic = hicstraw.HiCFile(path)
    chromnames = [chrom.name for chrom in hic.getChromosomes() if chrom.name != 'All']
    chromlengths = [chrom.length for chrom in hic.getChromosomes() if chrom.name != 'All']

    if len(chromlengths) == len(chromnames) == 1:  # prokaryotes
        name = chromnames[0]
        length = chromlengths[0]
        try:
            mat = hic.getMatrixZoomData(name, name, datatype, norm, 'BP', resolution)
            np_matrix = mat.getRecordsAsMatrix(0, length, 0, length)
            return np_matrix
        except Exception as e:
            raise e

    # TODO: rework and test next section

    # if len(chromnames) > 1:  # eukaryotes
    #     print("More than one chromosome found. Perhaps it's gonna take a while. Yet perhaps not.")
    #     np_matrices = {}
    #     for chr1, chr2 in itertools.combinations_with_replacement(chromnames, 2):
    #         try:

    #             mat = hic.getMatrixZoomData(chr1, chr2, datatype, norm, 'BP', resolution)

    #             chr1_len = chromlengths[chromnames.index(chr1)]
    #             chr2_len = chromlengths[chromnames.index(chr2)]

    #             np_mat = mat.getRecordsAsMatrix(0, chr1_len, 0, chr2_len)
    #             np_matrices[f"{chr1}_{chr2}"] = np_mat

    #             return np_matrices

    #         except Exception as e:
    #             print(e)
    #             continue

    elif len(chromnames) == 0:
        raise IndexError('No chromosomes found in provided .hic file.')


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

    stats = {}
    for line in lines:
        line_splitted = line.split(': ')
        stats[line_splitted[0]] = line_splitted[1].split(' ')
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


def cooler_poke(path: Union[str, os.PathLike],
                n_diags: int,
                output: Union[str, os.PathLike],
                resolution: int = 1000
                ) -> cooler.Cooler:
    """Poke out main diagonals of a contact matrix.
    :param path: path to a .cool file to operate with.
    :param n_diags: number of diagonals to poke out.
    :param resolution: specified resolution (if you are working with .mcool file).
    :param output: output file path.
    :returns: a cool file of specified resolution with n_diags diagonals nullified."""

    clr = read_cooler(cool_path=path,
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


def cooler_diff(clr1_path: Union[str, os.PathLike],
                clr2_path: Union[str, os.PathLike],
                out_dir: Union[str, os.PathLike],
                out_name: str,
                resolutions: Optional[List[int]] = None
                ) -> None:
    """Subtracts two matrices, then saves output as .mcool file. Equivalent of clr1_path - clr2_path.
    Creates a new UNbalanced cooler. Resolutions are inferred from the clr1_path by default.
    :param clr1_path:
    :param map2_path:
    :param out_dir: path to store output at.
    :param out_name: prefix name of resulting file without file extension.
    :returns: None
    """

    clr1 = read_cooler(clr1_path)
    clr2 = read_cooler(clr2_path)

    if resolutions is None:

        print(f'No resolutions provided; fetching resolutions list from{clr1_path}...')
        resolutions = [int(i.split('/')[-1]) for i in cooler.fileops.list_coolers(clr1_path)]
        resolutions.sort()
        resolutions = resolutions[1:]

    pixels_colnames = clr1.pixels()[:].columns.tolist()
    bins_colnames = ['chrom', 'start', 'end']  # to remove previously calculated weights

    mat1 = clr1.matrix(balance=False)[:]
    mat2 = clr2.matrix(balance=False)[:]

    print(f'Subtracting {clr2_path.split("/")[-1]} raw signal from {clr1_path.split("/")[-1]}...')

    diff_mat = np.subtract(mat1, mat2)
    triu_mat = np.triu(diff_mat)

    print('Generating coo matrix...')

    coo_mat = sparse.coo_matrix(triu_mat)
    coo_mat_data = np.transpose([coo_mat.row, coo_mat.col, coo_mat.data])
    pixels = pd.DataFrame(coo_mat_data, columns=pixels_colnames)
    pixels = pixels.sort_values(['bin1_id', 'bin2_id']).drop_duplicates()

    print('Writing temporary .cool...')

    bins = pd.DataFrame(clr1.bins()[:])[bins_colnames]
    cool_out = os.path.join(out_dir, out_name + '.cool')

    try:
        cooler.create_cooler(cool_uri=cool_out, bins=bins, pixels=pixels, symmetric_upper=True)
    except Exception as e:  # TODO: add human-readable specific cases.
        raise e

    outfile = os.path.join(out_dir, out_name + '.mcool')

    print(f'Zoomifying {cool_out} for resolutions {resolutions} to {outfile}...')

    try:
        cooler.zoomify_cooler(cool_out, outfile=outfile, resolutions=resolutions, chunksize=1000000, nproc=4)
    except ValueError as e:
        if 'cannot be derived from the base' in str(e):
            raise ValueError('Each item in resolutions must be a multiplier of the base resolution.')
    print('Cleaning up...')

    subprocess.run(['rm', cool_out])

    print('Completed.')

    return None


class DistillerStats:
    """Operates with different sections of distiller's .stats file.
    Can be initialized with custom tables (just in case).
    Different attributes of this class represent different sections of .stats file."""
    def __init__(self, path: Union[str, os.PathLike], **kwargs):
        self.path = path
        self._general = kwargs.get('general', None)
        self._pair_types = kwargs.get('pair_types', None)
        self._cis_trans = kwargs.get('cis_trans', None)
        self._chrom_freq = kwargs.get('chrom_freq', None)
        self._dist_freq = kwargs.get('dist_freq', None)

    def read_section(self, section_name: Union[str, tuple]) -> List[str]:
        try:
            with open(self.path, 'r') as inf:
                return [i.strip() for i in inf.readlines() if i.startswith(section_name)]
        except Exception as e:
            raise e

    @property
    def general(self) -> Union[pd.DataFrame, None]:
        if self._general is None:
            section_name = ('total', 'cis', 'trans')

            try:
                lines = self.read_section(section_name)
                lines_filtered = [line for line in lines if not line.startswith(('cis_', 'trans_'))]
                lines_filtered = [line.split('\t') for line in lines_filtered]

                general_df = pd.DataFrame(lines_filtered).set_index(0)
                general_df.index.name = None
                general_df.columns = ['Number of reads']

                return general_df

            except Exception as e:
                raise e

        return None

    @property
    def pair_types(self) -> Union[pd.DataFrame, None]:
        if self._pair_types is None:
            section_name = 'pair_types'

            try:
                lines = self.read_section(section_name)
                lines = [line.removeprefix('pair_types/').split('\t') for line in lines]

                pair_types_df = pd.DataFrame(lines, columns=['Type', 'Number of pairs']).set_index('Type')

                return pair_types_df

            except Exception as e:
                raise e

        return None

    @property
    def cis_trans(self) -> Union[pd.DataFrame, None]:
        if self._cis_trans is None:
            section_name = ('cis_', 'trans_')

            try:
                lines = self.read_section(section_name)
                lines = [line.split('\t') for line in lines]

                cis_trans_df = pd.DataFrame(lines)

                type_distance = cis_trans_df[0].str.split('_')

                cis_trans_df['Type'] = type_distance.str[0]
                cis_trans_df['Distance'] = type_distance.str[1]

                cis_trans_df = cis_trans_df.rename(columns={1: 'Number of contacts'})
                cis_trans_df = cis_trans_df.set_index(['Type', 'Distance']).drop(0, axis=1)

                return cis_trans_df

            except Exception as e:
                raise e

        return None

    @property
    def chrom_freq(self) -> Union[pd.DataFrame, None]:
        if self._chrom_freq is None:
            section_name = 'chrom_freq'

            try:
                lines = [line.removeprefix('chrom_freq/').split('\t') for line in self.read_section(section_name)]

                chrom_freq_df = pd.DataFrame(lines)

                chrom_freq_df[['chrA', 'chrB']] = chrom_freq_df[0].str.split('/', expand=True)

                chrom_freq_df = chrom_freq_df.rename(columns={1: 'Number of contacts'})
                chrom_freq_df = chrom_freq_df[['chrA', 'chrB', 'Number of contacts']]

                return chrom_freq_df

            except Exception as e:
                raise e

        return None

    @property
    def dist_freq(self) -> Union[pd.DataFrame, None]:
        if self._dist_freq is None:
            section_name = 'dist_freq'

            try:
                lines = [line.removeprefix('dist_freq/').split('\t') for line in self.read_section(section_name)]

                dist_freq_df = pd.DataFrame(lines)

                dist_freq_df[['Region', 'Contact_type']] = dist_freq_df[0].str.split('/', expand=True)

                dist_freq_df = dist_freq_df.rename(columns={1: 'Number of contacts'})
                dist_freq_df = dist_freq_df[['Region', 'Contact_type', 'Number of contacts']]

                return dist_freq_df

            except Exception as e:
                raise e

        return None
