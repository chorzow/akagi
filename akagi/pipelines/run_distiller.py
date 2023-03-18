import yaml
import click
import subprocess
import time
import glob
import os
import logging

LOG_FORMAT = "%(asctime)s %(levelname)s: %(message)s"
logging.basicConfig(filename='run_distiller.py.log',
                    level=logging.DEBUG,
                    format=LOG_FORMAT)
logger = logging.getLogger()


@click.command()
@click.argument('--sample', '-s', help='Sample name.')
@click.argument('-i', '--input_path', help='Input folder.')
@click.argument('-d', '--distiller_path', help='path to distiller.')
def run_distiller(input_path, distiller_path, sample):
    t = time.process_time()
    path = os.path.join(input_path, sample)
    config_template = {
        'input': {
            'raw_reads_paths': {
                'sample': {
                    'lane1':
                        glob.glob(f"{path}/input/*.fastq")

                }
            },
            'library_groups': {
                "all": [
                    sample
                ]
            },
            'genome': {
                'assembly_name': 'genome.fasta',
                'bwa_index_wildcard_path': f"{path}/input/genome.fasta*",
                'chrom_sizes_path': f"{path}/input/genome.chrom.sizes"
            }
        },
        'do_fastqc': False,
        'map': {
            'chunksize': 10000000,
            'mapping_options': '-SP5M',
            'trim_options': '',
            'use_custom_split': 'true'
        },
        'parse': {
            'make_pairsam': True,
            'drop_seq': False,
            'drop_readid': True,
            'keep_unparsed_bams': False,
            'parsing_options': '--add-columns mapq'
        },
        'dedup': {
            'max_mismatch_bp': 3
        },
        'bin': {
            'resolutions': [
                500000,
                250000,
                100000,
                50000,
                25000,
                10000,
                5000,
                4000,
                3000,
                2000,
                1000,
                500,
                100
            ],
            'balance': 'true',
            'balance_options': '',
            'filters': {
                'no_filter': "",
                'mapq_30': "(mapq1>=30) and (mapq2>=30)"
            }
        },
        'output': {
            'dirs': {
                'processed_fastqs': f'{path}/output/processed_fastqs/',
                'mapped_parsed_sorted_chunks': f'{path}/output/mapped_parsed_sorted_chunks/',
                'fastqc': f'{path}/output/fastqc/',
                'pairs_library': f'{path}/output/pairs_library/',
                'coolers_library': f'{path}/output/coolers_library/',
                'coolers_library_group': f'{path}/output/coolers_library_group/',
                'stats_library_group': f'{path}/output/stats_library_group/'
            }
        }
    }
    output_path = f"{path}/output/distiller_config.yml"
    with open(output_path, 'w') as yamlfile:
        logger.log(logging.INFO, f"writing config file to {output_path}")
        try:
            yaml.dump(config_template, yamlfile)
        except FileNotFoundError as e:
            logger.log(logging.ERROR, e)
            raise e

    distiller_command = ['nextflow', 'run', distiller_path, '-params-file', output_path]
    logger.log(logging.INFO,
               f"running distiller with the following command {distiller_command} and config file {config_template}")
    with open(f'{path}/output/distiller.log', 'w') as distiller_log:
        try:
            subprocess.Popen(distiller_command, stdout=distiller_log, stderr=distiller_log)
        except Exception as e:
            logger.log(logging.WARN, e)
        else:
            logger.log(logging.INFO, f'Process completed in {time.process_time() - t}')


if __name__ == "__main__":
    run_distiller()
