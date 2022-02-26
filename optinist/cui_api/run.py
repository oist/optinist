from snakemake import snakemake

from cui_api.const import OPTINIST_DIR
from cui_api.utils import join_file_path


def run_snakemake(snakemake_params):
    # run snakemake
    snakemake(
        join_file_path([OPTINIST_DIR, 'Snakefile']),
        **snakemake_params
    )