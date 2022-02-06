import os
import yaml
from snakemake import snakemake
from .const import OPTINIST_DIR


def write_snakemake_config(flow_config):
    with open(os.path.join(OPTINIST_DIR, 'config.yaml'), "w") as f:
        yaml.dump(flow_config, f)


def run_snakemake(snakemake_params):
    # run snakemake
    snakemake(
        os.path.join(OPTINIST_DIR, 'Snakefile'),
        **snakemake_params
    )