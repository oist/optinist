import os
import copy
import yaml
from cui_api.get_config_params import get_config_params


def get_snakemake_params(message_snakemake_params):
    # filepath = os.path.join(OPTINIST_DIR, 'config', f'snakemake.yaml')
    snakemake_params = get_config_params("snakemake")
    if message_snakemake_params != {}:
        snakemake_params = check_types(nest2dict(message_snakemake_params, snakemake_params))
    return snakemake_params

def get_nwbfile(message_nwb_params):
    # filepath = os.path.join(OPTINIST_DIR, 'config', f'nwb.yaml')
    nwb_params = get_config_params("nwb")
    if message_nwb_params != {}:
        nwb_pstsmd = check_types(nest2dict(message_nwb_params, nwb_params))
    return nwb_params


def get_algo_params(algo_name, message_params):
    # filepath = os.path.join(OPTINIST_DIR, 'config', f'{algo_name}.yaml')
    algo_params = get_config_params(algo_name)
    if message_params is not None:
        algo_params = check_types(nest2dict(message_params, algo_params))
    return algo_params