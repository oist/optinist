import os
import yaml
import sys
import copy
import inspect
import pickle
from snakemake import snakemake
from wrappers.data_wrapper import *
from wrappers import wrapper_dict
from .params import get_params
from .utils import nest2dict, check_types, algo_network


def create_snakemake_files(BASE_DIR, OPTINIST_DIR, message):

    # run snakemake
    snakemake(
        os.path.join(OPTINIST_DIR, 'Snakefile'),
        **snakemake_params
    )

    print("finish snakemake run")

    # Send message to the client
    for key, value in all_outputs.items():
        import time
        time.sleep(5)
        with open(key, "rb") as f:
            data = pickle.load(f)
            all_outputs[key]['info'] = data

    return all_outputs


def get_algo_params(OPTINIST_DIR, algo_path, node):
    filepath = os.path.join(OPTINIST_DIR, 'config', f'{algo_path.split("/")[-1]}.yaml')
    default_params = copy.deepcopy(get_params(filepath))
    params = default_params
    if 'param' in node['data'].keys() and node['data']['param'] is not None:
        params = nest2dict(node['data']['param'])
        params = check_types(params, default_params)
    return params
