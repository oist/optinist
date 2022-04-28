import os

_DEFAULT_DIR = '/tmp/optinist'
_ENV_DIR = os.environ.get('OPTINIST_DIR')

_DEFAULT_CLUSTER_INPUT_DIR = '/tmp/optinist_cluster/input'
_ENV_CLUSTER_INPUT_DIR = os.environ.get('OPTINIST_CLUSTER_INPUT_DIR')
_DEFAULT_CLUSTER_OUTPUT_DIR = '/tmp/optinist_cluster/output'
_ENV_CLUSTER_OUTPUT_DIR = os.environ.get('OPTINIST_CLUSTER_OUTPUT_DIR')

class DIRPATH:
    OPTINIST_DIR = _DEFAULT_DIR if _ENV_DIR is None else _ENV_DIR
    INPUT_DIR = f"{OPTINIST_DIR}/input"
    OUTPUT_DIR = f"{OPTINIST_DIR}/output"

    CLUSTER_INPUT_DIR = _DEFAULT_CLUSTER_INPUT_DIR\
        if _ENV_CLUSTER_INPUT_DIR is None else _ENV_CLUSTER_INPUT_DIR
    CLUSTER_OUTPUT_DIR = _DEFAULT_CLUSTER_OUTPUT_DIR\
        if _ENV_CLUSTER_OUTPUT_DIR is None else _ENV_CLUSTER_OUTPUT_DIR

    ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
    CONFIG_DIR = f"{ROOT_DIR}/config"

    if not os.path.exists(INPUT_DIR):
        os.makedirs(INPUT_DIR)
    assert os.path.exists(INPUT_DIR)

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    assert os.path.exists(OUTPUT_DIR)

    SNAKEMAKE_FILEPATH = f"{ROOT_DIR}/Snakefile"
    EXPERIMENT_YML = "experiment.yaml"
    SNAKEMAKE_CONFIG_YML = "config.yaml"
    SNAKEMAKE_CLUSTER_CONFIG_YML = "config_cluster.yaml"