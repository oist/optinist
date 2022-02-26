import os
import yaml
from datetime import datetime
from snakemake import snakemake

from cui_api.const import BASE_DIR, OPTINIST_DIR
from cui_api.utils import join_file_path


def write_snakemake_config(unique_id, flow_config):
    with open(join_file_path([OPTINIST_DIR, 'config.yaml']), "w") as f:
        yaml.dump(flow_config, f)

    save_path = join_file_path([BASE_DIR, unique_id])
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    with open(join_file_path([save_path, 'config.yaml']), "w") as f:
        yaml.dump(flow_config, f)


def write_experiment_config(unique_id, rule_config, runItem):
    experiment_config = {
        "timestamp": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        "name": "test",
        "unique_id": unique_id,
        "success": "running",
        "function": {},
    }

    for node in runItem.nodeList:
        name = node["data"]["label"]
        position = node["position"]

        experiment_config["function"][name] = {
            "unique_id": node["id"],
            "position": position,
            "success": "running",
        }
        
    save_path = join_file_path([BASE_DIR, unique_id])
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    with open(join_file_path([save_path, 'experiment.yaml']), "w") as f:
        yaml.dump(experiment_config, f)


def run_snakemake(snakemake_params):
    # run snakemake
    snakemake(
        join_file_path([OPTINIST_DIR, 'Snakefile']),
        **snakemake_params
    )