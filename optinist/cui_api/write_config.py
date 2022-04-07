import os
import yaml
from datetime import datetime

from optinist.cui_api.const import BASE_DIR, OPTINIST_DIR
from optinist.cui_api.utils import join_file_path


def write_snakemake_config(unique_id, flow_config):
    with open(join_file_path([OPTINIST_DIR, 'config.yaml']), "w") as f:
        yaml.dump(flow_config, f)

    save_path = join_file_path([BASE_DIR, unique_id])
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    with open(join_file_path([save_path, 'config.yaml']), "w") as f:
        yaml.dump(flow_config, f)


def write_experiment_config(unique_id, rule_config, runItem):

    save_path = join_file_path([BASE_DIR, unique_id])
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    exp_filepath = join_file_path([save_path, 'experiment.yaml'])

    if os.path.exists(exp_filepath):
        with open(exp_filepath, "r") as f:
            experiment_config = yaml.safe_load(f)
        experiment_config["timestamp"] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        experiment_config["nodeList"] = runItem.nodeList
        experiment_config["edgeList"] = runItem.edgeList
    else:
        experiment_config = {
            "timestamp": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "name": runItem.name,
            "unique_id": unique_id,
            "function": {},
            "nodeList": runItem.nodeList,
            "edgeList": runItem.edgeList,
        }

    for node in runItem.nodeList:
        name = node["data"]["label"]

        if node["data"]["type"] == "input":
            success = "success"
        else:
            success = "running"

        experiment_config["function"][node["id"]] = {
            "unique_id": node["id"],
            "name": name,
            "success": success,
        }

    with open(join_file_path([save_path, 'experiment.yaml']), "w") as f:
        yaml.dump(experiment_config, f)
