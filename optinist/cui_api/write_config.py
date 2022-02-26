import os
import yaml
from datetime import datetime

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
        "function": {},
        "nodeList": runItem.nodeList,
        "edgeList": runItem.edgeList,
    }

    for node in runItem.nodeList:
        name = node["data"]["label"]
        # position = node["position"]

        success = "running"
        if node["type"] == "ImageFileNode":
            success = "success"

        experiment_config["function"][name] = {
            "unique_id": node["id"],
            # "position": position,
            "success": success,
        }
        
    save_path = join_file_path([BASE_DIR, unique_id])
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    with open(join_file_path([save_path, 'experiment.yaml']), "w") as f:
        yaml.dump(experiment_config, f)
