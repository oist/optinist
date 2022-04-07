import pickle
import yaml
from collections import OrderedDict
from glob import glob
from optinist.wrappers.data_wrapper import *

from optinist.cui_api.utils import join_file_path
from optinist.cui_api.const import BASE_DIR


def get_results(unique_id, nodeIdList):
    runPaths = []
    for node_id in nodeIdList:
        print(join_file_path([BASE_DIR, unique_id, node_id, "*.pkl"]))
        for path in glob(join_file_path([BASE_DIR, unique_id, node_id, "*.pkl"])):
            runPaths.append(path)

    for i, path in enumerate(runPaths):
        runPaths[i] = path.replace("\\", "/")

    print(runPaths)

    results = {}
    for request_path in runPaths:
        node_id = request_path.split("/")[-2]
        algo_name = request_path.split("/")[-1].split(".")[0]

        results[node_id] = {}

        with open(request_path, "rb") as f:
            info = pickle.load(f)

        if isinstance(info, list) or isinstance(info, str):
            results[node_id] = get_error(info, node_id, unique_id)
        else:
            json_dir = "/".join(request_path.split("/")[:-1])
            results[node_id] = get_success(info, node_id, algo_name, json_dir, unique_id)

    return results


def get_error(info, node_id, unique_id):
    with open(join_file_path([BASE_DIR, unique_id, "experiment.yaml"]), "r") as f:
        config = yaml.safe_load(f)

    config["function"][node_id]["success"] = "error"

    with open(join_file_path([BASE_DIR, unique_id, "experiment.yaml"]), "w") as f:
        yaml.dump(config, f)

    if isinstance(info, str):
        error_message = info
    else:
        error_message = "\n".join(info)

    message = {
        "status": "error",
        "message": error_message,
    }

    return message


def get_success(info, node_id, algo_name, json_dir, unique_id):
    with open(join_file_path([BASE_DIR, unique_id, "experiment.yaml"]), "r") as f:
        config = yaml.safe_load(f)

    config["function"][node_id]["success"] = "success"

    with open(join_file_path([BASE_DIR, unique_id, "experiment.yaml"]), "w") as f:
        yaml.dump(config, f)

    message = {
        "status": "success",
        "message": f"{algo_name} success",
        "outputPaths": get_outputPaths(info, json_dir)
    }

    return message


def get_outputPaths(info, json_dir):
    outputPaths = {}
    for k, v in info.items():
        if isinstance(v, BaseData):
            v.save_json(json_dir)

        if isinstance(v, ImageData):
            outputPaths[k] = {}
            outputPaths[k]['path'] = v.json_path
            outputPaths[k]['type'] = 'images'
            if v.data.ndim == 3:
                outputPaths[k]['max_index'] = len(v.data)
            elif v.data.ndim == 2:
                outputPaths[k]['max_index'] = 1
        elif isinstance(v, TimeSeriesData):
            outputPaths[k] = {}
            outputPaths[k]['path'] = v.json_path
            outputPaths[k]['type'] = 'timeseries'
            outputPaths[k]['max_index'] = len(v.data)
        elif isinstance(v, CorrelationData):
            outputPaths[k] = {}
            outputPaths[k]['path'] = v.json_path
            outputPaths[k]['type'] = 'heatmap'
        elif isinstance(v, RoiData):
            outputPaths[k] = {}
            outputPaths[k]['path'] = v.json_path
            outputPaths[k]['type'] = 'roi'
        elif isinstance(v, ScatterData):
            outputPaths[k] = {}
            outputPaths[k]['path'] = v.json_path
            outputPaths[k]['type'] = 'scatter'
        elif isinstance(v, BarData):
            outputPaths[k] = {}
            outputPaths[k]['path'] = v.json_path
            outputPaths[k]['type'] = 'bar'
        elif isinstance(v, HTMLData):
            outputPaths[k] = {}
            outputPaths[k]['path'] = v.json_path
            outputPaths[k]['type'] = 'html'
        else:
            pass

    return outputPaths
