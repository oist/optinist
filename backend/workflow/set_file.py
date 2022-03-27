import pickle
import h5py

from wrappers.data_wrapper import *
from workflow.params import get_typecheck_params
from cui_api.const import BASE_DIR
from cui_api.utils import join_file_path
from wrappers.nwb_wrapper.const import NWBDATASET


def set_imagefile(unique_id, node, edgeList, nwbfile):
    for edge in edgeList:
        if node["id"] == edge["source"]:
            return_name = edge["sourceHandle"].split("--")[0]

    output_file = get_pickle_file(
        unique_id, node["id"], node["data"]["label"].split(".")[0])

    rule = {
        "rule_file": f"rules/smk/image.py",
        "input": node['data']['path'],
        "return_arg": return_name,
        "params": node["data"]["param"],
        "output": output_file,
        "type": "image",
        "nwbfile": nwbfile,
    }

    return rule


def set_csvfile(unique_id, node, edgeList, nwbfile):
    for edge in edgeList:
        if node["id"] == edge["source"]:
            return_name = edge["sourceHandle"].split("--")[0]

    output_file = get_pickle_file(
        unique_id, node["id"], node["data"]["label"].split(".")[0])

    rule = {
        "rule_file": f"rules/smk/csv.py",
        "input": node['data']['path'],
        "return_arg": return_name,
        "params": node["data"]["param"],
        "output": output_file,
        "type": "csv",
        "nwbfile": nwbfile,
    }

    return rule


def set_hdf5file(unique_id, node, edgeList, nwbfile):
    for edge in edgeList:
        if node["id"] == edge["source"]:
            return_name = edge["sourceHandle"].split("--")[0]

    output_file = get_pickle_file(
        unique_id, node["id"], node["data"]["label"].split(".")[0])

    rule = {
        "rule_file": f"rules/smk/hdf5.py",
        "input": node['data']['path'],
        "return_arg": return_name,
        "params": node["data"]["param"],
        "output": output_file,
        "type": "hdf5",
        "nwbfile": nwbfile,
        "hdf5Path": node["data"]["hdf5Path"],
    }

    return rule


def set_algofile(unique_id, node, edgeList, nodeDict):
    algo_input = []
    return_arg_names = {}
    for edge in edgeList:
        # inputとして入れる
        if node["id"] == edge["target"]:
            arg_name = edge["targetHandle"].split("--")[1]

            sourceNode = nodeDict[edge["source"]]
            if sourceNode["type"] == "AlgorithmNode":
                return_name = edge["sourceHandle"].split("--")[1]
                input_pickle_file = get_pickle_file(
                    unique_id, sourceNode["id"], sourceNode['data']['label'])
                algo_input.append(input_pickle_file)
            else:
                return_name = edge["sourceHandle"].split("--")[0]
                input_file = get_pickle_file(
                    unique_id, sourceNode["id"], sourceNode["data"]["label"].split(".")[0])
                algo_input.append(input_file)

            return_arg_names[return_name] = arg_name

    # parameter
    algo_path = node["data"]["path"]
    algo_name = node["data"]["label"]
    message_params = node["data"]["param"]
    params = get_typecheck_params(message_params, algo_name)

    algo_output = get_pickle_file(unique_id, node["id"], algo_name)

    rule = {
        "rule_file": f"rules/smk/{algo_path}.py",
        "input": algo_input,
        "return_arg": return_arg_names,
        "params": params,
        "output": algo_output,
        "path": algo_path,
        "type": algo_name,
    }

    return rule


def get_pickle_file(unique_id, node_id, algo_name):
    return join_file_path([BASE_DIR, unique_id, node_id, f"{algo_name}.pkl"])


def save_pickle(filepath, info):
    # output to pickle
    pickle_path = filepath.split(".")[0] + ".pkl"
    with open(pickle_path, 'wb') as f:
        pickle.dump(info, f)


def get_forcerun_list(unique_id, forceRunList):
    target_list = []
    for x in forceRunList:
        target_list.append(get_pickle_file(unique_id, x.nodeId, x.name))
    return target_list
