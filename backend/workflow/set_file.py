import pickle
import h5py

from wrappers.data_wrapper import *
from workflow.params import get_typecheck_params
from cui_api.const import BASE_DIR
from cui_api.utils import join_file_path


def set_imagefile(node, edgeList, nwbfile):
    for edge in edgeList:
        if node["id"] == edge["source"]:
            return_name = edge["sourceHandle"].split("--")[1]
            info = {return_name: ImageData(node['data']['path'], '')}
    
    # NWB file
    nwbfile['image_series']['external_file'] = info[return_name]
    info['nwbfile'] = nwbfile
    
    save_pickle(node["data"]["path"], info)

    return info


def set_csvfile(node, edgeList):
    for edge in edgeList:
        if node["id"] == edge["source"]:
            return_name = edge["sourceHandle"].split("--")[1]
            info = {return_name: TimeSeriesData(node['data']['path'], '')}

    info['nwbfile'] = None

    save_pickle(node["data"]["path"], info)

    return info


def set_hdf5file(node, edgeList):
    path = node["data"]["path"]
    h5path = node["data"]["hdf5Path"]

    with h5py.File(path, "r") as f:
        data = f[h5path][:]

    for edge in edgeList:
        if node["id"] == edge["source"]:
            return_name = edge["sourceHandle"].split("--")[1]
            if data.ndim == 3:
                info = {return_name: ImageData(data, '')}
            elif data.ndim == 2:
                info = {return_name: TimeSeriesData(data, '')}

    info['nwbfile'] = None

    save_pickle(node["data"]["path"], info)

    return info


def set_algofile(unique_id, node, edgeList, nodeDict):
    algo_input = []
    return_arg_names = {}
    for edge in edgeList:
        # inputとして入れる
        if node["id"] == edge["target"]:
            arg_name = edge["targetHandle"].split("--")[1]
            return_name = edge["sourceHandle"].split("--")[1]
            sourceNode = nodeDict[edge["source"]]
            if sourceNode["type"] == "AlgorithmNode":
                input_pickle_file = get_pickle_file(unique_id, sourceNode["id"], sourceNode['data']['label'])
                algo_input.append(input_pickle_file)
            else:
                algo_input.append(sourceNode["data"]["path"])

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
    }

    return rule


def get_pickle_file(unique_id, node_id, algo_name):
    return join_file_path([BASE_DIR, unique_id, node_id, f"{algo_name}_out.pkl"])


def save_pickle(filepath, info):
    # output to pickle
    pickle_path = filepath.split(".")[0] + ".pkl"
    with open(pickle_path, 'wb') as f:
        pickle.dump(info, f)
