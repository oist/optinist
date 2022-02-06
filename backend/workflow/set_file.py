import pickle

from wrappers.data_wrapper import *
from workflow.params import get_algo_params


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

    save_pickle(node["data"]["path"], info)

    return info


def set_algofile(node, edgeList, nodeDict):
    algo_input = []
    return_arg_names = {}
    for edge in edgeList:
        # inputとして入れる
        if node["id"] == edge["target"]:
            arg_name = edge["targetHandle"].split("--")[1]
            return_name = edge["sourceHandle"].split("--")[1]
            sourceNode = nodeDict[edge["source"]]
            if sourceNode["type"] == "AlgorithmNode":
                algo_input.append(os.path.join(
                    BASE_DIR, sourceNode["data"]["path"], f"{sourceNode['data']['label']}_out.pkl"))
            else:
                algo_input.append(sourceNode["data"]["path"])

            return_arg_names[return_name] = arg_name

    # parameter
    algo_path = node["data"]["path"]
    algo_name = node["data"]["label"]
    message_params = node["data"]["param"]
    params = get_algo_params(algo_name, message_params)

    algo_output = os.path.join(BASE_DIR, algo_path, f"{algo_name}_out.pkl")

    rule = {
        "rule_file": f"rules/smk/{algo_path}.py",
        "input": algo_input,
        "return_arg": return_arg_names,
        "params": params,
        "output": algo_output,
        "path": algo_path,
    }
    return rule


def save_pickle(filepath, info):
    # output to pickle
    pickle_path = filepath.split(".")[0] + ".pkl"
    with open(pickle_path, 'wb') as f:
        pickle.dump(info, f)
