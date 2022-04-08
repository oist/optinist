import yaml

from optinist.workflow.get_network import get_network
from optinist.workflow.params import get_typecheck_params
from optinist.workflow.set_file import set_imagefile, set_csvfile, set_algofile, set_hdf5file
from optinist.cui_api.experiment_config import exp_config_writer
from optinist.cui_api.snakemake_config import snakemake_config_writer


def set_workflow(unique_id, runItem):
    rules_to_execute, last_outputs, all_outputs = get_workflow(unique_id, runItem)

    flow_config = {
        "rules": rules_to_execute,
        "last_output": last_outputs,
    }

    snakemake_config_writer(unique_id, flow_config)
    exp_config_writer(unique_id, rules_to_execute, runItem)


def get_workflow(unique_id, runItem):
    # graph networkの解析
    nodeDict, edgeList, endNodeList = get_network(runItem)

    nwbfile = get_typecheck_params(runItem.nwbParam, "nwb")

    rules_to_execute = {}
    last_outputs = []
    all_outputs = {}

    for node in nodeDict.values():
        algo_label = node['data']['label']
        algo_path = node['data']['path']
        if node["type"] == "ImageFileNode":
            rule = set_imagefile(unique_id, node, edgeList, nwbfile)
            rules_to_execute[node["id"]] = rule
        elif node["type"] == "CsvFileNode":
            rule = set_csvfile(unique_id, node, edgeList, nwbfile)
            rules_to_execute[node["id"]] = rule
        elif node["type"] == "FluoFileNode":
            rule = set_csvfile(unique_id, node, edgeList, nwbfile)
            rules_to_execute[node["id"]] = rule
        elif node["type"] == "BehaviorFileNode":
            rule = set_csvfile(unique_id, node, edgeList, nwbfile, "behavior")
            rules_to_execute[node["id"]] = rule
        elif node["type"] == "HDF5FileNode":
            rule = set_hdf5file(unique_id, node, edgeList, nwbfile)
            rules_to_execute[node["id"]] = rule
        elif node["type"] == "AlgorithmNode":
            rule = set_algofile(unique_id, node, edgeList, nodeDict)
            rules_to_execute[node["id"]] = rule

            if node["id"] in endNodeList:
                last_outputs.append(rule["output"])

            all_outputs[rule["output"]] = {
                "label": algo_label,
                "path": algo_path,
            }
        else:
            assert False, "NodeType doesn't exists"

    return rules_to_execute, last_outputs, all_outputs
