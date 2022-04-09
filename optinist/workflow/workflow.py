from optinist.workflow.params import get_typecheck_params
from optinist.workflow.set_file import set_imagefile, set_csvfile, set_algofile, set_hdf5file
from optinist.cui_api.experiment_config import exp_config_writer
from optinist.cui_api.snakemake_config import snakemake_config_writer


def create_workflow(unique_id, runItem):
    rules_to_execute, last_outputs = create_rule_file(
        unique_id,
        runItem.nodeList,
        runItem.edgeList,
        runItem.nwbParam,
    )

    flow_config = {
        "rules": rules_to_execute,
        "last_output": last_outputs,
    }

    snakemake_config_writer(unique_id, flow_config)
    exp_config_writer(unique_id, rules_to_execute, runItem)


def create_rule_file(unique_id, nodeList, edgeList, nwbParam):
    nodeDict, endNodeList = get_nodeDict(nodeList)
    endNodeList = get_endNodeList(edgeList)

    nwbfile = get_typecheck_params(nwbParam, "nwb")

    rules_dict = {}
    last_outputs = []

    for node in nodeDict.values():
        if node["type"] == "ImageFileNode":
            rule = set_imagefile(unique_id, node, edgeList, nwbfile)
            rules_dict[node["id"]] = rule
        elif node["type"] == "CsvFileNode":
            rule = set_csvfile(unique_id, node, edgeList, nwbfile)
            rules_dict[node["id"]] = rule
        elif node["type"] == "FluoFileNode":
            rule = set_csvfile(unique_id, node, edgeList, nwbfile)
            rules_dict[node["id"]] = rule
        elif node["type"] == "BehaviorFileNode":
            rule = set_csvfile(unique_id, node, edgeList, nwbfile, "behavior")
            rules_dict[node["id"]] = rule
        elif node["type"] == "HDF5FileNode":
            rule = set_hdf5file(unique_id, node, edgeList, nwbfile)
            rules_dict[node["id"]] = rule
        elif node["type"] == "AlgorithmNode":
            rule = set_algofile(unique_id, node, edgeList, nodeDict)
            rules_dict[node["id"]] = rule

            if node["id"] in endNodeList:
                last_outputs.append(rule["output"])
        else:
            assert False, "NodeType doesn't exists"

    return rules_dict, last_outputs


def get_nodeDict(nodeList):
    # nodeを初期化
    nodeDict = {}
    for node in nodeList:
        nodeDict[node['id']] = node
    return nodeDict


def get_endNodeList(edgeList, nodeDict):
    returnCntDict = {key: 0 for key in nodeDict.keys()}
    for edge in edgeList:
        returnCntDict[edge["source"]] += 1

    endNodeList = []
    for key, value in returnCntDict.items():
        if value == 0:
            endNodeList.append(key)
    return endNodeList
