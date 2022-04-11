from typing import Dict, List
from dataclasses import asdict
from optinist.api.snakemake.snakemake import FlowConfig, Rule
from optinist.api.snakemake.snakemake_run import run_snakemake
from optinist.api.workflow.workflow import Node, NodeType

from optinist.api.workflow.workflow_params import get_typecheck_params
from optinist.api.snakemake.snakemake_setfile import (
    get_forcerun_list,
    set_imagefile,
    set_csvfile,
    set_algofile,
    set_hdf5file
)
from optinist.api.experiment.experiment_config_reader import ExpConfigReader
from optinist.api.experiment.experiment_config_writer import ExpConfigWriter
from optinist.api.snakemake.snakemake_config import snakemake_config_writer


def create_workflow(unique_id, name, nodeList, edgeList, nwbParam):
    rules, last_output = create_rule_file(
        unique_id,
        nodeList,
        edgeList,
        nwbParam,
    )

    flow_config = FlowConfig(
        rules=rules,
        last_output=last_output,
    )

    snakemake_config_writer(unique_id, asdict(flow_config))
    ExpConfigWriter.exp_config_writer(unique_id, name, nodeList, edgeList)


def create_rule_file(unique_id, nodeList, edgeList, nwbParam):
    nodeDict = get_nodeDict(nodeList)
    endNodeList = get_endNodeList(edgeList, nodeDict)

    nwbfile = get_typecheck_params(nwbParam, "nwb")

    rule_dict: Dict[str, Rule] = {}
    last_outputs = []

    for node in nodeDict.values():
        if node.type == NodeType.IMAGE:
            rule_dict[node.id] = set_imagefile(unique_id, node, edgeList, nwbfile)
        elif node.type == NodeType.CSV:
            rule_dict[node.id] = set_csvfile(unique_id, node, edgeList, nwbfile)
        elif node.type == NodeType.FLUO:
            rule_dict[node.id] = set_csvfile(unique_id, node, edgeList, nwbfile)
        elif node.type == NodeType.BEHAVIOR:
            rule_dict[node.id] = set_csvfile(unique_id, node, edgeList, nwbfile, "behavior")
        elif node.type == NodeType.HDF5:
            rule_dict[node.id] = set_hdf5file(unique_id, node, edgeList, nwbfile)
        elif node.type == NodeType.ALGO:
            rule = set_algofile(unique_id, node, edgeList, nodeDict)
            rule_dict[node.id] = rule

            if node.id in endNodeList:
                last_outputs.append(rule.output)
        else:
            assert False, "NodeType doesn't exists"

    return rule_dict, last_outputs


def get_nodeDict(nodeList: List[Node]) -> Dict[str, Node]:
    # nodeを初期化
    nodeDict = {}
    for node in nodeList:
        nodeDict[node.id] = node
    return nodeDict


def get_endNodeList(edgeList, nodeDict):
    returnCntDict = {key: 0 for key in nodeDict.keys()}
    for edge in edgeList:
        returnCntDict[edge.source] += 1

    endNodeList = []
    for key, value in returnCntDict.items():
        if value == 0:
            endNodeList.append(key)
    return endNodeList


def run_workflow(unique_id, background_tasks, runItem):
    create_workflow(
        unique_id,
        runItem.name,
        ExpConfigReader.nodeList_read(runItem.nodeList),
        ExpConfigReader.edgeList_read(runItem.edgeList),
        runItem.nwbParam
    )

    snakemake_params = get_typecheck_params(runItem.snakemakeParam, "snakemake")
    snakemake_params["forcerun"] = get_forcerun_list(unique_id, runItem.forceRunList)
    background_tasks.add_task(run_snakemake, snakemake_params)
