from typing import Dict, List
from dataclasses import asdict

from optinist.api.snakemake.smk import FlowConfig, ForceRun, Rule, SmkParam
from optinist.api.snakemake.snakemake_reader import SmkParamReader
from optinist.api.snakemake.snakemake_writer import SmkConfigWriter
from optinist.api.snakemake.snakemake_setfile import SmkSetfile
# from optinist.api.snakemake.snakemake_run import run_snakemake
from optinist.api.snakemake.snakemake_executor import snakemake_execute
from optinist.api.utils.filepath_creater import get_pickle_file
from optinist.api.workflow.workflow import Edge, Node, NodeType, RunItem
from optinist.api.workflow.workflow_params import get_typecheck_params
from optinist.api.experiment.experiment_reader import ExptConfigReader
from optinist.api.experiment.experiment_writer import ExptConfigWriter


class WorkflowRunner:
    @classmethod
    def run_workflow(cls, unique_id, background_tasks, runItem: RunItem):
        _create_workflow(
            unique_id,
            runItem.name,
            ExptConfigReader.nodeList_read(runItem.nodeList),
            ExptConfigReader.edgeList_read(runItem.edgeList),
            runItem.nwbParam
        )

        snakemake_params: SmkParam = get_typecheck_params(runItem.snakemakeParam, "snakemake")
        snakemake_params = SmkParamReader.read(snakemake_params)
        snakemake_params.forcerun = _get_forcerun_list(unique_id, runItem.forceRunList)
        # background_tasks.add_task(run_snakemake, snakemake_params)
        background_tasks.add_task(snakemake_execute, snakemake_params)


def _create_workflow(unique_id, name, nodeList, edgeList, nwbParam):
    rules, last_output = _rulefile(
        unique_id,
        nodeList,
        edgeList,
        nwbParam,
    )

    flow_config = FlowConfig(
        rules=rules,
        last_output=last_output,
    )

    SmkConfigWriter.write(unique_id, asdict(flow_config))
    ExptConfigWriter.write(unique_id, name, nodeList, edgeList)


def _rulefile(unique_id, nodeList, edgeList, nwbParam):
    nodeDict = _get_nodeDict(nodeList)
    endNodeList = _get_endNodeList(edgeList, nodeDict)

    nwbfile = get_typecheck_params(nwbParam, "nwb")

    rule_dict: Dict[str, Rule] = {}
    last_outputs = []

    for node in nodeDict.values():
        if node.type == NodeType.IMAGE:
            rule_dict[node.id] = SmkSetfile.image(unique_id, node, edgeList, nwbfile)
        elif node.type == NodeType.CSV:
            rule_dict[node.id] = SmkSetfile.csv(unique_id, node, edgeList, nwbfile)
        elif node.type == NodeType.FLUO:
            rule_dict[node.id] = SmkSetfile.csv(unique_id, node, edgeList, nwbfile)
        elif node.type == NodeType.BEHAVIOR:
            rule_dict[node.id] = SmkSetfile.csv(unique_id, node, edgeList, nwbfile, "behavior")
        elif node.type == NodeType.HDF5:
            rule_dict[node.id] = SmkSetfile.hdf5(unique_id, node, edgeList, nwbfile)
        elif node.type == NodeType.ALGO:
            rule = SmkSetfile.algo(unique_id, node, edgeList, nodeDict)
            rule_dict[node.id] = rule

            if node.id in endNodeList:
                last_outputs.append(rule.output)
        else:
            assert False, "NodeType doesn't exists"

    return rule_dict, last_outputs


def _get_nodeDict(nodeList: List[Node]) -> Dict[str, Node]:
    # nodeを初期化
    nodeDict = {}
    for node in nodeList:
        nodeDict[node.id] = node
    return nodeDict


def _get_endNodeList(edgeList: List[Edge], nodeDict):
    returnCntDict = {key: 0 for key in nodeDict.keys()}
    for edge in edgeList:
        returnCntDict[edge.source] += 1

    endNodeList = []
    for key, value in returnCntDict.items():
        if value == 0:
            endNodeList.append(key)
    return endNodeList


def _get_forcerun_list(unique_id, forceRunList: List[ForceRun]) -> List[str]:
    target_list = []
    for x in forceRunList:
        target_list.append(get_pickle_file(unique_id, x.nodeId, x.name))
    return target_list
