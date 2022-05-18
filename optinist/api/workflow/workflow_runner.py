from typing import Dict, List
from dataclasses import asdict

from optinist.api.snakemake.smk import FlowConfig, ForceRun, Rule, SmkParam
from optinist.api.snakemake.snakemake_reader import SmkParamReader
from optinist.api.snakemake.snakemake_writer import SmkConfigWriter
from optinist.api.snakemake.snakemake_setfile import SmkSetfile
from optinist.api.snakemake.snakemake_executor import delete_dependencies, snakemake_execute

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
            ExptConfigReader.nodeDict_read(runItem.nodeDict),
            ExptConfigReader.edgeDict_read(runItem.edgeDict),
            runItem.nwbParam
        )

        snakemake_params: SmkParam = get_typecheck_params(runItem.snakemakeParam, "snakemake")
        snakemake_params = SmkParamReader.read(snakemake_params)
        snakemake_params.forcerun = _get_forcerun_list(unique_id, runItem.forceRunList)
        if len(snakemake_params.forcerun) > 0:
            delete_dependencies(snakemake_params)
        background_tasks.add_task(snakemake_execute, snakemake_params)


def _create_workflow(unique_id, name, nodeDict, edgeDict, nwbParam):
    rules, last_output = _rulefile(
        unique_id,
        nodeDict,
        edgeDict,
        nwbParam,
    )

    flow_config = FlowConfig(
        rules=rules,
        last_output=last_output,
    )

    SmkConfigWriter.write(unique_id, asdict(flow_config))
    ExptConfigWriter.write(unique_id, name, nodeDict, edgeDict)


def _rulefile(unique_id: str, nodeDict: Dict[str, Node], edgeDict: Dict[str, Edge], nwbParam):
    endNodeList = _get_endNodeList(edgeDict, nodeDict)

    nwbfile = get_typecheck_params(nwbParam, "nwb")

    rule_dict: Dict[str, Rule] = {}
    last_outputs = []

    for node in nodeDict.values():
        if node.type == NodeType.IMAGE:
            rule_dict[node.id] = SmkSetfile.image(unique_id, node, edgeDict, nwbfile)
        elif node.type == NodeType.CSV:
            rule_dict[node.id] = SmkSetfile.csv(unique_id, node, edgeDict, nwbfile)
        elif node.type == NodeType.FLUO:
            rule_dict[node.id] = SmkSetfile.csv(unique_id, node, edgeDict, nwbfile)
        elif node.type == NodeType.BEHAVIOR:
            rule_dict[node.id] = SmkSetfile.csv(unique_id, node, edgeDict, nwbfile, "behavior")
        elif node.type == NodeType.HDF5:
            rule_dict[node.id] = SmkSetfile.hdf5(unique_id, node, edgeDict, nwbfile)
        elif node.type == NodeType.ALGO:
            rule = SmkSetfile.algo(unique_id, node, edgeDict, nodeDict)
            rule_dict[node.id] = rule

            if node.id in endNodeList:
                last_outputs.append(rule.output)
        else:
            assert False, "NodeType doesn't exists"

    return rule_dict, last_outputs


def _get_endNodeList(edgeDict: Dict[str, Edge], nodeDict: Dict[str, Node]):
    returnCntDict = {key: 0 for key in nodeDict.keys()}
    for edge in edgeDict.values():
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
