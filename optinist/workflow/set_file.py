
from dataclasses import dataclass
from typing import Dict, List
from pydantic import BaseModel
from optinist.cui_api.experiment_config import Edge, Node
from optinist.workflow.data_type import NodeType

from optinist.wrappers.data_wrapper import *
from optinist.workflow.params import get_typecheck_params
from optinist.cui_api.filepath_creater import get_pickle_file


class ForceRun(BaseModel):
    nodeId: str
    name: str


@dataclass
class Rule:
    rule_file: str
    input: str
    return_arg: str
    params: dict
    output: str
    type: str
    nwbfile: dict = None
    hdf5Path: str = None
    path: str = None


def set_imagefile(unique_id: str, node: Node, edgeList: List[Edge], nwbfile):
    for edge in edgeList:
        if node.id == edge.source:
            return_name = edge.sourceHandle.split("--")[0]

    output_file = get_pickle_file(
        unique_id, node.id, node.data.label.split(".")[0])

    return Rule(
        rule_file=f"rules/smk/image.smk",
        input=node.data.path,
        return_arg=return_name,
        params=node.data.param,
        output=output_file,
        type="image",
        nwbfile=nwbfile,
    )


def set_csvfile(unique_id, node: Node, edgeList: List[Edge], nwbfile, nodeType="csv"):
    for edge in edgeList:
        if node.id == edge.source:
            return_name = edge.sourceHandle.split("--")[0]

    output_file = get_pickle_file(
        unique_id, node.id, node.data.label.split(".")[0])

    return Rule(
        rule_file=f"rules/smk/{nodeType}.smk",
        input=node.data.path,
        return_arg=return_name,
        params=node.data.param,
        output=output_file,
        type=nodeType,
        nwbfile=nwbfile,
    )


def set_hdf5file(unique_id, node, edgeList, nwbfile):
    for edge in edgeList:
        if node.id == edge.source:
            return_name = edge.sourceHandle.split("--")[0]

    output_file = get_pickle_file(
        unique_id, node.id, node.data.label.split(".")[0])

    return Rule(
        rule_file=f"rules/smk/hdf5.smk",
        input=node.data.path,
        return_arg=return_name,
        params=node.data.param,
        output=output_file,
        type="hdf5",
        nwbfile=nwbfile,
        hdf5Path=node.data.hdf5Path,
    )


def set_algofile(unique_id, node: Node, edgeList: List[Edge], nodeDict: Dict[str, Node]):
    algo_input = []
    return_arg_names = {}
    for edge in edgeList:
        if node.id == edge.target:
            arg_name = edge.targetHandle.split("--")[1]

            sourceNode = nodeDict[edge.source]
            if sourceNode.type == NodeType.ALGO:
                return_name = edge.sourceHandle.split("--")[1]
                algo_input.append(get_pickle_file(
                    unique_id,
                    sourceNode.id,
                    sourceNode.data.label
                ))
            else:
                return_name = edge.sourceHandle.split("--")[0]
                algo_input.append(get_pickle_file(
                    unique_id,
                    sourceNode.id,
                    sourceNode.data.label.split(".")[0]
                ))

            return_arg_names[return_name] = arg_name

    params = get_typecheck_params(node.data.param, node.data.label)
    algo_output = get_pickle_file(unique_id, node.id, node.data.label)

    return Rule(
        rule_file=f"rules/smk/{node.data.path}.smk",
        input=algo_input,
        return_arg=return_arg_names,
        params=params,
        output=algo_output,
        path=node.data.path,
        type=node.data.label,
    )


def get_forcerun_list(unique_id, forceRunList: List[ForceRun]):
    target_list = []
    for x in forceRunList:
        target_list.append(get_pickle_file(unique_id, x.nodeId, x.name))
    return target_list
