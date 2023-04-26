from typing import Dict

from optinist.api.snakemake.smk import Rule
from optinist.api.snakemake.smk_builder import RuleBuilder
from optinist.api.utils.filepath_creater import get_pickle_file
from optinist.api.workflow.workflow import Edge, Node, NodeType
from optinist.api.workflow.workflow_params import get_typecheck_params


class SmkRule:
    def __init__(
        self, unique_id: str, node: Node, edgeDict: Dict[str, Edge], nwbfile=None
    ) -> None:
        self._unique_id = unique_id
        self._node = node
        self._edgeDict = edgeDict
        self._nwbfile = nwbfile

        _return_name = self.get_return_name()

        _output_file = get_pickle_file(
            self._unique_id, self._node.id, self._node.data.label.split(".")[0]
        )

        self.builder = RuleBuilder()
        (
            self.builder.set_input(self._node.data.path)
            .set_return_arg(_return_name)
            .set_params(self._node.data.param)
            .set_output(_output_file)
            .set_nwbfile(self._nwbfile)
        )

    def image(self) -> Rule:
        return self.builder.set_type("image").build()

    def csv(self, nodeType="csv") -> Rule:
        return self.builder.set_type(nodeType).build()

    def hdf5(self) -> Rule:
        return (
            self.builder.set_type("hdf5").set_hdf5Path(self._node.data.hdf5Path).build()
        )

    def algo(self, nodeDict: Dict[str, Node]) -> Rule:
        algo_input = []
        return_arg_names = {}
        for edge in self._edgeDict.values():
            if self._node.id == edge.target:
                arg_name = edge.targetHandle.split("--")[1]

                sourceNode = nodeDict[edge.source]
                if sourceNode.type == NodeType.ALGO:
                    return_name = edge.sourceHandle.split("--")[1]
                    funcname = sourceNode.data.label
                else:
                    return_name = edge.sourceHandle.split("--")[0]
                    funcname = sourceNode.data.label.split(".")[0]

                algo_input.append(
                    get_pickle_file(self._unique_id, sourceNode.id, funcname)
                )

                return_arg_names[return_name] = arg_name

        params = get_typecheck_params(self._node.data.param, self._node.data.label)
        algo_output = get_pickle_file(
            self._unique_id, self._node.id, self._node.data.label
        )

        return (
            self.builder.set_input(algo_input)
            .set_return_arg(return_arg_names)
            .set_params(params)
            .set_output(algo_output)
            .set_path(self._node.data.path)
            .set_type(self._node.data.label)
            .build()
        )

    def get_return_name(self) -> str or None:
        for edge in self._edgeDict.values():
            if self._node.id == edge.source:
                return edge.sourceHandle.split("--")[0]
        return None
