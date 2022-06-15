from typing import Dict, List
from dataclasses import asdict

from optinist.api.snakemake.smk import FlowConfig, Rule, SmkParam
from optinist.api.snakemake.snakemake_reader import SmkParamReader
from optinist.api.snakemake.snakemake_writer import SmkConfigWriter
from optinist.api.snakemake.snakemake_setfile import SmkSetfile
from optinist.api.snakemake.snakemake_executor import delete_dependencies, snakemake_execute

from optinist.api.utils.filepath_creater import get_pickle_file
from optinist.api.workflow.workflow_params import get_typecheck_params
from optinist.api.workflow.workflow import NodeType, RunItem
from optinist.api.experiment.experiment_reader import ExptConfigReader
from optinist.api.experiment.experiment_writer import ExptConfigWriter


class WorkflowRunner:

    def __init__(self, unique_id, runItem: RunItem) -> None:
        self.unique_id = unique_id
        self.runItem = runItem
        self.nodeDict = ExptConfigReader.read_nodeDict(self.runItem.nodeDict)
        self.edgeDict = ExptConfigReader.read_edgeDict(self.runItem.edgeDict)

        ExptConfigWriter.write(
            self.unique_id,
            self.runItem.name,
            self.nodeDict,
            self.edgeDict
        )

    def run_workflow(self, background_tasks):
        self.create_workflow()

        snakemake_params: SmkParam = get_typecheck_params(self.runItem.snakemakeParam, "snakemake")
        snakemake_params = SmkParamReader.read(snakemake_params)
        snakemake_params.forcerun = self.get_forceRunList()
        if len(snakemake_params.forcerun) > 0:
            delete_dependencies(snakemake_params)
        background_tasks.add_task(snakemake_execute, snakemake_params)

    def create_workflow(self):
        rules, last_output = self.rulefile()

        flow_config = FlowConfig(
            rules=rules,
            last_output=last_output,
        )

        SmkConfigWriter.write(self.unique_id, asdict(flow_config))

    def rulefile(self):
        endNodeList = self.get_endNodeList()

        nwbfile = get_typecheck_params(self.runItem.nwbParam, "nwb")

        rule_dict: Dict[str, Rule] = {}
        last_outputs = []

        for node in self.nodeDict.values():
            if node.type == NodeType.IMAGE:
                rule_dict[node.id] = SmkSetfile.image(self.unique_id, node, self.edgeDict, nwbfile)
            elif node.type == NodeType.CSV:
                rule_dict[node.id] = SmkSetfile.csv(self.unique_id, node, self.edgeDict, nwbfile)
            elif node.type == NodeType.FLUO:
                rule_dict[node.id] = SmkSetfile.csv(self.unique_id, node, self.edgeDict, nwbfile)
            elif node.type == NodeType.BEHAVIOR:
                rule_dict[node.id] = SmkSetfile.csv(self.unique_id, node, self.edgeDict, nwbfile, "behavior")
            elif node.type == NodeType.HDF5:
                rule_dict[node.id] = SmkSetfile.hdf5(self.unique_id, node, self.edgeDict, nwbfile)
            elif node.type == NodeType.ALGO:
                rule = SmkSetfile.algo(self.unique_id, node, self.edgeDict, self.nodeDict)
                rule_dict[node.id] = rule

                if node.id in endNodeList:
                    last_outputs.append(rule.output)
            else:
                assert False, "NodeType doesn't exists"

        return rule_dict, last_outputs

    def get_endNodeList(self):
        returnCntDict = {key: 0 for key in self.nodeDict.keys()}
        for edge in self.edgeDict.values():
            returnCntDict[edge.source] += 1

        endNodeList = []
        for key, value in returnCntDict.items():
            if value == 0:
                endNodeList.append(key)
        return endNodeList

    def get_forceRunList(self) -> List[str]:
        target_list = []
        for x in self.runItem.forceRunList:
            target_list.append(
                get_pickle_file(
                    self.unique_id,
                    x.nodeId,
                    x.name
                )
            )
        return target_list
