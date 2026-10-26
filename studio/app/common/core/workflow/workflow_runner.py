import uuid
from dataclasses import asdict
from datetime import datetime
from typing import Dict, List

from fastapi import BackgroundTasks

from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.experiment.experiment_record_services import (
    ExperimentRecordService,
)
from studio.app.common.core.experiment.experiment_writer import ExptConfigWriter
from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.snakemake.smk import FlowConfig, Rule, SmkParam
from studio.app.common.core.snakemake.snakemake_executor import (
    delete_dependencies,
    snakemake_execute,
)
from studio.app.common.core.snakemake.snakemake_reader import SmkParamReader
from studio.app.common.core.snakemake.snakemake_rule import SmkRule
from studio.app.common.core.snakemake.snakemake_writer import SmkConfigWriter
from studio.app.common.core.workflow.workflow import (
    NodeType,
    NodeTypeUtil,
    OutputPath,
    OutputType,
    RunItem,
    WorkflowRunStatus,
)
from studio.app.common.core.workflow.workflow_params import get_typecheck_params
from studio.app.common.core.workflow.workflow_writer import WorkflowConfigWriter
from studio.app.const import DATE_FORMAT


class WorkflowRunner:
    def __init__(self, workspace_id: str, unique_id: str, runItem: RunItem) -> None:
        self.workspace_id = workspace_id
        self.unique_id = unique_id
        self.runItem = runItem
        self.nodeDict = self.runItem.nodeDict
        self.edgeDict = self.runItem.edgeDict

        WorkflowConfigWriter(
            self.workspace_id,
            self.unique_id,
            self.nodeDict,
            self.edgeDict,
        ).write()

        ExptConfigWriter(
            self.workspace_id,
            self.unique_id,
            self.runItem.name,
            nwbfile=get_typecheck_params(self.runItem.nwbParam, "nwb"),
            snakemake=get_typecheck_params(self.runItem.snakemakeParam, "snakemake"),
        ).write()

        Runner.clear_pid_file(self.workspace_id, self.unique_id)

    @staticmethod
    def create_workflow_unique_id() -> str:
        new_unique_id = str(uuid.uuid4())[:8]
        return new_unique_id

    def run_workflow(self, background_tasks: BackgroundTasks):
        self.set_smk_config()

        snakemake_params: SmkParam = get_typecheck_params(
            self.runItem.snakemakeParam, "snakemake"
        )
        snakemake_params = SmkParamReader.read(snakemake_params)
        snakemake_params.forcerun = self.runItem.forceRunList
        if len(snakemake_params.forcerun) > 0:
            delete_dependencies(
                workspace_id=self.workspace_id,
                unique_id=self.unique_id,
                smk_params=snakemake_params,
                nodeDict=self.nodeDict,
                edgeDict=self.edgeDict,
            )

        background_tasks.add_task(
            snakemake_execute, self.workspace_id, self.unique_id, snakemake_params
        )

    def finish_workflow_without_run(
        self, status: WorkflowRunStatus = WorkflowRunStatus.SUCCESS
    ):
        """
        Saves the settings and finishes the workflow without actually running it.
        - Function solely for creating experiment record.
        """

        # Load current configs
        expt_config = ExptConfigReader.read(self.workspace_id, self.unique_id)

        # Construct update data (ExptConfig.*)
        update_expt_config = ExptConfigReader.create_empty_experiment_config()
        now = datetime.now().strftime(DATE_FORMAT)
        update_expt_config.success = status.value
        update_expt_config.finished_at = now
        update_expt_config.data_usage = 0

        # Construct update data (ExptConfig.function)
        update_expt_config.function = {}
        for node_id, function in expt_config.function.items():
            function.success = WorkflowRunStatus.SUCCESS.value
            function.outputPaths = {
                "empty": OutputPath(
                    path="empty",
                    type=OutputType.EMPTY,
                    max_index=1,
                )
            }

            update_expt_config.function[node_id] = function

        # Prepare data (dict variable) for overwriting the config file
        update_expt_config_dict = {
            k: v for k, v in asdict(update_expt_config).items() if v is not None
        }

        # Overwrite config file
        ExptConfigWriter(self.workspace_id, self.unique_id).overwrite(
            update_expt_config_dict
        )

        # Update experiment database record
        if ExperimentRecordService.is_available():
            ExperimentRecordService.regist_record_on_workflow_completed(
                self.workspace_id, self.unique_id
            )

    def set_smk_config(self):
        rules, last_output = self.rulefile()

        nwb_template = get_typecheck_params(self.runItem.nwbParam, "nwb")

        # Get client_id from current logging context to pass to snakemake workflow
        from studio.app.common.core.logger_context_helpers import (
            get_client_id_for_subprocess,
        )

        client_id = get_client_id_for_subprocess()

        flow_config = FlowConfig(
            rules=rules,
            last_output=last_output,
            nwb_template=nwb_template,
            client_id=client_id,
        )

        SmkConfigWriter.write_raw(
            self.workspace_id, self.unique_id, asdict(flow_config)
        )

    def rulefile(self) -> Dict[str, Rule]:
        endNodeList = self.get_endNodeList()

        nwbfile = get_typecheck_params(self.runItem.nwbParam, "nwb")

        rule_dict: Dict[str, Rule] = {}
        last_outputs = []

        for node in self.nodeDict.values():
            if NodeTypeUtil.check_nodetype(node.type) == NodeType.DATA:
                data_common_rule = SmkRule(
                    workspace_id=self.workspace_id,
                    unique_id=self.unique_id,
                    node=node,
                    edgeDict=self.edgeDict,
                    nwbfile=nwbfile,
                )
                data_rule = None

                if node.type == NodeType.IMAGE:
                    data_rule = data_common_rule.image()
                elif node.type == NodeType.CSV:
                    data_rule = data_common_rule.csv()
                elif node.type == NodeType.FLUO:
                    data_rule = data_common_rule.csv()
                elif node.type == NodeType.BEHAVIOR:
                    data_rule = data_common_rule.csv(nodeType="behavior")
                elif node.type == NodeType.HDF5:
                    data_rule = data_common_rule.hdf5()
                elif node.type == NodeType.MATLAB:
                    data_rule = data_common_rule.mat()
                elif node.type == NodeType.MICROSCOPE:
                    data_rule = data_common_rule.microscope()

                rule_dict[node.id] = data_rule

            elif NodeTypeUtil.check_nodetype(node.type) == NodeType.ALGO:
                algo_rule = SmkRule(
                    workspace_id=self.workspace_id,
                    unique_id=self.unique_id,
                    node=node,
                    edgeDict=self.edgeDict,
                ).algo(nodeDict=self.nodeDict)

                rule_dict[node.id] = algo_rule

                if node.id in endNodeList:
                    last_outputs.append(algo_rule.output)
            else:
                assert False, f"NodeType doesn't exists: {node.type}"

        return rule_dict, last_outputs

    def get_endNodeList(self) -> List[str]:
        returnCntDict = {key: 0 for key in self.nodeDict.keys()}
        for edge in self.edgeDict.values():
            returnCntDict[edge.source] += 1

        endNodeList = []
        for key, value in returnCntDict.items():
            if value == 0:
                endNodeList.append(key)
        return endNodeList
