import os
import signal
from dataclasses import asdict
from datetime import datetime
from glob import glob
from typing import Dict

from fastapi import HTTPException

from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.experiment.experiment_writer import ExptConfigWriter
from studio.app.common.core.snakemake.smk_status_logger import SmkStatusLogger
from studio.app.common.core.utils.file_reader import JsonReader
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.pickle_handler import PickleReader
from studio.app.common.core.workflow.workflow import Message, NodeRunStatus, OutputPath
from studio.app.common.dataclass import BaseData
from studio.app.const import DATE_FORMAT
from studio.app.dir_path import DIRPATH


class WorkflowResult:
    def __init__(self, workspace_id, unique_id):
        self.workspace_id = workspace_id
        self.unique_id = unique_id
        self.workflow_dirpath = join_filepath(
            [
                DIRPATH.OUTPUT_DIR,
                self.workspace_id,
                self.unique_id,
            ]
        )
        self.expt_filepath = join_filepath(
            [self.workflow_dirpath, DIRPATH.EXPERIMENT_YML]
        )

        self.pid_filepath = join_filepath([self.workflow_dirpath, "pid.json"])

    def observe(self, nodeIdList) -> Dict:
        """
        Perform the following operations for the specified workflow
          - Check and update the workflow execution status
          - Response with the confirmed workflow execution status
        """

        results: Dict[str, Message] = {}

        # check for workflow errors
        workflow_error = SmkStatusLogger.get_error_content(
            self.workspace_id, self.unique_id
        )

        # observe node list
        for node_id in nodeIdList:
            # Cases with errors in workflow
            if workflow_error["has_error"]:
                node_pickle_path = None
                node_result = NodeResult(
                    self.workspace_id,
                    self.unique_id,
                    node_id,
                    node_pickle_path,
                    workflow_error,
                )
                results[node_id] = node_result.observe()

            # Normal case
            else:
                # search node pickle files
                node_dirpath = join_filepath([self.workflow_dirpath, node_id])
                node_pickle_files = list(
                    set(glob(join_filepath([node_dirpath, "*.pkl"])))
                    - set(glob(join_filepath([node_dirpath, "tmp_*.pkl"])))
                )

                # process node pickle files
                for node_pickle_path in node_pickle_files:
                    node_result = NodeResult(
                        self.workspace_id,
                        self.unique_id,
                        node_id,
                        node_pickle_path,
                    )
                    results[node_id] = node_result.observe()

                    self.__check_has_nwb(node_id)

        # check workflow nwb
        self.__check_has_nwb()

        return results

    def __check_has_nwb(self, node_id=None):
        target_whole_nwb = node_id is None

        if target_whole_nwb:
            nwb_filepath_list = glob(join_filepath([self.workflow_dirpath, "*.nwb"]))
        else:
            node_dirpath = join_filepath([self.workflow_dirpath, node_id])
            nwb_filepath_list = glob(join_filepath([node_dirpath, "*.nwb"]))

        for nwb_filepath in nwb_filepath_list:
            if os.path.exists(nwb_filepath):
                config = ExptConfigReader.read(self.expt_filepath)

                if target_whole_nwb:
                    config.hasNWB = True
                else:
                    config.function[node_id].hasNWB = True

                # Update EXPERIMENT_YML
                ExptConfigWriter.write_raw(
                    self.workspace_id, self.unique_id, asdict(config)
                )

    def cancel(self):
        """
        The algorithm function of this workflow is being executed at the line:
        https://github.com/snakemake/snakemake/blob/27b224ed12448df8aebc7d1ff8f25e3bf7622232/snakemake/shell.py#L258
        ```
        proc = sp.Popen(
            cmd,
            bufsize=-1,
            shell=use_shell,
            stdout=stdout,
            universal_newlines=iterable or read or None,
            close_fds=close_fds,
            **cls._process_args,
            env=envvars,
        )
        ```
        The `cmd` argument has the following format:
        ```
        source ~/miniconda3/bin/activate
        '~/Documents/optinistfs/.snakemake/conda/491889952d2f07f3876bb801eea626e9_';
        set -euo pipefail;
        python ~/Documents/optinistfs/.snakemake/scripts/tmp03froqxo.func.py
        ```
        Interrupt the conda activate at the beginning of the process is impossible
        because it is only called when each algorithm function executes.
        This workflow is cancelled by killing process via PID of algorithm function
        saved in pid.json file
        Raises:
            HTTPException: if pid_filepath or last_script_file does not exist
        """
        if not os.path.exists(self.pid_filepath):
            raise HTTPException(status_code=404)

        pid_data = JsonReader.read(self.pid_filepath)

        if not os.path.exists(pid_data["last_script_file"]):
            raise HTTPException(status_code=404)

        os.remove(pid_data["last_script_file"])
        os.kill(pid_data["last_pid"], signal.SIGTERM)

        return True


class NodeResult:
    def __init__(
        self,
        workspace_id: str,
        unique_id: str,
        node_id: str,
        pickle_filepath: str,
        workflow_error: dict = None,
    ):
        self.workspace_id = workspace_id
        self.unique_id = unique_id
        self.workflow_dirpath = join_filepath(
            [
                DIRPATH.OUTPUT_DIR,
                self.workspace_id,
                self.unique_id,
            ]
        )
        self.node_id = node_id
        self.node_dirpath = join_filepath([self.workflow_dirpath, self.node_id])
        self.expt_filepath = join_filepath(
            [self.workflow_dirpath, DIRPATH.EXPERIMENT_YML]
        )
        self.workflow_has_error = (
            workflow_error["has_error"] if workflow_error else False
        )
        self.workflow_error_log = (
            workflow_error["error_log"] if workflow_error else None
        )

        if not self.workflow_has_error:
            pickle_filepath = pickle_filepath.replace("\\", "/")
            self.algo_name = os.path.splitext(os.path.basename(pickle_filepath))[0]
            try:
                self.info = PickleReader.read(pickle_filepath)
            except EOFError:
                self.info = None  # indicates error
        else:
            self.algo_name = None
            self.info = None

    def observe(self) -> Message:
        expt_config = ExptConfigReader.read(self.expt_filepath)

        # case) error throughout workflow
        if self.workflow_has_error:
            message = self.error(self.workflow_error_log)
        # case) error in node
        elif PickleReader.check_is_error_node_pickle(self.info):
            message = self.error()
        # case) success in node
        else:
            message = self.success()
            expt_config.function[self.node_id].outputPaths = message.outputPaths

        expt_config.function[self.node_id].success = message.status

        now = datetime.now().strftime(DATE_FORMAT)
        expt_config.function[self.node_id].finished_at = now
        expt_config.function[self.node_id].message = message.message

        statuses = list(map(lambda x: x.success, expt_config.function.values()))

        if NodeRunStatus.RUNNING.value not in statuses:
            expt_config.finished_at = now
            if NodeRunStatus.ERROR.value in statuses:
                expt_config.success = NodeRunStatus.ERROR.value
            else:
                expt_config.success = NodeRunStatus.SUCCESS.value

        # Update EXPERIMENT_YML
        ExptConfigWriter.write_raw(
            self.workspace_id, self.unique_id, asdict(expt_config)
        )
        return message

    def success(self) -> Message:
        return Message(
            status=NodeRunStatus.SUCCESS.value,
            message=f"{self.algo_name} success",
            outputPaths=self.output_paths(),
        )

    def error(self, message: str = None) -> Message:
        if message is None:
            if self.info is None:
                message = "Invalid node result info: None"
            else:
                message = (
                    "\n".join(self.info) if isinstance(self.info, list) else self.info
                )

        return Message(status=NodeRunStatus.ERROR.value, message=message)

    def output_paths(self) -> dict:
        output_paths: Dict[str, OutputPath] = {}
        for k, v in self.info.items():
            if isinstance(v, BaseData):
                v.save_json(self.node_dirpath)
                if v.output_path:
                    output_paths[k] = v.output_path

        return output_paths
