import os
import re
import signal
import time
from dataclasses import asdict
from datetime import datetime
from glob import glob
from typing import Dict, List

from fastapi import HTTPException, status
from psutil import AccessDenied, NoSuchProcess, Process, ZombieProcess, process_iter

from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.experiment.experiment_writer import ExptConfigWriter
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.snakemake.smk_status_logger import SmkStatusLogger
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.pickle_handler import PickleReader
from studio.app.common.core.workflow.workflow import Message, NodeRunStatus, OutputPath
from studio.app.common.dataclass import BaseData
from studio.app.common.schemas.workflow import (
    WorkflowErrorInfo,
    WorkflowPIDFileData,
    WorkflowProcessInfo,
)
from studio.app.const import DATE_FORMAT
from studio.app.dir_path import DIRPATH

logger = AppLogger.get_logger()


class WorkflowResult:
    def __init__(self, workspace_id: str, unique_id: str):
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
        self.monitor = WorkflowMonitor(workspace_id, unique_id)

    def observe(self, nodeIdList: List[str]) -> Dict:
        """
        Perform the following operations for the specified workflow
          - Check and update the workflow execution status
          - Response with the confirmed workflow execution status
        """
        # validate args
        if not nodeIdList:
            return {}

        # check for workflow errors
        workflow_error = SmkStatusLogger.get_error_content(
            self.workspace_id, self.unique_id
        )

        # observe node list
        results = self.__observe_node_list(nodeIdList, workflow_error)

        # check workflow status
        is_workflow_status_running = self.__is_workflow_status_running(
            nodeIdList, results
        )

        # If the workflow status is running (workflow is incomplete),
        # check whether the actual process exists.
        if is_workflow_status_running:
            # check workflow process exists
            current_process = self.monitor.search_process()

            # error handling for process not found
            if current_process is None:
                workflow_error = WorkflowErrorInfo(
                    has_error=True, error_log="No Snakemake process found."
                )

            # re-run observe node list (reflects workflow error)
            results = self.__observe_node_list(nodeIdList, workflow_error)

        return results

    def __observe_node_list(
        self, nodeIdList: List[str], workflow_error: WorkflowErrorInfo
    ) -> Dict[str, Message]:
        results: Dict[str, Message] = {}

        for node_id in nodeIdList:
            # Cases with errors in workflow
            if workflow_error.has_error:
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

    def __is_workflow_status_running(
        self, nodeIdList: List[str], messages: Dict[str, Message]
    ) -> bool:
        """
        By comparing the number of nodeIdList waiting for processing completion
        with the number of nodeIdList that has completed processing immediately before,
        is_running is determined.
        """
        is_running = len(nodeIdList) != len(messages.keys())

        logger.debug(
            "check wornflow running status "
            f"[{self.workspace_id}/{self.unique_id}] [is_running: {is_running}]"
        )

        return is_running

    def __check_has_nwb(self, node_id=None) -> None:
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


class NodeResult:
    def __init__(
        self,
        workspace_id: str,
        unique_id: str,
        node_id: str,
        pickle_filepath: str,
        workflow_error: WorkflowErrorInfo = None,
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
        self.workflow_has_error = workflow_error.has_error if workflow_error else False
        self.workflow_error_log = workflow_error.error_log if workflow_error else None

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


class WorkflowMonitor:
    PROCESS_SNAKEMAKE_CMDLINE = "python .*/\\.snakemake/scripts/"
    PROCESS_SNAKEMAKE_WAIT_TIMEOUT = 7200  # sec
    PROCESS_CONDA_CMDLINE = "conda env create .*/\\.snakemake/conda/"
    PROCESS_CONDA_WAIT_TIMEOUT = 3600  # sec

    def __init__(self, workspace_id: str, unique_id: str):
        self.workspace_id = workspace_id
        self.unique_id = unique_id
        self.expt_filepath = join_filepath(
            [
                DIRPATH.OUTPUT_DIR,
                self.workspace_id,
                self.unique_id,
                DIRPATH.EXPERIMENT_YML,
            ]
        )

    def search_process(self) -> WorkflowProcessInfo:
        pid_data = Runner.read_pid_file(self.workspace_id, self.unique_id)
        if pid_data is None:
            logger.warning(
                f"No workflow pid file found. [{self.workspace_id}/{self.unique_id}]"
            )

            # Refer experiment_data instead of pid_data
            expt_config = ExptConfigReader.read(self.expt_filepath)
            try:
                expt_started_time = datetime.strptime(
                    expt_config.started_at, DATE_FORMAT
                )
            except ValueError:
                expt_started_time = datetime.fromtimestamp(0)

            # Set dummy value to proceed to the next step.
            pid_data = WorkflowPIDFileData(
                last_pid=999999,
                last_script_file="__dummy_wrapper_function.py",
                create_time=expt_started_time.timestamp(),
            )

        process_data: WorkflowProcessInfo = None

        # Find the process corresponding to the pid in pid_data
        try:
            last_pid = pid_data.last_pid

            # get process
            process = Process(last_pid)
            logger.info(f"Found workflow process. {process}")

            # validate process name
            process_cmdline = " ".join(process.cmdline())
            if not re.search(self.PROCESS_SNAKEMAKE_CMDLINE, process_cmdline):
                logger.warning(
                    "Found another process with same PID:"
                    f" [{last_pid}] [{process_cmdline}]"
                )
                raise NoSuchProcess(last_pid)

            process_data = WorkflowProcessInfo(process=process, pid_data=pid_data)

        # If the target process does not exist,
        # check for the existence of the `conda env create` command process.
        except NoSuchProcess:
            logger.warning(f"No workflow process found. {pid_data}")

            # Search for the existence of a conda command process ("conda env create")
            conda_process = None
            for proc in process_iter(["pid", "name", "cmdline"]):
                try:
                    cmdline = proc.info.get("cmdline")
                    cmdline = " ".join(cmdline) if cmdline else ""

                    if re.search(self.PROCESS_CONDA_CMDLINE, cmdline):
                        conda_ps_create_elapsed = int(time.time() - proc.create_time())
                        logger.info(
                            f"Found conda process. [{proc}] [{cmdline}] "
                            f"[{conda_ps_create_elapsed} sec]",
                        )
                        conda_process = Process(proc.pid)

                        # Check elapsed time for process startup
                        #
                        # ATTENTION:
                        # The conda command process is a separate process from
                        #   the snakemake process (although it is a child process),
                        #   so it is difficult to identify the process with certainty.
                        # Therefore, the process start time is used here to determine
                        #   the process by estimation.
                        if conda_ps_create_elapsed < self.PROCESS_CONDA_WAIT_TIMEOUT:
                            process_data = WorkflowProcessInfo(
                                process=conda_process, pid_data=pid_data
                            )
                        else:
                            logger.warning(
                                "This conda command is "
                                "probably an irrelevant process.."
                                f"[{conda_process}] [{conda_ps_create_elapsed} sec]"
                            )
                    else:
                        continue  # skip that process

                except AccessDenied:
                    continue  # skip that process
                except ZombieProcess:
                    continue  # skip that process

        # Rescue action when process not found
        if process_data is None:
            # Check elapsed time for process startup
            # *Retry for a certain period of time even if process not found
            if pid_data.elapsed_time_float < self.PROCESS_SNAKEMAKE_WAIT_TIMEOUT:
                logger.debug(f"Set dummy workflow process tentatively. [{pid_data}]")
                process_data = WorkflowProcessInfo(process=None, pid_data=pid_data)
            else:
                logger.warning(f"No workflow process found at all. [{pid_data}]")

        return process_data

    def cancel_run(self):
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
        '/app/.snakemake/conda/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx_';
        set -euo pipefail;
        python /app/.snakemake/scripts/tmpxxxxxxxx.func.py
        ```
        Interrupt the conda activate at the beginning of the process is impossible
        because it is only called when each algorithm function executes.
        This workflow is cancelled by killing process via PID of algorithm function
        saved in `RUN_PROCESS_PID_FILE` file
        Raises:
            HTTPException: if pid_filepath or last_script_file does not exist
        """

        current_process = self.search_process()
        if current_process is None:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Current process not found",
            )
        elif current_process.process is None:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Current process not found",
            )
        pid_data = current_process.pid_data

        if os.path.exists(pid_data.last_script_file):
            # force remove run script file
            os.remove(pid_data.last_script_file)

        else:
            logger.warning(
                "The run script has not yet started. "
                f"script: {pid_data.last_script_file}"
            )

        # send kill to process
        current_process.process.send_signal(signal.SIGTERM)

        return True
