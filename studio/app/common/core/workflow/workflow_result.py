import copy
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

from studio.app.common.core.experiment.experiment import ExptConfig, ExptFunction
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
        self.monitor = WorkflowMonitor(workspace_id, unique_id)

    def observe(self, observe_node_ids: List[str]) -> Dict[str, Message]:
        """
        Perform the following operations for the specified workflow
          - Check and update the workflow execution status
          - Response with the confirmed workflow execution status
        """
        expt_config = ExptConfigReader.read(self.workspace_id, self.unique_id)

        # validate args
        if not observe_node_ids:
            return {}

        # check for workflow errors
        workflow_error = SmkStatusLogger.get_error_content(
            self.workspace_id, self.unique_id
        )

        # observe node list
        node_results = self.__observe_nodes(
            observe_node_ids, expt_config, workflow_error
        )

        # check workflow status
        is_workflow_observation_ongoing = self.__is_workflow_observation_ongoing(
            observe_node_ids, node_results
        )

        # If the workflow status observation is ongoing (maybe workflow is incomplete),
        # check whether the actual process exists.
        if is_workflow_observation_ongoing:
            # check workflow process exists
            current_process = self.monitor.search_process()

            # error handling for process not found
            if current_process is None:
                workflow_error = WorkflowErrorInfo(
                    has_error=True, error_log="No Snakemake process found."
                )

            # re-run observe node list (reflects workflow error)
            node_results = self.__observe_nodes(
                observe_node_ids, expt_config, workflow_error
            )

        return node_results

    def observe_overall(self) -> Dict[str, Message]:
        """
        Automatically observe all nodes
        """
        expt_config = ExptConfigReader.read(self.workspace_id, self.unique_id)

        # Observe all nodes
        observe_node_ids = expt_config.function.keys()

        return self.observe(observe_node_ids)

    def __observe_nodes(
        self,
        observe_node_ids: List[str],
        expt_config: ExptConfig,
        workflow_error: WorkflowErrorInfo,
    ) -> Dict[str, Message]:
        """
        Observe by Nodes
        """
        node_results: Dict[str, Message] = {}

        for node_id in observe_node_ids:
            # Check node processed status (inspect ExptConfig)
            expt_function = expt_config.function.get(node_id)
            if not expt_function:
                logger.warning(f"Invalid node_id [{node_id}]")
                continue

            # Perform observations of nodes
            node_result = NodeResult(
                self.workspace_id,
                self.unique_id,
                node_id,
                workflow_error=workflow_error,
            )
            if node_result.is_ready():
                node_results[node_id] = node_result.observe(expt_function)
            else:
                # # debug log.
                # logger.debug(
                #     "Target node not yet completed " f"[id: {node_id}]"
                # )
                pass

        # Check if whole.nwb file exists
        if not NodeResult.is_all_nodes_already_finished(expt_config):
            self.__check_has_whole_nwb()

        return node_results

    def __is_workflow_observation_ongoing(
        self, observe_node_ids: List[str], messages: Dict[str, Message]
    ) -> bool:
        """
        By comparing the number of observe_node_ids waiting for processing completion
        with the number of observe_node_ids that has completed processing immediately
        before, is_running is determined.
        """
        is_ongoing = len(observe_node_ids) != len(messages.keys())

        # # debug log.
        # logger.debug(
        #     "check wornflow ongoing status "
        #     f"[{self.workspace_id}/{self.unique_id}] [is_ongoing: {is_ongoing}]"
        # )

        return is_ongoing

    def __check_has_whole_nwb(self, node_id=None) -> None:
        nwb_filepath_list = glob(join_filepath([self.workflow_dirpath, "*.nwb"]))
        nwb_filepath = nwb_filepath_list[0] if nwb_filepath_list else None
        hasNWB = os.path.exists(nwb_filepath) if nwb_filepath else False

        # ------------------------------------------------------------
        # Update process of node processing status
        # ------------------------------------------------------------

        update_expt_config = ExptConfigReader.create_empty_experiment_config()

        update_expt_config.hasNWB = hasNWB

        # Make overwrite params
        update_expt_config_dict = {}
        if update_expt_config.hasNWB is not None:
            update_expt_config_dict["hasNWB"] = update_expt_config.hasNWB

        # Overwrite experiment config
        ExptConfigWriter(self.workspace_id, self.unique_id).overwrite(
            update_expt_config_dict
        )


class NodeResult:
    def __init__(
        self,
        workspace_id: str,
        unique_id: str,
        node_id: str,
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
        self.workflow_has_error = workflow_error.has_error if workflow_error else False
        self.workflow_error_log = workflow_error.error_log if workflow_error else None

        # Read node's output pickle file
        self.algo_name = None
        self.info = None
        if not self.workflow_has_error:
            node_pickle_path = PickleReader.search_node_pickle_path(self.node_dirpath)
            if node_pickle_path:
                node_pickle_path = node_pickle_path.replace("\\", "/")
                self.algo_name = os.path.splitext(os.path.basename(node_pickle_path))[0]
                try:
                    self.info = PickleReader.read(node_pickle_path)
                except Exception as e:
                    self.info = None  # indicates error
                    logger.error(e, exc_info=True)

    def observe(self, expt_function: ExptFunction) -> Message:
        # Generate node process status (workflow.Message)
        # case) error throughout workflow
        if self.workflow_has_error:
            message = self.error(self.workflow_error_log)
        # case) error in node
        elif PickleReader.check_is_error_node_pickle(self.info):
            message = self.error()
        # case) success in node
        else:
            message = self.success()

        # Determine if the node has already been processed
        # *If it has, skip subsequent ExptConfig update process.
        is_already_finished: bool = __class__.is_node_already_finished(expt_function)
        if is_already_finished:
            # # debug log.
            # logger.debug(
            #     "Target node has already been finished "
            #     f"[id: {self.node_id}][status: {expt_function.success}]"
            # )
            return message

        # ------------------------------------------------------------
        # Update process of node processing status
        # ------------------------------------------------------------

        # ----------------------------------------
        # Generate update parameters
        # ----------------------------------------

        original_expt_config = ExptConfigReader.read(
            self.workspace_id,
            self.unique_id,
        )
        update_expt_config = ExptConfigReader.create_empty_experiment_config()
        update_config_function: ExptFunction = copy.deepcopy(
            original_expt_config.function[self.node_id]
        )

        if message.status == NodeRunStatus.SUCCESS.value:
            update_config_function.outputPaths = message.outputPaths
        update_config_function.success = message.status

        # TODO: Set started_at on the node
        #   At the moment, there is insufficient data to obtain an accurate started_at,
        #  so set it to None.
        #   separate modification is required to record running info for each node.
        update_config_function.started_at = None

        now = datetime.now().strftime(DATE_FORMAT)
        update_config_function.finished_at = now
        update_config_function.message = message.message

        latest_statuses = [
            (
                update_config_function.success
                if k == update_config_function.unique_id
                else v.success
            )
            for k, v in original_expt_config.function.items()
        ]

        # Status check of each node
        if NodeRunStatus.RUNNING.value not in latest_statuses:
            update_expt_config.finished_at = now
            if NodeRunStatus.ERROR.value in latest_statuses:
                update_expt_config.success = NodeRunStatus.ERROR.value
            else:
                update_expt_config.success = NodeRunStatus.SUCCESS.value

        # Check if nwb file exists
        update_config_function.hasNWB = self.__check_has_nwb()

        # ----------------------------------------
        # ExptConfig update commit process
        # ----------------------------------------

        # Make overwrite params
        update_expt_config_dict = {}
        if update_config_function is not None:
            update_expt_config_dict["function"] = {
                update_config_function.unique_id: asdict(update_config_function)
            }
        if update_expt_config.success is not None:
            update_expt_config_dict["success"] = update_expt_config.success
        if update_expt_config.finished_at is not None:
            update_expt_config_dict["finished_at"] = update_expt_config.finished_at

        # Overwrite experiment config
        ExptConfigWriter(self.workspace_id, self.unique_id).overwrite(
            update_expt_config_dict
        )

        logger.debug(
            "Node observation completed "
            f"[id: {self.node_id}][status: {update_expt_config.success}]"
        )

        return message

    def is_ready(self) -> bool:
        is_ready = (self.info is not None) or self.workflow_has_error
        return is_ready

    @classmethod
    def is_node_already_finished(cls, expt_function: ExptFunction) -> bool:
        return expt_function.success != NodeRunStatus.RUNNING.value

    @classmethod
    def is_all_nodes_already_finished(cls, expt_config: ExptConfig) -> bool:
        expt_functions = expt_config.function.values()
        return all([cls.is_node_already_finished(v) for v in expt_functions])

    def __check_has_nwb(self) -> bool:
        nwb_filepath_list = glob(join_filepath([self.node_dirpath, "*.nwb"]))
        nwb_filepath = nwb_filepath_list[0] if nwb_filepath_list else None
        hasNWB = os.path.exists(nwb_filepath) if nwb_filepath else False

        return hasNWB

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
    PROCESS_SEARCH_NAMES = ["python", "conda"]
    PROCESS_SNAKEMAKE_CMDLINE = "\\b(?:python)\\b.*/\\.snakemake/scripts/"
    PROCESS_SNAKEMAKE_WAIT_TIMEOUT = 7200  # sec
    PROCESS_CONDA_CMDLINE = (
        "\\b(?:conda)\\b.*\\b(?:env)\\b.*\\b(?:create)\\b.*/\\.snakemake/conda/"
    )
    PROCESS_CONDA_WAIT_TIMEOUT = 3600  # sec

    def __init__(self, workspace_id: str, unique_id: str):
        self.workspace_id = workspace_id
        self.unique_id = unique_id

    def search_process(self) -> WorkflowProcessInfo:
        pid_data = Runner.read_pid_file(self.workspace_id, self.unique_id)
        if pid_data is None:
            # ATTENTION:
            # It should be a warning, but since it matches frequently,
            # it is temporarily set to debug.
            # logger.debug(
            #     f"No workflow pid file found. [{self.workspace_id}/{self.unique_id}]"
            # )

            # Refer experiment_data instead of pid_data
            expt_config = ExptConfigReader.read(self.workspace_id, self.unique_id)
            try:
                expt_started_time = datetime.strptime(
                    expt_config.started_at, DATE_FORMAT
                )
            except ValueError:
                expt_started_time = datetime.fromtimestamp(0)

            # Set dummy value to proceed to the next step.
            pid_data = WorkflowPIDFileData(
                last_pid=999999,
                func_name="__dummy_wrapper",
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
            process_cmdline = " ".join(process.cmdline()).replace("\\", "/")
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
            # ATTENTION:
            # It should be a warning, but since it matches frequently,
            # it is temporarily set to debug.
            # logger.debug(f"No workflow process found. {pid_data}")

            # Search for the existence of a conda command process ("conda env create")
            conda_process = None
            for proc in process_iter(["pid", "name"]):
                try:
                    # Targeting specific programs for backlog (for performance)
                    proc_name = proc.info.get("name")
                    if not any(v in proc_name for v in self.PROCESS_SEARCH_NAMES):
                        continue

                    # Get cmdline info
                    # Note:
                    #   Since "cmdline" is not specified in `process_iter`
                    #   (for performance), the cmdline is obtained using "proc.as_dict".
                    cmdline = proc.as_dict(attrs=["cmdline"]).get("cmdline")
                    cmdline = " ".join(cmdline) if cmdline else ""
                    cmdline = cmdline.replace("\\", "/")

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
                # logger.debug(f"Set dummy workflow process tentatively. [{pid_data}]")
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
