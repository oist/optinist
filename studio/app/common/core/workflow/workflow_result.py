import os
from dataclasses import asdict
from glob import glob
from typing import Dict

from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.utils.config_handler import ConfigWriter
from studio.app.common.core.utils.file_reader import Reader
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.pickle_handler import PickleReader
from studio.app.common.core.workflow.workflow import Message, OutputPath
from studio.app.common.dataclass import BaseData
from studio.app.dir_path import DIRPATH


class WorkflowResult:
    def __init__(self, unique_id):
        self.unique_id = unique_id
        self.workflow_dirpath = join_filepath(
            [
                DIRPATH.OUTPUT_DIR,
                self.unique_id,
            ]
        )
        self.expt_filepath = join_filepath(
            [self.workflow_dirpath, DIRPATH.EXPERIMENT_YML]
        )
        self.error_filepath = join_filepath([self.workflow_dirpath, "error.log"])

    def get(self, nodeIdList):
        results: Dict[str, Message] = {}
        for node_id in nodeIdList:
            if os.path.exists(self.error_filepath):
                error_message = Reader.read(self.error_filepath)
                if error_message != "":
                    results[node_id] = Message(
                        status="error",
                        message=error_message,
                    )

            glob_pickle_filepath = join_filepath(
                [self.workflow_dirpath, node_id, "*.pkl"]
            )
            for pickle_filepath in glob(glob_pickle_filepath):
                results[node_id] = NodeResult(
                    self.workflow_dirpath,
                    node_id,
                    pickle_filepath,
                ).get()
                self.has_nwb(node_id)

        self.has_nwb()

        return results

    def has_nwb(self, node_id=None):
        if node_id is None:
            nwb_filepath_list = glob(join_filepath([self.workflow_dirpath, "*.nwb"]))
        else:
            nwb_filepath_list = glob(
                join_filepath([self.workflow_dirpath, node_id, "*.nwb"])
            )

        for nwb_filepath in nwb_filepath_list:
            if os.path.exists(nwb_filepath):
                config = ExptConfigReader.read(self.expt_filepath)

                if node_id is None:
                    config.hasNWB = True
                else:
                    config.function[node_id].hasNWB = True

                ConfigWriter.write(
                    dirname=self.workflow_dirpath,
                    filename=DIRPATH.EXPERIMENT_YML,
                    config=asdict(config),
                )


class NodeResult:
    def __init__(self, workflow_dirpath, node_id, pickle_filepath):
        self.workflow_dirpath = workflow_dirpath
        self.node_id = node_id
        self.node_dirpath = join_filepath([self.workflow_dirpath, self.node_id])
        self.expt_filepath = join_filepath(
            [self.workflow_dirpath, DIRPATH.EXPERIMENT_YML]
        )

        pickle_filepath = pickle_filepath.replace("\\", "/")
        self.algo_name = os.path.splitext(os.path.basename(pickle_filepath))[0]
        self.info = PickleReader.read(pickle_filepath)

    def get(self):
        expt_config = ExptConfigReader.read(self.expt_filepath)
        if isinstance(self.info, (list, str)):
            expt_config.function[self.node_id].success = "error"
            message = self.error()
        else:
            expt_config.function[self.node_id].success = "success"
            message = self.success()

        ConfigWriter.write(
            dirname=self.workflow_dirpath,
            filename=DIRPATH.EXPERIMENT_YML,
            config=asdict(expt_config),
        )

        return message

    def success(self):
        return Message(
            status="success",
            message=f"{self.algo_name} success",
            outputPaths=self.outputPaths(),
        )

    def error(self):
        return Message(
            status="error",
            message=self.info if isinstance(self.info, str) else "\n".join(self.info),
        )

    def outputPaths(self):
        outputPaths: Dict[str, OutputPath] = {}
        for k, v in self.info.items():
            if isinstance(v, BaseData):
                v.save_json(self.node_dirpath)
                if v.output_path:
                    outputPaths[k] = v.output_path

        return outputPaths
