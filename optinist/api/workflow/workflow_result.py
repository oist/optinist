import os
from dataclasses import asdict
from glob import glob
from typing import Dict

from optinist.api.pickle.pickle_reader import PickleReader
from optinist.api.dataclass.dataclass import *
from optinist.api.workflow.workflow import Message, OutputPath, OutputType
from optinist.api.config.config_writer import ConfigWriter
from optinist.api.experiment.experiment_reader import ExptConfigReader
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.dir_path import DIRPATH


class WorkflowResult:
    @classmethod
    def get(cls, unique_id, nodeIdList):
        runPaths = []
        for node_id in nodeIdList:
            pickle_filepath = join_filepath([
                DIRPATH.OUTPUT_DIR,
                unique_id,
                node_id,
                "*.pkl"
            ])
            for path in glob(pickle_filepath):
                runPaths.append(path.replace("\\", "/"))

        results: Dict[str, Message] = {}
        for path in runPaths:
            node_id = path.split("/")[-2]
            algo_name = path.split("/")[-1].split(".")[0]

            info = PickleReader.read(path)

            ## success or error
            results[node_id] = {}
            if isinstance(info, (list, str)):
                results[node_id] = cls.error(info, node_id, unique_id)
            else:
                results[node_id] = cls.success(
                    info,
                    node_id,
                    algo_name,
                    unique_id,
                    join_filepath(path.split("/")[:-1])
                )

            # has nwb output or not
            cls.has_nwb(unique_id, node_id)

        cls.has_nwb(unique_id)

        return results

    @classmethod
    def has_nwb(cls, unique_id, node_id=None):
        if node_id is None:
            nwb_filepath_list = glob(join_filepath([
                DIRPATH.OUTPUT_DIR,
                unique_id,
                "*.nwb"
            ]))
        else:
            nwb_filepath_list = glob(join_filepath([
                DIRPATH.OUTPUT_DIR,
                unique_id,
                node_id,
                "*.nwb"
            ]))

        for nwb_filepath in nwb_filepath_list:
            if os.path.exists(nwb_filepath):
                expt_filepath = join_filepath([
                    DIRPATH.OUTPUT_DIR,
                    unique_id,
                    DIRPATH.EXPERIMENT_YML
                ])
                config = ExptConfigReader.read(expt_filepath)

                if node_id is None:
                    config.hasNWB = True
                else:
                    config.function[node_id].hasNWB = True

                ConfigWriter.write(
                    dirname=join_filepath([DIRPATH.OUTPUT_DIR, unique_id]),
                    filename=DIRPATH.EXPERIMENT_YML,
                    config=asdict(config),
                )

    @classmethod
    def success(cls, info, node_id, algo_name, unique_id, dirpath):
        expt_filepath = join_filepath([
            DIRPATH.OUTPUT_DIR,
            unique_id,
            DIRPATH.EXPERIMENT_YML
        ])
        config = ExptConfigReader.read(expt_filepath)
        config.function[node_id].success = "success"

        ConfigWriter.write(
            dirname=join_filepath([DIRPATH.OUTPUT_DIR, unique_id]),
            filename=DIRPATH.EXPERIMENT_YML,
            config=asdict(config),
        )

        return Message(
            status="success",
            message=f"{algo_name} success",
            outputPaths=cls.outputPaths(info, dirpath)
        )


    @classmethod
    def error(cls, info, node_id, unique_id):
        expt_filepath = join_filepath([
            DIRPATH.OUTPUT_DIR,
            unique_id,
            DIRPATH.EXPERIMENT_YML
        ])
        config = ExptConfigReader.read(expt_filepath)
        config.function[node_id].success = "error"

        ConfigWriter.write(
            dirname=join_filepath([DIRPATH.OUTPUT_DIR, unique_id]),
            filename=DIRPATH.EXPERIMENT_YML,
            config=asdict(config),
        )

        return Message(
            status="error",
            message=info if isinstance(info, str) else "\n".join(info),
        )

    @classmethod
    def outputPaths(cls, info, dirpath):
        outputPaths: Dict[str, OutputPath] = {}
        for k, v in info.items():
            if isinstance(v, BaseData):
                v.save_json(dirpath)

            if isinstance(v, ImageData):
                outputPaths[k] = OutputPath(
                    path=v.json_path,
                    type=OutputType.IMAGE,
                    max_index=len(v.data) if v.data.ndim == 3 else 1
                )
            elif isinstance(v, TimeSeriesData):
                outputPaths[k] = OutputPath(
                    path=v.json_path,
                    type=OutputType.TIMESERIES,
                    max_index=len(v.data)
                )
            elif isinstance(v, HeatMapData):
                outputPaths[k] = OutputPath(
                    path=v.json_path,
                    type=OutputType.HEATMAP,
                )
            elif isinstance(v, RoiData):
                outputPaths[k] = OutputPath(
                    path=v.json_path,
                    type=OutputType.ROI,
                )
            elif isinstance(v, ScatterData):
                outputPaths[k] = OutputPath(
                    path=v.json_path,
                    type=OutputType.SCATTER,
                )
            elif isinstance(v, BarData):
                outputPaths[k] = OutputPath(
                    path=v.json_path,
                    type=OutputType.BAR,
                )
            elif isinstance(v, HTMLData):
                outputPaths[k] = OutputPath(
                    path=v.json_path,
                    type=OutputType.HTML,
                )
            else:
                pass

        return outputPaths
