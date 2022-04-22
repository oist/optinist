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
            for path in glob(join_filepath([DIRPATH.OUTPUT_DIR, unique_id, node_id, "*.pkl"])):
                runPaths.append(path.replace("\\", "/"))

        print(runPaths)

        results: Dict[str, Message] = {}
        for path in runPaths:
            node_id = path.split("/")[-2]
            algo_name = path.split("/")[-1].split(".")[0]

            info = PickleReader.read(path)

            results[node_id] = {}
            if isinstance(info, (list, str)):
                results[node_id] = _error(info, node_id, unique_id)
            else:
                results[node_id] = _success(
                    info,
                    node_id,
                    algo_name,
                    unique_id,
                    join_filepath(path.split("/")[:-1])
                )

        return results


def _success(info, node_id, algo_name, unique_id, dirpath):
    config = ExptConfigReader.read(join_filepath(
        [DIRPATH.OUTPUT_DIR, unique_id, DIRPATH.EXPERIMENT_YML]))
    config.function[node_id].success = "success"

    ConfigWriter.write(
        dirname=join_filepath([DIRPATH.OUTPUT_DIR, unique_id]),
        filename=DIRPATH.EXPERIMENT_YML,
        config=asdict(config),
    )

    return Message(
        status="success",
        message=f"{algo_name} success",
        outputPaths=_outputPaths(info, dirpath)
    )


def _error(info, node_id, unique_id):
    config = ExptConfigReader.read(join_filepath([
        DIRPATH.OUTPUT_DIR, unique_id, DIRPATH.EXPERIMENT_YML]))
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


def _outputPaths(info, dirpath):
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
