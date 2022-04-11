import pickle
from dataclasses import asdict
from glob import glob
from typing import Dict

from optinist.api.workflow.workflow import Message
from optinist.api.config.config_writer import ConfigWriter
from optinist.api.experiment.experiment_config_reader import ExpConfigReader
from optinist.api.workflow.workflow import OutputPath
from optinist.wrappers.data_wrapper import *

from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.dir_path import DIRPATH


def get_results(unique_id, nodeIdList):
    runPaths = []
    for node_id in nodeIdList:
        for path in glob(join_filepath([DIRPATH.BASE_DIR, unique_id, node_id, "*.pkl"])):
            runPaths.append(path.replace("\\", "/"))

    print(runPaths)

    results: Dict[str, Message] = {}
    for path in runPaths:
        node_id = path.split("/")[-2]
        algo_name = path.split("/")[-1].split(".")[0]

        with open(path, "rb") as f:
            info = pickle.load(f)

        results[node_id] = {}
        if isinstance(info, (list, str)):
            results[node_id] = get_error(info, node_id, unique_id)
        else:
            json_dir = "/".join(path.split("/")[:-1])
            results[node_id] = get_success(info, node_id, algo_name, json_dir, unique_id)

    return results


def get_error(info, node_id, unique_id):
    config = ExpConfigReader.read(join_filepath([
        DIRPATH.BASE_DIR, unique_id, DIRPATH.EXPERIMENT_YML]))
    config.function[node_id].success = "error"

    ConfigWriter.write(
        dirname=join_filepath([DIRPATH.BASE_DIR, unique_id]),
        filename=DIRPATH.EXPERIMENT_YML,
        config=asdict(config),
    )

    return Message(
        status="error",
        message=info if isinstance(info, str) else "\n".join(info),
    )


def get_success(info, node_id, algo_name, json_dir, unique_id):
    config = ExpConfigReader.read(join_filepath(
        [DIRPATH.BASE_DIR, unique_id, DIRPATH.EXPERIMENT_YML]))
    config.function[node_id].success = "success"

    ConfigWriter.write(
        dirname=join_filepath([DIRPATH.BASE_DIR, unique_id]),
        filename=DIRPATH.EXPERIMENT_YML,
        config=asdict(config),
    )

    return Message(
        status="success",
        message=f"{algo_name} success",
        outputPaths=get_outputPaths(info, json_dir)
    )


def get_outputPaths(info, json_dir):
    outputPaths: Dict[str, OutputPath] = {}
    for k, v in info.items():
        if isinstance(v, BaseData):
            v.save_json(json_dir)

        if isinstance(v, ImageData):
            outputPaths[k] = OutputPath(
                path=v.json_path,
                type="images",
                max_index=len(v.data) if v.data.ndim == 3 else 1
            )
        elif isinstance(v, TimeSeriesData):
            outputPaths[k] = OutputPath(
                path=v.json_path,
                type="timeseries",
                max_index=len(v.data)
            )
        elif isinstance(v, CorrelationData):
            outputPaths[k] = OutputPath(
                path=v.json_path,
                type="heatmap",
            )
        elif isinstance(v, RoiData):
            outputPaths[k] = OutputPath(
                path=v.json_path,
                type="roi",
            )
        elif isinstance(v, ScatterData):
            outputPaths[k] = OutputPath(
                path=v.json_path,
                type="scatter",
            )
        elif isinstance(v, BarData):
            outputPaths[k] = OutputPath(
                path=v.json_path,
                type="bar",
            )
        elif isinstance(v, HTMLData):
            outputPaths[k] = OutputPath(
                path=v.json_path,
                type="html",
            )
        else:
            pass

    return outputPaths
