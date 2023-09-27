import copy
import gc
import os
from dataclasses import dataclass
from glob import glob

from fastapi import HTTPException, status
from snakemake import snakemake

from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.utils.config_handler import ConfigReader
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.filepath_finder import find_condaenv_filepath
from studio.app.common.core.utils.pickle_handler import PickleReader, PickleWriter
from studio.app.common.dataclass.base import BaseData
from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.edit_ROI.wrappers import edit_roi_wrapper_dict
from studio.app.optinist.core.nwb.nwb_creater import overwrite_nwb


@dataclass
class ACTION:
    EDIT_ROI: str = "edit_ROI"
    ADD: str = "add"
    MERGE: str = "merge"
    DELETE: str = "delete"


class EditROI:
    def __init__(self, action: str, filepath: str, params: dict) -> None:
        self.action = action
        self.filepath = filepath
        self.params = params
        self.algo = self.get_algo()
        self.node_dirpath = os.path.dirname(self.filepath)

    def get_algo(self):
        filepath = self.filepath

        if (
            not os.path.exists(filepath)
            and os.path.commonpath([DIRPATH.DATA_DIR, filepath]) != DIRPATH.DATA_DIR
        ):
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND)
        if os.path.basename(filepath) != "cell_roi.json":
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)

        algo_list = edit_roi_wrapper_dict.keys()

        algo = next((algo for algo in algo_list if algo in filepath), None)
        if not algo:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)
        return algo

    def excute(self):
        snakemake(
            DIRPATH.SNAKEMAKE_FILEPATH,
            use_conda=True,
            cores=2,
            workdir=f"{os.path.dirname(DIRPATH.STUDIO_DIR)}",
            config={
                "type": ACTION.EDIT_ROI,
                "action": self.action,
                "node_dirpath": self.node_dirpath,
                "algo": self.algo,
                "params": self.params,
            },
        )


class EditRoiUtils:
    @classmethod
    def conda(cls, config):
        algo = config["algo"]
        if "conda_name" in edit_roi_wrapper_dict[algo]:
            conda_name = edit_roi_wrapper_dict[algo]["conda_name"]
            return find_condaenv_filepath(conda_name) if conda_name else None

        return None

    @classmethod
    def excute(cls, config):
        algo = config["algo"]
        node_dirpath = config["node_dirpath"]
        action = config["action"]
        params = config["params"]

        func = copy.deepcopy(edit_roi_wrapper_dict[algo]["function"][action])
        output_info = func(node_dirpath, **params)
        del func
        gc.collect()

        workflow_dirpath = os.path.dirname(node_dirpath)
        smk_config_file = join_filepath(
            [workflow_dirpath, DIRPATH.SNAKEMAKE_CONFIG_YML]
        )
        smk_config = ConfigReader.read(smk_config_file)
        last_outputs = smk_config.get("last_output")
        pickle_files = glob(join_filepath([node_dirpath, "*.pkl"]))

        if len(pickle_files) > 0:
            func_name = os.path.splitext(os.path.basename(pickle_files[0]))[0]
            cls.__update_pickle_for_roi_edition(pickle_files[0], func_name, output_info)

            for last_output in last_outputs:
                last_output_path = join_filepath([DIRPATH.OUTPUT_DIR, last_output])
                last_output_info = cls.__update_pickle_for_roi_edition(
                    last_output_path, func_name, output_info
                )
                rule_type = os.path.splitext(os.path.basename(last_output_path))[0]
                whole_nwb_path = join_filepath(
                    [workflow_dirpath, f"whole_{rule_type}.nwb"]
                )

                Runner.save_all_nwb(whole_nwb_path, last_output_info["nwbfile"])

                del last_output_info
                gc.collect()

        for k, v in output_info.items():
            if isinstance(v, BaseData):
                v.save_json(node_dirpath)

            if k == "nwbfile":
                nwb_files = glob(join_filepath([node_dirpath, "*.nwb"]))
                if len(nwb_files) > 0:
                    overwrite_nwb(v, node_dirpath, os.path.basename(nwb_files[0]))

        del output_info
        gc.collect()

    @classmethod
    def __update_pickle_for_roi_edition(cls, filepath, func_name, new_output_info):
        output_info = PickleReader.read(filepath)
        for k, v in new_output_info.items():
            if k == "nwbfile":
                output_info[k][func_name] = v
            else:
                output_info[k] = v
        PickleWriter.overwrite(pickle_path=filepath, info=output_info)
        return output_info
