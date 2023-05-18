import copy
import gc
import os
from dataclasses import dataclass
from glob import glob

from fastapi import HTTPException, status
from snakemake import snakemake

from optinist.api.config.config_reader import ConfigReader
from optinist.api.config.config_writer import ConfigWriter
from optinist.api.dataclass.base import BaseData
from optinist.api.dir_path import DIRPATH
from optinist.api.nwb.nwb_creater import overwrite_nwb
from optinist.api.pickle.pickle_reader import PickleReader
from optinist.api.pickle.pickle_writer import PickleWriter
from optinist.api.rules.runner import Runner
from optinist.api.utils.filepath_creater import join_filepath


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
        from optinist.api.edit_ROI import edit_roi_wrapper_dict

        filepath = self.filepath

        if (
            not os.path.exists(filepath)
            and os.path.commonpath([DIRPATH.OPTINIST_DIR, filepath])
            != DIRPATH.OPTINIST_DIR
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
        self.set_smk_config()

        snakemake(
            DIRPATH.SNAKEMAKE_FILEPATH,
            use_conda=True,
            cores=2,
            workdir=f"{os.path.dirname(DIRPATH.ROOT_DIR)}",
        )

    def set_smk_config(self):
        config = {
            "type": ACTION.EDIT_ROI,
            "action": self.action,
            "node_dirpath": self.node_dirpath,
            "algo": self.algo,
            "params": self.params,
        }
        ConfigWriter.write(
            dirname=DIRPATH.ROOT_DIR,
            filename=DIRPATH.SNAKEMAKE_CONFIG_YML,
            config=config,
        )


class EditRoiUtils:
    @classmethod
    def conda(cls, config):
        from optinist.api.edit_ROI import edit_roi_wrapper_dict

        algo = config["algo"]
        if "conda_yaml" in edit_roi_wrapper_dict[algo]:
            conda_yaml = edit_roi_wrapper_dict[algo]["conda_yaml"]
            return f"{DIRPATH.CONDAYML_DIR}/{conda_yaml}" if conda_yaml else None

        return None

    @classmethod
    def excute(cls, config):
        from optinist.api.edit_ROI import edit_roi_wrapper_dict

        algo = config["algo"]
        node_dirpath = config["node_dirpath"]
        action = config["action"]
        params = config["params"]

        func = copy.deepcopy(edit_roi_wrapper_dict[algo]["function"][action])
        output_info = func(node_dirpath, **params)
        del func
        gc.collect()

        pickle_files = glob(join_filepath([node_dirpath, "*.pkl"]))
        if len(pickle_files) > 0:
            prev_output_info = PickleReader.read(pickle_files[0])
            func_name = os.path.splitext(os.path.basename(pickle_files[0]))[0]
            for k, v in output_info.items():
                if k == "nwbfile":
                    prev_output_info[k][func_name] = v
                else:
                    prev_output_info[k] = v
            PickleWriter.overwrite(pickle_path=pickle_files[0], info=prev_output_info)

        for k, v in output_info.items():
            if isinstance(v, BaseData):
                v.save_json(node_dirpath)

            if k == "nwbfile":
                nwb_files = glob(join_filepath([node_dirpath, "*.nwb"]))
                if len(nwb_files) > 0:
                    overwrite_nwb(v, node_dirpath, os.path.basename(nwb_files[0]))

        del prev_output_info, output_info
        gc.collect()

        workflow_dirpath = os.path.dirname(node_dirpath)
        smk_config_file = join_filepath(
            [workflow_dirpath, DIRPATH.SNAKEMAKE_CONFIG_YML]
        )
        smk_config = ConfigReader.read(smk_config_file)
        last_outputs = smk_config.get("last_output")

        for last_output in last_outputs:
            last_output_path = join_filepath([DIRPATH.OUTPUT_DIR, last_output])
            last_output_info = PickleReader.read(last_output_path)

            rule_type = os.path.splitext(os.path.basename(last_output_path))[0]
            whole_nwb_path = join_filepath([workflow_dirpath, f"whole_{rule_type}.nwb"])

            Runner.save_all_nwb(whole_nwb_path, last_output_info["nwbfile"])
            del last_output_info
            gc.collect()
