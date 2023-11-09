import os
from dataclasses import dataclass
from glob import glob
from typing import List

import numpy as np
from fastapi import HTTPException, status
from snakemake import snakemake

from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.utils.config_handler import ConfigReader
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.filepath_finder import find_condaenv_filepath
from studio.app.common.core.utils.pickle_handler import PickleReader, PickleWriter
from studio.app.common.dataclass.base import BaseData
from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.edit_ROI.utils import create_ellipse_mask
from studio.app.optinist.core.edit_ROI.wrappers import edit_roi_wrapper_dict
from studio.app.optinist.core.nwb.nwb_creater import overwrite_nwb
from studio.app.optinist.dataclass import *
from studio.app.optinist.schemas.roi import RoiStatus


@dataclass
class CellType:
    ROI = 1
    NON_ROI = 0
    TEMP_ADD = -1
    TEMP_DELETE = -2


class EditRoiUtils:
    @classmethod
    def conda(cls, config):
        algo = config["algo"]
        if "conda_name" in edit_roi_wrapper_dict[algo]:
            conda_name = edit_roi_wrapper_dict[algo]["conda_name"]
            return find_condaenv_filepath(conda_name) if conda_name else None

        return None

    @classmethod
    def get_algo(cls, filepath):
        algo_list = edit_roi_wrapper_dict.keys()

        algo = next((algo for algo in algo_list if algo in filepath), None)
        if not algo:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)
        return algo

    @classmethod
    def execute(cls, filepath):
        snakemake(
            DIRPATH.SNAKEMAKE_FILEPATH,
            use_conda=True,
            cores=2,
            workdir=f"{os.path.dirname(DIRPATH.STUDIO_DIR)}",
            config={
                "type": "EDIT_ROI",
                "algo": cls.get_algo(filepath),
                "file_path": filepath,
            },
        )


class EditROI:
    def __init__(self, file_path):
        self.node_dirpath = os.path.dirname(file_path)

        self.output_info = self.__read_output_info()
        self.data: EditRoiData = self.output_info.get("edit_roi_data")
        self.iscell = self.output_info.get("iscell").data
        self.fluorescence = self.output_info.get("fluorescence")

        print("start edit roi:", self.function_id)

    @property
    def function_id(self):
        return self.node_dirpath.split("/")[-1]

    @property
    def pickle_file_path(self):
        return glob(join_filepath([self.node_dirpath, "*.pkl"]))[0]

    @property
    def shape(self):
        return self.data.im.shape[1:]

    @property
    def num_cell(self):
        return self.data.im.shape[0]

    def get_status(self) -> RoiStatus:
        return self.data.status()

    def add(self, roi_pos):
        new_roi = create_ellipse_mask(self.shape, roi_pos)
        new_roi = new_roi[np.newaxis, :, :] * self.num_cell

        self.data.temp_add_roi[self.num_cell] = roi_pos
        self.iscell = np.append(self.iscell, CellType.TEMP_ADD)
        self.data.im = np.vstack((self.data.im, new_roi))

        info = {
            "cell_roi": RoiData(
                np.nanmax(self.data.im[self.iscell != CellType.NON_ROI], axis=0),
                output_dir=self.node_dirpath,
                file_name="cell_roi",
            ),
            "iscell": IscellData(self.iscell),
            "edit_roi_data": self.data,
        }
        self.__update_pickle_for_roi_edition(self.pickle_file_path, info)
        self.__save_json(info)

    def merge(self, ids: List[int]):
        merging_rois = self.data.im[ids, :, :]
        merging_rois[np.isnan(merging_rois)] = -np.inf
        merged_roi = np.maximum.reduce(merging_rois)
        merged_roi = np.where(merged_roi == -np.inf, np.nan, self.num_cell)

        self.data.temp_merge_roi[float(self.num_cell)] = ids
        self.data.im = np.vstack((self.data.im, merged_roi[np.newaxis, :, :]))

        self.iscell[ids] = CellType.TEMP_DELETE
        self.iscell = np.append(self.iscell, CellType.TEMP_ADD)

        info = {
            "cell_roi": RoiData(
                np.nanmax(self.data.im[self.iscell != CellType.NON_ROI], axis=0),
                output_dir=self.node_dirpath,
                file_name="cell_roi",
            ),
            "iscell": IscellData(self.iscell),
            "edit_roi_data": self.data,
        }

        self.__update_pickle_for_roi_edition(self.pickle_file_path, info)
        self.__save_json(info)

    def delele(self, ids: List[int]):
        self.iscell[ids] = CellType.TEMP_DELETE

        for id in ids:
            self.data.temp_delete_roi[id] = None

        info = {
            "iscell": IscellData(self.iscell),
            "edit_roi_data": self.data,
        }

        self.__update_pickle_for_roi_edition(self.pickle_file_path, info)
        self.__save_json(info)

    def commit(self):
        if "suite2p" in self.function_id:
            from studio.app.optinist.core.edit_ROI.wrappers.suite2p_edit_roi import (
                commit_edit as suite2p_commit,
            )

            info = suite2p_commit(
                self.data,
                self.output_info["ops"],
                self.iscell,
                self.node_dirpath,
                self.function_id,
            )
        elif "lccd" in self.function_id:
            from studio.app.optinist.core.edit_ROI.wrappers.lccd_edit_roi import (
                commit_edit as lccd_commit,
            )

            info = lccd_commit(
                self.data,
                self.fluorescence,
                self.iscell,
                self.node_dirpath,
                self.function_id,
            )

        elif "caiman" in self.function_id:
            from studio.app.optinist.core.edit_ROI.wrappers.caiman_edit_roi import (
                commit_edit as caiman_commit,
            )

            info = caiman_commit(
                self.data,
                self.fluorescence,
                self.iscell,
                self.node_dirpath,
                self.function_id,
            )

        self.__update_pickle_for_roi_edition(self.pickle_file_path, info)
        self.__save_json(info)
        self.__update_whole_nwb(info)

    def cancel(self):
        original_num_cell = len(self.fluorescence.data)
        self.data.im = self.data.im[:original_num_cell]
        self.iscell = self.iscell[:original_num_cell]
        self.data.cancel()

        info = {
            "cell_roi": RoiData(
                np.nanmax(self.data.im[self.iscell != CellType.NON_ROI], axis=0),
                output_dir=self.node_dirpath,
                file_name="cell_roi",
            ),
            "iscell": IscellData(self.iscell),
            "edit_roi_data": self.data,
        }
        self.__save_json(info)
        self.__update_pickle_for_roi_edition(self.pickle_file_path, info)

    def __update_whole_nwb(self, output_info):
        workflow_dirpath = os.path.dirname(self.node_dirpath)
        smk_config_file = join_filepath(
            [workflow_dirpath, DIRPATH.SNAKEMAKE_CONFIG_YML]
        )
        smk_config = ConfigReader.read(smk_config_file)
        last_outputs = smk_config.get("last_output")

        for last_output in last_outputs:
            last_output_path = join_filepath([DIRPATH.OUTPUT_DIR, last_output])
            last_output_info = self.__update_pickle_for_roi_edition(
                last_output_path, output_info
            )
            whole_nwb_path = join_filepath([workflow_dirpath, "whole.nwb"])

            Runner.save_all_nwb(whole_nwb_path, last_output_info["nwbfile"])

    def __save_json(self, output_info):
        for k, v in output_info.items():
            if isinstance(v, BaseData):
                v.save_json(self.node_dirpath)

            if k == "nwbfile":
                nwb_files = glob(join_filepath([self.node_dirpath, "[!tmp_]*.nwb"]))

                if len(nwb_files) > 0:
                    overwrite_nwb(v, self.node_dirpath, os.path.basename(nwb_files[0]))

    def __read_output_info(self) -> dict:
        return PickleReader.read(self.pickle_file_path)

    def __update_pickle_for_roi_edition(self, file_path, new_output_info):
        func_name = os.path.splitext(os.path.basename(self.pickle_file_path))[0]
        for k, v in new_output_info.items():
            if k == "nwbfile":
                self.output_info[k][func_name] = v
            else:
                self.output_info[k] = v
        PickleWriter.write(pickle_path=file_path, info=self.output_info)
        return self.output_info
