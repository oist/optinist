from dataclasses import dataclass
import os
from glob import glob
from typing import List, Tuple

import numpy as np

from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.utils.config_handler import ConfigReader
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.pickle_handler import PickleReader, PickleWriter
from studio.app.common.dataclass.base import BaseData
from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.nwb.nwb_creater import overwrite_nwb
from studio.app.optinist.dataclass import IscellData, RoiData, FluoData, EditRoiData
from studio.app.optinist.schemas.roi import RoiPos


@dataclass
class CellType:
    ROI = 1
    NON_ROI = 0
    TEMP_ROI = -1


class EditROI:
    def __init__(self, file_path):
        self.node_dirpath = os.path.dirname(file_path)

        self.output_info = self.__read_output_info()
        self.data: EditRoiData = self.output_info.get("edit_roi_data")
        self.iscell = self.output_info.get("iscell").data
        self.fluorescence = self.output_info.get("fluorescence").data

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

    def add(self, roi_pos):
        new_roi = self.create_ellipse_mask(self.shape, roi_pos)
        new_roi = new_roi[np.newaxis, :, :] * self.num_cell

        self.iscell = np.append(self.iscell, CellType.TEMP_ROI)
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

        self.data.im = np.vstack((self.data.im, merged_roi[np.newaxis, :, :]))

        self.iscell[ids] = CellType.NON_ROI
        self.iscell = np.append(self.iscell, CellType.TEMP_ROI)

        self.data.temp_merge_roi.append(float(self.num_cell))
        self.data.temp_merge_roi += ids
        self.data.temp_merge_roi.append((-1.0))

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
        self.iscell[ids] = CellType.NON_ROI

        info = {
            "cell_roi": RoiData(
                np.nanmax(self.data.im[self.iscell != CellType.NON_ROI], axis=0),
                output_dir=self.node_dirpath,
                file_name="cell_roi",
            ),
            "iscell": IscellData(self.iscell),
            "edit_roi_data": self.data,
        }
        self.data.temp_delete_roi += ids

        self.__update_pickle_for_roi_edition(self.pickle_file_path, info)
        self.__save_json(info)

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

    @classmethod
    def create_ellipse_mask(cls, shape: Tuple[int, int], roi_pos: RoiPos):
        import numpy as np

        x, y, width, height = (
            round(roi_pos.posx),
            round(roi_pos.posy),
            round(roi_pos.sizex),
            round(roi_pos.sizey),
        )

        x_coords = np.arange(0, shape[0])
        y_coords = np.arange(0, shape[1])
        xx, yy = np.meshgrid(x_coords, y_coords)

        # Calculate the distance of each pixel from the center of the ellipse
        a = width / 2
        b = height / 2
        distance = ((xx - x) / a) ** 2 + ((yy - y) / b) ** 2

        # Set the pixels within the ellipse to 1 and the pixels outside to NaN
        ellipse = np.empty(shape)
        ellipse[:] = np.nan
        ellipse[distance <= 1] = 1

        return ellipse
