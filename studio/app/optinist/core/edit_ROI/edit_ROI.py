import os
from dataclasses import dataclass
from glob import glob
from typing import Dict, List

import numpy as np
from fastapi import HTTPException, status

from studio.app.common.core.experiment.experiment import ExptOutputPathIds
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.snakemake.snakemake_reader import SmkConfigReader
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.pickle_handler import PickleReader, PickleWriter
from studio.app.common.dataclass.base import BaseData
from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.edit_ROI.utils import create_ellipse_mask
from studio.app.optinist.core.nwb.nwb_creater import overwrite_nwb
from studio.app.optinist.dataclass import EditRoiData, IscellData, RoiData
from studio.app.optinist.schemas.roi import RoiStatus

logger = AppLogger.get_logger()


@dataclass
class CellType:
    ROI = 1
    NON_ROI = 0
    TEMP_ADD = -1
    TEMP_DELETE = -2
    TEMP_PROMOTE = -3


class EditROI:
    def __init__(self, file_path):
        self.node_dirpath = os.path.dirname(file_path)
        self.workflow_dirpath = os.path.dirname(self.node_dirpath)
        self.workflow_ids = ExptOutputPathIds(self.node_dirpath)
        self.function_id = self.workflow_ids.function_id

        self.output_info: Dict = PickleReader.read(self.pickle_file_path)
        self.tmp_output_info: Dict = (
            PickleReader.read(self.tmp_pickle_file_path)
            if os.path.exists(self.tmp_pickle_file_path)
            else {}
        )

        self.data = self.output_info.get("edit_roi_data", {})
        self.tmp_data: EditRoiData = self.tmp_output_info.get(
            "edit_roi_data", self.data
        )

        if not isinstance(self.tmp_data, EditRoiData):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)

        self.tmp_data.images = None

        self.tmp_iscell = self.tmp_output_info.get(
            "iscell", self.output_info.get("iscell")
        ).data

        logger.info("start edit roi: %s", self.function_id)

    @property
    def pickle_file_path(self):
        files = list(
            set(glob(join_filepath([self.node_dirpath, "*.pkl"])))
            - set(glob(join_filepath([self.node_dirpath, "tmp_*.pkl"])))
        )
        if len(files) == 0:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)
        return files[0]

    @property
    def tmp_pickle_file_path(self):
        return join_filepath([self.node_dirpath, f"tmp_{self.function_id[:-11]}.pkl"])

    @property
    def shape(self):
        return self.tmp_data.im.shape[1:]

    @property
    def num_cell(self):
        return self.tmp_data.im.shape[0]

    @property
    def __cell_roi_im(self):
        # pending promotions and undone pending merges are not drawn as cells
        drawn = ~np.isin(self.tmp_iscell, [CellType.NON_ROI, CellType.TEMP_PROMOTE])
        for merged_id in self.tmp_data.temp_merge_roi:
            if self.tmp_iscell[int(merged_id)] == CellType.TEMP_DELETE:
                drawn[int(merged_id)] = False
        return np.fmax.reduce(
            self.tmp_data.im, axis=0, where=drawn[:, None, None], initial=np.nan
        )

    def get_status(self) -> RoiStatus:
        roi_status = self.tmp_data.status()
        roi_status.temp_promote_roi = np.where(
            self.tmp_iscell == CellType.TEMP_PROMOTE
        )[0].tolist()
        return roi_status

    def add(self, roi_pos):
        new_roi = create_ellipse_mask(self.shape, roi_pos)
        new_roi = new_roi[np.newaxis, :, :] * self.num_cell

        self.tmp_data.temp_add_roi[self.num_cell] = roi_pos
        self.tmp_iscell = np.append(self.tmp_iscell, CellType.TEMP_ADD)
        self.tmp_data.im = np.vstack((self.tmp_data.im, new_roi))

        info = {
            "cell_roi": RoiData(
                self.__cell_roi_im,
                output_dir=self.node_dirpath,
                file_name="cell_roi",
            ),
            "iscell": IscellData(self.tmp_iscell),
            "edit_roi_data": self.tmp_data,
        }
        self.__update_pickle_for_roi_edition(self.tmp_pickle_file_path, info)
        self.__save_json(info)

    def merge(self, ids: List[int]):
        merging_rois = self.tmp_data.im[ids, :, :]
        merging_rois[np.isnan(merging_rois)] = -np.inf
        merged_roi = np.maximum.reduce(merging_rois)
        merged_roi = np.where(merged_roi == -np.inf, np.nan, self.num_cell)

        self.tmp_data.temp_merge_roi[float(self.num_cell)] = ids
        self.tmp_data.im = np.vstack((self.tmp_data.im, merged_roi[np.newaxis, :, :]))

        self.tmp_iscell[ids] = CellType.TEMP_DELETE
        self.tmp_iscell = np.append(self.tmp_iscell, CellType.TEMP_ADD)

        info = {
            "cell_roi": RoiData(
                self.__cell_roi_im,
                output_dir=self.node_dirpath,
                file_name="cell_roi",
            ),
            "iscell": IscellData(self.tmp_iscell),
            "edit_roi_data": self.tmp_data,
        }

        self.__update_pickle_for_roi_edition(self.tmp_pickle_file_path, info)
        self.__save_json(info)

    def delete(self, ids: List[int]):
        # deleting a pending merge restores its sources; the entry stays so commit
        # still appends the merged ROI's trace and keeps fluorescence aligned
        for id in ids:
            for parent in self.tmp_data.temp_merge_roi.get(float(id), []):
                if parent in self.tmp_data.temp_delete_roi:
                    continue
                self.tmp_iscell[parent] = (
                    CellType.TEMP_ADD
                    if parent in self.tmp_data.temp_add_roi
                    else CellType.ROI
                )

        self.tmp_iscell[ids] = CellType.TEMP_DELETE

        for id in ids:
            self.tmp_data.temp_delete_roi[id] = None

        info = {
            "cell_roi": RoiData(
                self.__cell_roi_im,
                output_dir=self.node_dirpath,
                file_name="cell_roi",
            ),
            "iscell": IscellData(self.tmp_iscell),
            "edit_roi_data": self.tmp_data,
        }

        self.__update_pickle_for_roi_edition(self.tmp_pickle_file_path, info)
        self.__save_json(info)

    def promote(self, ids: List[int]):
        num_roi = len(self.tmp_iscell)
        not_promotable = [
            id
            for id in ids
            if not 0 <= id < num_roi or self.tmp_iscell[id] != CellType.NON_ROI
        ]
        if not_promotable:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"ROIs are not non-cell ROIs: {not_promotable}",
            )

        # a cell with no fluorescence record would have nothing to plot
        num_trace = len(self.output_info.get("fluorescence").data)
        without_trace = [id for id in ids if id >= num_trace]
        if without_trace:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"ROIs have no fluorescence record: {without_trace}",
            )

        self.tmp_iscell[ids] = CellType.TEMP_PROMOTE

        info = {
            "iscell": IscellData(self.tmp_iscell),
            "edit_roi_data": self.tmp_data,
        }

        self.__update_pickle_for_roi_edition(self.tmp_pickle_file_path, info)
        self.__save_json(info)

    def commit(self):
        self.tmp_iscell[self.tmp_iscell == CellType.TEMP_PROMOTE] = CellType.ROI

        if "suite2p" in self.function_id:
            from studio.app.optinist.core.edit_ROI.wrappers.suite2p_edit_roi import (
                commit_edit as suite2p_commit,
            )

            info = suite2p_commit(
                self.tmp_data,
                self.output_info["ops"],
                self.tmp_iscell,
                self.node_dirpath,
                self.function_id,
            )
        elif "lccd" in self.function_id:
            from studio.app.optinist.core.edit_ROI.wrappers.lccd_edit_roi import (
                commit_edit as lccd_commit,
            )

            info = lccd_commit(
                self.data.images,
                self.tmp_data,
                self.output_info.get("fluorescence"),
                self.tmp_iscell,
                self.node_dirpath,
                self.function_id,
            )

        elif "vacant_roi" in self.function_id:
            from studio.app.optinist.core.edit_ROI.wrappers.vacant_roi_edit_roi import (
                commit_edit as vacant_roi_commit,
            )

            info = vacant_roi_commit(
                self.data.images,
                self.tmp_data,
                self.output_info.get("fluorescence"),
                self.tmp_iscell,
                self.node_dirpath,
                self.function_id,
            )

        elif "caiman" in self.function_id:
            from studio.app.optinist.core.edit_ROI.wrappers.caiman_edit_roi import (
                commit_edit as caiman_commit,
            )

            info = caiman_commit(
                self.data.images,
                self.tmp_data,
                self.output_info.get("fluorescence"),
                self.tmp_iscell,
                self.node_dirpath,
                self.function_id,
            )

        iscell = info["iscell"].data
        non_cell_roi_file_name = self.__non_cell_roi_file_name()
        if non_cell_roi_file_name:
            im = info["edit_roi_data"].im
            # rows without a fluorescence record would offer a click that fails
            has_trace = np.arange(len(im)) < len(info["fluorescence"].data)
            non_cell = (iscell == CellType.NON_ROI) & has_trace
            # lowest id on top, so a deleted merge does not cover its sources
            info["non_cell_roi"] = RoiData(
                np.fmin.reduce(
                    im, axis=0, where=non_cell[:, None, None], initial=np.nan
                ),
                output_dir=self.node_dirpath,
                file_name=non_cell_roi_file_name,
            )

        info["edit_roi_data"].images = self.data.images

        self.__update_pickle_for_roi_edition(self.pickle_file_path, info)
        self.__save_json(info)
        self.__update_whole_nwb(info)

        (
            os.remove(self.tmp_pickle_file_path)
            if os.path.exists(self.tmp_pickle_file_path)
            else None
        )

    def cancel(self):
        original_num_cell = len(self.output_info.get("fluorescence").data)
        self.tmp_data.im = self.tmp_data.im[:original_num_cell]
        self.tmp_iscell = self.tmp_iscell[:original_num_cell]
        self.tmp_iscell[self.tmp_iscell == CellType.TEMP_PROMOTE] = CellType.NON_ROI
        self.tmp_data.cancel()

        info = {
            "cell_roi": RoiData(
                self.__cell_roi_im,
                output_dir=self.node_dirpath,
                file_name="cell_roi",
            ),
        }
        self.__save_json(info)
        (
            os.remove(self.tmp_pickle_file_path)
            if os.path.exists(self.tmp_pickle_file_path)
            else None
        )

    def __non_cell_roi_file_name(self):
        for file_name in ("non_cell_roi", "noncell_roi"):
            if os.path.exists(join_filepath([self.node_dirpath, f"{file_name}.json"])):
                return file_name
        return None

    def __update_whole_nwb(self, output_info):
        smk_config = SmkConfigReader.read(
            self.workflow_ids.workspace_id, self.workflow_ids.unique_id
        )
        last_outputs = smk_config.get("last_output")

        for last_output in last_outputs:
            last_output_path = join_filepath([DIRPATH.OUTPUT_DIR, last_output])
            last_output_info = self.__update_pickle_for_roi_edition(
                last_output_path, output_info
            )
            whole_nwb_path = join_filepath([self.workflow_dirpath, "whole.nwb"])

            Runner.save_all_nwb(whole_nwb_path, last_output_info["nwbfile"])

    def __save_json(self, output_info):
        for k, v in output_info.items():
            if isinstance(v, BaseData):
                v.save_json(self.node_dirpath)

            if k == "nwbfile":
                nwb_files = glob(join_filepath([self.node_dirpath, "[!tmp_]*.nwb"]))

                if len(nwb_files) > 0:
                    overwrite_nwb(v, self.node_dirpath, os.path.basename(nwb_files[0]))

    def __update_pickle_for_roi_edition(self, file_path, new_output_info):
        func_name = os.path.splitext(os.path.basename(self.pickle_file_path))[0]
        for k, v in new_output_info.items():
            if k == "nwbfile":
                self.output_info[k][func_name] = v
            else:
                self.output_info[k] = v
        PickleWriter.write(pickle_path=file_path, info=self.output_info)
        return self.output_info
