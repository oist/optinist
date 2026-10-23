import copy
import os
import shutil
from datetime import datetime
from glob import glob
from typing import Optional

import numpy as np

from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.utils.pickle_handler import PickleReader, PickleWriter
from studio.app.common.core.workflow.workflow import DataFilterParam, WorkflowRunStatus
from studio.app.common.core.workflow.workflow_reader import WorkflowConfigReader
from studio.app.common.core.workflow.workflow_writer import WorkflowConfigWriter
from studio.app.const import ORIGINAL_DATA_EXT
from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.core.nwb.nwb_creater import overwrite_nwb
from studio.app.optinist.dataclass import FluoData, IscellData, RoiData


class WorkflowNodeDataFilter:
    def __init__(self, workspace_id: str, unique_id: str, node_id: str) -> None:
        self.workspace_id = workspace_id
        self.unique_id = unique_id
        self.node_id = node_id

        self.workflow_dirpath = join_filepath(
            [DIRPATH.OUTPUT_DIR, workspace_id, unique_id]
        )
        self.workflow_config = WorkflowConfigReader.read(workspace_id, unique_id)
        self.node_dirpath = join_filepath([self.workflow_dirpath, node_id])

        # current output data path
        self.pkl_filepath = join_filepath(
            [
                self.node_dirpath,
                self.workflow_config.nodeDict[self.node_id].data.label.split(".")[0]
                + ".pkl",
            ]
        )
        self.cell_roi_filepath = join_filepath([self.node_dirpath, "cell_roi.json"])
        self.tiff_dirpath = join_filepath([self.node_dirpath, "tiff"])
        self.fluorescence_dirpath = join_filepath([self.node_dirpath, "fluorescence"])

        # original output data path
        self.original_pkl_filepath = self.pkl_filepath + ORIGINAL_DATA_EXT
        self.original_cell_roi_filepath = self.cell_roi_filepath + ORIGINAL_DATA_EXT
        self.original_tiff_dirpath = self.tiff_dirpath + ORIGINAL_DATA_EXT
        self.original_fluorescence_dirpath = (
            self.fluorescence_dirpath + ORIGINAL_DATA_EXT
        )

    def _check_data_filter(self):
        exp_config = ExptConfigReader.read(self.workspace_id, self.unique_id)

        assert (
            exp_config.function[self.node_id].success == WorkflowRunStatus.SUCCESS.value
        )
        assert os.path.exists(self.pkl_filepath)

    def filter_node_data(self, params: Optional[DataFilterParam]):
        self._check_data_filter()

        if params and not params.is_empty:
            if not os.path.exists(self.original_pkl_filepath):
                self._backup_original_data()

            original_output_info = PickleReader.read(self.original_pkl_filepath)
            filtered_output_info = copy.deepcopy(original_output_info)
            filtered_output_info = self.filter_data(
                filtered_output_info,
                params,
                type=self.workflow_config.nodeDict[self.node_id].data.label,
                output_dir=self.node_dirpath,
            )
            PickleWriter.write(self.pkl_filepath, filtered_output_info)
            self._save_json(filtered_output_info, self.node_dirpath)
        else:
            # reset filter
            if not os.path.exists(self.original_pkl_filepath):
                return
            self._recover_original_data()

        self._write_config(params)

    def _write_config(self, params):
        node_data = self.workflow_config.nodeDict[self.node_id].data
        node_data.draftDataFilterParam = params

        WorkflowConfigWriter(
            self.workspace_id,
            self.unique_id,
            self.workflow_config.nodeDict,
            self.workflow_config.edgeDict,
        ).write()

    def _backup_original_data(self):
        logger = AppLogger.get_logger()
        logger.info(f"Backing up data to {ORIGINAL_DATA_EXT} before applying filter")
        shutil.copyfile(self.pkl_filepath, self.original_pkl_filepath)

        # Back up NWB files in node directory
        nwb_files = glob(join_filepath([self.node_dirpath, "[!tmp_]*.nwb"]))
        for nwb_file in nwb_files:
            original_nwb_file = nwb_file + ORIGINAL_DATA_EXT
            logger.info(f"Backing up NWB file: {nwb_file} → {original_nwb_file}")
            shutil.copyfile(nwb_file, original_nwb_file)

        shutil.copyfile(self.cell_roi_filepath, self.original_cell_roi_filepath)
        shutil.copytree(
            self.tiff_dirpath,
            self.original_tiff_dirpath,
            dirs_exist_ok=True,
        )
        shutil.copytree(
            self.fluorescence_dirpath,
            self.original_fluorescence_dirpath,
            dirs_exist_ok=True,
        )

    def _recover_original_data(self):
        logger = AppLogger.get_logger()
        logger.info("Recovering original data after filter removed")

        # Restore original pickle file
        os.remove(self.pkl_filepath)
        shutil.move(self.original_pkl_filepath, self.pkl_filepath)

        # Trigger snakemake re-run next node by update modification time
        os.utime(
            self.pkl_filepath,
            (os.path.getctime(self.pkl_filepath), datetime.now().timestamp()),
        )

        # Restore node NWB files
        nwb_files = glob(join_filepath([self.node_dirpath, "[!tmp_]*.nwb"]))
        for nwb_file in nwb_files:
            original_nwb_file = nwb_file + ORIGINAL_DATA_EXT
            os.remove(nwb_file)
            shutil.move(original_nwb_file, nwb_file)
            logger.info(f"Restored NWB file: {original_nwb_file} → {nwb_file}")

        # Delete whole.nwb file to force regeneration without filter
        whole_nwb_path = join_filepath([self.workflow_dirpath, "whole.nwb"])
        if os.path.exists(whole_nwb_path):
            os.remove(whole_nwb_path)

        # Read the restored data to regenerate whole.nwb
        output_info = PickleReader.read(self.pkl_filepath)
        if "nwbfile" in output_info:
            Runner.save_all_nwb(whole_nwb_path, output_info["nwbfile"])
            logger.info(f"Regenerated whole.nwb: {whole_nwb_path}")

        os.remove(self.cell_roi_filepath)
        shutil.move(self.original_cell_roi_filepath, self.cell_roi_filepath)

        shutil.rmtree(self.tiff_dirpath)
        os.rename(self.original_tiff_dirpath, self.tiff_dirpath)

        shutil.rmtree(self.fluorescence_dirpath)
        os.rename(self.original_fluorescence_dirpath, self.fluorescence_dirpath)

    def _save_json(self, output_info, node_dirpath):
        for k, v in output_info.items():
            if isinstance(v, (FluoData, RoiData)):
                v.save_json(node_dirpath)

            if k == "nwbfile":
                # Update local node NWB file
                nwb_files = glob(join_filepath([node_dirpath, "[!tmp_]*.nwb"]))
                if len(nwb_files) > 0:
                    # Extract the node-specific data from nwbfile
                    type_key = self.workflow_config.nodeDict[self.node_id].data.label
                    if type_key in v:
                        # Pass the node-specific data to overwrite_nwb
                        overwrite_nwb(
                            v[type_key], node_dirpath, os.path.basename(nwb_files[0])
                        )
                    else:
                        # If type_key not in v, use the original method
                        overwrite_nwb(v, node_dirpath, os.path.basename(nwb_files[0]))

                # Update whole.nwb at workflow level
                whole_nwb_path = join_filepath([self.workflow_dirpath, "whole.nwb"])
                Runner.save_all_nwb(whole_nwb_path, v)

    @classmethod
    def filter_data(
        cls,
        output_info: dict,
        data_filter_param: DataFilterParam,
        type: str,
        output_dir,
    ) -> dict:
        logger = AppLogger.get_logger()

        # Deep copy all mutable data to avoid in-place modification
        filtered_output_info = copy.deepcopy(output_info)
        im = filtered_output_info["edit_roi_data"].im
        fluorescence = filtered_output_info["fluorescence"].data
        dff = (
            filtered_output_info["dff"].data
            if filtered_output_info.get("dff")
            else None
        )
        iscell = filtered_output_info["iscell"].data.copy()
        nwbfile = copy.deepcopy(filtered_output_info["nwbfile"])

        # Apply filters
        if data_filter_param.dim1:
            dim1_filter_mask = data_filter_param.dim1_mask(
                max_size=fluorescence.shape[1]
            )
            fluorescence = fluorescence[:, dim1_filter_mask]
            if dff is not None:
                dff = dff[:, dim1_filter_mask]

        if data_filter_param.roi:
            roi_filter_mask = data_filter_param.roi_mask(max_size=iscell.shape[0])
            iscell_filtered = iscell.copy()
            iscell_filtered[~roi_filter_mask] = False
        else:
            iscell_filtered = iscell

        function_id = list(nwbfile[type][NWBDATASET.POSTPROCESS].keys())[0]
        filtered_function_id = f"filtered_{function_id}"

        # 1. ROI section
        if NWBDATASET.ROI in nwbfile[type]:
            nwbfile[type][NWBDATASET.ROI][filtered_function_id] = copy.deepcopy(
                nwbfile[type][NWBDATASET.ROI][function_id]
            )

        # 2. COLUMN section
        if NWBDATASET.COLUMN in nwbfile[type]:
            nwbfile[type][NWBDATASET.COLUMN][filtered_function_id] = copy.deepcopy(
                nwbfile[type][NWBDATASET.COLUMN][function_id]
            )
            nwbfile[type][NWBDATASET.COLUMN][filtered_function_id][
                "data"
            ] = iscell_filtered

        # 3. FLUORESCENCE section
        if NWBDATASET.FLUORESCENCE in nwbfile[type]:
            nwbfile[type][NWBDATASET.FLUORESCENCE][
                filtered_function_id
            ] = copy.deepcopy(nwbfile[type][NWBDATASET.FLUORESCENCE][function_id])
            nwbfile[type][NWBDATASET.FLUORESCENCE][filtered_function_id][
                "Fluorescence"
            ]["data"] = fluorescence

        # 4. POSTPROCESS section updates
        logger.info(
            f"Saving filter ROI {data_filter_param.roi}, Time {data_filter_param.dim1}"
        )
        nwbfile[type][NWBDATASET.POSTPROCESS][filtered_function_id] = {}

        # Process ROI filter indices
        if data_filter_param.roi:
            filtered_roi_indices = []
            for range_param in data_filter_param.roi:
                if range_param.end:
                    filtered_roi_indices.extend(
                        range(range_param.start, range_param.end)
                    )
                else:
                    filtered_roi_indices.append(range_param.start)
            filtered_roi_indices = np.array(filtered_roi_indices, dtype="float")
            nwbfile[type][NWBDATASET.POSTPROCESS][filtered_function_id][
                "filter_roi_ind"
            ] = filtered_roi_indices

        # Process dim1 filter indices
        if data_filter_param.dim1:
            filtered_dim1_indices = []
            for range_param in data_filter_param.dim1:
                if range_param.end:
                    filtered_dim1_indices.extend(
                        range(range_param.start, range_param.end)
                    )
                else:
                    filtered_dim1_indices.append(range_param.start)
            filtered_dim1_indices = np.array(filtered_dim1_indices, dtype="float")
            nwbfile[type][NWBDATASET.POSTPROCESS][filtered_function_id][
                "filter_time_ind"
            ] = filtered_dim1_indices

        # Build return info with the copied and modified nwbfile
        info = {
            **filtered_output_info,
            "cell_roi": RoiData(
                np.nanmax(im[iscell_filtered != 0], axis=0, initial=np.nan),
                output_dir=output_dir,
                file_name="cell_roi",
            ),
            "fluorescence": FluoData(fluorescence, file_name="fluorescence"),
            "iscell": IscellData(iscell_filtered),
            "nwbfile": nwbfile,
        }

        if dff is not None:
            info["dff"] = FluoData(dff, file_name="dff")
        else:
            info["all_roi"] = RoiData(
                np.nanmax(im, axis=0, initial=np.nan),
                output_dir=output_dir,
                file_name="all_roi",
            )
            info["non_cell_roi"] = RoiData(
                np.nanmax(im[iscell_filtered == 0], axis=0, initial=np.nan),
                output_dir=output_dir,
                file_name="non_cell_roi",
            )

        return info
