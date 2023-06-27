import os
import shutil
from glob import glob

from optinist.api.dir_path import DIRPATH
from optinist.api.edit_ROI import ACTION, EditROI
from optinist.api.utils.filepath_creater import join_filepath

shutil.rmtree(f"{DIRPATH.OUTPUT_DIR}/caiman", ignore_errors=True)
shutil.copytree(
    f"{DIRPATH.OPTINIST_DIR}/output_test/caiman",
    f"{DIRPATH.OUTPUT_DIR}/caiman",
    dirs_exist_ok=True,
)
file_path = f"{DIRPATH.OUTPUT_DIR}/caiman/caiman_cnmf/cell_roi.json"
node_dirpath = os.path.dirname(file_path)

cur_num_rois = lambda: len(  # noqa: E731
    glob(join_filepath(f'{node_dirpath}/fluorescence/*.json'))
)


def test_suit2p_add_roi():
    cur_roi = cur_num_rois()
    EditROI(
        action=ACTION.ADD,
        filepath=file_path,
        params=dict(
            posx=40,
            posy=40,
            sizex=10,
            sizey=10,
        ),
    ).excute()
    # breakpoint()
    assert cur_num_rois() == cur_roi + 1


def test_suit2p_merge_roi():
    cur_roi = cur_num_rois()
    EditROI(
        action=ACTION.MERGE,
        filepath=file_path,
        params=dict(ids=[1, 2, 3]),
    ).excute()

    assert cur_num_rois() == cur_roi + 1


def test_suit2p_delete_roi():
    cur_roi = cur_num_rois()
    EditROI(
        action=ACTION.DELETE,
        filepath=file_path,
        params=dict(ids=[4, 5, 6]),
    ).excute()

    assert cur_num_rois() == cur_roi
