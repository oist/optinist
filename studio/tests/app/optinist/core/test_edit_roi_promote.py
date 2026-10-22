import json
import os
import shutil

import numpy as np
import pytest
import yaml
from fastapi import HTTPException

from studio.app.common.core.utils.pickle_handler import PickleWriter
from studio.app.dir_path import DIRPATH
from studio.app.optinist.core.edit_ROI.edit_ROI import CellType, EditROI
from studio.app.optinist.dataclass import EditRoiData, FluoData, IscellData

NUM_ROI = 3
SHAPE = (4, 4)
NUM_FRAME = 5
WORKSPACE = f"{DIRPATH.OUTPUT_DIR}/edit_roi_promote_ws"


@pytest.fixture
def node_dirpath(request):
    # under OUTPUT_DIR: commit() dispatches on the function_id read off the path
    yield f"{WORKSPACE}/{request.node.name}/vacant_roi_00000000000"
    shutil.rmtree(WORKSPACE, ignore_errors=True)


def build_node(node_dirpath, iscell):
    os.makedirs(node_dirpath)
    with open(f"{os.path.dirname(node_dirpath)}/snakemake.yaml", "w") as f:
        yaml.dump({"last_output": []}, f)

    im = np.full((NUM_ROI, *SHAPE), np.nan)
    for i in range(NUM_ROI):
        im[i, i, :] = i
    PickleWriter.write(
        pickle_path=f"{node_dirpath}/vacant_roi.pkl",
        info={
            "edit_roi_data": EditRoiData(images=np.ones((NUM_FRAME, *SHAPE)), im=im),
            "iscell": IscellData(np.array(iscell)),
            "fluorescence": FluoData(np.zeros((NUM_ROI, NUM_FRAME))),
            "nwbfile": {},
        },
    )
    rows = list(range(SHAPE[0]))
    with open(f"{node_dirpath}/non_cell_roi.json", "w") as f:
        json.dump({"columns": rows, "index": rows, "data": [[None] * 4] * 4}, f)
    return f"{node_dirpath}/cell_roi.json"


def drawn_ids(node_dirpath, name):
    with open(f"{node_dirpath}/{name}.json") as f:
        projection = np.array(json.load(f)["data"], dtype=float)
    return sorted({int(v) for v in projection[~np.isnan(projection)]})


def test_a_promoted_roi_moves_from_non_cell_roi_to_cell_roi(node_dirpath):
    file_path = build_node(
        node_dirpath, [CellType.ROI, CellType.NON_ROI, CellType.NON_ROI]
    )

    EditROI(file_path=file_path).promote([2])
    assert EditROI(file_path=file_path).get_status().temp_promote_roi == [2]
    EditROI(file_path=file_path).commit()

    assert drawn_ids(node_dirpath, "cell_roi") == [0, 2]
    assert drawn_ids(node_dirpath, "non_cell_roi") == [1]


def test_a_deleted_roi_appears_in_non_cell_roi_after_commit(node_dirpath):
    file_path = build_node(node_dirpath, [CellType.ROI, CellType.ROI, CellType.NON_ROI])

    EditROI(file_path=file_path).delete([1])
    EditROI(file_path=file_path).commit()

    assert drawn_ids(node_dirpath, "cell_roi") == [0]
    assert drawn_ids(node_dirpath, "non_cell_roi") == [1, 2]


def test_promote_refuses_an_roi_that_is_already_a_cell(node_dirpath):
    file_path = build_node(
        node_dirpath, [CellType.ROI, CellType.NON_ROI, CellType.NON_ROI]
    )

    with pytest.raises(HTTPException) as excinfo:
        EditROI(file_path=file_path).promote([0])
    assert excinfo.value.status_code == 400


def test_deleting_a_pending_merge_restores_its_sources(node_dirpath):
    file_path = build_node(node_dirpath, [CellType.ROI, CellType.ROI, CellType.NON_ROI])

    EditROI(file_path=file_path).merge([0, 1])
    EditROI(file_path=file_path).delete([NUM_ROI])

    edit_roi = EditROI(file_path=file_path)
    assert list(edit_roi.tmp_iscell[:2]) == [CellType.ROI, CellType.ROI]
    assert not {0, 1} & set(edit_roi.get_status().temp_merge_roi)
    # the undone merge no longer covers its sources, so they can be clicked again
    assert drawn_ids(node_dirpath, "cell_roi") == [0, 1]

    edit_roi.commit()
    assert drawn_ids(node_dirpath, "cell_roi") == [0, 1]


def test_undoing_a_merge_keeps_a_deleted_source_deleted(node_dirpath):
    file_path = build_node(node_dirpath, [CellType.ROI, CellType.ROI, CellType.NON_ROI])

    EditROI(file_path=file_path).delete([0])
    EditROI(file_path=file_path).merge([0, 1])
    EditROI(file_path=file_path).delete([NUM_ROI])
    EditROI(file_path=file_path).commit()

    assert drawn_ids(node_dirpath, "cell_roi") == [1]


def test_a_deleted_committed_merge_leaves_its_sources_selectable(node_dirpath):
    file_path = build_node(node_dirpath, [CellType.ROI, CellType.ROI, CellType.NON_ROI])

    EditROI(file_path=file_path).merge([0, 1])
    EditROI(file_path=file_path).commit()
    EditROI(file_path=file_path).delete([NUM_ROI])
    EditROI(file_path=file_path).commit()

    assert drawn_ids(node_dirpath, "cell_roi") == []
    assert {0, 1} <= set(drawn_ids(node_dirpath, "non_cell_roi"))


def test_cancel_drops_a_pending_promotion(node_dirpath):
    file_path = build_node(
        node_dirpath, [CellType.ROI, CellType.NON_ROI, CellType.NON_ROI]
    )

    EditROI(file_path=file_path).promote([2])
    EditROI(file_path=file_path).cancel()

    assert EditROI(file_path=file_path).get_status().temp_promote_roi == []
    assert drawn_ids(node_dirpath, "cell_roi") == [0]


def test_promote_refuses_an_roi_without_a_trace(node_dirpath):
    file_path = build_node(node_dirpath, [CellType.ROI] * NUM_ROI)

    # deleting every ROI empties the fluorescence, so none can be promoted back
    EditROI(file_path=file_path).delete(list(range(NUM_ROI)))
    EditROI(file_path=file_path).commit()
    assert drawn_ids(node_dirpath, "non_cell_roi") == []

    with pytest.raises(HTTPException) as excinfo:
        EditROI(file_path=file_path).promote([0])
    assert excinfo.value.status_code == 400
