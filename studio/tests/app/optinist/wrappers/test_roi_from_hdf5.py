import numpy as np
import pytest

from studio.app.common.dataclass import ImageData
from studio.app.optinist.dataclass import FluoData, IscellData
from studio.app.optinist.wrappers.optinist.visualize_utils.roi_fluo_from_hdf5 import (
    roi_fluo_from_hdf5,
)
from studio.app.optinist.wrappers.optinist.visualize_utils.roi_from_hdf5 import (
    roi_from_hdf5,
)

LOADERS = {
    "roi_from_hdf5": lambda img, out, iscell: roi_from_hdf5(img, out, iscell=iscell),
    "roi_fluo_from_hdf5": lambda img, out, iscell: roi_fluo_from_hdf5(
        img, FluoData(np.ones((2, 5))), out, iscell=iscell, params={"transpose": False}
    ),
}


@pytest.mark.parametrize("loader", LOADERS)
@pytest.mark.parametrize("label", [True, False], ids=["all_cells", "no_cells"])
def test_an_empty_cell_or_non_cell_set_gives_an_all_nan_image(tmp_path, loader, label):
    stack = np.zeros((2, 6, 6), np.uint8)
    stack[0, 1, 1], stack[1, 4, 4] = 1, 2
    cell_img = ImageData(stack, output_dir=str(tmp_path), file_name="cell_img")

    info = LOADERS[loader](cell_img, str(tmp_path), IscellData(np.full(2, label)))

    filled, empty = ("cell_roi", "non_cell_roi")[:: 1 if label else -1]
    assert np.isnan(info[empty].data).all()
    np.testing.assert_array_equal(info[filled].data, stack.max(axis=0))
