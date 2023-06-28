import os

import numpy as np

from optinist.api.dataclass.dataclass import (
    CsvData,
    HTMLData,
    ImageData,
    IscellData,
    RoiData,
    TimeSeriesData,
)
from optinist.api.dir_path import DIRPATH


def test_image_data():
    data = ImageData(np.random.rand(10, 10), file_name="test_image_data")

    assert os.path.exists(data.path[0])


def test_roi_data():
    data = RoiData(np.random.rand(10, 10), file_name="rois")

    assert os.path.exists(data.path)


def test_timeseries_data():
    data = TimeSeriesData(np.random.rand(10, 10), file_name="timeseries")
    data.save_json(DIRPATH.OUTPUT_DIR)
    assert os.path.exists(data.json_path)


def test_iscell_data():
    IscellData(np.random.rand(10))


def test_html_data():
    HTMLData(np.random.rand(10))


def test_csv_data():
    data = CsvData(data=np.random.rand(10, 10), params={})
    data.save_json(DIRPATH.OUTPUT_DIR)
    assert os.path.exists(data.json_path)
