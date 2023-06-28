import os

from optinist.api.dir_path import DIRPATH
from optinist.api.nwb.nwb_creater import save_nwb

input_config = {
    "device": {
        "description": "Microscope Information",
        "manufacturer": "Microscope Manufacture",
        "name": "Microscope device",
    },
    "experiment_description": "None",
    "identifier": "optinist",
    "image_series": {
        "starting_frame": [0],
        "starting_time": 0,
    },
    "imaging_plane": {
        "description": "standard",
        "excitation_lambda": 600.0,
        "imaging_rate": 30.0,
        "indicator": "GCaMap",
        "location": "V1",
        "name": "ImagingPlane",
    },
    "ophys": {"plane_segmentation": {"description": "", "name": "PlaneSegmentation"}},
    "optical_channel": {
        "description": "optical channel",
        "emission_lambda": 500.0,
        "name": "OpticalChannel",
    },
    "session_description": "optinist",
}
config = {}


def test_save_nwb():
    save_nwb(f"{DIRPATH.OUTPUT_DIR}/test.nwb", input_config, config)
    assert os.path.exists(f"{DIRPATH.OUTPUT_DIR}/test.nwb")
