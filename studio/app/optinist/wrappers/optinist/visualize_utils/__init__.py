from studio.app.optinist.wrappers.optinist.visualize_utils.fluo_from_hdf5 import (
    fluo_from_hdf5,
)
from studio.app.optinist.wrappers.optinist.visualize_utils.roi_from_hdf5 import (
    roi_from_hdf5,
)

utils_wrapper_dict = {
    "fluo_from_hdf5": {
        "function": fluo_from_hdf5,
        "conda_name": "optinist",
    },
    "roi_from_hdf5": {
        "function": roi_from_hdf5,
        "conda_name": "optinist",
    },
}
