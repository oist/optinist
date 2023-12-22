from dataclasses import dataclass


@dataclass
class FILETYPE:
    IMAGE: str = "image"
    CSV: str = "csv"
    HDF5: str = "hdf5"
    BEHAVIOR: str = "behavior"
    MATLAB: str = "mat"


ACCEPT_TIFF_EXT = [".tif", ".tiff", ".TIF", ".TIFF"]
ACCEPT_CSV_EXT = [".csv"]
ACCEPT_HDF5_EXT = [".hdf5", ".nwb", ".HDF5", ".NWB"]
ACCEPR_MATLAB_EXT = [".mat"]

NOT_DISPLAY_ARGS_LIST = ["params", "output_dir", "nwbfile", "kwargs"]

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
