from dataclasses import dataclass
from enum import Enum


@dataclass
class FILETYPE:
    IMAGE: str = "image"
    CSV: str = "csv"
    HDF5: str = "hdf5"
    BEHAVIOR: str = "behavior"
    MATLAB: str = "matlab"
    MICROSCOPE: str = "microscope"


class ACCEPT_FILE_EXT(Enum):
    TIFF_EXT = [".tif", ".tiff", ".TIF", ".TIFF"]
    CSV_EXT = [".csv"]
    HDF5_EXT = [".hdf5", ".h5", ".nwb", ".HDF5", ".NWB"]
    MATLAB_EXT = [".mat"]
    MICROSCOPE_EXT = [".nd2", ".oir", ".isxd", ".thor.zip"]

    ALL_EXT = TIFF_EXT + CSV_EXT + HDF5_EXT + MATLAB_EXT + MICROSCOPE_EXT


ORIGINAL_DATA_EXT = ".orig"

NOT_DISPLAY_ARGS_LIST = ["params", "output_dir", "nwbfile", "kwargs"]

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
