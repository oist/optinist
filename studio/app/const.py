from dataclasses import dataclass


@dataclass
class FILETYPE:
    IMAGE: str = "image"
    CSV: str = "csv"
    HDF5: str = "hdf5"
    BEHAVIOR: str = "behavior"
    MATLAB: str = "matlab"
    MICROSCOPE: str = "microscope"


ACCEPT_TIFF_EXT = [".tif", ".tiff", ".TIF", ".TIFF"]
ACCEPT_CSV_EXT = [".csv"]
ACCEPT_HDF5_EXT = [".hdf5", ".nwb", ".HDF5", ".NWB"]
ACCEPT_MATLAB_EXT = [".mat"]
ACCEPT_MICROSCOPE_EXT = [".nd2", ".oir", ".isxd", ".thor.zip"]

NOT_DISPLAY_ARGS_LIST = ["params", "output_dir", "nwbfile", "kwargs"]

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
