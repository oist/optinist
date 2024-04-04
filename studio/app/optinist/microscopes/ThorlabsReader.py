import os
import re
import shutil
import sys
import zipfile
from datetime import datetime
from glob import glob

import tifffile
import xmltodict

from studio.app.optinist.microscopes.MicroscopeDataReaderBase import (
    MicroscopeDataReaderBase,
    OMEDataModel,
)
from studio.app.optinist.microscopes.MicroscopeDataReaderUtils import (
    MicroscopeDataFileExt,
)


class ThorlabsReader(MicroscopeDataReaderBase):
    """Thorlabs(ThorImageLS) TIFF(OME-TIFF) reader"""

    # Note: Thorlabs image data format is OME-TIFF,
    #   so any module that supports it (e.g., tifffile) can be used.
    SDK_MODULE_NAME = "tifffile"

    RAW_ARCHIVE_FILE_PATTERN = MicroscopeDataFileExt.THOR_FILE_EXT.value

    OME_TIFF_EXT = ".tif"

    METADATA_EXPERIMENT_FILENAME = "Experiment.xml"

    @staticmethod
    def get_library_path() -> str:
        """Returns the path of the library (dll) file"""
        # Note: Thorlabs data format does not use library (ddl).
        return None  # do nothing.

    @staticmethod
    def is_available() -> bool:
        """Determine if library is available"""
        return __class__.SDK_MODULE_NAME in sys.modules

    def __init__(self):
        super().__init__()

        self.__data_extracted_paths: dict = None

    def _init_library(self):
        # Note: Thorlabs data format does not use library (ddl).
        pass  # do nothing.

    def _load_file(self, data_file_path: str) -> object:
        handle_tiff = None

        # If path is a directory, it is assumed to be the directory
        #   containing OME-TIFF files and processing is performed.
        if os.path.isdir(data_file_path):
            data_extracted_path = data_file_path
            self.__data_extracted_paths = {
                "path": data_extracted_path,
                "is_temporary": False,
            }

            handle_tiff = self.__load_ome_tiff(data_extracted_path)

        # If path is a file, it is assumed to be an archive file in the directory
        #   containing OME-TIFF files, and processing is performed.
        elif os.path.isfile(data_file_path) and re.match(
            __class__.RAW_ARCHIVE_FILE_PATTERN, data_file_path
        ):
            data_extracted_path = self.__get_data_extracted_path(data_file_path)
            self.__extract_raw_archive(data_file_path, data_extracted_path)

            handle_tiff = self.__load_ome_tiff(data_extracted_path)
            self.__data_extracted_paths = {
                "path": data_extracted_path,
                "is_temporary": True,
            }

        else:
            raise FileNotFoundError(data_file_path)

        return (handle_tiff,)

    def __get_data_extracted_path(self, data_file_path: str) -> str:
        if os.path.isdir(data_file_path):
            raise AssertionError("Unexpected case")
        elif os.path.isfile(data_file_path):
            path = "{}/{}_extracted-{}".format(
                os.path.dirname(data_file_path),
                os.path.splitext(os.path.basename(data_file_path))[0],
                datetime.now().strftime("%Y%m%d%H%M%S"),
            )
            return path
        else:
            raise FileNotFoundError(self.data_file_path)

    def __extract_raw_archive(self, data_file_path: str, data_extracted_path: str):
        if not os.path.isfile(data_file_path) or not re.match(
            __class__.RAW_ARCHIVE_FILE_PATTERN, data_file_path
        ):
            raise FileNotFoundError(data_file_path)

        with zipfile.ZipFile(data_file_path) as zf:
            # Prepare required file names in the archive
            # *Inspect file contents before extraction
            in_arvhice_subdir_name = os.path.splitext(os.path.basename(data_file_path))[
                0
            ]
            bare_experiment_filename = __class__.METADATA_EXPERIMENT_FILENAME
            subdir_experiment_filename = (
                f"{in_arvhice_subdir_name}/{bare_experiment_filename}"
            )

            # Case #1) Case without subdirectories (stored in root of archive)
            if bare_experiment_filename in zf.namelist():
                zf.extractall(data_extracted_path)

            # Case #2) Cases containing subdirectory (same dirname as archive)
            elif subdir_experiment_filename in zf.namelist():
                zf.extractall(data_extracted_path)

                # Move files under subdirectory to directly
                # under the expansion directory
                for f in glob(f"{data_extracted_path}/{in_arvhice_subdir_name}/*"):
                    shutil.move(f, data_extracted_path)

            else:
                raise AssertionError(f"Invalid raw archive file. [{data_file_path}]")

    def __cleanup_raw_extracted_path(self):
        """
        Clean up temporary directories extracted from Raw archive file.
        """

        # Note: To ensure safety,
        #   the contents of the target directory are checked and then deleted.

        data_extracted_path = self.__data_extracted_paths["path"]
        is_temporary = self.__data_extracted_paths["is_temporary"]
        metadata_experiment_path = os.path.join(
            data_extracted_path, __class__.METADATA_EXPERIMENT_FILENAME
        )

        # Perform cleanup only on temporary folder
        if not is_temporary:
            return

        if os.path.isdir(data_extracted_path) and os.path.isfile(
            metadata_experiment_path
        ):
            shutil.rmtree(data_extracted_path)
        else:
            raise FileNotFoundError(data_extracted_path)

    def __load_ome_tiff(self, data_extracted_path: str) -> object:
        # read one tiffflie at the beginning.
        #   (in OME-TIFF, if one file is read, other tiff files
        #   are automatically read by the tiffflie module as series)
        tiff_files = sorted(
            glob(os.path.join(data_extracted_path, f"*{__class__.OME_TIFF_EXT}"))
        )
        if len(tiff_files) == 0:
            raise FileNotFoundError(data_extracted_path)

        # Note:
        # In the tiffflie module, if a split tiff file
        #   for a series of OME-TIFF is missing, a warning is output,
        #   but processing continues without interruption.
        #   -> The missing tiff file is processed as a black image.
        #   -> The caller does not seem to be able to determine
        #      whether the above warning exists or not.

        tiff_file = tiff_files[0]
        handle_tiff = tifffile.TiffFile(tiff_file)

        return handle_tiff

    def _build_original_metadata(self, data_name: str) -> dict:
        handle_tiff = None
        (handle_tiff,) = self.resource_handles

        # Read extended metadata files in ThorImageLS
        metadata_experiment_path = os.path.join(
            self.__data_extracted_paths["path"], __class__.METADATA_EXPERIMENT_FILENAME
        )
        if os.path.isfile(metadata_experiment_path):
            with open(metadata_experiment_path) as f:
                experiments = xmltodict.parse(f.read())

        else:
            experiments = None

        original_metadata = {
            "data_name": data_name,
            "experiments": experiments,
        }

        return original_metadata

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        Note: Thorlabs format is not supported in OME Bio-Formats
        """

        handle_tiff = None
        (handle_tiff,) = self.resource_handles

        ome_metadata = xmltodict.parse(handle_tiff.ome_metadata)
        ome_metadata = ome_metadata["OME"]
        ome_metadata_image_pixcels = ome_metadata["Image"]["Pixels"]
        experiments = original_metadata["experiments"]["ThorImageExperiment"]

        imaging_rate = round(float(experiments["LSM"]["@frameRate"]), 2)

        omeData = OMEDataModel(
            image_name=original_metadata["data_name"],
            size_x=int(ome_metadata_image_pixcels["@SizeX"]),
            size_y=int(ome_metadata_image_pixcels["@SizeY"]),
            size_t=int(ome_metadata_image_pixcels["@SizeT"]),
            size_z=int(ome_metadata_image_pixcels["@SizeZ"]),
            size_c=int(ome_metadata_image_pixcels["@SizeC"]),
            physical_sizex=float(ome_metadata_image_pixcels["@PhysicalSizeX"]),
            physical_sizey=float(ome_metadata_image_pixcels["@PhysicalSizeY"]),
            depth=OMEDataModel.get_depth_from_pixel_type(
                ome_metadata_image_pixcels["@Type"]
            ),
            significant_bits=0,  # Note: currently unsettled
            acquisition_date=ome_metadata["Image"]["@AcquiredDate"],
            objective_model=None,  # TODO: 取得方法 要検討
            imaging_rate=imaging_rate,
        )

        return omeData

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        # Note: Not currently supported
        return None

    def _release_resources(self) -> None:
        self.__cleanup_raw_extracted_path()

    def _get_image_stacks(self) -> list:
        handle_tiff = None
        (handle_tiff,) = self.resource_handles

        # get return value buffer (all channel's stack)
        # *numpy.ndarray
        result_channels_stacks = handle_tiff.asarray()

        # reshape/transpose operation.
        raw_result_channels_stacks = result_channels_stacks
        result_channels_stacks = raw_result_channels_stacks.transpose(
            1, 0, 2, 3
        )  # transpose to XYCT -> XYTC

        return result_channels_stacks
