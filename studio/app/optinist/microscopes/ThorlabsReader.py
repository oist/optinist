import os
import sys
from glob import glob

import tifffile
import xmltodict
from MicroscopeDataReaderBase import MicroscopeDataReaderBase, OMEDataModel


class ThorlabsReader(MicroscopeDataReaderBase):
    """Thorlabs(ThorImageLS) TIFF(OME-TIFF) reader"""

    # Note: Thorlabs image data format is OME-TIFF,
    #   so any module that supports it (e.g., tifffile) can be used.
    SDK_MODULE_NAME = "tifffile"

    OME_TIFF_EXT = ".tif"

    @staticmethod
    def get_library_path() -> str:
        """Returns the path of the library (dll) file"""
        # Note: Thorlabs data format does not use library (ddl).
        return None  # do nothing.

    @staticmethod
    def is_available() -> bool:
        """Determine if library is available"""
        return __class__.SDK_MODULE_NAME in sys.modules

    def _init_library(self):
        # Note: Thorlabs data format does not use library (ddl).
        pass  # do nothing.

    def _load_file(self, data_file_path: str) -> object:
        handle_tiff = None

        # If path is a directory, it is assumed to be the directory
        #   containing OME-TIFF files and processing is performed.
        if os.path.isdir(data_file_path):
            # read one tiffflie at the beginning.
            #   (in OME-TIFF, if one file is read, other tiff files
            #   are automatically read by the tiffflie module as series)
            tiff_files = sorted(
                glob(os.path.join(data_file_path, f"*{__class__.OME_TIFF_EXT}"))
            )
            if len(tiff_files) == 0:
                raise FileNotFoundError(data_file_path)

            # Note:
            # In the tiffflie module, if a split tiff file
            #   for a series of OME-TIFF is missing, a warning is output,
            #   but processing continues without interruption.
            #   -> The missing tiff file is processed as a black image.
            #   -> The caller does not seem to be able to determine
            #      whether the above warning exists or not.

            tiff_file = tiff_files[0]
            handle_tiff = tifffile.TiffFile(tiff_file)

        # If path is a file, it is assumed to be an archive file in the directory
        #   containing OME-TIFF files, and processing is performed.
        elif os.path.isfile(data_file_path):
            # TODO: 実装前
            # # 想定処理内容
            # - zipアーカイブの展開（tmpdirへ）
            # - data_file_path の更新（展開後のディレクトリパスを保持）
            # - 展開フォルダの後処理（_release_resources() での実施想定）
            # など
            raise RuntimeError("Under construction")

        else:
            raise FileNotFoundError(data_file_path)

        return (handle_tiff,)

    def _build_original_metadata(self, data_name: str) -> dict:
        handle_tiff = None
        (handle_tiff,) = self.resource_handles

        # Read extended metadata files in ThorImageLS
        metadata_experiment_path = os.path.join(self.data_file_path, "Experiment.xml")
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

        omeData = OMEDataModel(
            image_name=original_metadata["data_name"],
            size_x=int(ome_metadata["Image"]["Pixels"]["@SizeX"]),
            size_y=int(ome_metadata["Image"]["Pixels"]["@SizeY"]),
            size_t=int(ome_metadata["Image"]["Pixels"]["@SizeT"]),
            size_z=int(ome_metadata["Image"]["Pixels"]["@SizeZ"]),
            size_c=int(ome_metadata["Image"]["Pixels"]["@SizeC"]),
            acquisition_date=ome_metadata["Image"]["@AcquiredDate"],
            objective_model=None,  # TODO: 取得方法調査
            fps=0,  # TODO: 取得方法調査
        )

        return omeData

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        # Note: Not currently supported
        return None

    def _release_resources(self) -> None:
        # Note: Thorlabs data format does not use library (ddl).
        pass  # do nothing.

    def _get_image_stacks(self) -> list:
        handle_tiff = None
        (handle_tiff,) = self.resource_handles

        return handle_tiff.asarray()
