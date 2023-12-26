import ctypes

# import json
import os
import platform

from MicroscopeDataReaderBase import MicroscopeDataReaderBase, OMEDataModel


class OIRReader(MicroscopeDataReaderBase):
    """Olympus OIR data reader"""

    SDK_LIBRARY_FILES = {
        "Windows": {
            "main": "/olympus/windows/Idaldll.dll",
        },
        "Linux": {
            "main": "/olympus/linux/libIdaldll.so",
            "dependencies": (),
        },
    }

    @staticmethod
    def get_library_path() -> str:
        """Returns the path of the library (dll) file"""
        platform_name = platform.system()

        if __class__.LIBRARY_DIR_KEY not in os.environ:
            return None

        if platform_name not in __class__.SDK_LIBRARY_FILES:
            return None

        return (
            os.environ.get(__class__.LIBRARY_DIR_KEY)
            + __class__.SDK_LIBRARY_FILES[platform_name]["main"]
        )

    @staticmethod
    def is_available() -> bool:
        """Determine if library is available"""
        return (__class__.LIBRARY_DIR_KEY in os.environ) and os.path.isfile(
            __class__.get_library_path()
        )

    def _init_library(self):
        # load sdk libraries (dependencies)
        if "dependencies" in __class__.SDK_LIBRARY_FILES[platform.system()]:
            platform_library_dir = os.path.dirname(__class__.get_library_path())
            dependencies = __class__.SDK_LIBRARY_FILES[platform.system()][
                "dependencies"
            ]

            for dependency in dependencies:
                dependency_path = f"{platform_library_dir}/{dependency}"
                self.__dll = ctypes.cdll.LoadLibrary(dependency_path)

        self.__dll = ctypes.cdll.LoadLibrary(__class__.get_library_path())

        # load sdk library
        self.__dll.Initialize()

    def _load_data_file(self, data_file_path: str) -> object:
        # TODO: Under construction
        return None

    def _build_original_metadata(self, handle: object, data_name: str) -> dict:
        # TODO: Under construction
        return None

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        @link OME/NativeND2Reader
        """

        # TODO: Under construction
        return None

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        # TODO: Under construction
        return None

    def _release_resources(self, handle: object) -> None:
        # TODO: Under construction
        pass

    def get_images_stack(self) -> list:
        # TODO: under construction
        return None
