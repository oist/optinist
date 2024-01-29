import ctypes

# import json
import os
import platform

from MicroscopeDataReaderBase import MicroscopeDataReaderBase, OMEDataModel

import studio.app.optinist.microscopes.modules.olympus.lib as lib
from studio.app.optinist.microscopes.modules.olympus.h_ida import IDA_OpenMode


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
                ctypes.cdll.LoadLibrary(dependency_path)

        # load sdk library
        self.__dll = lib.load_library(__class__.get_library_path())

        # initialize sdk library
        self.__dll.Initialize()

    def _load_data_file(self, data_file_path: str) -> object:
        ida = self.__dll
        hAccessor = self.__hAccessor = ctypes.c_void_p()
        hFile = self.__hFile = ctypes.c_void_p()

        # Get Accessor
        ida.GetAccessor(data_file_path, ctypes.byref(hAccessor))
        if not hAccessor:
            # TODO: raise exception
            print("Please check the File path")
            return

        # Connect
        ida.Connect(hAccessor)

        # Open file
        # TODO: process exception
        ida.Open(
            # result = ida.Open(
            hAccessor,
            data_file_path,
            IDA_OpenMode.IDA_OM_READ,
            ctypes.byref(hFile),
        )

        # GetNumberOfGroup
        num_of_group = ctypes.c_int()
        ida.GetNumOfGroups(hAccessor, hFile, ctypes.byref(num_of_group))

        # Get Group Handle
        hGroup = self.__hGroup = ctypes.c_void_p()
        specify_group = (
            0  # OIR Data has only 1 group, omp2info file may have more groups
        )
        ida.GetGroup(hAccessor, hFile, specify_group, ctypes.byref(hGroup))

        # GetNumberOfLevels
        num_of_layer = ctypes.c_int()
        ida.GetNumOfLevels(hAccessor, hGroup, ctypes.byref(num_of_layer))
        print("-- GetNumOfLevels - num_of_layer:", num_of_layer)

        return (hAccessor, hFile, hGroup)

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
