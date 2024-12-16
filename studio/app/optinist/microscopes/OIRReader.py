import ctypes
import os
import platform
import shutil

import numpy as np
import requests

import studio.app.optinist.microscopes.modules.olympus.lib as lib
from studio.app.dir_path import DIRPATH
from studio.app.optinist.microscopes.MicroscopeDataReaderBase import (
    MicroscopeDataReaderBase,
    OMEDataModel,
)
from studio.app.optinist.microscopes.modules.olympus.area_image_size import (
    AreaImageSize,
)
from studio.app.optinist.microscopes.modules.olympus.axis_info import AxisInfo
from studio.app.optinist.microscopes.modules.olympus.channel_info import ChannelInfo
from studio.app.optinist.microscopes.modules.olympus.file_creation_time import (
    FileCreationTime,
)
from studio.app.optinist.microscopes.modules.olympus.frame_manager import FrameManager
from studio.app.optinist.microscopes.modules.olympus.h_ida import (
    CMN_RECT,
    IDA_AXIS_INFO,
    IDA_OpenMode,
    IDA_Result,
)
from studio.app.optinist.microscopes.modules.olympus.objective_lens_info import (
    ObjectiveLensInfo,
)
from studio.app.optinist.microscopes.modules.olympus.pixel_length import PixelLength
from studio.app.optinist.microscopes.modules.olympus.roi_collection import RoiCollection
from studio.app.optinist.microscopes.modules.olympus.scanner_settings import (
    ScannerSettings,
)
from studio.app.optinist.microscopes.modules.olympus.system_info import SystemInfo
from studio.app.optinist.microscopes.modules.olympus.user_comment import UserComment


class OIRReader(MicroscopeDataReaderBase):
    """Olympus OIR data reader

    * IDAL SDK usage method is based on IDA_Sample/IDA_Sample.cpp
    """

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
    def unpack_libs():
        """Unpack library files"""
        if not os.path.isdir(DIRPATH.MICROSCOPE_LIB_DIR):
            if not os.path.exists(DIRPATH.MICROSCOPE_LIB_ZIP):
                response = requests.get(
                    "https://github.com/oist/optinist/raw/v1.2.1/studio/app/optinist/microscopes/libs.zip"  # noqa: E501
                )
                with open(DIRPATH.MICROSCOPE_LIB_ZIP, "wb") as f:
                    f.write(response.content)

            shutil.unpack_archive(
                DIRPATH.MICROSCOPE_LIB_ZIP, DIRPATH.MICROSCOPE_LIB_DIR
            )

    @staticmethod
    def get_library_path() -> str:
        """Returns the path of the library (dll) file"""
        platform_name = platform.system()

        if not os.path.isdir(DIRPATH.MICROSCOPE_LIB_DIR):
            return None

        if platform_name not in __class__.SDK_LIBRARY_FILES:
            return None

        return (
            DIRPATH.MICROSCOPE_LIB_DIR
            + __class__.SDK_LIBRARY_FILES[platform_name]["main"]
        )

    @staticmethod
    def is_available() -> bool:
        """Determine if library is available"""
        __class__.unpack_libs()
        return os.path.isfile(__class__.get_library_path())

    def _init_library(self):
        # load sdk libraries (dependencies)
        __class__.unpack_libs()

        if "dependencies" in __class__.SDK_LIBRARY_FILES[platform.system()]:
            platform_library_dir = os.path.dirname(__class__.get_library_path())
            dependencies = __class__.SDK_LIBRARY_FILES[platform.system()][
                "dependencies"
            ]

            for dependency in dependencies:
                dependency_path = f"{platform_library_dir}/{dependency}"
                ctypes.cdll.LoadLibrary(dependency_path)

        # load sdk library
        ida = self.__dll = lib.load_library(__class__.get_library_path())

        # initialize sdk library
        ida_result = ida.Initialize()
        if ida_result != IDA_Result.IDA_RESULT_SUCCESS:
            raise Exception("IDA Initialize Error")

    def _load_file(self, data_file_path: str) -> object:
        ida = self.__dll

        # Get Accessor
        hAccessor = self.__hAccessor = ctypes.c_void_p()
        ida_result = ida.GetAccessor(data_file_path, ctypes.byref(hAccessor))
        if ida_result != IDA_Result.IDA_RESULT_SUCCESS or not hAccessor:
            raise FileNotFoundError(f"IDA GetAccessor Error: {data_file_path}")

        # Connect
        ida.Connect(hAccessor)

        # Open file
        hFile = ctypes.c_void_p()
        ida_result = ida.Open(
            hAccessor,
            data_file_path,
            IDA_OpenMode.IDA_OM_READ,
            ctypes.byref(hFile),
        )
        if ida_result != IDA_Result.IDA_RESULT_SUCCESS or not hFile:
            raise FileNotFoundError(f"IDA Open Error: {data_file_path}")

        # Get Group Handle
        hGroup = ctypes.c_void_p()
        specify_group = (
            0  # OIR Data has only 1 group, omp2info file may have more groups
        )
        ida.GetGroup(hAccessor, hFile, specify_group, ctypes.byref(hGroup))

        # GetArea
        hArea = ctypes.c_void_p()
        specify_layer = 0  # OIR and omp2info file has only 1 layer
        specify_area = ctypes.c_int()
        ida.GetArea(hAccessor, hGroup, specify_layer, specify_area, ctypes.byref(hArea))

        return (hAccessor, hFile, hGroup, hArea)

    def _build_original_metadata(self, data_name: str) -> dict:
        ida = self.__dll

        (hAccessor, hFile, hGroup, hArea) = self.resource_handles

        # --------------------------------------------------
        # Get data from API
        # --------------------------------------------------

        # GetNumberOfGroup
        num_of_groups = ctypes.c_int()
        ida.GetNumOfGroups(hAccessor, hFile, ctypes.byref(num_of_groups))

        # GetNumberOfLevels
        num_of_levels = ctypes.c_int()
        ida.GetNumOfLevels(hAccessor, hGroup, ctypes.byref(num_of_levels))

        # GetLevelImageSize
        rect = CMN_RECT()
        specify_layer = 0  # OIR and omp2info file has only 1 layer
        ida.GetLevelImageSize(hAccessor, hGroup, specify_layer, ctypes.byref(rect))

        # GetNumOfArea
        num_of_area = ctypes.c_int()
        ida.GetNumOfArea(hAccessor, hGroup, specify_layer, ctypes.byref(num_of_area))

        # Channel Information
        channel_info = ChannelInfo(hAccessor, hArea)

        # Axes Information
        axis_info = AxisInfo(hAccessor, hArea)

        # Pixel Length
        pixel_length = PixelLength(hAccessor, hArea)

        # Objective Lens Info
        objective_lens_info = ObjectiveLensInfo(hAccessor, hArea)

        # Scanner Settings
        scanner_settings = ScannerSettings(hAccessor, hArea)

        # File Creation Time
        file_creation_time = FileCreationTime(hAccessor, hArea)

        # System Information
        system_info = SystemInfo(hAccessor, hArea)

        # User Comment
        user_comment = UserComment(hAccessor, hArea)

        # --------------------------------------------------
        # Construct metadata
        # --------------------------------------------------

        original_metadata = {
            "data_name": data_name,
            "num_of_groups": num_of_groups.value,
            "num_of_levels": num_of_levels.value,
            "num_of_area": num_of_area.value,
            "rect": {
                "x": rect.x,
                "y": rect.y,
                "width": rect.width,
                "height": rect.height,
            },
            "axis_info": axis_info.get_values(),
            "pixel_length": pixel_length.get_values(),
            "channel_info": channel_info.get_values(),
            "objective_lens_info": objective_lens_info.get_values(),
            "scanner_settings": scanner_settings.get_values(),
            "file_creation_time": file_creation_time.get_values(),
            "system_info": system_info.get_values(),
            "user_comment": user_comment.get_values(),
        }

        return original_metadata

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        @link OME/OIRReader
        """

        rect = original_metadata["rect"]
        axis_info = original_metadata["axis_info"]
        pixel_length = original_metadata["pixel_length"]
        channel_info = original_metadata["channel_info"]
        objective_lens_info = original_metadata["objective_lens_info"]
        scanner_settings = original_metadata["scanner_settings"]
        file_creation_time = original_metadata["file_creation_time"]

        # get sequence counts
        # Note: Use the largest sequence count for each axis.
        axes_sequence_counts = []
        for axis_name, axis in axis_info.items():
            axes_sequence_counts.append(axis["max"])
        sequence_count = max(axes_sequence_counts)

        # get axis z frame count
        nZLoop = axis_info["ZSTACK"]["max"] if "ZSTACK" in axis_info else 1

        imaging_rate = round(1000 / scanner_settings["frame_speed"], 2)

        omeData = OMEDataModel(
            image_name=original_metadata["data_name"],
            size_x=rect["width"],
            size_y=rect["height"],
            size_t=sequence_count,
            size_z=nZLoop,
            size_c=len(channel_info),
            physical_sizex=pixel_length["x"],
            physical_sizey=pixel_length["y"],
            depth=int(channel_info[0]["depth"]) * 8,
            significant_bits=channel_info[0]["bit_count"],
            acquisition_date=file_creation_time["creation_time"],
            objective_model=objective_lens_info["name"],
            imaging_rate=imaging_rate,
        )

        return omeData

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        # --------------------------------------------------
        # Get parameters from original_metadata
        # --------------------------------------------------

        num_of_groups = original_metadata["num_of_groups"]
        num_of_levels = original_metadata["num_of_levels"]
        num_of_area = original_metadata["num_of_area"]
        rect = original_metadata["rect"]
        pixel_length = original_metadata["pixel_length"]
        channel_info = original_metadata["channel_info"]
        objective_lens_info = original_metadata["objective_lens_info"]
        file_creation_time = original_metadata["file_creation_time"]
        system_info = original_metadata["system_info"]
        user_comment = original_metadata["user_comment"]

        # --------------------------------------------------
        # Make parameters
        # --------------------------------------------------

        (hAccessor, hFile, hGroup, hArea) = self.resource_handles
        del hFile, hGroup

        axis_info = AxisInfo(hAccessor, hArea)

        nLLoop = nTLoop = nZLoop = 0
        Zstep = Zstart = Zend = Tstep = 0.0

        if axis_info.exist("LAMBDA"):
            axis = axis_info.get_axis("LAMBDA")
            nLLoop = axis.get_max()
        nLLoop = nLLoop or 1

        if axis_info.exist("ZSTACK"):
            axis = axis_info.get_axis("ZSTACK")
            nZLoop = axis.get_max()
            Zstep = axis.get_step()
            Zstart = axis.get_start()
            Zend = axis.get_end()
        nZLoop = nZLoop or 1

        if axis_info.exist("TIMELAPSE"):
            axis = axis_info.get_axis("TIMELAPSE")
            nTLoop = axis.get_max()
            Tstep = axis.get_step()
        nTLoop = nTLoop or 1

        # --------------------------------------------------
        # Construct metadata
        # --------------------------------------------------

        lab_specific_metadata = {
            "uiWidth": rect["width"],
            "uiHeight": rect["height"],
            "Loops": nTLoop,
            "ZSlicenum": nZLoop,
            "nChannel": len(channel_info),
            "PixelLengthX": pixel_length["x"],
            "PixelLengthY": pixel_length["y"],
            "ZInterval": Zstep,
            "TInterval": Tstep,
            "ZStart": Zstart,
            "ZEnd": Zend,
            "ObjectiveName": objective_lens_info["name"],
            "ObjectiveMag": objective_lens_info["magnification"],
            "ObjectiveNA": objective_lens_info["na"],
            "ReflectiveIndex": objective_lens_info["reflective_index"],
            "Immersion": objective_lens_info["immersion"],
            "Date": file_creation_time["creation_time"],
            "NumberOfGroup": num_of_groups,
            "NumberOfLevel": num_of_levels,
            "NumberOfArea": num_of_area,
            "ByteDepthCh0": channel_info[0]["depth"] if len(channel_info) > 0 else None,
            "SystemName": system_info["system_name"],
            "SystemVersion": system_info["system_version"],
            "DeviceName": system_info["device_name"],
            "UserName": system_info["user_name"],
            "CommentByUser": user_comment["comment"],
        }

        return lab_specific_metadata

    def _release_resources(self) -> None:
        # ----------------------------------------
        # Release each resource
        # ----------------------------------------

        ida = self.__dll

        (hAccessor, hFile, hGroup, hArea) = self.resource_handles

        ida.ReleaseArea(hAccessor, hArea)
        ida.ReleaseGroup(hAccessor, hGroup)
        ida.Close(hAccessor, hFile)
        ida.Disconnect(hAccessor)
        ida.ReleaseAccessor(ctypes.byref(hAccessor))
        ida.Terminate()

    def _get_image_stacks(self) -> list:
        """Return microscope image stacks"""

        (hAccessor, hFile, hGroup, hArea) = self.resource_handles

        rect = CMN_RECT()

        # Imaging ROI Information
        imaging_roi = RoiCollection(
            hAccessor, hArea, "ImagingROIList", "ImagingROIInfo"
        )

        # Channel Information
        channel_info = ChannelInfo(hAccessor, hArea)

        # Image Size
        area_image_size = AreaImageSize(hAccessor, hArea)

        # Axes Information
        axis_info = AxisInfo(hAccessor, hArea)

        pAxes = (IDA_AXIS_INFO * 3)()

        nLLoop = nTLoop = nZLoop = 0

        # For Max Loop Values for lambda, z, t
        if axis_info.exist("LAMBDA"):
            nLLoop = axis_info.get_axis("LAMBDA").get_max()
        if axis_info.exist("ZSTACK"):
            nZLoop = axis_info.get_axis("ZSTACK").get_max()
        if axis_info.exist("TIMELAPSE"):
            nTLoop = axis_info.get_axis("TIMELAPSE").get_max()

        nLLoop = nLLoop or 1
        nTLoop = nTLoop or 1
        nZLoop = nZLoop or 1

        # Retrieve all imaged area
        rect.width = area_image_size.get_x()
        rect.height = area_image_size.get_y()

        # Get the number of channels
        channels_count = channel_info.get_num_of_channel()

        # allocate return value buffer (all channel's stack)
        # *using numpy.ndarray
        result_channels_stacks = np.empty(
            [channels_count, (nLLoop * nTLoop * nZLoop), rect.height, rect.width],
            dtype=self.ome_metadata.pixel_np_dtype,
        )

        # Retrieve Image data and TimeStamp frame-by-frame
        for channel_no in range(channels_count):
            # Sequential index through nLLoop/nZLoop/nTLoop
            serial_loops_index = 0

            for i in range(nLLoop):
                for j in range(nZLoop):
                    for k in range(nTLoop):
                        nAxisCount = lib.set_frame_axis_index(
                            i, j, k, imaging_roi, axis_info, pAxes, 0
                        )

                        # Create Frame Manager
                        frame_manager = FrameManager(
                            hAccessor,
                            hArea,
                            channel_info.get_channel_id(channel_no),
                            pAxes,
                            nAxisCount,
                        )

                        # Get Image Body
                        buffer_pointer = frame_manager.get_image_body(rect)
                        ctypes_buffer_ptr = buffer_pointer[1]

                        # Obtain image data in ndarray format
                        single_plane_buffer = np.ctypeslib.as_array(ctypes_buffer_ptr)

                        # construct return value (each channel's stack)
                        result_channels_stacks[
                            channel_no, serial_loops_index
                        ] = single_plane_buffer
                        serial_loops_index += 1

                        frame_manager.release_image_body()

        # reshape/transpose operation.
        # Note: To be performed for 4D data
        # TODO: Need to test
        if self.ome_metadata.size_z > 1 and self.ome_metadata.size_t > 1:
            # 1. For OIR 4D data(XYTZ), reshape 3D format(XY(T*Z)) to 4D format(XYTZ)
            # 2. For OIR 4D data(XYTZ), transpose to 4D format(XYZT)
            raw_result_channels_stacks = result_channels_stacks
            result_channels_stacks = raw_result_channels_stacks.reshape(
                self.ome_metadata.size_c,
                self.ome_metadata.size_z,
                self.ome_metadata.size_t,
                self.ome_metadata.size_y,
                self.ome_metadata.size_x,
            ).transpose(
                0, 2, 1, 3, 4
            )  # transpose Z<->T

        return result_channels_stacks
