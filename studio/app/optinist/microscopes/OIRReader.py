import ctypes
import os
import platform

import numpy as np
from MicroscopeDataReaderBase import MicroscopeDataReaderBase, OMEDataModel

import studio.app.optinist.microscopes.modules.olympus.lib as lib
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
)
from studio.app.optinist.microscopes.modules.olympus.objective_lens_info import (
    ObjectiveLensInfo,
)
from studio.app.optinist.microscopes.modules.olympus.pixel_length import PixelLength
from studio.app.optinist.microscopes.modules.olympus.roi_collection import RoiCollection
from studio.app.optinist.microscopes.modules.olympus.system_info import SystemInfo
from studio.app.optinist.microscopes.modules.olympus.user_comment import UserComment


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

    def _load_data_file(self, data_file_path: str) -> object:
        ida = self.__dll

        hAccessor = self.__hAccessor = ctypes.c_void_p()
        hFile = self.__hFile = ctypes.c_void_p()

        # initialize sdk library
        # TODO: 要リファクタリング：_load_data_file の外部でのコールとする
        ida.Initialize()

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

        # Get Group Handle
        hGroup = self.__hGroup = ctypes.c_void_p()
        specify_group = (
            0  # OIR Data has only 1 group, omp2info file may have more groups
        )
        ida.GetGroup(hAccessor, hFile, specify_group, ctypes.byref(hGroup))

        # GetArea
        hArea = self.__hArea = ctypes.c_void_p()
        # TODO: `specify_layer = 0` is valid?
        specify_layer = 0  # OIR and omp2info file has only 1 layer
        specify_area = ctypes.c_int()
        ida.GetArea(hAccessor, hGroup, specify_layer, specify_area, ctypes.byref(hArea))

        # TODO: return type を要リファクタリング（他のReaderのI/Fとも要整合）
        return (hAccessor, hFile, hGroup, hArea)

    def _build_original_metadata(self, handle: object, data_name: str) -> dict:
        # TODO: 各データ取得クラス（ChannelInfo, etc）に get_metadata を定義、情報を取得・構築する
        return {}

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        @link OME/NativeND2Reader
        """

        # TODO: Under construction
        return None

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        # TODO: 以下の各値は original_metadata より取得する形式とする

        ida = self.__dll

        hAccessor = self.__hAccessor
        hFile = self.__hFile
        hGroup = self.__hGroup
        hArea = self.__hArea

        # GetNumberOfGroup
        num_of_group = ctypes.c_int()
        ida.GetNumOfGroups(hAccessor, hFile, ctypes.byref(num_of_group))

        # GetNumberOfLevels
        num_of_layer = ctypes.c_int()
        ida.GetNumOfLevels(hAccessor, hGroup, ctypes.byref(num_of_layer))

        # GetLevelImageSize
        rect = CMN_RECT()
        specify_layer = 0  # OIR and omp2info file has only 1 layer
        ida.GetLevelImageSize(hAccessor, hGroup, specify_layer, ctypes.byref(rect))

        # GetNumOfArea
        num_of_area = ctypes.c_int()
        ida.GetNumOfArea(hAccessor, hGroup, specify_layer, ctypes.byref(num_of_area))

        # Channel Information
        channel_info = ChannelInfo(hAccessor, hArea)
        channel_info.print()

        # Image Size
        area_image_size = AreaImageSize(hAccessor, hArea)
        area_image_size.print()

        # Axes Information
        axis_info = AxisInfo(hAccessor, hArea)
        axis_info.print()

        # Pixel Length
        pixel_length = PixelLength(hAccessor, hArea)
        pixel_length.print()

        # Objective Lens Info
        objective_lens_info = ObjectiveLensInfo(hAccessor, hArea)
        objective_lens_info.print()

        # File Creation Time
        file_creation_time = FileCreationTime(hAccessor, hArea)
        file_creation_time.print()

        # System Information
        system_info = SystemInfo(hAccessor, hArea)
        system_info.print()

        # User Comment
        user_comment = UserComment(hAccessor, hArea)
        user_comment.print()

        # ====================
        # Something here
        # ====================

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

        # ====================
        # Output
        # ====================

        result_data = {
            "uiWidth": rect.width,
            "uiHeight": rect.height,
            "Loops": nTLoop,
            "ZSlicenum": nZLoop,
            "nChannel": channel_info.get_num_of_channel(),
            "PixelLengthX": pixel_length.get_pixel_length_x(),
            "PixelLengthY": pixel_length.get_pixel_length_y(),
            "ZInterval": Zstep,
            "TInterval": Tstep,
            "ZStart": Zstart,
            "ZEnd": Zend,
            "ObjectiveName": objective_lens_info.name,
            "ObjectiveMag": objective_lens_info.magnification,
            "ObjectiveNA": objective_lens_info.na,
            "ReflectiveIndex": objective_lens_info.reflective_index,
            "Immersion": objective_lens_info.immersion,
            "Date": file_creation_time.creation_time,
            "NumberOfGroup": num_of_group.value,
            "NumberOfLevel": num_of_layer.value,
            "NumberOfArea": num_of_area.value,
            "ByteDepthCh0": channel_info.depth_of_ch0,
            "SystemName": system_info.system_name,
            "SystemVersion": system_info.system_version,
            "DeviceName": system_info.device_name,
            "UserName": system_info.user_name,
            "CommentByUser": user_comment.comment,
        }

        return result_data

    def _release_resources(self, handle: object) -> None:
        # ----------------------------------------
        # Release each resource
        # ----------------------------------------

        ida = self.__dll

        hAccessor = self.__hAccessor
        hFile = self.__hFile
        hGroup = self.__hGroup
        hArea = self.__hArea

        ida.ReleaseArea(hAccessor, hArea)
        ida.ReleaseGroup(hAccessor, hGroup)
        ida.Close(hAccessor, hFile)
        ida.Disconnect(hAccessor)
        ida.ReleaseAccessor(ctypes.byref(hAccessor))
        ida.Terminate()

    def get_images_stack(self) -> list:
        """Return microscope image stacks"""

        # initialization
        (hAccessor, hFile, hGroup, hArea) = self._load_data_file(self.data_path)

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

        # initialize return value (each channel's stack)
        result_channels_stacks = []

        # Retrieve Image data and TimeStamp frame-by-frame
        for channel_no in range(channel_info.get_num_of_channel()):
            # Variable for storing results (image stack)
            result_stack = []

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
                        # m_pucImageBuffer = frame_manager.get_image_body(rect)
                        frame_manager.get_image_body(rect)

                        # NOTE:
                        #   Since there are concerns about the efficiency
                        #     of this process (acquiring pixel data one dot at a time),
                        #     another process (using ndarray) is used.
                        # # Store Image Data Pixel by Pixel
                        # frame_manager.pucBuffer_to_WORD_TM()
                        # for nDataCnt in range(rect.width * rect.height):
                        #     result = frame_manager.get_pixel_value_tm(nDataCnt)
                        #     result += 1

                        # Obtain image data in ndarray format
                        pucBuffer_to_WORD_TM = frame_manager.pucBuffer_to_WORD_TM(
                            area_image_size.get_x(),
                            area_image_size.get_y(),
                        )
                        pucBuffer_ndarray = np.ctypeslib.as_array(pucBuffer_to_WORD_TM)
                        result_stack.append(pucBuffer_ndarray)

                        frame_manager.release_image_body()

                        # TODO: 要否確認
                        # frame_manager.get_frame_position()
                        # frame_manager.write_frame_position()

            # construct return value (each channel's stack)
            result_channels_stacks.append(result_stack)

        # do release resources
        # TODO: 要リファクタリング：引数仕様整理
        self._release_resources((hAccessor, hFile, hGroup, hArea))

        return result_channels_stacks
