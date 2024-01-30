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
    IDA_Result,
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
            raise Exception("GetAccessor Error: Please check the File path.")

        # Connect
        ida.Connect(hAccessor)

        # Open file
        result = ida.Open(
            hAccessor,
            data_file_path,
            IDA_OpenMode.IDA_OM_READ,
            ctypes.byref(hFile),
        )
        if result != IDA_Result.IDA_RESULT_SUCCESS:
            raise Exception("Open Error")

        # Get Group Handle
        hGroup = self.__hGroup = ctypes.c_void_p()
        specify_group = (
            0  # OIR Data has only 1 group, omp2info file may have more groups
        )
        ida.GetGroup(hAccessor, hFile, specify_group, ctypes.byref(hGroup))

        # GetArea
        hArea = self.__hArea = ctypes.c_void_p()
        specify_layer = 0  # OIR and omp2info file has only 1 layer
        specify_area = ctypes.c_int()
        ida.GetArea(hAccessor, hGroup, specify_layer, specify_area, ctypes.byref(hArea))

        # TODO: return type を要リファクタリング（他のReaderのI/Fとも要整合）
        return (hAccessor, hFile, hGroup, hArea)

    def _build_original_metadata(self, handle: object, data_name: str) -> dict:
        ida = self.__dll

        hAccessor = self.__hAccessor
        hFile = self.__hFile
        hGroup = self.__hGroup
        hArea = self.__hArea

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
            "file_creation_time": file_creation_time.get_values(),
            "system_info": system_info.get_values(),
            "user_comment": user_comment.get_values(),
        }

        return original_metadata

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        @link OME/NativeND2Reader
        """

        rect = original_metadata["rect"]
        axis_info = original_metadata["axis_info"]

        # get sequence counts
        # Note: Use the largest sequence count for each axis.
        axes_sequence_counts = []
        for axis_name, axis in axis_info.items():
            axes_sequence_counts.append(axis["max"])
        sequence_count = max(axes_sequence_counts)

        fps = 0  # TODO: 計算対象予定

        omeData = OMEDataModel(
            image_name=original_metadata["data_name"],
            size_x=rect["width"],
            size_y=rect["height"],
            size_t=sequence_count,
            size_c=0,
            fps=fps,
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

        hAccessor = self.__hAccessor
        hArea = self.__hArea

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
