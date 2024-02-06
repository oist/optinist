import ctypes
import json
import os
import platform
import re
from enum import Enum, IntEnum

import numpy as np
from MicroscopeDataReaderBase import MicroscopeDataReaderBase, OMEDataModel


class LimCode(IntEnum):
    LIM_OK = 0
    LIM_ERR_UNEXPECTED = -1


class ExperimentLoopType(Enum):
    Unknown = "Unknown"
    TimeLoop = "TimeLoop"
    XYPosLoop = "XYPosLoop"
    ZStackLoop = "ZStackLoop"
    NETimeLoop = "NETimeLoop"


class DimensionType(Enum):
    D_2D_SINGLE = "2Ds"
    D_2D_MULTI = "2Dm"
    D_3D_SINGLE = "3Ds"
    D_3D_MULTI = "3Dm"


class LIMPICTURE(ctypes.Structure):
    _fields_ = [
        ("uiWidth", ctypes.c_uint),
        ("uiHeight", ctypes.c_uint),
        ("uiBitsPerComp", ctypes.c_uint),
        ("uiComponents", ctypes.c_uint),
        ("uiWidthBytes", ctypes.c_size_t),
        ("uiSize", ctypes.c_size_t),
        ("pImageData", ctypes.c_void_p),
    ]


class ND2Reader(MicroscopeDataReaderBase):
    """Nikon ND2 data reader"""

    SDK_LIBRARY_FILES = {
        "Windows": {
            "main": "/nikon/windows/nd2readsdk-shared.dll",
        },
        "Linux": {
            "main": "/nikon/linux/libnd2readsdk-shared.so",
            "dependencies": ("libjpeg.so.8", "libtiff.so.5"),
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
        self.__dll = ctypes.cdll.LoadLibrary(__class__.get_library_path())

        # define ctypes interfaces
        self.__dll.Lim_FileOpenForReadUtf8.argtypes = (ctypes.c_char_p,)
        self.__dll.Lim_FileOpenForReadUtf8.restype = ctypes.c_void_p
        self.__dll.Lim_FileGetMetadata.argtypes = (ctypes.c_void_p,)
        self.__dll.Lim_FileGetMetadata.restype = ctypes.c_char_p
        self.__dll.Lim_FileGetAttributes.argtypes = (ctypes.c_void_p,)
        self.__dll.Lim_FileGetAttributes.restype = ctypes.c_char_p
        self.__dll.Lim_FileGetTextinfo.argtypes = (ctypes.c_void_p,)
        self.__dll.Lim_FileGetTextinfo.restype = ctypes.c_char_p
        self.__dll.Lim_FileGetExperiment.argtypes = (ctypes.c_void_p,)
        self.__dll.Lim_FileGetExperiment.restype = ctypes.c_char_p
        self.__dll.Lim_FileGetFrameMetadata.argtypes = (ctypes.c_void_p, ctypes.c_uint)
        self.__dll.Lim_FileGetFrameMetadata.restype = ctypes.c_char_p
        self.__dll.Lim_FileGetSeqCount.argtypes = (ctypes.c_void_p,)
        self.__dll.Lim_FileGetSeqCount.restype = ctypes.c_uint
        self.__dll.Lim_FileGetCoordSize.argtypes = (ctypes.c_void_p,)
        self.__dll.Lim_FileGetCoordSize.restype = ctypes.c_size_t
        self.__dll.Lim_FileGetCoordsFromSeqIndex.argtypes = (
            ctypes.c_void_p,
            ctypes.c_uint,
        )
        self.__dll.Lim_FileGetCoordsFromSeqIndex.restype = ctypes.c_size_t
        self.__dll.Lim_FileGetImageData.argtypes = (
            ctypes.c_void_p,
            ctypes.c_uint,
            ctypes.c_void_p,
        )
        self.__dll.Lim_FileGetImageData.restype = ctypes.c_int

        self.__dll.Lim_FileClose.argtypes = (ctypes.c_void_p,)

    def _load_file(self, data_file_path: str) -> object:
        handle = self.__dll.Lim_FileOpenForReadUtf8(data_file_path.encode("utf-8"))

        if handle is None:
            raise FileNotFoundError(f"Open Error: {data_file_path}")

        return (handle,)

    def _build_original_metadata(self, data_name: str) -> dict:
        (handle,) = self.resource_handles

        attributes = self.__dll.Lim_FileGetAttributes(handle)
        metadata = self.__dll.Lim_FileGetMetadata(handle)
        textinfo = self.__dll.Lim_FileGetTextinfo(handle)
        experiments = self.__dll.Lim_FileGetExperiment(handle)
        frame_metadata = self.__dll.Lim_FileGetFrameMetadata(
            handle, 0
        )  # get 1st frame info

        attributes = json.loads(attributes)
        metadata = json.loads(metadata)
        textinfo = json.loads(textinfo) if textinfo is not None else None
        experiments = json.loads(experiments) if experiments is not None else None
        frame_metadata = json.loads(frame_metadata)

        # frame_metadata は、必要な項目のみに絞る
        frame_metadata_single = frame_metadata["channels"][0]

        original_metadata = {
            "data_name": data_name,
            "attributes": attributes,
            "metadata": metadata,
            "textinfo": textinfo,
            "experiments": experiments,
            "frame_metadata": frame_metadata_single,
        }

        return original_metadata

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        @link OME/NativeND2Reader
        """

        attributes = original_metadata["attributes"]
        metadata = original_metadata["metadata"]
        textinfo = original_metadata["textinfo"]
        experiments = original_metadata["experiments"]
        first_experiment_parameters = (
            experiments[0]["parameters"] if experiments else {}
        )
        metadata_ch0_microscope = (
            metadata["channels"][0]["microscope"] if experiments else {}
        )

        # experiment, periods, の参照は先頭データの内容から取得
        if "periods" in first_experiment_parameters:
            try:
                fps = (
                    1000
                    / first_experiment_parameters["periods"][0]["periodDiff"]["avg"]
                )
            except:  # noqa: E722
                fps = 1000 / first_experiment_parameters["periodDiff"]["avg"]
        else:
            fps = 0

        omeData = OMEDataModel(
            image_name=original_metadata["data_name"],
            size_x=attributes["widthPx"],
            size_y=attributes["heightPx"],
            size_t=0,  # size_z は後続処理で計算・設定する
            size_z=0,  # size_z は後続処理で計算・設定する
            size_c=len(metadata["channels"]),
            acquisition_date=re.sub(" +", " ", textinfo["date"]),
            objective_model=metadata_ch0_microscope.get("objectiveName", None),
            fps=fps,
        )

        return omeData

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        """Build metadata in lab-specific format"""

        attributes = original_metadata["attributes"]
        metadata = original_metadata["metadata"]
        textinfo = original_metadata["textinfo"]
        experiments = original_metadata["experiments"]
        frame_metadata = original_metadata["frame_metadata"]

        # ----------------------------------------
        # データの次元タイプ別でのパラメータ構築
        # ----------------------------------------

        experiments_count = len(experiments) if experiments is not None else 0
        dimension_type = None
        loop_count = 0
        z_slicenum = 0
        z_interval = 0

        if experiments_count == 0:
            dimension_type = DimensionType.D_2D_SINGLE.value
        elif experiments_count == 2:
            dimension_type = DimensionType.D_3D_MULTI.value
            loop_count = experiments[0]["count"]
            z_slicenum = experiments[1]["count"]
            z_interval = experiments[1]["parameters"]["stepUm"]
        else:
            if experiments[0]["type"] == ExperimentLoopType.NETimeLoop.value:
                dimension_type = DimensionType.D_2D_MULTI.value
                loop_count = experiments[0]["count"]
            elif experiments[0]["type"] == ExperimentLoopType.ZStackLoop.value:
                dimension_type = DimensionType.D_3D_SINGLE.value
                z_slicenum = experiments[0]["count"]
                z_interval = experiments[0]["parameters"]["stepUm"]

        # ※ome_metadata の一部項目をアップデート
        self.ome_metadata.size_t = loop_count
        self.ome_metadata.size_z = z_slicenum

        # ----------------------------------------
        # Lab固有metadata変数構築
        # ----------------------------------------

        # ※一部のデータ項目は ch0 より取得
        metadata_ch0_microscope = metadata["channels"][0]["microscope"]
        metadata_ch0_volume = metadata["channels"][0]["volume"]
        axes_calibration_x = metadata_ch0_volume["axesCalibration"][0]
        axes_calibration_y = metadata_ch0_volume["axesCalibration"][1]

        lab_specific_metadata = {
            # /* 11 cinoibebts (From Lab MEX) */
            "uiWidth": attributes["widthPx"],
            "uiHeight": attributes["heightPx"],
            "uiWidthBytes": attributes["widthBytes"],
            "uiColor": attributes["componentCount"],
            "uiBpcInMemory": attributes["bitsPerComponentInMemory"],
            "uiBpcSignificant": attributes["bitsPerComponentSignificant"],
            "uiSequenceCount": attributes["sequenceCount"],
            # Note: There are minor differences in values between SDK versions.
            # (uiTileWidth)
            "uiTileWidth": attributes.get("tileWidthPx", None),  # Optional
            # Note: There are minor differences in values between SDK versions.
            # (uiTileHeight)
            "uiTileHeight": attributes.get("tileHeightPx", None),  # Optional
            # Note: There are minor differences in values between SDK versions.
            # (uiCompression)
            "uiCompression": attributes.get("compressionType", None),  # Optional
            # Note: There are minor differences in values between SDK versions.
            # (uiQuality)
            "uiQuality": attributes.get("compressionLevel", None),  # Optional
            # /* 4 components (From Lab MEX) */
            "Type": dimension_type,
            "Loops": loop_count,
            "ZSlicenum": z_slicenum,
            "ZInterval": z_interval,
            # /* 7 components (From Lab MEX) */
            "dTimeStart": frame_metadata["time"]["absoluteJulianDayNumber"],
            # Note: There are minor differences in values between SDK versions.
            # (MicroPerPixel)
            "MicroPerPixel": axes_calibration_x,
            "PixelAspect": axes_calibration_x / axes_calibration_y,
            "ObjectiveName": metadata_ch0_microscope.get(
                "objectiveName", None
            ),  # Optional, channels[0] からの取得
            # Note: There are minor differences in values between SDK versions.
            # (dObjectiveMag)
            "dObjectiveMag": metadata_ch0_microscope.get(
                "objectiveMagnification", None
            ),  # Optional, channels[0] からの取得
            "dObjectiveNA": metadata_ch0_microscope.get(
                "objectiveNumericalAperture", None
            ),  # Optional, channels[0] からの取得
            "dZoom": metadata_ch0_microscope.get(
                "zoomMagnification", None
            ),  # Optional, channels[0] からの取得
            # /* 3 components (From Lab MEX) */
            # Note: All elements of textinfo are optional.
            "wszDescription": textinfo.get("description", None),
            "wszCapturing": textinfo.get("capturing", None),
            "Date": textinfo.get("date", None),
        }

        return lab_specific_metadata

    def _release_resources(self) -> None:
        (handle,) = self.resource_handles

        self.__dll.Lim_FileClose(handle)

    def _get_image_stacks(self) -> list:
        """Return microscope image stacks"""

        (handle,) = self.resource_handles

        # initialization
        pic: LIMPICTURE = LIMPICTURE()
        seq_count = self.__dll.Lim_FileGetSeqCount(handle)

        # read image attributes
        attributes = self.__dll.Lim_FileGetAttributes(handle)
        attributes = json.loads(attributes)
        image_component_count = int(attributes["componentCount"])

        # initialize return value (each channel's stack)
        result_channels_stacks = [[] for i in range(image_component_count)]

        # loop for each sequence
        for seq_idx in range(seq_count):
            # read image attributes
            if LimCode.LIM_OK != self.__dll.Lim_FileGetImageData(
                handle, seq_idx, ctypes.byref(pic)
            ):
                break

            # calculate pixel byte size
            # Note: pixcel_bytes は [1/2/4/6/8] を取りうる想定
            pixcel_bytes = int(pic.uiComponents * int((pic.uiBitsPerComp + 7) / 8))
            if pixcel_bytes not in (1, 2, 4, 6, 8):
                raise AttributeError(f"Invalid pixcel_bytes: {pixcel_bytes}")

            # # debug print
            # print(
            #     f"Picture info: #{seq_idx + 1}",
            #     pic.uiWidth, pic.uiHeight, bytes, pic.uiComponents,
            #     pic.uiBitsPerComp, pixcel_bytes,
            # )

            # scan image lines
            lines_buffer = []
            for line_idx in range(pic.uiHeight):
                # Calculate the number of bytes per line (ctypes._CData format)
                # * Probably the same value as pic.uiWidthBytes,
                #     but ported based on the logic of the SDK sample code.
                # * It appears that the unit should be ctypes.c_ushort.
                line_bytes = ctypes.c_ushort * int(pixcel_bytes * pic.uiWidth / 2)

                # Data acquisition for line and conversion to np.ndarray format
                line_buffer = line_bytes.from_address(
                    pic.pImageData + (line_idx * pic.uiWidthBytes)
                )
                line_buffer_array = np.ctypeslib.as_array(line_buffer)

                # stored in line buffer stack
                lines_buffer.append(line_buffer_array)

            # allocate planar image buffer (np.ndarray)
            raw_plane_buffer = np.concatenate([[v] for v in lines_buffer])

            # extract image data for each channel(component)
            for component_idx in range(pic.uiComponents):
                # In nikon nd2 format, "component" in effect indicates "channel".
                channel_idx = component_idx

                # A frame image is cut out from a raw frame image
                #     in units of channel(component).
                # Note: The pixel values of each component are adjacent to each other,
                #     one pixel at a time.
                #   Image: [px1: [c1][c2]..[cN]]..[pxN: [c1][c2]..[cN]]
                component_pixel_indexes = [
                    (i * image_component_count) + component_idx
                    for i in range(pic.uiWidth)
                ]
                channel_plane_buffer = raw_plane_buffer[:, component_pixel_indexes]

                # construct return value (each channel's stack)
                result_channels_stacks[channel_idx].append(channel_plane_buffer)

        return result_channels_stacks
