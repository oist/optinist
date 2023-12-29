import ctypes
import json
import os
import platform
from enum import IntEnum

import numpy as np
from MicroscopeDataReaderBase import MicroscopeDataReaderBase, OMEDataModel


class LimCode(IntEnum):
    LIM_OK = 0
    LIM_ERR_UNEXPECTED = -1


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
                self.__dll = ctypes.cdll.LoadLibrary(dependency_path)

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

    def _load_data_file(self, data_path: str) -> object:
        handle = self.__dll.Lim_FileOpenForReadUtf8(data_path.encode("utf-8"))

        if handle is None:
            raise FileNotFoundError(data_path)

        return handle

    def _build_original_metadata(self, handle: object, data_name: str) -> dict:
        attributes = self.__dll.Lim_FileGetAttributes(handle)
        metadata = self.__dll.Lim_FileGetMetadata(handle)
        textinfo = self.__dll.Lim_FileGetTextinfo(handle)
        experiment = self.__dll.Lim_FileGetExperiment(handle)

        attributes = json.loads(attributes)
        metadata = json.loads(metadata)
        textinfo = json.loads(textinfo) if textinfo is not None else None
        experiment = json.loads(experiment) if experiment is not None else None

        original_metadata = {
            "data_name": data_name,
            "attributes": attributes,
            "metadata": metadata,
            "textinfo": textinfo,
            "experiment": experiment,
        }

        return original_metadata

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        @link OME/NativeND2Reader
        """

        attributes = original_metadata["attributes"]
        experiment = original_metadata["experiment"]

        # TODO: experiment, periods, の参照は先頭データに固定でOK？
        if (experiment is not None) and ("periods" in experiment[0]["parameters"]):
            period_ms = float(experiment[0]["parameters"]["periods"][0]["periodMs"])
            fps = 1000 / period_ms
        else:
            period_ms = 0
            fps = 0

        omeData = OMEDataModel(
            image_name=original_metadata["data_name"],
            size_x=attributes["widthPx"],
            size_y=attributes["heightPx"],
            size_t=attributes["sequenceCount"],
            size_c=attributes["componentCount"],  # TODO: この内容が正しいか要確認
            fps=fps,
        )

        return omeData

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        # ----------------------------------------
        # Lab固有仕様のmetadata作成
        # ----------------------------------------

        # TODO: 既存のMEXコードより、以下のコード移植が必要だが、最新版ライブラリでの各値の取得方法が不明となっている。
        #     - m_Experiment.uiLevelCount
        #     - m_Experiment.pAllocatedLevels
        #     - m_Experiment.pAllocatedLevels[n].uiLoopSize
        #     - m_Experiment.pAllocatedLevels[n].dInterval
        #
        # if(m_Experiment.uiLevelCount==0){
        #       *type=(char *)"2Ds";
        # } else if (m_Experiment.uiLevelCount==2){
        #       *type=(char *)"3Dm";
        #       *mxGetPr(Loops)= (double)m_Experiment.pAllocatedLevels[0].
        #           uiLoopSize;
        #       *mxGetPr(ZSlicenum)=(double)m_Experiment.pAllocatedLevels[1].uiLoopSize;
        #       *mxGetPr(ZInterval)=(double)m_Experiment.pAllocatedLevels[1].dInterval;
        # } else {
        #     if (m_Experiment.pAllocatedLevels[0].uiExpType==2){
        #           *type=(char *)"3Ds";
        #           *mxGetPr(ZSlicenum)=(double)m_Experiment.pAllocatedLevels[0].uiLoopSize;
        #           *mxGetPr(ZInterval)=(double)m_Experiment.pAllocatedLevels[0].dInterval;
        #     }
        #     if (m_Experiment.pAllocatedLevels[0].uiExpType==0){
        #           *type=(char *)"2Dm";
        #           *mxGetPr(Loops)=(double)m_Experiment.pAllocatedLevels[0].uiLoopSize;
        #     }
        # }

        attributes = original_metadata["attributes"]
        metadata = original_metadata["metadata"]
        textinfo = original_metadata["textinfo"]
        # experiment = original_metadata["experiment"]

        # ※一部のデータ項目は ch0 より取得
        metadata_ch0_microscope = metadata["channels"][0]["microscope"]

        lab_specific_metadata = {
            # /* 11 cinoibebts (From Lab MEX) */
            "uiWidth": attributes["widthPx"],
            "uiHeight": attributes["heightPx"],
            "uiWidthBytes": attributes["widthBytes"],
            "uiColor": attributes["componentCount"],
            "uiBpcInMemory": attributes["bitsPerComponentInMemory"],
            "uiBpcSignificant": attributes["bitsPerComponentSignificant"],
            "uiSequenceCount": attributes["sequenceCount"],
            "uiTileWidth": attributes.get("tileWidthPx", None),  # Optional
            "uiTileHeight": attributes.get("tileHeightPx", None),  # Optional
            # TODO: 最新版ライブラリでは該当するパラメータがない？
            # (compressionLevel or compressionType)
            "uiCompression": None,  # TODO: Optional ?
            # TODO: 最新版ライブラリでは該当するパラメータがない？
            "uiQuality": None,
            # /* 4 components (From Lab MEX) */
            "Type": None,  # TODO: 要設定
            "Loops": None,  # TODO: 要設定
            "ZSlicenum": None,  # TODO: 要設定
            "ZInterval": None,  # TODO: 要設定
            # /* 7 components (From Lab MEX) */
            "dTimeStart": None,  # TODO: 最新版ライブラリでは該当するパラメータがない？
            "MicroPerPixel": None,  # TODO: 最新版ライブラリでは該当するパラメータがない？
            # metadata.volume.axesCalibrated が該当？
            "PixelAspect": None,  # TODO: 最新版ライブラリでは該当するパラメータがない？
            "ObjectiveName": metadata_ch0_microscope.get(
                "objectiveName", None
            ),  # Optional, channels[0] からの取得
            # TODO: ライブラリバージョンにより、取得値が異なる可能性がある？
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

    def _release_resources(self, handle: object) -> None:
        self.__dll.Lim_FileClose(handle)

    def get_images_stack(self) -> list:
        """Return microscope image stacks"""

        # initialization
        handle = self._load_data_file(self.data_path)
        pic: LIMPICTURE = LIMPICTURE()
        seq_count = self.__dll.Lim_FileGetSeqCount(handle)
        result_stack = []

        # loop for each sequence
        for seq_idx in range(seq_count):
            # read image attributes
            if LimCode.LIM_OK != self.__dll.Lim_FileGetImageData(
                handle, seq_idx, ctypes.byref(pic)
            ):
                break

            # calculate pixel byte size
            # TODO: 以下の計算式で適切であるか画像フォーマットの正確な仕様を確認必要
            #       ※一旦、pixcel_bytes=2 (uiComponents=1, pic.uiBitsPerComp) の動作OKは確認
            #       ・pixcel_bytes は最大9バイト(3*3)が想定される？
            #       ・移植元コードでは、uiComponents へのアクセスに、ビットシフトを利用している模様？
            # TODO: またその他のパラメーターの考慮も必要？
            #       ・"Channel" を考慮した画像取得（アドレス計算）が、要件への対応に必要
            #       ・3D（"ZSlicenum" ?）画像を考慮した画像取得（アドレス計算）が必要？
            pixcel_bytes = int(pic.uiComponents * ((pic.uiBitsPerComp + 7) / 8))
            if pixcel_bytes == 1:
                pixcel_data_type = ctypes.c_byte
            elif pixcel_bytes > 1 and pixcel_bytes <= 2:
                pixcel_data_type = ctypes.c_ushort
            elif pixcel_bytes > 2:
                pixcel_data_type = ctypes.c_uint
            else:
                raise AttributeError(f"Invalid pixcel_bytes: {pixcel_bytes}")

            # # debug print
            # print(
            #     f"Picture info: #{seq_idx + 1}",
            #     pic.uiWidth, pic.uiHeight, bytes, pic.uiComponents,
            #     pic.uiBitsPerComp, pixcel_bytes, pixcel_data_type,
            # )

            # scan each lines
            lines_buffer = []
            for line_idx in range(pic.uiHeight):
                # Data acquisition for line and conversion to np.ndarray format
                line_buffer = (pixcel_data_type * pic.uiWidth).from_address(
                    pic.pImageData + (line_idx * pic.uiWidthBytes)
                )
                line_buffer_array = np.ctypeslib.as_array(line_buffer)

                # stored in line data stack
                lines_buffer.append(line_buffer_array)

            # allocate planar image buffer (nb.ndarray format)
            plane_buffer = np.concatenate([[v] for v in lines_buffer])
            result_stack.append(plane_buffer)

        # release resources
        self._release_resources(handle)

        return result_stack
