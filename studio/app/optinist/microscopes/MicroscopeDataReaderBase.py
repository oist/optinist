import os
from abc import ABCMeta, abstractmethod
from dataclasses import dataclass


@dataclass
class OMEDataModel:
    """OME(Open Microscopy Environment) aware metadata
    @link https://ome-model.readthedocs.io/en/stable/ome-xml/
    """

    image_name: str
    size_x: int
    size_y: int
    size_t: int
    size_c: int
    fps: int  # frames_per_second  # TODO: OME標準の類似項目に合わせる

    # TODO: 以下今後追加想定
    # SizeZ: int
    # AcquisitionDate: date
    # Instrument/(Laser|Detector): str


class MicroscopeDataReaderBase(metaclass=ABCMeta):
    """Microscope data reader base class"""

    LIBRARY_DIR_KEY = "MICROSCOPES_LIBRARY_DIR"

    def __init__(self):
        """
        Initialization
        """
        self._init_library()

        # init members
        self.__data_path = None
        self.__data_handle = None
        self.__original_metadata = None
        self.__ome_metadata = None
        self.__lab_specific_metadata = None

    def load(self, data_path: str):
        """
        Reset data
        """
        self.__original_metadata = None
        self.__ome_metadata = None
        self.__lab_specific_metadata = None

        """
        Load data
        """
        handle = self._load_data_file(data_path)
        self.__data_handle = handle
        self.__data_path = data_path
        data_name = os.path.basename(data_path)

        """
        Read metadata
        """
        self.__original_metadata = self._build_original_metadata(handle, data_name)
        self.__ome_metadata = self._build_ome_metadata(self.__original_metadata)
        self.__lab_specific_metadata = self._build_lab_specific_metadata(
            self.__original_metadata
        )

        """
        Release resources
        """
        self._release_resources(handle)

    @abstractmethod
    def _init_library(self) -> dict:
        """Initialize microscope library"""
        pass

    @abstractmethod
    def _load_data_file(self, data_path: str) -> object:
        """Return metadata specific to microscope instruments"""
        pass

    @abstractmethod
    def _build_original_metadata(self, handle: object, data_name: str) -> dict:
        """Build metadata specific to microscope instruments"""
        pass

    @abstractmethod
    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """Build OME(Open Microscopy Environment) aware metadata"""
        pass

    @abstractmethod
    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        """Build metadata in lab-specific format"""
        pass

    @abstractmethod
    def _release_resources(self, handle: object) -> None:
        """Release microscope library resources"""
        pass

    @abstractmethod
    def get_images_stack(self) -> list:
        """Return microscope image stacks"""
        pass

    @property
    def data_path(self) -> str:
        return self.__data_path

    @property
    def data_handle(self) -> object:
        return self.__data_handle

    @property
    def original_metadata(self) -> dict:
        return self.__original_metadata

    @property
    def ome_metadata(self) -> dict:
        return self.__ome_metadata

    @property
    def lab_specific_metadata(self) -> dict:
        return self.__lab_specific_metadata
