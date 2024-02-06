import os
from abc import ABCMeta, abstractmethod
from dataclasses import dataclass


@dataclass
class OMEDataModel:
    """OME(Open Microscopy Environment) aware metadata
    @link https://ome-model.readthedocs.io/en/stable/ome-xml/
    """

    image_name: str
    size_x: int  # width
    size_y: int  # height
    size_t: int  # time frames
    size_z: int  # axis z frames
    size_c: int  # channels
    acquisition_date: str
    objective_model: str  # objective lens model
    fps: int  # frames_per_second  # Note: extended from OME format

    def get_ome_values(self):
        return {
            "Image": {
                "Name": self.image_name,
                "AcquisitionDate": self.acquisition_date,
                "Fps": self.fps,
            },
            "Pixels": {
                "SizeC": self.size_c,
                "SizeT": self.size_t,
                "SizeX": self.size_x,
                "SizeY": self.size_y,
                "SizeZ": self.size_z,
            },
            "Objective": {
                "Model": self.objective_model,
            },
        }


class MicroscopeDataReaderBase(metaclass=ABCMeta):
    """Microscope data reader base class"""

    LIBRARY_DIR_KEY = "MICROSCOPES_LIBRARY_DIR"

    def __init__(self):
        """
        Initialization
        """
        self._init_library()

        # init members
        self.__data_file_path = None
        self.__resource_handles = None
        self.__original_metadata = None
        self.__ome_metadata = None
        self.__lab_specific_metadata = None

    def __del__(self):
        """
        Destructor
        """
        if self.__resource_handles is not None:
            self._release_resources()
            self.__resource_handles = None

    def load(self, data_file_path: str):
        """
        Release resources
        """
        if self.__resource_handles is not None:
            raise Exception("Reader module already initialized.")

        """
        Reset data
        """
        self.__original_metadata = None
        self.__ome_metadata = None
        self.__lab_specific_metadata = None

        """
        Load data file
        """
        handles = self._load_file(data_file_path)
        self.__resource_handles = handles
        self.__data_file_path = data_file_path
        data_name = os.path.basename(data_file_path)

        """
        Read metadata
        """
        self.__original_metadata = self._build_original_metadata(data_name)
        self.__ome_metadata = self._build_ome_metadata(self.__original_metadata)
        self.__lab_specific_metadata = self._build_lab_specific_metadata(
            self.__original_metadata
        )

    def get_image_stacks(self) -> list:
        """Return microscope image stacks"""
        return self._get_image_stacks()

    @abstractmethod
    def _init_library(self) -> dict:
        """Initialize microscope library"""
        pass

    @abstractmethod
    def _load_file(self, data_file_path: str) -> object:
        """Return metadata specific to microscope instruments"""
        pass

    @abstractmethod
    def _build_original_metadata(self, data_name: str) -> dict:
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
    def _release_resources() -> None:
        """Release microscope library resources"""
        pass

    @abstractmethod
    def _get_image_stacks(self) -> list:
        """Return microscope image stacks"""
        pass

    @property
    def data_file_path(self) -> str:
        return self.__data_file_path

    @property
    def resource_handles(self) -> list:
        return self.__resource_handles

    @property
    def original_metadata(self) -> dict:
        return self.__original_metadata

    @property
    def ome_metadata(self) -> dict:
        return self.__ome_metadata

    @property
    def lab_specific_metadata(self) -> dict:
        return self.__lab_specific_metadata
