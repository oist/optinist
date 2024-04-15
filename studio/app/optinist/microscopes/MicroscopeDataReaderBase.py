import os
from abc import ABCMeta, abstractmethod
from ctypes import c_uint8, c_uint16, c_uint32
from dataclasses import dataclass

from numpy import uint8 as np_uint8
from numpy import uint16 as np_uint16
from numpy import uint32 as np_uint32
from numpy.typing import DTypeLike


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
    depth: int  # pixel depth (bits per pixel)
    significant_bits: int  # pixel significant-bits
    acquisition_date: str
    objective_model: str  # objective lens model
    imaging_rate: int  # rate images are acquired (Hz)  # *extended from OME format

    physical_sizex: float = None  # physical size of a pixel
    physical_sizex_unit: str = "µm"  # units of the physical_sizex.
    physical_sizey: float = None  # physical size of a pixel
    physical_sizey_unit: str = "µm"  # units of the physical_sizey.

    def get_ome_values(self):
        return {
            "Image": {
                "Name": self.image_name,
                "AcquisitionDate": self.acquisition_date,
                "ImagingRate": self.imaging_rate,
            },
            "Pixels": {
                "SizeC": self.size_c,
                "SizeT": self.size_t,
                "SizeX": self.size_x,
                "SizeY": self.size_y,
                "SizeZ": self.size_z,
                "PhysicalSizeX": self.physical_sizex,
                "PhysicalSizeXUnit": self.physical_sizex_unit,
                "PhysicalSizeY": self.physical_sizey,
                "PhysicalSizeYUnit": self.physical_sizey_unit,
                "SignificantBits": self.significant_bits,
                "Type": self.pixel_type,
            },
            "Objective": {
                "Model": self.objective_model,
            },
        }

    @property
    def pixel_type(self) -> str:
        """
        @see https://www.openmicroscopy.org/Schemas \
            /Documentation/Generated/OME-2016-06/ome_xsd.html#Pixels_Type
        """

        depth = self.depth
        if depth > 0 and depth <= 8:
            return "uint8"
        elif depth > 8 and depth <= 16:
            return "uint16"
        elif depth > 16 and depth <= 32:
            return "uint32"
        else:
            return "other"

    @property
    def pixel_np_dtype(self) -> DTypeLike:
        depth = self.depth
        if depth > 0 and depth <= 8:
            return np_uint8
        elif depth > 8 and depth <= 16:
            return np_uint16
        elif depth > 16 and depth <= 32:
            return np_uint32
        else:
            return np_uint8

    @property
    def pixel_ct_type(self):
        depth = self.depth
        if depth > 0 and depth <= 8:
            return c_uint8
        elif depth > 8 and depth <= 16:
            return c_uint16
        elif depth > 16 and depth <= 32:
            return c_uint32
        else:
            return c_uint8

    @staticmethod
    def get_depth_from_pixel_type(type_label: str):
        if type_label == "unit8":
            return 8
        elif type_label == "uint16":
            return 16
        elif type_label == "uint32":
            return 32
        else:
            return 0


class MicroscopeDataReaderBase(metaclass=ABCMeta):
    """Microscope data reader base class"""

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

        """
        Loading completion handler
        """
        self._load_complete_handler()

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

    # Note: not @abstractmethod
    def _load_complete_handler(self):
        """Loading completion handler"""

        # Overwrite if necessary processing.
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
