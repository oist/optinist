import os
import re
from enum import Enum


class MicroscopeDataFileExt(Enum):
    ND2_FILE_EXT = ".*\\.nd2$"
    OIR_FILE_EXT = ".*\\.oir$"
    ISX_FILE_EXT = ".*\\.isx$"
    THOR_FILE_EXT = ".*\\.thor.zip$"


class MicroscopeDataReaderUtils:
    """Microscope data reader utilities class"""

    @staticmethod
    def get_reader(data_file_path: str):
        """
        Automatic generation of Reader from file type information
        """
        from studio.app.optinist.microscopes.IsxdReader import IsxdReader
        from studio.app.optinist.microscopes.ND2Reader import ND2Reader
        from studio.app.optinist.microscopes.OIRReader import OIRReader
        from studio.app.optinist.microscopes.ThorlabsReader import ThorlabsReader

        if re.match(MicroscopeDataFileExt.ND2_FILE_EXT.value, data_file_path):
            assert ND2Reader.is_available(), "ND2Reader is not available."
            reader = ND2Reader()
        elif re.match(MicroscopeDataFileExt.OIR_FILE_EXT.value, data_file_path):
            assert OIRReader.is_available(), "OIRReader is not available."
            reader = OIRReader()
        elif re.match(MicroscopeDataFileExt.ISX_FILE_EXT.value, data_file_path):
            assert IsxdReader.is_available(), "IsxdReader is not available."
            reader = IsxdReader()
        elif re.match(MicroscopeDataFileExt.THOR_FILE_EXT.value, data_file_path):
            assert ThorlabsReader.is_available(), "ThorlabsReader is not available."
            reader = ThorlabsReader()
        else:
            filename = os.path.basename(data_file_path)
            raise Exception(f"Unsupported file type: {filename}")

        reader.load(data_file_path)

        return reader
