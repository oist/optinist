import logging
import os

from ND2Reader import ND2Reader

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
LIBRARY_DIR = CURRENT_DIR_PATH + "/dll"
TEST_DATA_PATH = CURRENT_DIR_PATH + "/testdata/nikon/nikon-oriOD001.nd2"

os.environ[ND2Reader.LIBRARY_DIR_KEY] = LIBRARY_DIR


def test_nd2_reader():
    if not ND2Reader.is_available():
        # Note: To output the logging contents to the console,
        #       specify the following options to pytest
        #   > pytest.exe --log-cli-level=DEBUG
        logging.warning("ND2Reader is not available.")
        return

    data_reader = ND2Reader()
    data_reader.load(TEST_DATA_PATH)

    # # debug print.
    import json

    print(json.dumps(data_reader.original_metadata))
    # print(data_reader.ome_metadata)
    # print(json.dumps(data_reader.lab_specific_metadata))

    assert data_reader.original_metadata["attributes"]["widthPx"] > 0
    assert data_reader.ome_metadata.size_x > 0
    assert data_reader.lab_specific_metadata["uiWidth"] > 0


if __name__ == "__main__":
    test_nd2_reader()
