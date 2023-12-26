import logging
import os
from pprint import pprint

from IsxdReader import IsxdReader

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_PATH = (
    CURRENT_DIR_PATH + "/testdata/inscopix/fixed-oist-sample_data_short.isxd"
)


def test_isxd_reader():
    if not IsxdReader.is_available():
        # Note: To output the logging contents to the console,
        #       specify the following options to pytest
        #   > pytest --log-cli-level=DEBUG
        logging.warning("IsxdReader is not available.")
        return

    # initialize
    data_reader = IsxdReader()
    data_reader.load(TEST_DATA_PATH)

    # debug print.
    import json

    # dump attributes
    pprint(json.dumps(data_reader.original_metadata))
    pprint(data_reader.ome_metadata)
    pprint(json.dumps(data_reader.lab_specific_metadata))

    # dump image stack
    # images_stack = data_reader.get_images_stack()
    # pprint(len(images_stack))

    # asserts
    assert data_reader.original_metadata["spacing"]["width"] > 0
    assert data_reader.ome_metadata.size_x > 0


if __name__ == "__main__":
    test_isxd_reader()
