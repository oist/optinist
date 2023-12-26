import logging
import os
from pprint import pprint

from OIRReader import OIRReader

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
LIBRARY_DIR = CURRENT_DIR_PATH + "/dll"
TEST_DATA_PATH = CURRENT_DIR_PATH + "/testdata/olympus/olympus-xyt005_0001.oir"

os.environ[OIRReader.LIBRARY_DIR_KEY] = LIBRARY_DIR


def test_oir_reader():
    if not OIRReader.is_available():
        # Note: To output the logging contents to the console,
        #       specify the following options to pytest
        #   > pytest --log-cli-level=DEBUG
        logging.warning("OIRReader is not available.")
        return

    # initialize
    data_reader = OIRReader()

    pprint(data_reader)

    # TODO: Under construction


if __name__ == "__main__":
    test_oir_reader()
