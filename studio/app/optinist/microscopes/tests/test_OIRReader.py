import json
import logging
import os
import sys
from pprint import pprint

# add sys.path for conda env
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../../../../")

from studio.app.optinist.microscopes.OIRReader import OIRReader  # NOQA

TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_PATH = TEST_DIR_PATH + "/test_data/olympus/olympus-xyt005_0001.oir"


def test_oir_reader(dump_metadata=True, dump_stack=True):
    if not OIRReader.is_available():
        # Note: To output the logging contents to the console,
        #       specify the following options to pytest
        #   > pytest --log-cli-level=DEBUG
        logging.warning("OIRReader is not available.")
        return

    # initialize
    data_reader = OIRReader()
    data_reader.load(TEST_DATA_PATH)

    # dump attributes
    if dump_metadata:
        print(
            "[original_metadata]", json.dumps(data_reader.original_metadata, indent=2)
        )
        pprint(data_reader.ome_metadata)
        pprint(data_reader.ome_metadata.get_ome_values())
        print(
            "[lab_specific_metadata]",
            json.dumps(data_reader.lab_specific_metadata, indent=2),
        )

    # get & dump image stack
    if dump_stack:
        # get image stacks (for all channels)
        channels_stacks = data_reader.get_image_stacks()

        # save tiff image (multi page) test
        if (channels_stacks.shape[0] > 0) and (channels_stacks.shape[1] > 0):
            import tifffile

            for channel_idx, image_stack in enumerate(channels_stacks):
                save_path = "{}/{}.out.ch{}.tiff".format(
                    TEST_DIR_PATH, os.path.basename(TEST_DATA_PATH), channel_idx + 1
                )
                print(f"save image: {save_path}")

                tifffile.imwrite(save_path, image_stack)

    # asserts
    assert data_reader.original_metadata["rect"]["width"] > 0
    assert data_reader.ome_metadata.size_x > 0
    assert data_reader.lab_specific_metadata["uiWidth"] > 0


if __name__ == "__main__":
    test_oir_reader(dump_stack=True)
