import json
import logging
import os
from pprint import pprint

from ND2Reader import ND2Reader

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
LIBRARY_DIR = CURRENT_DIR_PATH + "/dll"
TEST_DATA_PATH = CURRENT_DIR_PATH + "/testdata/nikon/pia_volume_area1.nd2"

os.environ[ND2Reader.LIBRARY_DIR_KEY] = LIBRARY_DIR


def test_nd2_reader():
    if not ND2Reader.is_available():
        # Note: To output the logging contents to the console,
        #       specify the following options to pytest
        #   > pytest --log-cli-level=DEBUG
        logging.warning("ND2Reader is not available.")
        return

    # initialize
    data_reader = ND2Reader()
    data_reader.load(TEST_DATA_PATH)

    # dump attributes
    print("[original_metadata]", json.dumps(data_reader.original_metadata, indent=2))
    pprint(data_reader.ome_metadata)
    pprint(data_reader.ome_metadata.get_ome_values())
    print(
        "[lab_specific_metadata]",
        json.dumps(data_reader.lab_specific_metadata, indent=2),
    )

    # get image stacks (for all channels)
    channels_stacks = data_reader.get_image_stacks()

    # save tiff image (multi page) test
    if (len(channels_stacks) > 0) and (len(channels_stacks[0]) > 0):
        from PIL import Image

        # save stacks for all channels
        for channel_idx, image_stack in enumerate(channels_stacks):
            save_stack = [Image.fromarray(frame) for frame in image_stack]
            save_path = (
                os.path.basename(TEST_DATA_PATH) + f".out.ch{channel_idx+1}.tiff"
            )
            print(f"save image: {save_path}")

            save_stack[0].save(
                save_path,
                compression="tiff_deflate",
                save_all=True,
                append_images=save_stack[1:],
            )

    # asserts
    assert data_reader.original_metadata["attributes"]["widthPx"] > 0
    assert data_reader.ome_metadata.size_x > 0
    assert data_reader.lab_specific_metadata["uiWidth"] > 0


if __name__ == "__main__":
    test_nd2_reader()
