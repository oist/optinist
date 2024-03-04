import json
import logging
import os
from pprint import pprint

from studio.app.optinist.microscopes.ThorlabsReader import ThorlabsReader

TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_PATH = TEST_DIR_PATH + "/test_data/thorlabs/timelapse600"


def test_thorlabs_reader():
    if not ThorlabsReader.is_available():
        # Note: To output the logging contents to the console,
        #       specify the following options to pytest
        #   > pytest --log-cli-level=DEBUG
        logging.warning("ThorlabsReader is not available.")
        return

    # initialize
    data_reader = ThorlabsReader()
    data_reader.load(TEST_DATA_PATH)

    # dump attributes
    print("[original_metadata]", json.dumps(data_reader.original_metadata, indent=2))
    pprint(data_reader.ome_metadata)
    pprint(data_reader.ome_metadata.get_ome_values())
    print(
        "[lab_specific_metadata]",
        json.dumps(data_reader.lab_specific_metadata, indent=2),
    )

    # get image stacks
    channels_stacks = data_reader.get_image_stacks()
    channel_len = channels_stacks.shape[1]

    # save tiff image (multi page) test
    if channel_len > 0:
        import tifffile

        for channel_idx in range(channel_len):
            save_path = "{}/{}.out.tiff".format(
                TEST_DIR_PATH, os.path.basename(TEST_DATA_PATH)
            )
            print(f"save image: {save_path}")

            tifffile.imwrite(save_path, channels_stacks[:, channel_idx, :, :])

    # asserts
    assert data_reader.ome_metadata.size_x > 0


if __name__ == "__main__":
    test_thorlabs_reader()
