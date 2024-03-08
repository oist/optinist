import json
import logging
import os
import sys
from pprint import pprint

# add sys.path for conda env
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../../../../")

from studio.app.optinist.microscopes.IsxdReader import IsxdReader  # NOQA

TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_PATH = (
    TEST_DIR_PATH + "/test_data/inscopix/oist_short_example_preprocessed.isxd"
)


def test_isxd_reader(dump_metadata=True, dump_stack=True):
    if not IsxdReader.is_available():
        # Note: To output the logging contents to the console,
        #       specify the following options to pytest
        #   > pytest --log-cli-level=DEBUG
        logging.warning("IsxdReader is not available.")
        return

    # initialize
    data_reader = IsxdReader()
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
        # get image stacks
        image_stack = data_reader.get_image_stacks()

        # save tiff image (multi page) test
        if len(image_stack) > 0:
            from PIL import Image

            save_stack = [Image.fromarray(frame) for frame in image_stack]
            save_path = "{}/{}.out.tiff".format(
                TEST_DIR_PATH, os.path.basename(TEST_DATA_PATH)
            )
            print(f"save image: {save_path}")

            save_stack[0].save(
                save_path,
                compression="tiff_deflate",
                save_all=True,
                append_images=save_stack[1:],
            )

    # asserts
    assert data_reader.original_metadata["spacing"]["width"] > 0
    assert data_reader.ome_metadata.size_x > 0


if __name__ == "__main__":
    test_isxd_reader(dump_stack=True)
