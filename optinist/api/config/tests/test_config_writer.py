import pytest
import os

from optinist.api.config.config_writer import ConfigWriter
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath


def test_config_writer():
    dirpath = join_filepath([DIRPATH.ROOT_DIR, "test_data"])
    filename = "test.yaml"
    filepath = join_filepath([dirpath, filename])
    if os.path.exists(filepath):
        os.remove(filepath)

    ConfigWriter.write(
        dirpath,
        filename,
        {"a": "a"}
    )

    assert os.path.exists(filepath)