from genericpath import exists
import pytest
import os

from studio.api.config.config_writer import ConfigWriter
from studio.api.utils.filepath_creater import join_filepath

dirpath = "/tmp/studio/output"
filename = "test.yaml"


def test_config_writer():
    filepath = join_filepath([dirpath, filename])

    if os.path.exists(filepath):
        os.remove(filepath)

    ConfigWriter.write(
        dirpath,
        filename,
        {"test": "test"}
    )

    assert os.path.exists(filepath)
