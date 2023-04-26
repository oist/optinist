import os

import pytest
from genericpath import exists

from optinist.api.config.config_writer import ConfigWriter
from optinist.api.utils.filepath_creater import join_filepath

dirpath = "/tmp/optinist/output"
filename = "test.yaml"


def test_config_writer():
    filepath = join_filepath([dirpath, filename])

    if os.path.exists(filepath):
        os.remove(filepath)

    ConfigWriter.write(dirpath, filename, {"test": "test"})

    assert os.path.exists(filepath)
