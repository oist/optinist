import os

from studio.core.dir_path import DIRPATH
from studio.core.utils.config_handler import ConfigWriter
from studio.core.utils.filepath_creater import join_filepath

dirpath = f"{DIRPATH.OPTINIST_DIR}/output"
filename = "test.yaml"


def test_config_writer():
    filepath = join_filepath([dirpath, filename])

    if os.path.exists(filepath):
        os.remove(filepath)

    ConfigWriter.write(dirpath, filename, {"test": "test"})

    assert os.path.exists(filepath)
