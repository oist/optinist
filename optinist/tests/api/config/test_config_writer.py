import os

from optinist.services.config.config_writer import ConfigWriter
from optinist.services.dir_path import DIRPATH
from optinist.services.utils.filepath_creater import join_filepath

dirpath = f"{DIRPATH.OPTINIST_DIR}/output"
filename = "test.yaml"


def test_config_writer():
    filepath = join_filepath([dirpath, filename])

    if os.path.exists(filepath):
        os.remove(filepath)

    ConfigWriter.write(dirpath, filename, {"test": "test"})

    assert os.path.exists(filepath)
