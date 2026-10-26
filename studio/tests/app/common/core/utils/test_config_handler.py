import os

from studio.app.common.core.utils.config_handler import ConfigWriter
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.core.workflow.workflow_params import read_default_params
from studio.app.dir_path import DIRPATH

dirpath = DIRPATH.OUTPUT_DIR
filename = "test.yaml"


def test_config_reader():
    node_name = "eta"
    config = read_default_params(node_name)

    assert isinstance(config, dict)
    assert len(config) > 0

    node_name = "not_exist_config"
    config = read_default_params(node_name)

    assert isinstance(config, dict)
    assert len(config) == 0


def test_config_writer():
    filepath = join_filepath([dirpath, filename])

    if os.path.exists(filepath):
        os.remove(filepath)

    ConfigWriter.write(dirpath, filename, {"test": "test"})

    assert os.path.exists(filepath)
