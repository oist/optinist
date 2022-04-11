import pytest

from optinist.api.dir_path import DIRPATH
from optinist.api.config.config_reader import ConfigReader
from optinist.api.utils.filepath_creater import join_filepath

def test_config_reader():
    filename = "eta"
    filepath = join_filepath([DIRPATH.ROOT_DIR, 'config', f'{filename}.yaml'])
    config = ConfigReader.read(filepath)

    assert isinstance(config, dict)
    assert len(config) > 0

    filename = "a"
    filepath = join_filepath([DIRPATH.ROOT_DIR, 'config', f'{filename}.yaml'])
    config = ConfigReader.read(filepath)

    assert isinstance(config, dict)
    assert len(config) == 0