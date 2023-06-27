from studio.services.config.config_reader import ConfigReader
from studio.services.utils.filepath_finder import find_param_filepath


def test_reader():
    filename = "eta"
    filepath = find_param_filepath(filename)
    config = ConfigReader.read(filepath)

    assert isinstance(config, dict)
    assert len(config) > 0

    filename = "not_exist_config"
    filepath = find_param_filepath(filename)
    config = ConfigReader.read(filepath)

    assert isinstance(config, dict)
    assert len(config) == 0
