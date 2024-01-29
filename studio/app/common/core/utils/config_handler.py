import os

import yaml

from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)


class ConfigReader:
    @classmethod
    def read(cls, filepath):
        config = {}
        if filepath is not None and os.path.exists(filepath):
            with open(filepath) as f:
                config = yaml.safe_load(f)
        return config


class ConfigWriter:
    @classmethod
    def write(cls, dirname, filename, config):
        create_directory(dirname)

        with open(join_filepath([dirname, filename]), "w") as f:
            yaml.dump(config, f, sort_keys=False)
