import os

import yaml


class ConfigReader:
    @classmethod
    def read(cls, filepath):
        config = {}
        if filepath is not None and os.path.exists(filepath):
            with open(filepath) as f:
                config = yaml.safe_load(f)
        return config
