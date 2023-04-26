import os

import yaml


class ConfigReader:
    @classmethod
    def read(cls, filepath):
        config = {}
        if os.path.exists(filepath):
            with open(filepath) as f:
                config = yaml.safe_load(f)
        return config
