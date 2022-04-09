import os
import yaml
from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import join_filepath


class ConfigReader:
    @classmethod
    def read(cls, filepath):
        config = {}
        if os.path.exists(filepath):
            with open(filepath) as f:
                config = yaml.safe_load(f)
        return config
