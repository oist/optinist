import os
import yaml

from optinist.api.utils.filepath_creater import join_filepath


class ConfigWriter:
    @classmethod
    def write(cls, dirname, filename, config):
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        with open(join_filepath([dirname, filename]), "w") as f:
            yaml.dump(config, f)
