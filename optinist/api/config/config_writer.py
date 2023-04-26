import yaml

from optinist.api.utils.filepath_creater import create_directory, join_filepath


class ConfigWriter:
    @classmethod
    def write(cls, dirname, filename, config):
        create_directory(dirname)

        with open(join_filepath([dirname, filename]), "w") as f:
            yaml.dump(config, f)
