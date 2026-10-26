import os

import yaml
from filelock import FileLock

from studio.app.common.core.utils.filelock_handler import FileLockUtils
from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)


def differential_deep_merge(d1: dict, d2: dict) -> dict:
    """
    Deep merge only the differences to avoid destroying existing elements
    """
    result = d1.copy()
    for key, value in d2.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = differential_deep_merge(result[key], value)
        else:
            result[key] = value
    return result


class ConfigReader:
    @classmethod
    def read(cls, filepath: str) -> dict:
        config = {}

        if filepath is not None and os.path.exists(filepath):
            with open(filepath) as f:
                loaded_config = yaml.safe_load(f)
                if loaded_config is not None:
                    config = loaded_config

        return config

    @classmethod
    def read_from_bytes(cls, content: bytes) -> dict:
        config = yaml.safe_load(content)
        return config


class ConfigWriter:
    FILE_LOCK_TIMEOUT = 60

    @classmethod
    def write(cls, dirname: str, filename: str, config: dict, auto_file_lock=True):
        create_directory(dirname)

        config_path = join_filepath([dirname, filename])

        if auto_file_lock:
            # Exclusive control for parallel updates from multiple processes.
            lock_path = FileLockUtils.get_lockfile_path(config_path)
            with FileLock(lock_path, cls.FILE_LOCK_TIMEOUT):
                cls.__write(config_path, config)
        else:
            cls.__write(config_path, config)

    @classmethod
    def __write(cls, config_path: str, config: dict):
        config_tmp_path = f"{config_path}.tmp"

        # First write to a temporary file
        # (a measure to avoid read conflicts due to write delays)
        with open(config_tmp_path, "w") as f:
            yaml.dump(config, f, sort_keys=False)

        # Write to the original file path
        # (write atomically by using os.replace)
        os.replace(config_tmp_path, config_path)
