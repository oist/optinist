import os
import yaml


def get_params(filepath):
    config = {}
    if os.path.exists(filepath):
        with open(filepath) as f:
            config = yaml.safe_load(f)
    return config

