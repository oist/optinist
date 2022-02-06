import os
import yaml
from .const import OPTINIST_DIR


def get_config_params(filename):
    filepath = os.path.join(OPTINIST_DIR, 'config', f'{filename}.yaml')
    print(filepath)
    if os.path.exists(filepath):
        with open(filepath) as f:
            config = yaml.safe_load(f)
        return config
    else:
        return {}
