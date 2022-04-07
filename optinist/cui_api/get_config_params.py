import os
import yaml
from optinist.cui_api.const import OPTINIST_DIR
from optinist.cui_api.utils import join_file_path


def get_config_params(filename):
    filepath = join_file_path([OPTINIST_DIR, 'config', f'{filename}.yaml'])
    print(filepath)
    if os.path.exists(filepath):
        with open(filepath) as f:
            config = yaml.safe_load(f)
        return config
    else:
        return {}
