import os

_DEFAULT_DIR = '/tmp/optinist'
_ENV_DIR = os.environ.get('OPTINIST_DIR')

class DIRPATH:
    BASE_DIR = _DEFAULT_DIR if _ENV_DIR is None else _ENV_DIR
    ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
    CONFIG_DIR = f"{ROOT_DIR}/config"

    if not os.path.exists(BASE_DIR):
        os.makedirs(BASE_DIR, exist_ok=True)
    assert os.path.exists(BASE_DIR)

    SNAKEMAKE_FILEPATH = f"{ROOT_DIR}/Snakefile"