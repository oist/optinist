import os
from .utils import join_file_path

_OPTINIST_DIR = os.environ.get('OPTINIST_DIR')
if _OPTINIST_DIR is None:
    BASE_DIR = '/tmp/optinist'
else:
    BASE_DIR = _OPTINIST_DIR

OPTINIST_DIR = join_file_path([os.path.dirname(__file__), ".."])

if not os.path.exists(BASE_DIR):
    os.makedirs(BASE_DIR, exist_ok=True)