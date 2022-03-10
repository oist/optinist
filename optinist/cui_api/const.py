import os
from .utils import join_file_path

TMP_DIR = os.environ.get('OPTINIST_DIR')
print("OPTINIST_DIR: ", TMP_DIR)
if TMP_DIR is None:
    BASE_DIR = '/tmp/optinist'
else:
    BASE_DIR = TMP_DIR

OPTINIST_DIR = join_file_path([os.path.dirname(__file__), ".."])

if not os.path.exists(BASE_DIR):
    os.makedirs(BASE_DIR, exist_ok=True)

assert os.path.exists(BASE_DIR)