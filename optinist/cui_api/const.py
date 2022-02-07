import os
from .utils import join_file_path

BASE_DIR = '/tmp/optinist'
OPTINIST_DIR = join_file_path([os.path.dirname(__file__), ".."])

if not os.path.exists(BASE_DIR):
    os.makedirs(BASE_DIR, exist_ok=True)