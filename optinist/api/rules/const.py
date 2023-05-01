import os

optinist_dirname = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
OPTINIST_DIRPATH = f"{optinist_dirname}"
