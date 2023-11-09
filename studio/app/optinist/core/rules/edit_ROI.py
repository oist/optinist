# flake8: noqa
# Exclude from lint for the following reason
# This file is executed by snakemake and cause the following lint errors
# - E402: sys.path.append is required to import optinist modules
# - F821: do not import snakemake
import sys
from os.path import abspath, dirname

ROOT_DIRPATH = dirname(dirname(dirname(dirname(dirname(dirname(abspath(__file__)))))))
sys.path.append(ROOT_DIRPATH)

from studio.app.optinist.core.edit_ROI import EditROI

if __name__ == "__main__":
    config = snakemake.config
    EditROI(file_path=config["file_path"]).commit()
