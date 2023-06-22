# flake8: noqa
# Exclude from lint for the following reason
# This file is executed by snakemake and cause the following lint errors
# - E402: sys.path.append is required to import optinist modules
# - F821: do not import snakemake
import sys

from const import OPTINIST_DIRPATH

sys.path.append(OPTINIST_DIRPATH)

from optinist.api.edit_ROI import EditRoiUtils

if __name__ == "__main__":
    config = snakemake.config
    EditRoiUtils.excute(config)
