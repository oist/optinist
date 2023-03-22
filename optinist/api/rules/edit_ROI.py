import sys
from const import OPTINIST_DIRPATH
sys.path.append(OPTINIST_DIRPATH)

from optinist.wrappers.edit_roi import EditRoiUtils

if __name__ == '__main__':
    config = snakemake.config
    EditRoiUtils.excute(config)