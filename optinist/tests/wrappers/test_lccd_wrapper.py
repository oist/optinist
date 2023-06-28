import os
import shutil

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import SmkParam
from optinist.api.snakemake.snakemake_executor import snakemake_execute


def test_lccd_cell_detection():
    unique_id = "57d9349b"
    smk_param = SmkParam(
        use_conda=True,
        cores=2,
        forceall=True,
        forcetargets=True,
        lock=False,
    )
    shutil.copyfile(
        f"{DIRPATH.OPTINIST_DIR}/config/lccd_config.yaml",
        f"{DIRPATH.ROOT_DIR}/config.yaml",
    )
    snakemake_execute(unique_id, smk_param)

    OUTPUT_DIRPATH = f"{DIRPATH.OUTPUT_DIR}/{unique_id}"
    assert os.path.exists(
        f"{OUTPUT_DIRPATH}/lccd_cell_detection/lccd_cell_detection.pkl"
    )
    assert os.path.exists(f"{OUTPUT_DIRPATH}/lccd_cell_detection/tiff")
