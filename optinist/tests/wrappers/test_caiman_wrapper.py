import os
import shutil

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import SmkParam
from optinist.api.snakemake.snakemake_executor import snakemake_execute


def test_caiman_mc_cnmf():
    unique_id = "b9d1513e"
    smk_param = SmkParam(
        use_conda=True,
        cores=2,
        forceall=True,
        forcetargets=True,
        lock=False,
    )

    shutil.copyfile(
        f"{DIRPATH.OPTINIST_DIR}/caiman_config.yaml",
        f"{DIRPATH.ROOT_DIR}/config.yaml",
    )
    snakemake_execute(unique_id, smk_param)

    OUTPUT_DIRPATH = f"{DIRPATH.OUTPUT_DIR}/{unique_id}"
    assert os.path.exists(f"{OUTPUT_DIRPATH}/caiman_mc/caiman_mc.pkl")
    assert os.path.exists(f"{OUTPUT_DIRPATH}/caiman_cnmf/caiman_cnmf.pkl")
    assert os.path.exists(f"{OUTPUT_DIRPATH}/caiman_cnmf/tiff")
