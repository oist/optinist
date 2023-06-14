import os
import shutil

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import SmkParam
from optinist.api.snakemake.snakemake_executor import snakemake_execute


def test_dimension_reduction_config():
    """
    cca
    pca
    tsne
    """
    unique_id = "281f3b24"
    smk_param = SmkParam(
        use_conda=True,
        cores=2,
        forceall=True,
        forcetargets=True,
        lock=False,
    )

    shutil.copyfile(
        f"{DIRPATH.OPTINIST_DIR}/dimension_reduction_config.yaml",
        f"{DIRPATH.ROOT_DIR}/config.yaml",
    )
    snakemake_execute(unique_id, smk_param)

    OUTPUT_DIRPATH = f"{DIRPATH.OUTPUT_DIR}/{unique_id}"
    assert os.path.exists(f"{OUTPUT_DIRPATH}/cca/cca.pkl")
    assert os.path.exists(f"{OUTPUT_DIRPATH}/pca/pca.pkl")
    assert os.path.exists(f"{OUTPUT_DIRPATH}/tsne/tsne.pkl")
