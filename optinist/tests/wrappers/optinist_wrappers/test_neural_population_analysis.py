import os
import shutil

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import SmkParam
from optinist.api.snakemake.snakemake_executor import snakemake_execute


def test_neural_population_analysis_config():
    """
    cross_correlation
    correlation
    granger
    """
    unique_id = "ae974587"
    smk_param = SmkParam(
        use_conda=True,
        cores=2,
        forceall=True,
        forcetargets=True,
        lock=False,
    )

    shutil.copyfile(
        f"{DIRPATH.OPTINIST_DIR}/neural_population_analysis_config.yaml",
        f"{DIRPATH.ROOT_DIR}/config.yaml",
    )
    snakemake_execute(unique_id, smk_param)

    OUTPUT_DIRPATH = f"{DIRPATH.OUTPUT_DIR}/{unique_id}"
    assert os.path.exists(f"{OUTPUT_DIRPATH}/granger/granger.pkl")
    assert os.path.exists(f"{OUTPUT_DIRPATH}/correlation/correlation.pkl")
    assert os.path.exists(f"{OUTPUT_DIRPATH}/cross_correlation/cross_correlation.pkl")
