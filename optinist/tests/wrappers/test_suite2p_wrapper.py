import os
import shutil

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import SmkParam
from optinist.api.snakemake.snakemake_executor import snakemake_execute


def test_suite2p():
    """
    suite2p_file_convert -> suite2p_registration -> suite2p_roi -> suite2p_spike_deconv
    """
    unique_id = "fd536a4f"
    smk_param = SmkParam(
        use_conda=True,
        cores=2,
        forceall=True,
        forcetargets=True,
        lock=False,
    )
    shutil.copyfile(
        f"{DIRPATH.OPTINIST_DIR}/config/suite2p_config.yaml",
        f"{DIRPATH.ROOT_DIR}/config.yaml",
    )
    snakemake_execute(unique_id, smk_param)

    OUTPUT_DIRPATH = f"{DIRPATH.OUTPUT_DIR}/{unique_id}"
    assert os.path.exists(
        f"{OUTPUT_DIRPATH}/suite2p_file_convert/suite2p_file_convert.pkl"
    )
    assert os.path.exists(
        f"{OUTPUT_DIRPATH}/suite2p_registration/suite2p_registration.pkl"
    )
    assert os.path.exists(f"{OUTPUT_DIRPATH}/suite2p_roi/suite2p_roi.pkl")
    assert os.path.exists(
        f"{OUTPUT_DIRPATH}/suite2p_spike_deconv/suite2p_spike_deconv.pkl"
    )
