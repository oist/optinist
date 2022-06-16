import pytest
import shutil

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import SmkParam
from optinist.api.snakemake.snakemake_executor import snakemake_execute


tiff_filename = "test.tif"

smk_param = SmkParam(
    use_conda=False,
    cores=2,
    forceall=True,
    forcetargets=True,
    lock=False,
)

shutil.copyfile(
    "/tmp/optinist/config.yaml",
    f"{DIRPATH.ROOT_DIR}/config.yaml",
)


def test_snakemake_execute():
    snakemake_execute(smk_param)
