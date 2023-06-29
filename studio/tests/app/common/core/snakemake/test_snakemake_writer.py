import os

from studio.app.common.core.snakemake.snakemake_writer import SmkConfigWriter
from studio.config.dir_path import DIRPATH

unique_id = "smk_test"
output_fileapth = f"{DIRPATH.DATA_DIR}/output/{unique_id}/config.yaml"


def test():
    if os.path.exists(output_fileapth):
        os.remove(output_fileapth)

    SmkConfigWriter.write(
        unique_id,
        {"test"},
    )

    assert os.path.exists(output_fileapth)
