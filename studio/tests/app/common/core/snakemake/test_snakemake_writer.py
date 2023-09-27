import os

from studio.app.common.core.snakemake.snakemake_writer import SmkConfigWriter
from studio.app.dir_path import DIRPATH

workspace_id = "default"
unique_id = "smk_test"
output_fileapth = f"{DIRPATH.DATA_DIR}/output/{workspace_id}/{unique_id}/snakemake.yaml"


def test():
    if os.path.exists(output_fileapth):
        os.remove(output_fileapth)

    SmkConfigWriter.write(
        workspace_id,
        unique_id,
        {"test"},
    )

    assert os.path.exists(output_fileapth)
