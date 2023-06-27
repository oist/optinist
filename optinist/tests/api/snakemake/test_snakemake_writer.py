import os

from optinist.services.dir_path import DIRPATH
from optinist.services.snakemake.snakemake_writer import SmkConfigWriter

unique_id = "smk_test"
output_fileapth = f"{DIRPATH.OPTINIST_DIR}/output/{unique_id}/config.yaml"


def test():
    if os.path.exists(output_fileapth):
        os.remove(output_fileapth)

    SmkConfigWriter.write(
        unique_id,
        {"test"},
    )

    assert os.path.exists(output_fileapth)
