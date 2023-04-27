import os

from optinist.api.snakemake.snakemake_writer import SmkConfigWriter

unique_id = "smk_test"
output_fileapth = f"/tmp/optinist/output/{unique_id}/config.yaml"


def test():
    if os.path.exists(output_fileapth):
        os.remove(output_fileapth)

    SmkConfigWriter.write(
        unique_id,
        {"test"},
    )

    assert os.path.exists(output_fileapth)
