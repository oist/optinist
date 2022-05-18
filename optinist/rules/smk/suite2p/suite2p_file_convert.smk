from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk_dir import smk_input, smk_output

name = "suite2p_file_convert"

rule:
    input:
        smk_input(config, name)
    output:
        smk_output(config, name)
    params:
        name = name
    conda:
        f'{DIRPATH.ROOT_DIR}/rules/envs/suite2p_env.yaml'
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/func.py'