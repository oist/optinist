from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk_dir import smk_input, smk_output

name = "glm"

rule:
    input:
        smk_input(config, name)
    output:
        smk_output(config, name)
    conda:
        f'{DIRPATH.ROOT_DIR}/rules/envs/optinist_env.yaml'
    params:
        name = name
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/func.py'