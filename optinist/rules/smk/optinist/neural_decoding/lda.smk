from optinist.api.dir_path import DIRPATH

name = "lda"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    conda:
        f'{DIRPATH.ROOT_DIR}/rules/envs/optinist_env.yaml'
    params:
        name = name
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/func.py'