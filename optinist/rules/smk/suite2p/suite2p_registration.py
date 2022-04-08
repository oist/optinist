from optinist.cui_api.dir_path import DIRPATH

name = "suite2p_registration"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    params:
        name = name
    conda:
        f'{DIRPATH.ROOT_DIR}/rules/envs/suite2p_env.yaml'
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/func.py'