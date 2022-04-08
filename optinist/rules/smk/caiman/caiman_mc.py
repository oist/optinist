from optinist.cui_api.dir_path import DIRPATH

name = "caiman_mc"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    params:
        name = name
    # conda:
    #     "../../envs/caiman_env.yaml"
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/func.py'