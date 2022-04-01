from cui_api.const import OPTINIST_DIR

name = "caiman_cnmf"

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
        f'{OPTINIST_DIR}/rules/scripts/func.py'