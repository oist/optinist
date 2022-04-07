from optinist.cui_api.const import OPTINIST_DIR

name = "cross_correlation"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    conda:
        f'{OPTINIST_DIR}/rules/envs/optinist_env.yaml'
    params:
        name = name
    script:
        f'{OPTINIST_DIR}/rules/scripts/func.py'