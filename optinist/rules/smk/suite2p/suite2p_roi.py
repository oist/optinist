name = "suite2p_roi"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    params:
        name = name
    conda:
        "../../envs/suite2p_env.yaml"
    script:
        '../../scripts/func.py'
