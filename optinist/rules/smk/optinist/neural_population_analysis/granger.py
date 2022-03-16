name = "granger"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    conda:
        "../../../envs/optinist_env.yaml"
    params:
        name = name
    script:
        '../../../scripts/func.py'