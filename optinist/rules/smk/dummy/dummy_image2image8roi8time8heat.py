from cui_api.const import OPTINIST_DIR

name = "dummy_image2image8roi8time8heat"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    params:
        name = name
    script:
        f'{OPTINIST_DIR}/rules/scripts/func.py'