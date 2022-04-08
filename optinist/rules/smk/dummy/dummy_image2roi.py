from cui_api.const import ROOT_DIR

name = "dummy_image2roi"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    params:
        name = name
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/func.py'