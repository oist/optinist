from optinist.cui_api.dir_path import DIRPATH
name="csv"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    params:
        name=name
    script:
        f'{ROOT_DIR/rules/scripts/data.py'