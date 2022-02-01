# from rules.utils import run_script

rule:
    input:
        config["rules"]["dummy_image2image8roi8time8heat"]["input"]
    output:
        touch(config["rules"]["dummy_image2image8roi8time8heat"]["output"])
    # run:
    #     __func_config = config["rules"]["dummy_image2image8roi8time8heat"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_image2image8roi8time8heat.py"