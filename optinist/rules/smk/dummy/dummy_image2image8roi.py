# from rules.utils import run_script

rule:
    input:
        config["rules"]["dummy_image2image8roi"]["input"]
    output:
        touch(config["rules"]["dummy_image2image8roi"]["output"])
    # run:
    #     __func_config = config["rules"]["dummy_image2image8roi"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_image2image8roi.py"