# from rules.utils import run_script

rule:
    input:
        config["rules"]["dummy_image2roi"]["input"]
    output:
        touch(config["rules"]["dummy_image2roi"]["output"])
    # run:
    #     __func_config = config["rules"]["dummy_image2roi"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_image2roi.py"