# from rules.utils import run_script

rule:
    input:
        config["rules"]["dummy_image2time"]["input"]
    output:
        touch(config["rules"]["dummy_image2time"]["output"])
    # run:
    #     __func_config = config["rules"]["dummy_image2time"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_image2time.py"