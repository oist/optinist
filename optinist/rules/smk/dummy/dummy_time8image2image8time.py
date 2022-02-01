# from rules.utils import run_script

rule:
    input:
        config["rules"]["dummy_time8image2image8time"]["input"]
    output:
        touch(config["rules"]["dummy_time8image2image8time"]["output"])
    # run:
    #     __func_config = config["rules"]["dummy_time8image2image8time"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_time8image2image8time.py"