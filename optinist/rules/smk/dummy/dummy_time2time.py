# from rules.utils import run_script

rule:
    input:
        config["rules"]["dummy_time2time"]["input"]
    output:
        touch(config["rules"]["dummy_time2time"]["output"])
    # run:
    #     __func_config = config["rules"]["dummy_time2time"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_time2time.py"