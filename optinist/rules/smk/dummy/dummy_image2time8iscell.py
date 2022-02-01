rule:
    input:
        config["rules"]["dummy_image2time8iscell"]["input"]
    output:
        config["rules"]["dummy_image2time8iscell"]["output"]
    # run:
    #     __func_config = config["rules"]["dummy_image2time8iscell"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_image2time8iscell.py"