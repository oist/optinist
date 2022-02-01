rule:
    input:
        config["rules"]["dummy_image8image2image8time"]["input"]
    output:
        config["rules"]["dummy_image8image2image8time"]["output"]
    # run:
    #     __func_config = config["rules"]["dummy_image8image2image8time"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_image8image2image8time.py"