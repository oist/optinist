rule:
    input:
        config["rules"]["dummy_image2image8roi"]["input"]
    output:
        config["rules"]["dummy_image2image8roi"]["output"]
    # run:
    #     __func_config = config["rules"]["dummy_image2image8roi"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_image2image8roi.py"