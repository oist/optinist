rule:
    input:
        config["rules"]["dummy_image2image"]["input"]
    output:
        config["rules"]["dummy_image2image"]["output"]
    # run:
    #     __func_config = config["rules"]["dummy_image2image"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_image2image.py"