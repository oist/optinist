rule:
    input:
        config["rules"]["pca"]["input"]
    output:
        touch(config["rules"]["pca"]["output"])
    # run:
    #     __func_config = config["rules"]["pca"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/dimension_reduction/pca.py"