rule:
    input:
        config["rules"]["cca"]["input"]
    output:
        touch(config["rules"]["cca"]["output"])
    # run:
    #     __func_config = config["rules"]["cca"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/dimension_reduction/cca.py"