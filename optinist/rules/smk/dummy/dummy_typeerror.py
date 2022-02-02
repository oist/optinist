rule:
    input:
        config["rules"]["dummy_typeerror"]["input"]
    output:
        config["rules"]["dummy_typeerror"]["output"]
    # run:
    #     __func_config = config["rules"]["dummy_typeerror"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_typeerror.py"