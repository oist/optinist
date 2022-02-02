rule:
    input:
        config["rules"]["dummy_keyerror"]["input"]
    output:
        config["rules"]["dummy_keyerror"]["output"]
    # run:
    #     __func_config = config["rules"]["dummy_keyerror"]
    #     run_script(__func_config)
    script:
        "../../scripts/dummy/dummy_keyerror.py"