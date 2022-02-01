rule:
    input:
        config["rules"]["glm"]["input"]
    output:
        config["rules"]["glm"]["output"]
    # run:
    #     __func_config = config["rules"]["glm"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/neural_decoding/glm.py"