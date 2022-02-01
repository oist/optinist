rule:
    input:
        config["rules"]["svm"]["input"]
    output:
        config["rules"]["svm"]["output"]
    # run:
    #     __func_config = config["rules"]["svm"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/neural_decoding/svm.py"