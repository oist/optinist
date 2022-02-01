rule:
    input:
        config["rules"]["granger"]["input"]
    output:
        config["rules"]["granger"]["output"]
    # run:
    #     __func_config = config["rules"]["svm"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/neural_population_analysis/granger.py"