rule:
    input:
        config["rules"]["cross_correlation"]["input"]
    output:
        touch(config["rules"]["cross_correlation"]["output"])
    # run:
    #     __func_config = config["rules"]["svm"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/neural_population_analysis/cross_correlation.py"