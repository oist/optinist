rule:
    input:
        config["rules"]["correlation"]["input"]
    output:
        touch(config["rules"]["correlation"]["output"])
    # run:
    #     __func_config = config["rules"]["svm"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/neural_population_analysis/correlation.py"