rule:
    input:
        config["rules"]["lda"]["input"]
    output:
        touch(config["rules"]["lda"]["output"])
    # run:
    #     __func_config = config["rules"]["lda"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/neural_decoding/lda.py"