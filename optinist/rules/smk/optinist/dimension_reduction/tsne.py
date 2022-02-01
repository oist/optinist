rule:
    input:
        config["rules"]["tsne"]["input"]
    output:
        touch(config["rules"]["tsne"]["output"])
    # run:
    #     __func_config = config["rules"]["tsne"]
    #     run_script(__func_config)
    script:
        "../../../scripts/optinist/dimension_reduction/tsne.py"