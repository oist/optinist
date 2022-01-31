from rules.utils import run_script

rule:
    input:
        config["rules"]["dummy_image2image"]["input"]
    output:
        touch(config["rules"]["dummy_image2image"]["output"])
    run:
        __func_config = config["rules"]["dummy_image2image"]
        run_script(__func_config)
