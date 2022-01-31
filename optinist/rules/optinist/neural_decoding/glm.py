from rules.utils import run_script

rule:
    input:
        config["rules"]["glm"]["input"]
    output:
        touch(config["rules"]["glm"]["output"])
    run:
        __func_config = config["rules"]["glm"]
        run_script(__func_config)
