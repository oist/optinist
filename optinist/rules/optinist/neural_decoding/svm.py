from rules.utils import run_script

rule:
    input:
        config["rules"]["svm"]["input"]
    output:
        touch(config["rules"]["svm"]["output"])
    run:
        __func_config = config["rules"]["svm"]
        run_script(__func_config)
