from rules.utils import run_script

rule:
    input:
        config["rules"]["lda"]["input"]
    output:
        touch(config["rules"]["lda"]["output"])
    run:
        __func_config = config["rules"]["lda"]
        run_script(__func_config)
