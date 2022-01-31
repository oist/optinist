from rules.utils import run_script

rule:
    input:
        config["rules"]["suite2p_registration"]["input"]
    output:
        touch(config["rules"]["suite2p_registration"]["output"])
    run:
        __func_config = config["rules"]["suite2p_registration"]
        run_script(__func_config)
