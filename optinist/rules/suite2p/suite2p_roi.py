from rules.utils import run_script

rule:
    input:
        config["rules"]["suite2p_roi"]["input"]
    output:
        touch(config["rules"]["suite2p_roi"]["output"])
    run:
        __func_config = config["rules"]["suite2p_roi"]
        run_script(__func_config)
