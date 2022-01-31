from rules.utils import run_script

rule:
    input:
        config["rules"]["suite2p_file_convert"]["input"]
    output:
        touch(config["rules"]["suite2p_file_convert"]["output"])
    run:
        __func_config = config["rules"]["suite2p_file_convert"]
        run_script(__func_config)
