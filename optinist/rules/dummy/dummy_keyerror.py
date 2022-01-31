from rules.utils import run_script

rule:
    input:
        config["rules"]["dummy_keyerror"]["input"]
    output:
        touch(config["rules"]["dummy_keyerror"]["output"])
    run:
        __func_config = config["rules"]["dummy_keyerror"]
        run_script(__func_config)
