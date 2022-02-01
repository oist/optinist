# from rules.utils import run_script

rule:
    input:
        config["rules"]["suite2p_registration"]["input"]
    output:
        touch(config["rules"]["suite2p_registration"]["output"])
    # run:
    #     __func_config = config["rules"]["suite2p_registration"]
    #     run_script(__func_config)
    conda:
        "../envs/suite2p_env.yaml"
    script:
        "../../scripts/suite2p/suite2p_registration.py"