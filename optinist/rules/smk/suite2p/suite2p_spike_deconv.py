# from rules.utils import run_script

rule:
    input:
        config["rules"]["suite2p_spike_deconv"]["input"]
    output:
        touch(config["rules"]["suite2p_spike_deconv"]["output"])
    # run:
    #     __func_config = config["rules"]["suite2p_spike_deconv"]
    #     run_script(__func_config)
    conda:
        "../envs/suite2p_env.yaml"
    script:
        "../../scripts/suite2p/suite2p_spike_deconv.py"