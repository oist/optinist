rule suite2p_spike_deconv:
    input: 
        config["rules"]["suite2p_spike_deconv"]["input"]
    output: 
        config["rules"]["suite2p_spike_deconv"]["output"]
    script:
        "../scripts/suite2p_spike_deconv.py"