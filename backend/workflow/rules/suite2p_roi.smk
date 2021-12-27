rule suite2p_roi:
    input: 
        config["rules"]["suite2p_roi"]["input"]
    output: 
        config["rules"]["suite2p_roi"]["output"]
    script:
        "../scripts/suite2p_roi.py"