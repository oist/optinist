rule suite2p_registration:
    input: 
        config["rules"]["suite2p_registration"]["input"]
    output: 
        config["rules"]["suite2p_registration"]["output"]
    script:
        "../scripts/suite2p_registration.py"