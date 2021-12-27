rule suite2p_fileconvert:
    input: 
        config["rules"]["suite2p_file_convert"]["input"]
    output: 
        config["rules"]["suite2p_file_convert"]["output"]
    script:
        "../scripts/suite2p_file_convert.py"