rule:
    input:
        config["rules"]["eta"]["input"]
    output:
        config["rules"]["eta"]["output"]
    script:
        "../../../scripts/optinist/basic_neural_analysis/eta.py"