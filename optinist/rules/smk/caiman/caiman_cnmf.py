rule:
    input:
        config["rules"]["caiman_cnmf"]["input"]
    output:
        config["rules"]["caiman_cnmf"]["output"]
    # conda:
    #     "../../envs/caiman_env.yaml"
    script:
        "../../scripts/caiman/caiman_cnmf.py"