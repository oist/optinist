import sys
sys.path.append('../optinist')

if __name__ == '__main__':
    from rules.utils import run_script
    __func_config = snakemake.config["rules"]["caiman_mc"]
    last_output = snakemake.config["last_output"]
    run_script(__func_config, last_output)
