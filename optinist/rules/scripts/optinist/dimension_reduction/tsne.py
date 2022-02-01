import sys
sys.path.append('../optinist')

if __name__ == '__main__':
    from rules.utils import run_script
    __func_config = snakemake.config["rules"]["tsne"]
    run_script(__func_config)
