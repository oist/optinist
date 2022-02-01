import sys
sys.path.append('../optinist')

if __name__ == '__main__':
    from rules.utils import run_script
    __func_config = snakemake.config["rules"]["dummy_image2image8roi8time8heat"]
    run_script(__func_config)
