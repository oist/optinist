import sys
sys.path.append('../optinist')

if __name__ == '__main__':
    from rules.utils import run_script
    __func_config = snakemake.config["rules"]["dummy_time2time"]
    run_script(__func_config)
