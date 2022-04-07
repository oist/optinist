import sys
sys.path.append('../optinist')
from optinist.rules.utils import run_script

if __name__ == '__main__':
    last_output = snakemake.config["last_output"]

    for rule in snakemake.config["rules"].values():
        if rule['type'] == snakemake.params.name:       
            run_script(rule, last_output)
