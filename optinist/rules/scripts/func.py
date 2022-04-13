import sys
sys.path.append('../optinist')

from optinist.api.snakemake.snakemake_reader import RuleConfigReader
from optinist.rules.scripts.runner import Runner

if __name__ == '__main__':
    last_output = snakemake.config["last_output"]

    for rule_config in snakemake.config["rules"].values():
        rule_config = RuleConfigReader.read(rule_config)

        if rule_config.type == snakemake.params.name:  
            Runner.run(rule_config, last_output)
