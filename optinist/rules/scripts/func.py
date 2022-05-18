import sys
sys.path.append('../optinist')

from optinist.api.snakemake.snakemake_reader import RuleConfigReader
from optinist.rules.scripts.runner import Runner
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath


if __name__ == '__main__':
    last_output = snakemake.config["last_output"]

    for rule_config in snakemake.config["rules"].values():
        rule_config = RuleConfigReader.read(rule_config)

        rule_config.input = [
            join_filepath([DIRPATH.OUTPUT_DIR, x])
            for x in rule_config.input
        ]
        rule_config.output = join_filepath([DIRPATH.OUTPUT_DIR, rule_config.output])

        if rule_config.type == snakemake.params.name:  
            Runner.run(rule_config, last_output)
