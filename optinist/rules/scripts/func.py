import sys
sys.path.append('../optinist')

from optinist.api.snakemake.snakemake_reader import RuleConfigReader
from optinist.rules.scripts.runner import Runner
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath


if __name__ == '__main__':
    last_output = [
        join_filepath([DIRPATH.OUTPUT_DIR, x])
        for x in snakemake.config["last_output"]
    ]

    rule_config = RuleConfigReader.read(snakemake.params.name)

    rule_config.input = [
        join_filepath([DIRPATH.OUTPUT_DIR, x])
        for x in rule_config.input
    ]
    rule_config.output = join_filepath([DIRPATH.OUTPUT_DIR, rule_config.output])

    Runner.run(rule_config, last_output)
