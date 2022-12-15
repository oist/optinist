import sys
from const import studio_DIRPATH
sys.path.append(studio_DIRPATH)

from studio.api.snakemake.snakemake_reader import RuleConfigReader
from studio.api.rules.runner import Runner
from studio.api.dir_path import DIRPATH
from studio.api.utils.filepath_creater import join_filepath


if __name__ == '__main__':
    last_output = [
        join_filepath([DIRPATH.OUTPUT_DIR, x])
        for x in snakemake.config["last_output"]
    ]

    rule_config = RuleConfigReader.read(snakemake.params.name)

    rule_config.input = snakemake.input
    rule_config.output = snakemake.output[0]

    Runner.run(rule_config, last_output)
