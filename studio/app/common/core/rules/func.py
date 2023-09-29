# flake8: noqa
# Exclude from lint for the following reason
# This file is executed by snakemake and cause the following lint errors
# - E402: sys.path.append is required to import optinist modules
# - F821: do not import snakemake
import json
import os
import sys
from os.path import abspath, dirname
from pathlib import Path

ROOT_DIRPATH = dirname(dirname(dirname(dirname(dirname(dirname(abspath(__file__)))))))

sys.path.append(ROOT_DIRPATH)

from studio.app.common.core.rules.runner import Runner
from studio.app.common.core.snakemake.snakemake_reader import RuleConfigReader
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.dir_path import DIRPATH

if __name__ == "__main__":
    last_output = [
        join_filepath([DIRPATH.OUTPUT_DIR, x]) for x in snakemake.config["last_output"]
    ]

    rule_config = RuleConfigReader.read(snakemake.params.name)

    rule_config.input = snakemake.input
    rule_config.output = snakemake.output[0]

    # save snakemake script file path and PID of current running algo function
    pid_data = {"last_pid": os.getpid(), "last_script_file": sys.argv[0]}
    with open(Path(rule_config.output).parent.parent / "pid.json", "w") as f:
        json.dump(pid_data, f)

    Runner.run(rule_config, last_output)
