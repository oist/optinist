import sys
sys.path.append('../optinist')

from optinist.api.pickle.pickle_writer import PickleWriter
from optinist.rules.scripts.file_writer import FileWriter
from optinist.api.snakemake.snakemake_reader import SmkConfigReader
from optinist.routers.model import FILETYPE
from optinist.rules.scripts.runner import Runner

if __name__ == '__main__':
    last_output = snakemake.config["last_output"]

    for rule_config in snakemake.config["rules"].values():
        rule_config = SmkConfigReader.read(rule_config)

        if rule_config.type == snakemake.params.name:  
            Runner.run(rule_config, last_output)
