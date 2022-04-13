import sys
sys.path.append('../optinist')

from optinist.api.pickle.pickle_writer import PickleWriter
from optinist.rules.scripts.file_writer import FileWriter
from optinist.api.snakemake.snakemake_reader import RuleConfigReader
from optinist.routers.model import FILETYPE


if __name__ == '__main__':
    last_output = snakemake.config["last_output"]

    for rule_config in snakemake.config["rules"].values():
        rule_config = RuleConfigReader.read(rule_config)
        
        if rule_config.type in [FILETYPE.CSV, FILETYPE.BEHAVIOR]:
            outputfile = FileWriter.csv(rule_config, rule_config.type)
            PickleWriter.write(rule_config.output, outputfile)
        elif rule_config.type == FILETYPE.IMAGE:
            outputfile = FileWriter.image(rule_config)
            PickleWriter.write(rule_config.output, outputfile)
        elif rule_config.type == FILETYPE.HDF5:
            outputfile = FileWriter.hdf5(rule_config)
            PickleWriter.write(rule_config.output, outputfile)
