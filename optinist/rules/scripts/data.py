import sys

sys.path.append('../optinist')

from optinist.api.pickle.pickle_writer import PickleWriter
from optinist.rules.scripts.file_writer import FileWriter
from optinist.api.snakemake.snakemake_reader import RuleConfigReader
from optinist.routers.model import FILETYPE
from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath


if __name__ == '__main__':
    last_output = snakemake.config["last_output"]

    rule_config = RuleConfigReader.read(snakemake.params.name)
    if rule_config.type in [FILETYPE.IMAGE]:
        rule_config.input = [
            join_filepath([DIRPATH.INPUT_DIR, x])
            for x in rule_config.input
        ]
    elif rule_config.type in [FILETYPE.CSV, FILETYPE.BEHAVIOR, FILETYPE.HDF5]:
        rule_config.input = join_filepath([DIRPATH.INPUT_DIR, rule_config.input])

    rule_config.output = join_filepath([DIRPATH.OUTPUT_DIR, rule_config.output])

    if rule_config.type in [FILETYPE.CSV, FILETYPE.BEHAVIOR]:
        outputfile = FileWriter.csv(rule_config, rule_config.type)
        PickleWriter.write(rule_config.output, outputfile)
    elif rule_config.type == FILETYPE.IMAGE:
        outputfile = FileWriter.image(rule_config)
        PickleWriter.write(rule_config.output, outputfile)
    elif rule_config.type == FILETYPE.HDF5:
        outputfile = FileWriter.hdf5(rule_config)
        PickleWriter.write(rule_config.output, outputfile)
