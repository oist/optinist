
from dataclasses import dataclass

from optinist.cui_api.dir_path import DIRPATH
from optinist.cui_api.filepath_creater import join_filepath
from optinist.cui_api.config_writer import ConfigWriter


@dataclass
class Rule:
    input: list
    output: str
    type: str
    rule_file: str
    return_arg: str
    params: list
    nwbfile: list


@dataclass
class SnakemakeConfig:
    last_output: list
    rules: Rule


def snakemake_config_writer(unique_id, flow_config):
    ConfigWriter.write(
        dirname=DIRPATH.ROOT_DIR,
        filename='config.yaml',
        config=flow_config
    )

    ConfigWriter.write(
        dirname=join_filepath([DIRPATH.BASE_DIR, unique_id]),
        filename='config.yaml',
        config=flow_config
    )
