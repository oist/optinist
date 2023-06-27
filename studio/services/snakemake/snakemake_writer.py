from studio.services.config.config_writer import ConfigWriter
from studio.services.dir_path import DIRPATH
from studio.services.utils.filepath_creater import join_filepath


class SmkConfigWriter:
    @classmethod
    def write(cls, unique_id, flow_config):
        ConfigWriter.write(
            dirname=DIRPATH.ROOT_DIR,
            filename=DIRPATH.SNAKEMAKE_CONFIG_YML,
            config=flow_config,
        )

        ConfigWriter.write(
            dirname=join_filepath([DIRPATH.OUTPUT_DIR, unique_id]),
            filename=DIRPATH.SNAKEMAKE_CONFIG_YML,
            config=flow_config,
        )
