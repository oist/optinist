from studio.app.common.core.utils.config_handler import ConfigWriter
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.dir_path import DIRPATH


class SmkConfigWriter:
    @classmethod
    def write(cls, workspace_id, unique_id, flow_config):
        ConfigWriter.write(
            dirname=join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, unique_id]),
            filename=DIRPATH.SNAKEMAKE_CONFIG_YML,
            config=flow_config,
        )
