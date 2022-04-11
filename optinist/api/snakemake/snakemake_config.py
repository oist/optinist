from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath
from optinist.api.config.config_writer import ConfigWriter


def snakemake_config_writer(unique_id, flow_config):
    ConfigWriter.write(
        dirname=DIRPATH.ROOT_DIR,
        filename=DIRPATH.SNAKEMAKE_CONFIG_YML,
        config=flow_config
    )

    ConfigWriter.write(
        dirname=join_filepath([DIRPATH.BASE_DIR, unique_id]),
        filename=DIRPATH.SNAKEMAKE_CONFIG_YML,
        config=flow_config
    )
