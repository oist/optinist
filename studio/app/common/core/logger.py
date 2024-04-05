import logging
import logging.config

import yaml

from studio.app.common.core.mode import MODE
from studio.app.dir_path import DIRPATH


class AppLogger:
    """
    Generic Application Logger
    """

    LOGGER_NAME = None  # Note: use root logger (empty name)

    @staticmethod
    def init_logger():
        """
        Note #1.
            At the time of starting to use this Logger,
            the logging initialization process has already been performed
            at the following location,
            so no explicit initialization process is required.

            - logger initialization location
              - Web App ... studio.__main_unit__
              - Batch App ... studio.app.optinist.core.expdb.batch_runner

        Note #2.
            However, only in the case of the snakemake process,
            the initialization process is required because it is a separate process.
        """

        log_config_file = (
            f"{DIRPATH.CONFIG_DIR}/logging.yaml"
            if MODE.IS_STANDALONE
            else f"{DIRPATH.CONFIG_DIR}/logging.multiuser.yaml"
        )

        with open(log_config_file) as file:
            logging.config.dictConfig(yaml.load(file.read(), yaml.FullLoader))

    @staticmethod
    def get_logger():
        logger = logging.getLogger(__class__.LOGGER_NAME)

        # If before initialization, call init
        if not logger.handlers:
            __class__.init_logger()

        return logger
