from typing import Dict
import logging

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import create_directory, join_filepath


def get_logger(unique_id: str) -> logging.Logger:
    output_dirpath = join_filepath([
        DIRPATH.OUTPUT_DIR,
        unique_id,
    ])
    create_directory(output_dirpath)

    logging.basicConfig(
        filename=f"{output_dirpath}/error.log",
        level=logging.DEBUG,
        format="%(asctime)s : %(levelname)s - %(filename)s - %(message)s",
    )
    logger = logging.getLogger(unique_id)

    return logger


class Logger:
    def __init__(self, unique_id):
        self.unique_id = unique_id
        self.logger: logging.Logger = get_logger(unique_id)

    def smk_logger(self, msg: Dict[str, str] = None):
        """
        msg:
            level:
                "job_info": jobの始まりを通知。inputやoutputがある
                "job_finished": jobの終了を通知。
        """

        # エラーした
        if "exception" in msg:
            print("Error: ", msg)
            self.logger.debug("Error")

        if "level" in msg and "debug" in msg["level"]:
            level = msg["level"]
            if "debug" in level and "msg" in msg:
                self.logger.debug("error message: ", msg)
