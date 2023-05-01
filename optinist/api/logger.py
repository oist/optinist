import logging
import os
from typing import Dict

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import create_directory, join_filepath


def get_logger(unique_id: str) -> logging.Logger:
    output_dirpath = join_filepath(
        [
            DIRPATH.OUTPUT_DIR,
            unique_id,
        ]
    )
    create_directory(output_dirpath)

    filepath = f"{output_dirpath}/error.log"
    if os.path.exists(filepath):
        try:
            os.remove(filepath)
        except Exception as e:
            print("[Exception][Logger]", e)

    logger = logging.getLogger(unique_id)

    # FileHandlerの設定
    fh = logging.FileHandler(filepath)
    fh.setLevel(logging.ERROR)
    fmt = logging.Formatter("%(asctime)s : %(levelname)s - %(filename)s - %(message)s")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

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
        # pass
        # # エラーした
        # if "exception" in msg:
        #     self.logger.error(msg)

        if "level" in msg and "debug" in msg["level"]:
            level = msg["level"]
            if "debug" in level and "msg" in msg:
                if "Traceback" in msg["msg"]:
                    self.logger.error(msg)
