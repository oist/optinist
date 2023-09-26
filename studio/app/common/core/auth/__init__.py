import json
import logging

import pyrebase

from studio.app.common.core.mode import MODE
from studio.app.dir_path import DIRPATH

try:
    pyrebase_app = pyrebase.initialize_app(
        json.load(open(DIRPATH.FIREBASE_CONFIG_PATH))
    )
except FileNotFoundError as e:
    if MODE.IS_STANDALONE:
        pyrebase_app = None
    else:
        logging.getLogger().error("Firebase config file not found.")
        raise e
except Exception as e:
    logging.getLogger().error(e)
    raise e
