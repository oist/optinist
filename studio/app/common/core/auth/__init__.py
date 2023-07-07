import json

import pyrebase

from studio.app.common.core.auth.auth_config import AUTH_CONFIG
from studio.app.dir_path import DIRPATH

pyrebase_app = (
    pyrebase.initialize_app(json.load(open(DIRPATH.FIREBASE_CONFIG_PATH)))
    if not AUTH_CONFIG.IS_STANDALONE
    else None
)
