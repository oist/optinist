from firebase_admin import credentials, initialize_app

from studio.app.common.core.auth.auth_config import AUTH_CONFIG
from studio.app.dir_path import DIRPATH

if not AUTH_CONFIG.IS_STANDALONE:
    initialize_app(credentials.Certificate(DIRPATH.FIREBASE_PRIVATE_PATH))
