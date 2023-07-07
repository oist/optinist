from firebase_admin import credentials, initialize_app

from studio.app.dir_path import DIRPATH

try:
    initialize_app(credentials.Certificate(DIRPATH.FIREBASE_PRIVATE_PATH))
except Exception:
    pass
