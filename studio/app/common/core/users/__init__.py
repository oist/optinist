from firebase_admin import credentials, initialize_app

from studio.app.dir_path import DIRPATH

initialize_app(credentials.Certificate(DIRPATH.FIREBASE_PRIVATE_PATH))
