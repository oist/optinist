import json

import pyrebase

from studio.app.dir_path import DIRPATH

pyrebase_app = pyrebase.initialize_app(json.load(open(DIRPATH.FIREBASE_CONFIG_PATH)))
