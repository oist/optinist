BASE_DIR = '/tmp/optinist'

import os
if not os.path.exists(BASE_DIR):
    os.makedirs(BASE_DIR, exist_ok=True)
