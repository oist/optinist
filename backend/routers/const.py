BASE_DIR = '/tmp/optinist'
OPTINIST_DIR = '../optinist'

import os
if not os.path.exists(BASE_DIR):
	os.makedirs(BASE_DIR, exist_ok=True)
