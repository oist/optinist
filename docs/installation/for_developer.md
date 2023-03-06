For developer
=================

- [Installation](#installation)
  - [1. Make backend environment](#2-make-backend-environment)
    - [Clone repository](#clone-repository)
    - [Create anaconda environment](#create-anaconda-environment)
    - [Install requirements](#install-requirements)
    - [Set saving directory](#set-saving-directory)
  - [2. Create virtualenv](#3-create-virtualenv)
  - [3. Run backend](#4-run-backend)
    - [Launch browser.  http://localhost:8000](#launch-browser--httplocalhost8000)

## Installation
We introduce how to install optinist for developer.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

## 1. Make backend environment
### Clone repository
```
git clone https://github.com/oist/optinist.git
cd ./optinist
```
### Install [Anaconda](https://www.anaconda.com/products/individual)
### Create anaconda environment
```
conda create -n optinist python=3.8
conda activate optinist
```

### Install requirements

```bash
pip install -r requirements.txt
```
### Set saving directory
Optinist default saving directory is `/tmp/optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```bash
export OPTINIST_DIR="your_saving_dir"
```

## 2. Create virtualenv
Under maintenance...
## 3. Run backend
```
python main.py
```
- `python main.py` log is as blow:
```
$ run_optinist
INFO:     Will watch for changes in these directories: ['/home/oist/optinist/backend']
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:     Started reloader process [6520] using statreload
INFO:     Started server process [6557]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```
### Launch browser.  http://localhost:8000
It opens correctly!

Done!
