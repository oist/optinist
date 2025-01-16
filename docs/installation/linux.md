Linux
=================

```{contents}
:depth: 4
```

## Installation

We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

```{eval-rst}
.. caution::
    We confirmed them on Ubuntu 18.04/20.04/22.04.
```

## 1. Make Backend Environment

### Install Tools

#### Install gcc, g++

- For install CaImAn, you need to install gcc and g++.

```bash
sudo apt install gcc g++
```

#### Install Anaconda

```bash
# *The latest version of the module is ok.
ANACONDA_VERSION=2022.10
wget https://repo.anaconda.com/archive/Anaconda3-${ANACONDA_VERSION}-Linux-x86_64.sh
bash Anaconda3-${ANACONDA_VERSION}-Linux-x86_64.sh
```

### Create Conda Environment

```bash
conda create -n optinist python=3.8
conda activate optinist
```

### Install Library

```bash
pip install optinist
```

### Set Saving Directory

Optinist default saving directory is `/tmp/studio`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.

```bash
export OPTINIST_DIR="your_saving_dir"
```

## 2. Run Backend

```bash
run_optinist
```

- `run_optinist` log is as blow:

```bash
$ run_optinist
INFO:     Will watch for changes in these directories: ['/home/oist/optinist/backend']
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:     Started reloader process [6520] using statreload
INFO:     Started server process [6557]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```

- Launch browser, and go to http://localhost:8000

Done!
