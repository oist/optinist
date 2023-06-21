Linux
=================

* [Installation](#installation)
* [1. Make backend environment](#2-make-backend-environment)
   * [Install gcc, g++](#install-gcc-g)
   * [Install Anaconda](#install-anaconda)
   * [Create anaconda environment](#create-anaconda-environment)
   * [Install mamba](#install-mamba)
   * [Install library](#install-library)
   * [Set saving directory](#set-saving-directory)
* [2. Create virtualenv](#3-create-virtualenv)
* [3. Run backend](#4-run-backend)
   * [Launch browser.  <a href="http://localhost:8000" rel="nofollow">http://localhost:8000</a>](#launch-browser--httplocalhost8000)

## Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

**CAUTION**: We confirmed them on Ubuntu 18.04 or 20.04.

## 1. Make backend environment
### Install gcc, g++
- For install CaImAn, you need to install gcc and g++.
```
sudo apt install gcc g++
```
### Install Anaconda
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.11-Linux-x86_64.sh
```
### Create anaconda environment
```
conda create -n optinist python=3.8
conda activate optinist
```

```
conda config --set channel_priority strict
```

### Install mamba
We use snakemake library, and it requires mamba.
```
conda install -n base -c conda-forge mamba
```

### Install library

```bash
pip install optinist
```
### Set saving directory
Optinist default saving directory is `/tmp/optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```bash
export OPTINIST_DIR="your_saving_dir"
```

## 2. Create virtualenv
Under maintenance...
<!-- In snakemake used by optinist, a virtual environment is created and executed for each function.
The procedure for first creating a virtual environment for processing suite2p, caiman, pca, etc. is described in the following link.

*It is possible to run snakemake without creating a virtual environment in advance, but it is recommended to create a virtual environment in advance because of the higher possibility of errors during execution.

Follow this [link](create_virtualenv.md). -->

## 3. Run backend
```
run_optinist
```
- `run_optinist` log is as blow:
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
