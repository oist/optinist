Mac
=================

* [Installation](#installation)
* [0. GitHub SSH access settings](#0-github-ssh-access-settings)
* [1. Clone optinist repository](#1-clone-optinist-repository)
* [2. Make backend environment](#2-make-backend-environment)
   * [Install Anaconda](#install-anaconda)
   * [Create anaconda environment](#create-anaconda-environment)
   * [Install library](#install-library)
   * [Set saving directory](#set-saving-directory)
* [3. Create virtualenv](#3-create-virtualenv)
* [4. Run backend](#4-run-backend)
   * [Launch browser.  <a href="http://localhost:8000" rel="nofollow">http://localhost:8000</a>](#launch-browser--httplocalhost8000)
* [FAQ](#faq)

## Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

**CAUTION for M1 Mac User**
We use tensorflow in caiman code. We know that M1 mac doesn't install tensorflow easily, so if there is a problem, skip install caiman. (Release in progress…)

## 0. GitHub SSH access settings
**You only need to do the following once.**

Follow this [link](settings.md).

## 1. Clone optinist repository

First, you get optinist code from github repository.
```
cd "your working repository"
git clone git@github.com:oist/optinist.git
```
<br />

## 2. Make backend environment

### Install Anaconda
Download https://repo.anaconda.com/archive/Anaconda3-2021.11-MacOSX-x86_64.pkg

or Miniconda https://docs.conda.io/en/latest/miniconda.html

Install it.

### Create anaconda environment
```
conda create -n optinist python=3.8
conda activate optinist
```

```
conda config --set channel_priority strict
```

### Install library
```bash
pip install -r requirements.txt
```

In case an error occuered when you install ``` requirements.txt ```, pip upgrade command below may solve the error.
```
pip install --upgrade pip
```

### Set saving directory
Optinist default saving directory is `/tmp/optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```bash
export OPTINIST_DIR="your_saving_dir"
```

## 3. Create virtualenv
In snakemake used by optinist, a virtual environment is created and executed for each function.
The procedure for first creating a virtual environment for processing suite2p, caiman, pca, etc. is described in the following link.

*It is possible to run snakemake without creating a virtual environment in advance, but it is recommended to create a virtual environment in advance because of the higher possibility of errors during execution.

Follow this [link](create_virtualenv.md).


## 4. Run backend
```
python main.py
```
- `python main.py` log is as blow:
```
$ python main.py
INFO:   Will watch for changes in these directories: [‘/Users/oist/optinist/backend’]
INFO:   Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:   Started reloader process [5811] using statreload
INFO:   Started server process [5820]
INFO:   Waiting for application startup.
INFO:   Application startup complete.
```
### Launch browser.  http://localhost:8000
It opens correctly!

Done!


## FAQ
If you get the warning message shown below, we recommend `rm -rf /Users/usename/opt/anaconda3/envs/optinist ` and recreate conda environment
or add `--force` option when you create optinist conda environment.
> WARNING: A directory already exists at the target location '/Users/usename/opt/anaconda3/envs/optinist' but it is not a conda environment.
