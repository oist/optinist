# Mac

```{contents}
:depth: 4
```

## Installation

We introduce how to install optinist.
Please follow instructions below.

<br />

## 1. Make Backend Environment

### Install Tools

(mac-install-conda)=

#### Install Conda Tool

- Install [Miniforge](https://github.com/conda-forge/miniforge)

```bash
# install (*latest version is fine)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-x86_64.sh"
bash Miniforge3-$(uname)-x86_64.sh

# initial setting (re-login required)
$HOME/miniforge3/bin/conda init
logout

# setting configs
conda config --set channel_priority flexible
```

```{eval-rst}
.. caution::
   Even if you're using arm64 (Apple Silicon, M1, M2...) architecture's Mac, the x86_64 version is required.
   Some modules cannot be installed by conda install or pip install in the arm64 version.
   Installing the x86_64 version of conda can be done using `rosetta`.

   1. Install Rosetta using the terminal:

        .. code-block:: bash

          /usr/sbin/softwareupdate --install-rosetta --agree-to-license

   2. Open a Terminal Session in Rosetta:
      - Open your existing Terminal (which is running natively on ARM).
      - Start a new Terminal session that emulates the x86_64 architecture using the following command:

        .. code-block:: bash

           arch -x86_64 /usr/bin/env bash

   Now continue creating the optinist environment using conda
```

```{eval-rst}
.. note::
   Some packages no longer publish prebuilt wheels for macOS x86_64
   (e.g. ``cryptography`` 49 and later), so pip builds them from source,
   which can produce a broken binary. OptiNiSt pins such packages.
   If a workflow fails with ``Unknown OpenSSL error`` or the process crashes,
   update OptiNiSt and recreate the conda environment.
```

### Create Conda Environment

```bash
conda create -n optinist python=3.9
conda activate optinist
```

### Install Library

```bash
pip install optinist
```

### Set Saving Directory

Optinist default saving directory is `/tmp/studio`. If you reboot your PC, this repository content is deleted. And setting the saving directory in environment path.

```bash
export OPTINIST_DIR="your_saving_dir"
```

## 2. Run Backend

```bash
run_optinist
```

- `run_optinist` log is as below:

```bash
$ run_optinist
INFO:   Will watch for changes in these directories: [‘/Users/oist/optinist’]
INFO:   Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:   Started reloader process [xxxx] using statreload
INFO:   Started server process [xxxx]
INFO:   Waiting for application startup.
INFO:   Application startup complete.
```

- Launch browser, and go to `http://localhost:8000`

Done!
