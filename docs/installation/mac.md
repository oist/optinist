Mac
=================

```{contents}
:depth: 4
```

## Installation

We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

## 1. Make backend environment

### Install Tools

(mac-install-anaconda)=

#### Install Anaconda

- Download and install the package:
  - [Anaconda Archive](https://repo.anaconda.com/archive/)
    - Download the latest version: `Anaconda3-*.MacOSX-x86_64.pkg`
      - *The latest version of the module is fine.*

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

   3. Download and install Miniconda:

        .. code-block:: bash

           curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

           bash Miniconda3-latest-MacOSX-x86_64.sh

           /Users/MYUSERNAME/miniconda3/bin/conda init

           conda activate

   Now continue creating the optinist environment using conda
```

### Create anaconda environment

```
conda create -n optinist python=3.8
conda activate optinist
```


### Install library

```
pip install optinist
```

### Set saving directory

Optinist default saving directory is `/tmp/studio`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```
export OPTINIST_DIR="your_saving_dir"
```

## 2. Run backend

```
run_optinist
```
- `run_optinist` log is as blow:
```
$ run_optinist
INFO:   Will watch for changes in these directories: [‘/Users/oist/optinist/backend’]
INFO:   Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:   Started reloader process [5811] using statreload
INFO:   Started server process [5820]
INFO:   Waiting for application startup.
INFO:   Application startup complete.
```
- Launch browser, and go to http://localhost:8000

It opens correctly!

Done!
