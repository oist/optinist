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

- Download and install package.
  - https://repo.anaconda.com/archive/
    - Anaconda3-\*.\*-MacOSX-x86_64.pkg
      - *The latest version of the module is ok.

```{eval-rst}
.. caution::
   Even if you're using arm64 (Apple Sillicon, M1, M2...) architecture's Mac, x86_64 version is required.
   Some modules cannot be installed by conda install or pip install in arm64 version.
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
