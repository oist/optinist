Linux
=================

```{contents}
:depth: 4
```

## Installation

We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

**CAUTION**: We confirmed them on Ubuntu 18.04/20.04/22.04.

## 1. Make backend environment

### Install Tools

#### Install gcc, g++

- For install CaImAn, you need to install gcc and g++.
```
sudo apt install gcc g++
```

#### Install Anaconda

```
# *The latest version of the module is ok.
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
bash Anaconda3-2022.10-Linux-x86_64.sh
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
INFO:     Will watch for changes in these directories: ['/home/oist/optinist/backend']
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:     Started reloader process [6520] using statreload
INFO:     Started server process [6557]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```
- Launch browser, and go to http://localhost:8000

It opens correctly!

Done!
