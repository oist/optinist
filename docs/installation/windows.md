Windows
=================

```{contents}
:depth: 4
```

## Installation

We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

```{eval-rst}
.. caution::
   For WSL2, we confirmed them on Ubuntu 20.04/22.04.
```

## For Windows (PowerShell)

### 1. Make Backend Environment

#### Install Tools

(windows-install-visual-studio-build-tools)=

##### Install Visual Studio Build Tools

- For install CaImAn, you need to install Visual Studio Build Tools.
- Download `Build Tools for Visual Studio 2022` from https://visualstudio.microsoft.com/ja/downloads/
- In insteraller, select `Desktop Application for C++`

(windows-install-anaconda)=

##### Install Anaconda

- Install Anaconda for Windows
  - https://www.anaconda.com/products/individual

#### Create Conda Environment

On the Anaconda PowerShell Prompt(anaconda3),

```bash
conda create -n optinist python=3.8
conda activate optinist
```

#### Install Library

```bash
pip install optinist
```

#### Set Saving Directory

Optinist default saving directory is `C:\tmp\optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.

```PowerShell
$ENV:OPTINIST_DIR="your_saving_dir"
# for example
# $ENV:OPTINIST_DIR="C:\optinist_data"
```

### 2. Run Backend

On the Anaconda PowerShell Prompt(anaconda3),

```bash
run_optinist
```

- `run_optinist` log is as blow:

```bash
(optinist) PS C:\optinist\backend> run_optinist
INFO:     Will watch for changes in these directories: ['C:\\optinist\\backend']
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:     Started reloader process 16312 using statreload
INFO:     Started server process 34084
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```

- Launch browser, and go to http://localhost:8000

Done!

## For Windows (WSL2)

### 1. Make Backend Environment

#### Install Tools

##### Install gcc, g++

- For install CaImAn, you need to install gcc and g++.

```bash
sudo apt update
sudo apt install gcc g++
```

##### Install Anaconda

```bash
# *The latest version of the module is ok.
ANACONDA_VERSION=2022.10
wget https://repo.anaconda.com/archive/Anaconda3-${ANACONDA_VERSION}-Linux-x86_64.sh
bash Anaconda3-${ANACONDA_VERSION}-Linux-x86_64.sh
```

#### Create Conda Environment

```bash
conda create -n optinist python=3.8
conda activate optinist
```

#### Install Library

```bash
pip install optinist
```

#### Set Saving Directory

Optinist default saving directory is `/tmp/studio`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.

```bash
export OPTINIST_DIR="your_saving_dir"
```

### 2. Run Backend

```bash
run_optinist
```

- `run_optinist` log is as blow:

```bash
$ run_optinist
INFO:     Will watch for changes in these directories: ['/home/oist/optinist/backend']
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:     Started reloader process [3268] using statreload
INFO:     Started server process [3311]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```

- Launch browser, and go to http://localhost:8000

Done!
