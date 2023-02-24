Windows
=================

* [Installation](#installation)
* [1. Make backend environment](#2-make-backend-environment)
   * [For Windows(PowerShell) Users](#for-windowspowershell-users)
      * [Install Visutal Studio Build Tools](#install-visutal-studio-build-tools)
      * [Install Anaconda](#install-anaconda)
      * [Create anaconda environment](#create-anaconda-environment)
      <!-- * [Install mamba](#install-mamba) -->
      * [Install library](#install-library)
      * [Set saving directory](#set-saving-directory)
      * [Run backend](#run-backend)
      * [Launch browser.  <a href="http://localhost:8000" rel="nofollow">http://localhost:8000</a>](#launch-browser--httplocalhost8000)
   * [For Windows (WSL2) Users](#for-windows-wsl2-users)
      * [Install gcc, g++](#install-gcc-g)
      * [Install Anaconda](#install-anaconda-1)
      * [Create anaconda environment](#create-anaconda-environment-1)
      <!-- * [Install mamba](#install-mamba-1) -->
      * [Install library](#install-library-1)
      * [Set saving directory](#set-saving-directory-1)
* [2. Create virtualenv](#3-create-virtualenv)
* [3. Run backend](#4-run-backend)
   * [Launch browser.  <a href="http://localhost:8000" rel="nofollow">http://localhost:8000</a>](#launch-browser--httplocalhost8000-1)

## Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

**CAUTION**: For WSL2, we confirmed them on [Ubuntu 20.04](https://www.microsoft.com/ja-jp/p/ubuntu-2004-lts/9n6svws3rx71).


## 1. Make backend environment

### For Windows(PowerShell) Users
#### Install Visutal Studio Build Tools
- For install CaImAn, you need to install Visual Studio Build Tools.
- Download `Build Tools for Visual Studio 2022` from https://visualstudio.microsoft.com/ja/downloads/
- In insteraller, select `Desktop Application for C++`

#### Install Anaconda
Install [Anaconda for Windows](https://www.anaconda.com/products/individual)
#### Create anaconda environment
On the Anaconda PowerShell Prompt(anaconda3),
```
conda create -n optinist python=3.8
conda activate optinist
```

<!-- ```
conda config --set channel_priority strict
``` -->

<!-- #### Install mamba
We use snakemake library, and it requires mamba.
On the Anaconda PowerShell Prompt(anaconda3),
```
conda install -n base -c conda-forge mamba
``` -->
<!-- #### Install library
On the Anaconda PowerShell Prompt(anaconda3),
```bash
pip install optinist
# for suite2p
pip install "PyQt5<=5.15.1" "PyQt5-sip<=12.8.1" "pyqtgraph<=0.11.0" "pandas" "suite2p<=0.10.3" "tifffile<=v2022.3.25"
# for GLM
pip install "sklearn" "statsmodels<=0.13.1" "pynwb"
# for CaImAn
pip install "cython" "opencv-python" "matplotlib" "scikit-image==0.18.0" "scikit-learn" "ipyparallel" "holoviews" "watershed" "tensorflow"
git clone https://github.com/flatironinstitute/CaImAn.git
cd CaImAn
pip install -e .
cd ..
``` -->

#### Set saving directory
Optinist default saving directory is `C:\tmp\optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```PowerShell
$ENV:OPTINIST_DIR="your_saving_dir"
# for example
# $ENV:OPTINIST_DIR="C:\optinist_data"
```

#### Run backend
On the Anaconda PowerShell Prompt(anaconda3),
```
run_optinist
```
- `run_optinist` log is as blow:
```
(optinist) PS C:\optinist\backend> run_optinist
[32mINFO[0m:     Will watch for changes in these directories: ['C:\\optinist\\backend']
[32mINFO[0m:     Uvicorn running on [1mhttp://0.0.0.0:8000[0m (Press CTRL+C to quit)
[32mINFO[0m:     Started reloader process [[36m[1m16312[0m] using [36m[1mstatreload[0m
[32mINFO[0m:     Started server process [[36m34084[0m]
[32mINFO[0m:     Waiting for application startup.
[32mINFO[0m:     Application startup complete.
```
#### Launch browser.  http://localhost:8000
It opens correctly!

Done!

### For Windows (WSL2) Users
#### Install gcc, g++
- For install CaImAn, you need to install gcc and g++.
```
sudo apt update
sudo apt install gcc g++
```
#### Install Anaconda
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.11-Linux-x86_64.sh
```
#### Create anaconda environment
```
conda create -n optinist python=3.8
conda activate optinist
```

<!-- ```
conda config --set channel_priority strict
``` -->

<!-- #### Install mamba
We use snakemake library, and it requires mamba.
```
conda install -n base -c conda-forge mamba
``` -->
#### Install library
```bash
pip install optinist
```
#### Set saving directory
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
INFO:     Started reloader process [3268] using statreload
INFO:     Started server process [3311]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```
#### Launch browser.  http://localhost:8000
It opens correctly!

Done!
