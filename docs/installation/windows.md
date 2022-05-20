Windows
=================

* [Installation](#installation)
* [0. GitHub SSH access settings](#0-github-ssh-access-settings)
* [1. Clone optinist repository](#1-clone-optinist-repository)
* [2. Make backend environment](#2-make-backend-environment)
   * [For Windows(PowerShell) Users](#for-windowspowershell-users)
      * [Install Visutal Studio Build Tools](#install-visutal-studio-build-tools)
      * [Install Anaconda](#install-anaconda)
      * [Create anaconda environment](#create-anaconda-environment)
      * [Install mamba](#install-mamba)
      * [Install library](#install-library)
      * [Set saving directory](#set-saving-directory)
      * [Run backend](#run-backend)
      * [Launch browser.  <a href="http://localhost:8000" rel="nofollow">http://localhost:8000</a>](#launch-browser--httplocalhost8000)
   * [For Windows (WSL2) Users](#for-windows-wsl2-users)
      * [Install gcc, g++](#install-gcc-g)
      * [Install Anaconda](#install-anaconda-1)
      * [Create anaconda environment](#create-anaconda-environment-1)
      * [Install mamba](#install-mamba-1)
      * [Install library](#install-library-1)
      * [Set saving directory](#set-saving-directory-1)
      * [Run backend](#run-backend-1)
      * [Launch browser.  <a href="http://localhost:8000" rel="nofollow">http://localhost:8000</a>](#launch-browser--httplocalhost8000-1)

## Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

**CAUTION**: For WSL2, we confirmed them on [Ubuntu 20.04](https://www.microsoft.com/ja-jp/p/ubuntu-2004-lts/9n6svws3rx71).


## 0. GitHub SSH access settings
**You only need to do the following once.**

Follow this [link](settings.md).

## 1. Clone optinist repository
First, you get optinist code from github repository.
```
cd "your working repository"
git clone git@github.com:oist/optinist.git
```
Clone CaImAn repository.
```
cd optinist
git clone https://github.com/flatironinstitute/CaImAn -b v1.9.7
```
<br />

## 2. Make backend environment

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
#### Install mamba
We use snakemake library, and it requires mamba.
On the Anaconda PowerShell Prompt(anaconda3),
```
conda install -n base -c conda-forge mamba
```
#### Install library
On the Anaconda PowerShell Prompt(anaconda3),
```bash
pip install -r requirements.txt
# for suite2p
pip install PyQt5<=5.15.1 PyQt5-sip<=12.8.1 pyqtgraph<=0.11.0 pandas suite2p<=0.10.3
# for GLM
pip install sklearn statsmodels<=0.13.1 pynwb
# for CaImAn
pip install cython opencv-python matplotlib scikit-image==0.18.0 scikit-learn ipyparallel holoviews watershed tensorflow
cd CaImAn
pip install -e .
cd ..
```

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
python main.py
```
- `python main.py` log is as blow:
```
(optinist) PS C:\optinist\backend> python main.py
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
#### Install mamba
We use snakemake library, and it requires mamba.
```
conda install -n base -c conda-forge mamba
```
#### Install library
```bash
pip install -r requirements.txt
# for CaImAn
git clone https://github.com/flatironinstitute/CaImAn -b v1.9.7
pip install cython opencv-python matplotlib scikit-image==0.18.0 scikit-learn ipyparallel holoviews watershed tensorflow
cd CaImAn
pip install -e .
cd ..
```
#### Set saving directory
Optinist default saving directory is `/tmp/optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```bash
export OPTINIST_DIR="your_saving_dir"
```
## 3. Run backend
```
python main.py
```
- `python main.py` log is as blow:
```
$ python main.py
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
