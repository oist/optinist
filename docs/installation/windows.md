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

**CAUTION**: For WSL2, we confirmed them on [Ubuntu 20.04](https://www.microsoft.com/ja-jp/p/ubuntu-2004-lts/9n6svws3rx71).


## For Windows (PowerShell)

### 1. Make backend environment

#### Install Tools

(windows-install-visual-studio-build-tools)=

##### Install Visual Studio Build Tools

- For install CaImAn, you need to install Visual Studio Build Tools.
- Download `Build Tools for Visual Studio 2022` from https://visualstudio.microsoft.com/ja/downloads/
- In insteraller, select `Desktop Application for C++`

(windows-install-anaconda)=

##### Install Anaconda

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

<!--
#### Install mamba

We use snakemake library, and it requires mamba.
On the Anaconda PowerShell Prompt(anaconda3),
```
conda install -n base -c conda-forge mamba
```
-->

#### Install library

```
pip install optinist
```

<!--
On the Anaconda PowerShell Prompt(anaconda3),
```
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
```
-->

#### Set saving directory

Optinist default saving directory is `C:\tmp\optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```PowerShell
$ENV:OPTINIST_DIR="your_saving_dir"
# for example
# $ENV:OPTINIST_DIR="C:\optinist_data"
```

### 2. Run backend

On the Anaconda PowerShell Prompt(anaconda3),
```
run_optinist
```
- `run_optinist` log is as blow:
```
(optinist) PS C:\optinist\backend> run_optinist
INFO:     Will watch for changes in these directories: ['C:\\optinist\\backend']
INFO:     Uvicorn running on [1mhttp://0.0.0.0:8000[0m (Press CTRL+C to quit)
INFO:     Started reloader process [[36m[1m16312[0m] using [36m[1mstatreload[0m
INFO:     Started server process [[36m34084[0m]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```
- Launch browser, and go to http://localhost:8000

It opens correctly!

Done!

## For Windows (WSL2)

### 1. Make backend environment

#### Install Tools

##### Install gcc, g++

- For install CaImAn, you need to install gcc and g++.
```
sudo apt update
sudo apt install gcc g++
```

##### Install Anaconda

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

<!--
#### Install mamba

We use snakemake library, and it requires mamba.
```
conda install -n base -c conda-forge mamba
```
-->

#### Install library

```
pip install optinist
```

#### Set saving directory

Optinist default saving directory is `/tmp/optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```bash
export OPTINIST_DIR="your_saving_dir"
```

<!--
## 2. Create virtualenv

Under maintenance...
-->
<!--
In snakemake used by optinist, a virtual environment is created and executed for each function.
The procedure for first creating a virtual environment for processing suite2p, caiman, pca, etc. is described in the following link.

*It is possible to run snakemake without creating a virtual environment in advance, but it is recommended to create a virtual environment in advance because of the higher possibility of errors during execution.

Follow this [link](create_virtualenv.md).
-->

### 2. Run backend

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
- Launch browser, and go to http://localhost:8000

It opens correctly!

Done!
