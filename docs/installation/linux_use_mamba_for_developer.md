Linux for Developer (Use Mamba)
=================

```{contents}
:depth: 4
```

## Installation

We introduce how to install optinist for developer.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

## Unix-like platforms (Linux, Windows WSL, Mac)

### 1. Make backend environment

#### Install Tools

##### Install Mamba

- https://mamba.readthedocs.io/en/latest/installation.html
  - Use Mambaforge
    - https://github.com/conda-forge/miniforge#mambaforge
  - Install
    ```
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
    bash Mambaforge-$(uname)-$(uname -m).sh

    # restart shell, and check mamba installed. (show version)
    mamba --version
    ```

#### Clone repository

```
git clone https://github.com/oist/optinist.git
cd ./optinist
```

#### Create mamba(anaconda) environment

```
mamba create -c conda-forge -c bioconda -n optinist_dev python=3.8 snakemake
conda activate optinist_dev
```

- Note:
  - If mamba is already installed, snakamake prefers to use mamba over conda.


#### Install requirements

```
pip install -r requirements.txt
```

#### Set saving directory

Optinist default saving directory is `/tmp/optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```bash
export OPTINIST_DIR="your_saving_dir"
```

<!--
### 2. Create virtualenv

Under maintenance...
-->

### 2. Run backend

```
python main.py
```
- `python main.py` log is as blow:
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
