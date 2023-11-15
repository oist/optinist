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

## Open Terminal

<img width="674" alt="Terminal photo" src="https://github.com/josh-oloro/optinist/assets/72240796/651fcb9d-8ade-4be8-bda2-a53ff83ef083">

## 1. Make backend environment

### Install Tools

(mac-install-anaconda)=

#### Install Anaconda

- Download and install package.
  - https://repo.anaconda.com/archive/
    - Anaconda3-\*.\*-MacOSX-x86_64.pkg
      - *The latest version of the module is ok.

- **CAUTION**
  - On mac, you must use the x86_64 version.
    - The arm64 version has some modules that cannot be installed by conda install or pip install.
      
  - For new Mac Silicon users, check if Anaconda and Python package are in Intel using Activity Monitor (while the program is running).
    <img width="924" alt="Activity monitor check" src="https://github.com/josh-oloro/optinist/assets/72240796/ac0ea102-8478-4d97-ab1e-a429cd11052e">


### Create anaconda environment

```
conda create -n optinist python=3.8
conda activate optinist
```

<!-- ```
conda config --set channel_priority strict
``` -->

### Install library

```
pip install optinist
```

### Set saving directory

Optinist default saving directory is `/tmp/studio`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```
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

**It opens correctly!**

<img width="1726" alt="Sample photo of OptiNiSt" src="https://github.com/josh-oloro/optinist/assets/72240796/fc686064-b6b7-491d-b8b8-b80ab43af31d">


Done!

## When restarting

Open Terminal

```
run_optinist
```
<img width="906" alt="Restart prompt OptiNiSt" src="https://github.com/josh-oloro/optinist/assets/72240796/b5a1d1ee-7921-416c-afd0-4f5f087cf107">



