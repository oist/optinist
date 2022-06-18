CUI execution
=================
This section describes how to run a workflow created in GUI on a cluster or in CUI.

* [1. Config file settings](#1-config-file-settings)
* [2. Set environment variables for save paths](#2-set-environment-variables-for-save-paths)
* [3. Input File Settings](#3-input-file-settings)
* [4. Execution](#4-execution)
* [5.Output result](#5output-result)

## 1. Config file settings
Place the config file required by Snakemake at the appropriate location.
The config file can be downloaded from **Record** on the GUI.

## 2. Set environment variables for save paths
Change environment variables. Change the environment variable as follows, because the default setting refers to the directory under `/tmp/optinist`.
```bash
export OPTINIST_DIR="your_saving_dir"
```

With environment variables set, input refers to `{OPTINIST_DIR}/input/` and results are output to `{OPTINIST_DIR}/output`.

## 3. Input File Settings
Input files are stored under `/{OPTINIST_DIR}/input/`.
For example, if there is `mouse2p_2_donotouse.tiff`, it is stored as `/{OPTINIST_DIR}/input/mouse2p_2_donotouse.tiff`.

## 4. Execution
It can be executed in CUI by running `run_cluster.py`.

The parameter arguments are as follows.
- config(string): path of config.yaml file set in step 1
- cores(int): Specifies the number of CPU cores. (defalt: 2, cores >= 2)
- forceall(bool): Whether to overwrite existing results or not. (default: false)
- use_conda(bool): Whether to use the conda virtual environment or not. If not, it is necessary to have caiman, suite2p, etc. INSTALL the current environment into the execution environment.


The command executes the following
```bash
python run_cluster.py --config="{config.yaml file path}"
```

## 5.Output result
The results are stored in `{OPTINIST_DIR}/output/{unique_id}`.
If you want to visualize the results, you can move this directory to local, and you can check the results in GUI.
