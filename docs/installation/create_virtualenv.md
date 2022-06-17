Run workflow from an existing conda environment
=================

snakemake can be executed using an existing virtual environment.
The following is the procedure to execute a function in a virtual environment created in advance.

* [create a virtual environment](#create-a-virtual-environment)
  * [create suite2p environment](#create-suite2p-environment)
  * [create optinist postprocessing(PCA, ETA, etc.) environment](#create-optinist-postprocessingpca-eta-etc-environment)
  * [create caiman environment](#create-caiman-environment)
* [FAQ](#faq)


## create a virtual environment
Create a virtual environment for suite2p. (Make sure you are on the optinist root directory.)

### create suite2p environment
```
conda env create --prefix ./conda/envs/suite2p -f ./conda/yaml/suite2p_env.yaml --force
```

### create optinist postprocessing(PCA, ETA, etc.) environment
```
conda env create --prefix ./conda/envs/optinist -f ./conda/yaml/optinist_env.yaml --force
```

### create caiman environment
````
conda env create --prefix ./conda/envs/caiman -f ./conda/yaml/caiman_env.yaml --force
```

For M1 mac, re-install tensorflow.
```
conda activate ./conda/envs/caiman && pip install ./conda/yaml/tensorflow-2.4.1-py3-none-any.whl --no-dependencies --force-reinstall
````

## FAQ

A virtual environment will be created in `optinist/conda/envs`.

*If you have errors here, it would be helpful if you could check conda's install errors or ask questions on issue, etc.
