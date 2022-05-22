How to add original algorithm

## Algorithm list directory

[current algorithm list](https://github.com/oist/optinist/tree/develop/optinist/wrappers) status
- caiman
- suite2p
- optinist
    - basic neural analysis
    - dimension reduction
    - neural population analysis
    - neural decoding

It hierarcy structure is writen in `__init__.py` file.
So, if you add new algorithm, regist in dictionary of `__init__.py` file.


## Add new algorithm

### 1. Description algorithm

#### 1.1 import
Import dataclass which use for input, output datatype.
And NWB Dataset to save variable as nwb format.
```python
from optinist.api.dataclass.dataclass import *
from optinist.api.nwb.nwb import NWBDATASET
```

#### 1.2 Argument Datatype
Describe function receving datatype.
Exlain `correlation` function as example.  
- Correlation function get `Fluorescence` and `iscell` as input.（Each datatypes are `FluoData`, `IscellData`).
- Iscell is whether cells are 0 or 1.
- If default value is None, it running not to connect in GUI.
```python
def correlation(
        neural_data: FluoData,
        iscell: IscellData=None,
        params: dict=None
    ):
    neural_data = neural_data.data
    iscell = iscell.data
```

#### 1.3 Return DataType
- Function return as dictionary format.
- correlation function output HeatMapData datatype. so output wrap as heatmap datatype. 
- Describe `->{'corr': Correlation}` and argument handle.

```python
def correlation(
        ・・・・・・
    ) -> {'corr': HeatMapData}:
　　　　　　 　・・・・・・
    info = {
        'corr': HeatMapData(corr, file_name='corr'),
        'nwbfile': nwbfile,
    }
    return info
```

#### 1.4 NWB Register
- Postprocess function register `corr`.

```python
def correlation(
        ・・・・・・
    ) -> {'corr': HeatMapData}:
　　　　　　 　・・・・・・
    # NWB
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        'corr': corr,
    }
```


### 1.5 Final result correlation
Final result

```python
from optinist.api.dataclass.dataclass import *
from optinist.api.nwb.nwb import NWBDATASET

def correlation(
        neural_data: FluoData,
        iscell: IscellData=None,
        params: dict=None
    ) -> dict():

    neural_data = neural_data.data

    # data shold be time x component matrix
    if params['transpose']:
        X = neural_data.transpose()
    else:
        X = neural_data

    if iscell is not None:
        iscell = iscell.data
        ind  = np.where(iscell > 0)[0]
        X = X[ind, :]

    num_cell = X.shape[0]

    # calculate correlation
    corr = np.corrcoef(X)
    for i in range(num_cell):
        corr[i, i] = np.nan

    # NWB
    nwbfile = {}
    nwbfile[NWBDATASET.POSTPROCESS] = {
        'corr': corr,
    }

    info = {
        'corr': HeatMapData(corr, file_name='corr'),
        'nwbfile': nwbfile,
    }

    return info
```


## 2. Algorithm register on GUI
Need to Register algorithm [list file](https://github.com/oist/optinist/blob/develop/optinist/wrappers/optinist_wrapper/neural_population_analysis/__init__.py).
- function name as key, function as value.

```python
from .correlation import correlation

original_wrapper_dict = {
    'correlation': correlation
}
```

## 3. Parameter register
- Algorithm parameters are saved in [config directory](https://github.com/oist/optinist/tree/main/optinist/config).
- File is written in yaml format.
- It is related registed function name, so correlation parameter file is named `correlation.yaml` file.

ex)
lda.yaml
```
# whether standardize the data or not
standard_x_mean: True
standard_x_std: True

transpose_x: False
transpose_y: False

target_index: 1

# n_splits = int, default=5
# Number of folds. Must be at least 2.
CV:
  n_splits: 5
  shuffle: False

LDA:
  solver: 'svd'
  shrinkage:
  priors:
  n_components: 2
  store_covariance: False
  tol: 0.0001
  covariance_estimator:
```


## 4. Snakemake register
To use function in workflow, register in snakefile too.
https://github.com/oist/optinist/blob/develop/optinist/rules/smk/optinist/neural_population_analysis/correlation.py

- File contents are almost same as other file, copy and change file name.
- conda environment `conda:`.

```
from cui_api.const import ROOT_DIR

name = "correlation"

rule:
    input:
        smk_input(config, name)
    output:
        smk_output(config, name)
    conda:
        f'{DIRPATH.ROOT_DIR}/rules/envs/optinist_env.yaml'
    params:
        name = name
    conda:
        f'{DIRPATH.ROOT_DIR}/rules/envs/optinist_env.yaml'
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/func.py'
```
