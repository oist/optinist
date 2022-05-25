Add algorithm
=================

How to add original algorithm

## Create algorithm file
First, create python file in algorithm list directory.
Optinist support algorithm list is below link.

[current algorithm list](https://github.com/oist/optinist/tree/develop/optinist/wrappers) 


- `__init__.py` - ①
- caiman_wrapper
    - `__init__.py` - ②
    - caiman_mc.py
    - caiman_cnmf.py
- suite2p_wrapper
    - `__init__.py`
    - suite2p_file_convert.py
    - suite2p_registration.py
    - suite2p_roi.py
    - suite2p_cnmf.py
- optinist_wrapper
    - basic neural analysis
    - dimension reduction
    - neural population analysis
    - neural decoding
    - `new_algorithm.py` - *

*: create file and, describe sample code.

new_algorithm.py
```python
def new_algorithm():
    return
```



## Register algorithm
①: For example, base `__init__.py` of ① is written like below.
```python
from .caiman_wrapper import caiman_wrapper_dict
from .suite2p_wrapper import suite2p_wrapper_dict
from .optinist_wrapper import optinist_wrapper_dict

wrapper_dict = {}
wrapper_dict.update(**caiman_wrapper_dict)
wrapper_dict.update(**suite2p_wrapper_dict)
wrapper_dict.update(**optinist_wrapper_dict)
```

Getting each function dictionary and it's correspond to GUI algorithm tree view.


②: And `caiman_wrapper/__init__py` is written like below.
```python
from .motion_correction import caiman_mc
from .cnmf import caiman_cnmf

caiman_wrapper_dict = {
    'caiman': {
        'caiman_mc': { 
            'function': caiman_mc,
        },
        'caiman_cnmf': {
            'function': caiman_cnmf,
        },
    }
}

```

So, whole dictionary structure are
```python
wrapper_dict = {
    'caiman': {
        'caiman_mc': {
            'function': caiman_mc,
        },
        'caiman_cnmf': {
            'function': caiman_cnmf,
        },
    },
    'suite2p': {
        'suite2p_file_convert': {
            'function': suite2p_file_convert
        },
        ...
    },
    'optinist': {
        ...
    }
}
```

After creating python file, register your file in `optinist_wrapper/__init__.py` like below.
```python
from .basic_neural_analysis import basic_neural_analysis_wrapper_dict
from .dimension_reduction import dimension_reduction_wrapper_dict
from .neural_population_analysis import neural_population_analysis_wrapper_dict
from .neural_decoding import neural_decoding_wrapper_dict

# ↓↓new add↓↓
from .new_algorithm import new_algorithm
# ↑↑new add↑↑

optinist_wrapper_dict = {
    'optinist': {
        'basic_neural_analysis': basic_neural_analysis_wrapper_dict,
        'dimension_reduction': dimension_reduction_wrapper_dict ,
        'neural_population_analysis': neural_population_analysis_wrapper_dict,
        'neural_decoding': neural_decoding_wrapper_dict,
        # ↓↓new add↓↓
        'new_algorithm': {
            'function': new_algorithm
        }
        # ↑↑new add↑↑
    }
}
```


Re-launch GUI, your creating algorithm node appear in GUI tree view.
<p align="center">
<img width="300px" src="../_static/add_algorithm/new_algorithm.png" alt="new_algorithm" />
</p>


## Describe Algorithm

### 1. Input & Output

#### 1.1 import
Import dataclass which use for input, output datatype.
```python
from optinist.api.dataclass.dataclass import *
```

Optinist support datatype.
- ImageData
- TimeSeriesData
- FluoData
- BehaviorData
- IscellData
- Suite2pData
- ScatterData
- BarData

It's correspond to GUI handle color.

#### 1.2 Input & Output handle
If your algorithm get `ImageData` datatype, arugument is `ImageData`.
And return `FluoData` dataype, return is `FluoData`.
```python
def new_algorithm(
        image_data: ImageData,
        params: dict=None
    ) -> dict(fluo=FluoData):
```

Re-launch GUI, and algorithm node input & output handle change datatype.
<p align="center">
<img width="200px" src="../_static/add_algorithm/input_output.png" alt="input_output" />
</p>


#### 1.3 Visualize output result
- 上では、nodeのinputとoutputのhandleについて記述した、ここでは、結果の可視化について説明する。
- 関数の出力はdictionaryを指定する。
- まず、`new_algorithm`関数は`fluo`変数を`FluoData`でWrapして出力する。
- それ以外に、可視化したい変数については、そのデータ型でWrapし出力する。

```python
def new_algorithm(
        image_data: ImageData,
        params: dict=None
    ) -> dict(fluo=FluoData):
    import numpy as np
    info = {
        "fluo": FluoData(np.random.rand(100, 20), file_name="fluo"),
        "image": ImageData(np.random.rand(10, 100, 100), file_name="image"),
        "heatmap": HeatMapData(np.random.rand(20, 20), file_name="heatmap")
    }
    return info
```

#### 1.4 Snakemakeの登録
関数を実行するためにSnakemakeファイルを記述する。
Snakemakeファイルは、関数と同じdirectory構造で以下のように記述されている。[snakemake list](https://github.com/oist/optinist/tree/develop/optinist/rules/smk) 

したがって、ここでは`smk/optinist/new_algorithm.smk`というファイルを作成する。
中身は他のファイルをコピペし、`name`変数を`new_algorithm`にする。

```
from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk_dir import smk_input, smk_output

name = "new_algorithm"

rule:
    input:
        smk_input(config, name)
    output:
        smk_output(config, name)
    params:
        name = name
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/func.py'
```


Re-launch GUI, and algorithm node input & output handle change datatype.
<p align="center">
<img width="300px" src="../_static/add_algorithm/run.png" alt="run" />
</p>

<p align="center">
<img width="240px" src="../_static/add_algorithm/visualize_output.png" alt="output" />
</p>


** 2 ~ 3秒で終わる処理なので、処理が終わらない場合にはエラーをしている可能性がある。console画面で赤文字のエラー部分を確認して頂きたい。エラーが解決できない場合には、slackやissueに貼って貰えると解決できる場合もある。