Procedure for adding a new algorithm
=================

How to add algorithms

* [1. Create algorithm function file](#1-create-algorithm-function-file)
* [2. Register your algorithm](#2-register-your-algorithm)
* [3. Describe function processing](#3-describe-function-processing)
    * [3.1 import](#31-import)
    * [3.2 Input &amp; Output handle](#32-input--output-handle)
    * [3.3 Drawing output results](#33-drawing-output-results)
    * [4. How to add parameters](#4-how-to-add-parameters)

## 1. Create algorithm function file
First, create a python file at the appropriate location in the following directory.
Here is an example of how to create a function named `new_algorithm`. The file should be created at the following location [optinist/wrappers/optinist_wrapper](https://github.com/oist/optinist/tree/main/optinist/wrappers/optinist_wrapper).

Create it under the name **new_algorith.py**.

```python:new_algorithm.py
def new_algorithm():
    return
```

[algorithm list](https://github.com/oist/optinist/tree/develop/optinist/wrappers) 

- \_\_init__.py - ①
- caiman_wrapper
    - \_\_init__.py - ②
    - caiman_mc.py
    - caiman_cnmf.py
- suite2p_wrapper
    - \_\_init__.py
    - suite2p_file_convert.py
    - suite2p_registration.py
    - suite2p_roi.py
    - suite2p_cnmf.py
- optinist_wrapper
    - \_\_init__.py
    - basic neural analysis
    - dimension reduction
    - neural population analysis
    - neural decoding
    - `new_algorithm.py` - ＊

<br />

## 2. Register your algorithm
①: To be able to use the created **new_algorithm** on the GUI, it is necessary to register it in the \_\_init__.py. For example, [optinist/wrappers/\_\_init__.py](https://github.com/oist/optinist/blob/main/optinist/wrappers/__init__.py) in ① is as follows.
These are reading one level down, \_\_init__.py.


```python
from .caiman_wrapper import caiman_wrapper_dict
from .suite2p_wrapper import suite2p_wrapper_dict
from .optinist_wrapper import optinist_wrapper_dict

wrapper_dict = {}
wrapper_dict.update(**caiman_wrapper_dict)
wrapper_dict.update(**suite2p_wrapper_dict)
wrapper_dict.update(**optinist_wrapper_dict)
```


②: In [optinist/wrappers/caiman_wrapper/\_\_init__py](https://github.com/oist/optinist/blob/main/optinist/wrappers/caiman_wrapper/__init__.py) function is defined concretely and is written as follows. It can be registered as a function by writing `function name: {'function': function name}`.

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


Actually register the **new_algorithm** function in [optinist/wrappers/optinist_wrapper/\_init__.py](https://github.com/oist/optinist/blob/main/optinist/wrappers/optinist_wrapper/__init__.py).

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


Restart the GUI and check TreeView, you can actually see the **new_algorithm**.
<p align="center">
<img width="300px" src="https://github.com/oist/optinist/blob/main/docs/_static/add_algorithm/new_algorithm.png" alt="new_algorithm" />
</p>


## 3. Describe function processing
### 3.1 import
Next, data inputs and outputs are defined.
Optinist defines several DataClasses to ensure consistency between Input and Output types. The main data types are as follows. These correspond to the color of each Node's handle.

Optinist support datatype.
- ImageData
- TimeSeriesData
- FluoData
- BehaviorData
- IscellData
- Suite2pData
- ScatterData
- BarData


### 3.2 Input & Output handle
As an example, the **new_algorithm** function takes **ImageData** and returns **FluoData**.
The `from optinist.api.dataclass.dataclass import *` statement is the file where the dataclass is defined. params is necessary because it contains parameters for this function.
```python
from optinist.api.dataclass.dataclass import *

def new_algorithm(
        image_data: ImageData,
        params: dict=None
    ) -> dict(fluo=FluoData):
    return
```

Restart the GUI and put **new_algorithm**, and you will see that the handle color has changed.
<p align="center">
<img width="200px" src="https://github.com/oist/optinist/blob/main/docs/_static/add_algorithm/input_output.png" alt="input_output" />
</p>


### 3.3 Drawing output results
- Above we described the node input and output handle, here we describe the visualization of the result.
- The output of the function is a dictionary. (Here we use the variable name **info**.)
- First, the **fluo** variable that is the return value of the **new_algorithm function** is output by Wrap with **FluoData**. The name of the key in this case must match the **fluo** of the return value when declaring the function.
- In addition, variables to be visualized are wrapped with their data types and output. In this example, **ImageData** and **HeatMap** are output.

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

Restart the GUI, connect imageNode and run it, and you will see the output as follows.

** Note: This process takes only 2~3 seconds, so if the process does not finish, there may be an error. If the error cannot be resolved, please post a message on slack or an issue.

<p align="center">
<img width="300px" src="https://github.com/oist/optinist/blob/main/docs/_static/add_algorithm/run.png" alt="run" />
</p>

<p align="center">
<img width="240px" src="https://github.com/oist/optinist/blob/main/docs/_static/add_algorithm/visualize_output.png" alt="output" />
</p>


### 4. How to add parameters
The parameters are stored under [optinist/config](https://github.com/oist/optinist/tree/main/optinist/config) with **the same file name as the function name**.
The file name is `new_algorithm.yaml` and can be registered by creating it.

