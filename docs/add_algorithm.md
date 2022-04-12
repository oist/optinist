# アルゴリズム一覧のディレクトリ
アルゴリズム一覧は以下のディレクトリに記録されている。
https://github.com/oist/optinist/tree/develop/optinist/wrappers

以下ではcorrelation関数を登録する手順を例に説明する。
githubのブランチは最新版に更新しておく。
```
git pull origin develop
```

# 新規アルゴリズムの登録

## 1. アルゴリズムの記述
correlationアルゴリズムの記述手順を説明していく。

### 1.1 import文の記述
optinistでは型を定義することで、inputとoutputの整合性を取っている。
まずは、optinistで用意されている型をimportする。
```python
from wrappers.data_wrapper import *
```

### 1.2 引数の設定
関数が受け取る型を定義する。
ここでは、`correlation`関数を例に説明する。  
- correlationは変数`timeseries`、`iscell`を受け取る。（それぞれの型は`TimeSeriesData`, `IscellData`である。）  
- timeseriesデータは時系列データである。iscellはcellかどうかを0,1で定義している。
- default値をNoneとすると、GUI上でedgeのconnectionがなくても動作する。
- nwbfileとparamsは全ての関数に共通であるため、以下のようにdefault値=Noneとして記述。
- 以上により`correlation`関数の引数の定義は次のようになる。
```python
def correlation(
        timeseries: TimeSeriesData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ):
    timeseries = timeseries.data
    iscell = iscell.data
　　　　　　　 ・・・・・・・
```

### 1.3 返り値の設定
- optinistで定義されている関数の返り値は辞書型で返す。nfoという辞書型の変数を返り値とする。
- correlation関数で出力した変数をheatmapで描画したい場合には、CorrelationDataクラスとして、Wrapする。
- GUI上で引数ハンドルを作成したい場合には `->{'corr': Correlation}`と記述すれば、GUI上での返り値に加えられる。

```python
def correlation(
        timeseries: TimeSeriesData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {'corr': CorrelationData}:
　　　　　　 　・・・・・・
    info = {}
    info['corr'] = CorrelationData(
        corr,
        file_name='corr'
    )
    return info
```

### 1.4 NWB登録
- POSTPROCESSの場合、その変数の名前と変数を以下のようにnwbfileに登録すれば、nwbfileとして登録される。

```python
def correlation(
        timeseries: TimeSeriesData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {'corr': CorrelationData}:
　　　　　　 　・・・・・・
    # NWB追加
    if nwbfile is not None:
        nwbfile[NWBDATASET.POSTPROCESS] = {
            'corr': corr,
        }

    info['nwbfile'] = nwbfile

    return info
```


### 1.5 最終的なcorrelation関数
以上の手順により出来上がったcorrelation関数は以下の通りである。

```python
from wrappers.data_wrapper import *

def correlation(
        neural_data: TimeSeriesData,
        iscell: IscellData=None,
        nwbfile: NWBFile=None,
        params: dict=None
    ) -> {}:

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

    info = {}
    info['corr'] = CorrelationData(
        corr,
        file_name='corr'
    )

    # NWB追加
    if nwbfile is not None:
        nwbfile[NWBDATASET.POSTPROCESS] = {
            'corr': corr,
        }

    info['nwbfile'] = nwbfile

    return info

```


## 2. アルゴリズムの登録
作成したアルゴリズムをGUI側で使いたい場合に、登録が必要である。
- correlationの登録は以下に記述されている。
- 以下のように、関数名をkeyに、その関数をvalueにして登録すればGUI上で使用できる。
https://github.com/oist/optinist/blob/develop/optinist/wrappers/optinist_wrapper/neural_population_analysis/__init__.py

- それぞれのディレクトリに__init__ファイルがあるため、それぞれで名前をつければ、TreeView上で階層的に定義される。

```python
from .correlation import correlation

original_wrapper_dict = {
    'correlation': correlation
}
```

## 3. パラメータの登録
- アルゴリズムのパラメータはconfigディレクトリに保存されている。
- ファイルはyamlファイルで書かれている。
- 2の手順で関数名をcorrelationと登録した場合には、correlation.yamlとして作ると自動でファイルを参照する。
https://github.com/oist/optinist/tree/develop/optinist/config  
*現在の仕様ではアルゴリズム名と一対一で対応するようにしているため、ファイル名は関数名と同じでないとエラーする。

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


## 4. Snakemakeの登録
Snakemakeへの登録は上で追加した関数と同じディレクトリ構造になるようにファイルを作成する。
https://github.com/oist/optinist/blob/develop/optinist/rules/smk/optinist/neural_population_analysis/correlation.py

- ファイルの中身は下のテンプレをコピーして、`name`を関数名にする。
- conda環境を作成する場合は、`conda:`にinstallパッケージを書く。

```
from cui_api.const import ROOT_DIR

name = "correlation"

rule:
    input:
        [x["input"] for x in config["rules"].values() if x["type"] == name]
    output:
        [x["output"] for x in config["rules"].values() if x["type"] == name]
    conda:
        f'{DIRPATH.ROOT_DIR}/rules/envs/optinist_env.yaml'
    params:
        name = name
    conda:
        f'{DIRPATH.ROOT_DIR}/rules/envs/optinist_env.yaml'
    script:
        f'{DIRPATH.ROOT_DIR}/rules/scripts/main.py'
```