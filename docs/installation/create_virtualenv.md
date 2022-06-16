# 既存のconda環境からworkflowを実行する

snakemakeでは、既にある仮想環境を用いた実行が可能である。
以下ではあらかじめ作った仮想環境内で関数を実行する手順を説明する。

## 1. 仮想環境を作成する
以下ではsuite2pの仮想環境を作成する。（optinist上にいることを確認）
```
conda env create --prefix ./conda/envs/caiman -f ./conda/yaml/caiman_env.yaml --force
conda env create --prefix ./conda/envs/suite2p -f ./conda/yaml/suite2p_env.yaml --force
conda env create --prefix ./conda/envs/optinist -f ./conda/yaml/optinist_env.yaml --force
```

M1 macの場合は、tensorflowを入れ直す。
```
conda activate ./conda/envs/caiman && pip install ./conda/yaml/tensorflow-2.4.1-py3-none-any.whl --no-dependencies --force-reinstall
```

仮想環境は、`optinist/conda/envs`に作成される。

*ここでエラーする場合には、condaのinstallエラーを調べてもらうか、issuなどで質問してもらえると助かる。


## 2. snakemakeのinstall
snakemakeを以下のコマンドでインストールする。（20220年6月15日現在バグがあるため、修正したrepositoryからinstallする必要がある。）

```
pip uninstall snakemake
```

```
pip install git+https://github.com/ShogoAkiyama/snakemake@main#egg=snakemake
```


3. 実行
GUIを立ち上げ、suite2pを実行する。
