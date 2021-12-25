# optinist

## install

### CaImAnのインストール
あらかじめpipのupgradeとCythonのインストールが必要
```
pip install --upgrade pip
pip install -r requierements.txt
```

[CaImAn](https://github.com/flatironinstitute/CaImAn)のリポジトリの指示に従いCaImAnをインストールする。condaとcloneでのインストール方法がある。

CaImAnをインストールディレクトリからリポジトリをクローンしCaImAnをインストール。
```
cd CaImAn
pip install .
python setup.py build_ext -i
python caimanmanager.py install --inplace
```
