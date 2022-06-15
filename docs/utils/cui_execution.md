# CUI execution
ここでは、GUIで作成したworkflowをcluster上やCUIで実行する方法について説明する。


## 1. configファイルの設定
Snakemakeで必要となるconfigファイルを適当な位置に設置する。
configファイルはGUI上の**Record**からdownloadができる。

## 2. 保存パスの環境変数の設定
環境変数を変更する。デフォルトの設定では、`/tmp/optinist`以下を参照するため、以下のように環境変数を変更する。
```bash
export OPTINIST_DIR="your_saving_dir"
```

環境変数を設定すると、入力は`{OPTINIST_DIR}/input/`を参照し、結果は`{OPTINIST_DIR}/output`に出力される。

## 3. 入力ファイルの設定
入力ファイルを`/{OPTINIST_DIR}/input/`以下に格納する。
例えば、`mouse2p_2_donotouse.tiff`あれば、`/{OPTINIST_DIR}/input/mouse2p_2_donotouse.tiff`というように格納する。

## 4. 実行
`run_cluster.py`を実行することで、CUIでの実行が可能である。

パラメータ引数は以下の通りである。
- config(string)：step1で設定したconfig.yamlファイルのpath
- cores(int)：CPUコア数の指定。（defalt: 2, cores >= 2）
- forceall(bool)：既存の結果を上書きするかどうか。（default: False）
- use_conda(bool)：condaの仮想環境を使うかどうか。使わない場合には、現在の環境にcaimanやsuite2pなどを実行環境にinstallしてもらう必要がある。


コマンドは以下を実行する
```bash
python run_cluster.py --config="{config.yaml file path}"
```

## 5.結果
結果は、`{OPTINIST_DIR}/output/{unique_id}`に保存されている。
可視化などする場合には、このdirectoryをlocalに移動すれば、GUIで結果を確認することができる。
