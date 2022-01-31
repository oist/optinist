# optinist
Hot to install and use optinist
- [Installation](https://github.com/oist/optinist/blob/develop/docs/installation.md)
- [How to use gui flowchart](https://github.com/oist/optinist/blob/develop/docs/gui_flowchart.md)
- [How to use gui visualize](https://github.com/oist/optinist/blob/develop/docs/gui_visualize.md)


<br>

# 1.ディレクトリ構成

```
optinist
 - frontend
 - backend
 - optinist
```

### frontend

backend ディレクトリは主にフロント側の処理を記述する。
React で記述されている。

### backend

backend ディレクトリは主にサーバー側の処理を記述する。
fastapi で記述されている。

### optinist

optinist ディレクトリは Python での CUI 処理を記述する。
NeuroScience の解析部分の一連のフローなどはこちらに記述する。


# 2. 環境構築
2022/1月時点では、frontendとbackendを別々に環境構築するのを推奨する。
frontendはdockerで構築する。
backendはメモリを多く使用するためローカル環境で構築する。

## 2.1 frontendの環境構築

### docker build
optinistディレクトリに移動し、dockerをbuildする。
```
cd optinist
docker-compose build frontend
```
### frontendをdockerで起動
dockerでfront側を起動する。
```
docker-compose up frontend
```

browserを起動。

[http://localhost:3000](http://localhost:3000)

## 2.2 backendの環境構築
手元の環境によって、使いやすいものを選択するのをおすすめする。
ここでは、anacondaやvirutualenvで作った仮想環境にインストールする方法と、pipenvで作った仮想環境での方法を提示する。

## anacondaやvirtualenvで起動
#### 仮想環境を作成
```
conda create -n optinist python=3.9.7
conda activate optinist
```

#### 必要なライブラリをインストール
```
cd optinist/backend
pip install -r requirements.txt
```

- For m1 mac user
```
cd optinist/backend
pip install -r requirements_m1mac.txt
```

#### caimanを個別にインストール
* m1 macの場合、tensorflowがうまくインストールできず、caimanを動かすことができないため、ここはスキップ。

optinistの仮想環境に入る。（既に入っている場合には省略可）
```
conda activate optinist
```

caimanをインストールする。インストールする場所はユーザに任せる。
```
git clone https://github.com/flatironinstitute/CaImAn -b v1.9.7
cd /CaImAn && pip install -e .
```

* ```pip list```などでcaimanやsuite2pなどが正しくインストールできていることを確認すると良い。
* caimanをインストール後にnumpyなどを再インストールすると、numpyのエラーが起きることが分かっているが、これはnumpyの仕様が変更されたためである。これに対する明確な解決策が提示されていない。
versionを上げるとcaimanのエラーは解決されるが、versionを挙げたことによって逆にsuite2pが使っているnumbaライブラリでのエラーがでる。
[https://stackoverflow.com/questions/66060487/valueerror-numpy-ndarray-size-changed-may-indicate-binary-incompatibility-exp]

#### 実行
backendディレクトリに移動
```
cd optinist/backend
```

実行
gunicornで起動
```
> gunicorn
```


### ファイルの保存場所
`/tmp/optinist`に保存されるため、直接tifファイルをこちらに入れてもらうと、ファイルをアップロードする必要がなくなる。
* `/tmp`フォルダはPCを再起動するとデータが消えるため注意が必要。
