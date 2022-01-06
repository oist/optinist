# optinist

optinist のディレクトリ。

## ディレクトリ構成

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


# 環境構築
2022/1月時点では、frontendとbackendを別々に環境構築するのを推奨する。
frontendはdockerで構築する。
backendはメモリを多く使用するためローカル環境で構築する。

## frontendの環境構築

### docker build
｀optinist｀ディレクトリに移動し、dockerをbuildする。
```
cd optinist
docker-compose build
```
### frontendをdockerで起動
dockerでfront側を起動する。
```
docker-compose up frontend
```

## backendの環境構築

### pipenvで起動
### pipenvをインストール
ローカル環境にpipenvがない場合はインストール。
```
pip install pipenv
```

pyenvもインストールする。
macの場合
```
brew install pyenv
```

windowsの場合のpyenvのインストール[https://github.com/pyenv/pyenv#windows]

<br />

### 仮想環境の構築
backendディレクトリに移動し、packageをインストールする。
```
cd optinist/backend
pipenv install --skip-lock
```

### backendの起動
backendをpipenv環境で実行する。
```
pipenv run python main.py
```

### ファイルの保存場所
`/tmp/optinist`に保存されるため、直接tifファイルをこちらに入れてもらうと、ファイルをアップロードする必要がなくなる。


<!-- # 開発環境

## backend

Docker コンテナ上で API サーバーを動かします。pip 等で依存ライブラリをインストールする必要はありません。

### [Docker](https://docs.docker.com/)のインストール

[https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)からダウンロードする。

### Docker イメージ作成

```
$ docker-compose build backend
```

### Docker コンテナ起動

```
$ docker-compose up backend
```

## frontend

### [yarn](https://yarnpkg.com/) のインストール

1. [https://nodejs.org/](https://nodejs.org/)から Node.js をインストール

2. Node.js がインストールできたら、yarn をインストールする
   ```
   $ npm install -g yarn
   ```

### プロジェクトの依存パッケージをインストール

```
$ cd ./frontend
$ yarn
```

### アプリの実行

```
$ cd ./frontend
$ yarn start
```

- [http://localhost:3000](http://localhost:3000)にアクセス -->
