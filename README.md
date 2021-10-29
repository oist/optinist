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

# 開発環境

## backend

Docker コンテナ上で API サーバーを動かします。pip 等で依存ライブラリをインストールする必要はありません。

### [Docker](https://docs.docker.com/)のインストール

[https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)からダウンロードする。

### Docker イメージ作成

```
$ docker-compose build
```

### Docker コンテナ起動

```
$ docker-compose up
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

- [http://localhost:3000](http://localhost:3000)にアクセス
