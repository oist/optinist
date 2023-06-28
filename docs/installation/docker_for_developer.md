Docker for Developer
=================

```{contents}
:depth: 4
```

## Installation

We introduce how to install optinist for developer.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

## 1. Make backend environment

### Install Tools

- Install [Docker](https://www.docker.com/products/docker-desktop/)

### Clone repository

```
git clone https://github.com/oist/optinist.git
cd ./optinist
```

### Start docker container

```
docker compose -f docker-compose.dev.yml up
```

- add `-d` option to make container run in background

<!--
## 2. Create virtualenv

Under maintenance...
-->

## 2. Access to backend

- Launch browser, and go to http://localhost:3000
- Your local code change will be applied on save.

```{eval-rst}
.. note::
    dev container uses port 3000,
    while production docker image uses 8000.
```

Done!
