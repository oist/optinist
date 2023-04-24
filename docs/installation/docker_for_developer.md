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

### Build docker image

```
docker build --no-cache -t optinist-dev -f Dockerfile .
```

### Start docker container

```
docker run -it --shm-size=2G --name optinist_dev_container -d -p 8000:8000 --restart=unless-stopped optinist-dev

```

<!--
## 2. Create virtualenv

Under maintenance...
-->

## 2. Access to backend

- Launch browser, and go to http://localhost:8000

It opens correctly!

Done!
