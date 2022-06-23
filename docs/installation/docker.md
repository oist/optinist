Docker
=================

* [Installation](#installation)
* [1. Make docker image container](#1-make-docker-image-container)
   * [Make docker image](#make-docker-image)
   * [Set saving directory](#set-saving-directory)

## Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

## 1. Make docker image container

### Make docker image
Pull the latest docker image from docker hub.
```
docker pull oistncu/optinist
```
Start docker container.
```
docker run -it --shm-size=2G --name optinist_container -d -p 8000:8000 --add-host="localhost:127.0.0.1" --restart=unless-stopped oistncu/optinist:latest
```

Execute in terminal
```
docker exec -it optinist_container /bin/bash
```

### Set saving directory
Optinist default saving directory is `/tmp/optinist`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```bash
export OPTINIST_DIR="your_saving_dir"
```

Done!
