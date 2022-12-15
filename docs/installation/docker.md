Docker
=================

* [Installation](#installation)
* [1. Make docker image container](#1-make-docker-image-container)
   * [Make docker image](#make-docker-image)
   * [Set saving directory](#set-saving-directory)

## Installation
We introduce how to install studio.
We have developed studio python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

## 1. Make docker image container

### Make docker image
Pull the latest docker image from docker hub.
```
docker pull oistncu/studio
```
Start docker container.
```
docker run -it --shm-size=2G --name studio_container -d -p 8000:8000 --restart=unless-stopped oistncu/studio:latest
```

Execute in terminal
```
docker exec -it studio_container /bin/bash
```

### Set saving directory
studio default saving directory is `/tmp/studio`. If you reboot your PC, this repogitory content is deleted. And setting the saving directory in environment path.
```bash
export studio_DIR="your_saving_dir"
```

Done!
