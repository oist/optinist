Table of Contents
=================

* [Installation](#installation)
* [0. GitHub SSH access settings](#0-github-ssh-access-settings)
* [1. Clone optinist repository](#1-clone-optinist-repository)
* [2. Make docker image container](#2-Make-docker-image-container)
   * [Make docker image](#make-docker-image)
   * [run backend](#run-backend)
   * [Launch browser.  <a href="http://localhost:8000" rel="nofollow">http://localhost:8000</a>](#launch-browser--httplocalhost8000)
# Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

**CAUTION for M1 Mac User**
We use tensorflow in caiman code. We know that M1 mac doesn't install tensorflow easily, so if there is a problem, skip install caiman. (Release in progress…)

# 0. GitHub SSH access settings
**You only need to do the following once.**

Follow this [link](installation_github_settings.md).

# 1. Clone optinist repository

First, you get optinist code from github repository.
```
cd "your working repository"
git clone git@github.com:oist/optinist.git
```
<br />

# 2. Make docker image container

## Make docker image
Buid docker
```
docker build --rm -t optinist:latest .
```

**CAUTION for M1 Mac User**

For M1 mac user, please add option `--platform=linux/amd64`
```
docker build --rm -t optinist:latest . --platform=linux/amd64
```

Start docker container
```
docker run -it --shm-size=2G --name optinist_container -d -p 8000:8000 --restart=unless-stopped optinist:latest
```

Execute in terminal
```
docker exec -it optinist_container /bin/bash
```



## run backend
Change backend directory
```
cd /app/backend
```


```
python main.py
```
- `python main.py` log is as blow:
```
$ python main.py
INFO:   Will watch for changes in these directories: [‘/Users/oist/optinist/backend’]
INFO:   Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:   Started reloader process [5811] using statreload
INFO:   Started server process [5820]
INFO:   Waiting for application startup.
INFO:   Application startup complete.
```
## Launch browser.  http://localhost:8000
It opens correctly!

Done!
