# Installation
We introduce how to install optinist.
We develop optinist python(backend) and typescript(frontend), so you need to make both of environment.
Please follow below instrunction.

<br />

# 1. Clone optinist repository

First, you get optinist code from github repository.
```
cd "your working repository"
git clone git@github.com:oist/optinist.git
```

If you don't have access authority, please contact with in charge of person.

<br />

# 2. Make frontend environment

## Install docker

Install docker from docker site. [https://docs.docker.com/get-docker/]

**** If you get memory error, increase memory to use. ****

Launching docker and increase memory capacity.

How to open setting display.
- [For Windows](https://docs.docker.com/desktop/windows/)
- [For Mac](https://docs.docker.com/desktop/mac/)
- [For Linux](https://docs.docker.com/desktop/linux/)

**** memory setting ****

Please, Check whether docker open or not in the command.
```
docker ps
```

If you get below message, docker doesn't launch. so click docker icon in your LaunchPad(mac).

> Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running?


## docker build
check your working repository in `optinist`.
```
cd optinist
```

Build frontend code.
```
docker-compose build frontend
```

## docker up
Launch frontend.
```
docker-compose up frontend
```

Launch browser.  http://localhost:3000  
It open correctly!

<br />

# 3. make backend environment

<br />

## Create anaconda or virtualenv
Prepare virtual environment, anaconda or virtualenv.
You can install in local environment directly, but I recommend using virtual environment.


We introduce how to setup in anaconda environment.

Create virtual environment, named `optinist`
```
conda create -n optinist python=3.9.7
conda activate optinist
```

<br />

## Install library
Check you're working directory in `optinist/backend`.
```
cd optinist/backend
```

Install library from requirements.txt
```
pip install -r requirements.txt
```

We know that M1 mac doesn't install tensorflow easily, but we use tensorflow in caiman code, so comment out tensorflow from requirements.txt and skip caiman. (Release in progress…)


<br />

## Install caiman
** If you use m1 mac, skip this.

Check you're working in `optinist` environment.
```
conda activate optinist
```

Install caiman, donwload directory is up to you.
```
git clone https://github.com/flatironinstitute/CaImAn -b v1.9.7
cd /CaImAn && pip install -e .
```

<br />

## Check library
* You can check install library from ```pip list``` command.
* If you install numpy after install caiman, you get numpy error. Because of numpy and Cython version conflict. 
[https://stackoverflow.com/questions/66060487/valueerror-numpy-ndarray-size-changed-may-indicate-binary-incompatibility-exp]

## Run backend
Check you're working directory in `optinist/backend`
```
cd optinist/backend
```

run with gunicorn
```
> gunicorn
```

It open correctly!
