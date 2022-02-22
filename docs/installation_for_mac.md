Table of Contents
=================

* [Installation](#installation)
* [0. GitHub SSH access settings](#0-github-ssh-access-settings)
* [1. Clone optinist repository](#1-clone-optinist-repository)
* [2. Make frontend environment](#2-make-frontend-environment)
      * [Install Docker desktop](#install-docker-desktop)
      * [docker-compose build frontend](#docker-compose-build-frontend)
      * [docker-compose up frontend](#docker-compose-up-frontend)
      * [Launch browser.  <a href="http://localhost:3000" rel="nofollow">http://localhost:3000</a>](#launch-browser--httplocalhost3000)
* [3. Make backend environment](#3-make-backend-environment)
      * [Install Anaconda](#install-anaconda)
      * [Create anaconda environment](#create-anaconda-environment)
      * [Install mamba](#install-mamba)
      * [Install library](#install-library)
      * [Install CaImAn](#install-caiman)
      * [run backend](#run-backend)
      * [Reload your browser. If you can select the "Algorithm" tab, it is working correctly!](#reload-your-browser-if-you-can-select-the-algorithm-tab-it-is-working-correctly)
      * [FAQ](#faq)
# Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

**CAUTION for M1 Mac User**
We use tensorflow in caiman code. We know that M1 mac doesn't install tensorflow easily, so if there is a problem, skip install caiman. (Release in progress…)

# 0. GitHub SSH access settings
**You only need to do the following once.**

Open terminal App,
```
git config --global user.name "user.name"
git config --global user.email "user@oist.jp"
git config --global core.quotepath false
```
make your ssh key
```
mkdir ~/.ssh
cd ~/.ssh
ssh-keygen -t rsa -C 'user@oist.jp'
#Enter
#Input passphrase (empty for no passphrase):
#Input passphrase (again)
```
Open the generated public key (rsa.pub) with a text editor and copy all the contents.
Go to GitHub and follow the steps below to register your public key.

1. Log in to GitHub and select `Settings` from the menu on the top right
2. Select `SSH and GPG keys`
3. Press `New SSH Key`
4. Enter the `Title` (free) and `Key` (paste the copied content) and press `Add SSH` key.

This completes the SSH connection settings!

# 1. Clone optinist repository

First, you get optinist code from github repository.
```
cd "your working repository"
git clone git@github.com:oist/optinist.git
```
<br />

# 2. Make frontend environment
### Install Docker desktop
- Install [Docker Desktop](https://www.docker.com/products/docker-desktop) for Mac.

### docker-compose build frontend
```
docker-compose build frontend
```
- `docker-compose build frontend` log is as blow:
```
$ docker-compose build frontend
[+] Building 301.5s (9/11)
[+] Building 329.7s (9/11)
[+] Building 455.9s (12/12) FINISHED
 => [internal] load build definition from Dockerfile                           0.0s
 => => transferring dockerfile: 224B                                           0.0s
 => [internal] load .dockerignore                                              0.0s
 => => transferring context: 70B                                               0.0s
 => [internal] load metadata for docker.io/library/nginx:1.20.1                4.3s
 => [internal] load metadata for docker.io/library/node:14                     4.3s
 => [internal] load build context                                              0.2s
 => => transferring context: 1.08MB                                            0.2s
 => CACHED [node 1/4] FROM docker.io/library/node:14@sha256:b2c75df8c9706156c3 0.0s
 => CACHED [stage-1 1/2] FROM docker.io/library/nginx:1.20.1@sha256:a98c2360dc 0.0s
 => [node 2/4] COPY . ./                                                       0.1s
 => [node 3/4] RUN yarn install                                              123.0s
 => [node 4/4] RUN yarn build                                                316.1s
 => [stage-1 2/2] COPY --from=node ./build /usr/share/nginx/html               0.1s
 => exporting to image                                                         0.1s
 => => exporting layers                                                        0.1s
 => => writing image sha256:4953a0a63a7842e7a108d4d32c567e04de5fbae87d83d9f792 0.0s
 => => naming to docker.io/library/optinist_frontend                           0.0s
Use ‘docker scan’ to run Snyk tests against images to find vulnerabilities and learn how to fix them
$
```

If you haven't installed docker-compose, you need to install docker-compose.
```
wget https://github.com/docker/compose/releases/download/v2.2.3/docker-compose-darwin-x86_64
mv docker-compose-darwin-x86_64 docker-compose
sudo mv docker-compose /usr/local/bin/docker-compose
```

### docker-compose up frontend
```
docker-compose up frontend
```
- `docker-compose up frontend` log is as blow:
```
$ docker-compose up frontend
+] Running 1/1
⠿ Container optinist-frontend-1 Recreated                   0.2s
Attaching to optinist-frontend-1
optinist-frontend-1 | 2022/02/14 06:38:16 [notice] 7#7: using the “epoll” event method
optinist-frontend-1 | 2022/02/14 06:38:16 [notice] 7#7: nginx/1.20.1
optinist-frontend-1 | 2022/02/14 06:38:16 [notice] 7#7: built by gcc 8.3.0 (Debian 8.3.0-6)
optinist-frontend-1 | 2022/02/14 06:38:16 [notice] 7#7: OS: Linux 5.10.76-linuxkit
optinist-frontend-1 | 2022/02/14 06:38:16 [notice] 7#7: getrlimit(RLIMIT_NOFILE): 1048576:1048576
optinist-frontend-1 | 2022/02/14 06:38:16 [notice] 7#7: start worker processes
optinist-frontend-1 | 2022/02/14 06:38:16 [notice] 7#7: start worker process 8
optinist-frontend-1 | 2022/02/14 06:38:16 [notice] 7#7: start worker process 9
```
### Launch browser.  http://localhost:3000
It opens correctly!

<br />

# 3. Make backend environment

### Install Anaconda
Download https://repo.anaconda.com/archive/Anaconda3-2021.11-MacOSX-x86_64.pkg

Install

### Create anaconda environment
```
conda create -n optinist python=3.8
conda activate optinist
```
### Install mamba
We use snakemake library, and it requires mamba.
```
conda install -n base -c conda-forge mamba
```
### Install library
```bash
cd backend
pip install -r requirements.txt
cd ..
```

In case an error occuered when you install ``` requirements.txt ```, pip upgrade command below may solve the error.
```
pip install --upgrade pip
```

### Install CaImAn
**CAUTION for M1 Mac User**
We use tensorflow in caiman code. We know that M1 mac doesn't install tensorflow easily, so if there is a problem, skip install caiman. (Release in progress…)

If you are m1 mac user, skip this.

```bash
# for CaImAn
git clone https://github.com/flatironinstitute/CaImAn -b v1.9.7
pip install cython opencv-python matplotlib scikit-image==0.18.0 scikit-learn ipyparallel holoviews watershed tensorflow
cd CaImAn
pip install -e .
cd ..
```
### run backend
```
python main.py
```
- `python main.py` log is as blow:
```
$ python main.py
INFO:   Will watch for changes in these directories: [‘/Users/oist/work/optinist/backend’]
INFO:   Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:   Started reloader process [2678] using statreload
INFO:   Started server process [2690]
INFO:   Waiting for application startup.
INFO:   Application startup complete.
```
### Reload your browser. If you can select the "Algorithm" tab, it is working correctly!
Done!


### FAQ
If you get the warning message shown below, we recommend `rm -rf /Users/usename/opt/anaconda3/envs/optinist ` and recreate conda environment.
> WARNING: A directory already exists at the target location '/Users/usename/opt/anaconda3/envs/optinist' but it is not a conda environment.
