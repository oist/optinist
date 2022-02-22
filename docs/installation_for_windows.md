Table of Contents
=================

* [Installation](#installation)
* [0. GitHub SSH access settings](#0-github-ssh-access-settings)
* [1. Clone optinist repository](#1-clone-optinist-repository)
* [2. Make frontend environment](#2-make-frontend-environment)
   * [For Windows(PowerShell) Users](#for-windowspowershell-users)
      * [Install nodejs and yarn](#install-nodejs-and-yarn)
      * [Launch browser.  <a href="http://localhost:3000" rel="nofollow">http://localhost:3000</a>](#launch-browser--httplocalhost3000)
   * [For Windows(WSL2) Users](#for-windowswsl2-users)
      * [Install Docker desktop for Windows and Settings](#install-docker-desktop-for-windows-and-settings)
      * [Install docker-compose](#install-docker-compose)
      * [docker-compose build](#docker-compose-build)
      * [docker-compose up frontend](#docker-compose-up-frontend)
      * [Launch browser.  <a href="http://localhost:3000" rel="nofollow">http://localhost:3000</a>](#launch-browser--httplocalhost3000-1)
* [3. Make backend environment](#3-make-backend-environment)
   * [For Windows(PowerShell) Users](#for-windowspowershell-users-1)
      * [Install Visutal Studio Build Tools](#install-visutal-studio-build-tools)
      * [Install Anaconda](#install-anaconda)
      * [Create anaconda environment](#create-anaconda-environment)
      * [Install mamba](#install-mamba)
      * [Install library](#install-library)
      * [run backend](#run-backend)
      * [Reload your browser. If you can select the "Algorithm" tab, it is working correctly!](#reload-your-browser-if-you-can-select-the-algorithm-tab-it-is-working-correctly)
   * [For Windows (WSL2) Users](#for-windows-wsl2-users)
      * [Install gcc, g++](#install-gcc-g)
      * [Install Anaconda](#install-anaconda-1)
      * [Create anaconda environment](#create-anaconda-environment-1)
      * [Install mamba](#install-mamba-1)
      * [Install library](#install-library-1)
      * [run backend](#run-backend-1)
      * [Reload your browser. If you can select the "Algorithm" tab, it is working correctly!](#reload-your-browser-if-you-can-select-the-algorithm-tab-it-is-working-correctly-1)
# Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

**CAUTION**: We strongly recommend **Windows 10** (not 11).

**CAUTION**: For WSL2, we confirmed them on [Ubuntu 20.04](https://www.microsoft.com/ja-jp/p/ubuntu-2004-lts/9n6svws3rx71).


# 0. GitHub SSH access settings
**You only need to do the following once.**

If you don't have git installed, you need to install  [git-for-windows](https://git-scm.com/download/win) or, git command in WSL2.

On Git Bash or WSL2
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
Clone CaImAn repository.
```
cd optinist/backend
git clone https://github.com/flatironinstitute/CaImAn -b v1.9.7
```
<br />

# 2. Make frontend environment
## For Windows(PowerShell) Users
### Install nodejs and yarn
You need to install nodejs and yarn as below:
- Install [nodejs](https://nodejs.org/ja/download/) from .msi installer.
- After installation, start PowerShell with administrator, run the following command:
```
npm install --global yarn
Set-ExecutionPolicy RemoteSigned
```
On PowerShell, start frontend.
```
cd optinist\frontend
yarn
yarn start
```
- `yarn` and `yarn start` command log is as blow:
```
PS C:\optinist\frontend> yarn
yarn install v1.22.17
[1/4] Resolving packages...
[2/4] Fetching packages...
[3/4] Linking dependencies...
warning "@emotion/styled > @emotion/babel-plugin@11.7.2" has unmet peer dependency "@babel/core@^7.0.0".
warning "@emotion/styled > @emotion/babel-plugin > @babel/plugin-syntax-jsx@7.16.7" has unmet peer dependency "@babel/core@^7.0.0-0".
warning " > @mui/x-data-grid@5.5.1" has unmet peer dependency "@mui/system@^5.2.8".
warning " > @testing-library/user-event@12.8.3" has unmet peer dependency "@testing-library/dom@>=7.21.4".
warning " > react-plotlyjs-ts@2.2.2" has incorrect peer dependency "plotly.js@^1.31.2".
warning " > react-plotlyjs-ts@2.2.2" has incorrect peer dependency "react@^16.0.0".
warning " > react-plotlyjs-ts@2.2.2" has incorrect peer dependency "react-dom@^16.0.0".
warning " > react-plotlyjs-ts@2.2.2" has incorrect peer dependency "typescript@^2.5.3".
warning " > styled-components@5.3.3" has unmet peer dependency "react-is@>= 16.8.0".
[4/4] Building fresh packages...
$ cd .. && husky install frontend/.husky
Done in 29.67s.
PS C:\optinist\frontend> yarn start
yarn run v1.22.17
$ react-scripts start
i ｢wds｣: Project is running at http://192.168.1.xx/
i ｢wds｣: webpack output is served from
i ｢wds｣: Content not from webpack is served from C:\Users\oist\Documents\optinist\frontend\public
i ｢wds｣: 404s will fallback to /
Starting the development server...

=============

WARNING: You are currently running a version of TypeScript which is not officially supported by @typescript-eslint/typescript-estree.

You may find that it works just fine, or you may not.

SUPPORTED TYPESCRIPT VERSIONS: >=3.3.1 <4.5.0

YOUR TYPESCRIPT VERSION: 4.5.5

Please only submit bug reports when using the officially supported version.

=============
Compiled successfully!

You can now view optinist in the browser.

  Local:            http://localhost:3000
  On Your Network:  http://192.168.1.xx:3000

Note that the development build is not optimized.
To create a production build, use yarn build.
```

### Launch browser.  http://localhost:3000
It opens correctly!

## For Windows(WSL2) Users
### Install Docker desktop for Windows and Settings
- Install [Docker Desktop](https://www.docker.com/products/docker-desktop) for Windows.
  - During installation, mark checkbox for `Install required Windows components for WSL2`
- Reboot Windows
- From Docker Desktop Settings, enable the following items,
  - General -> `Use the WSL2 based engine`
  - Resources -> WSL INTEGRATION -> `Enables integration with additional distros:` and `Ubuntu 20.04`
  - Apply and Restart

### Install docker-compose
```
wget https://github.com/docker/compose/releases/download/v2.2.3/docker-compose-linux-x86_64
mv docker-compose-linux-x86_64 docker-compose
sudo mv docker-compose /usr/local/bin/docker-compose
```

### docker-compose build
- To avoid using docker commands with sudo, you need to set $USER to docker user group.
```
sudo usermod -aG docker $USER
```
- Download and start docker image
```
cd optinist
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

### docker-compose up frontend
```
docker-compose up frontend
```
- `docker-compose up frontend` log is as blow:
```
$ docker-compose up frontend
Attaching to optinist-frontend-1
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: using the "epoll" event method
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: nginx/1.20.1
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: built by gcc 8.3.0 (Debian 8.3.0-6)
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: OS: Linux 5.10.16.3-microsoft-standard-WSL2
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: getrlimit(RLIMIT_NOFILE): 1048576:1048576
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker processes
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 8
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 9
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 10
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 11
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 12
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 13
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 14
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 15
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 16
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 17
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 18
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 19
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 20
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 21
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 22
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 23
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 24
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 25
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 26
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 27
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 28
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 29
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 30
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 31
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 32
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 33
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 34
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 35
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 36
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 37
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 38
optinist-frontend-1  | 2022/02/14 04:06:31 [notice] 7#7: start worker process 39
```
### Launch browser.  http://localhost:3000
It opens correctly!

<br />

# 3. Make backend environment

## For Windows(PowerShell) Users
### Install Visutal Studio Build Tools
- For install CaImAn, you need to install Visual Studio Build Tools.
  - Download `Build Tools for Visual Studio 2022` from https://visualstudio.microsoft.com/ja/downloads/
  - In insteraller, select `Desktop Application for C++`
### Install Anaconda
Install [Anaconda for Windows](https://www.anaconda.com/products/individual)
### Create anaconda environment
On the Anaconda PowerShell Prompt(anaconda3),
```
conda create -n optinist python=3.8
conda activate optinist
```
### Install mamba
We use snakemake library, and it requires mamba.
On the Anaconda PowerShell Prompt(anaconda3),
```
conda install -n base -c conda-forge mamba
```
### Install library
On the Anaconda PowerShell Prompt(anaconda3),
```bash
cd backend
pip install -r requirements.txt
# for suite2p
pip install PyQt5<=5.15.1 PyQt5-sip<=12.8.1 pyqtgraph<=0.11.0 pandas suite2p<=0.10.3
# for GLM
pip install sklearn statsmodels<=0.13.1 pynwb
# for CaImAn
pip install cython opencv-python matplotlib scikit-image==0.18.0 scikit-learn ipyparallel holoviews watershed tensorflow
cd CaImAn
pip install -e .
cd ..
```
### run backend
On the Anaconda PowerShell Prompt(anaconda3),
```
python main.py
```
- `python main.py` log is as blow:
```
(optinist) PS C:\optinist\backend> python main.py
[32mINFO[0m:     Will watch for changes in these directories: ['C:\\optinist\\backend']
[32mINFO[0m:     Uvicorn running on [1mhttp://0.0.0.0:8000[0m (Press CTRL+C to quit)
[32mINFO[0m:     Started reloader process [[36m[1m22996[0m] using [36m[1mstatreload[0m
[32mINFO[0m:     Started server process [[36m19940[0m]
[32mINFO[0m:     Waiting for application startup.
[32mINFO[0m:     Application startup complete.
```
### Reload your browser. If you can select the "Algorithm" tab, it is working correctly!
Done!

## For Windows (WSL2) Users
### Install gcc, g++
- For install CaImAn, you need to install gcc and g++.
```
sudo apt install gcc g++
```
### Install Anaconda
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
bash Anaconda3-2021.11-Linux-x86_64.sh
```
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
INFO:     Will watch for changes in these directories: ['/home/oist/optinist/backend']
INFO:     Uvicorn running on http://0.0.0.0:8000 (Press CTRL+C to quit)
INFO:     Started reloader process [13700] using statreload
INFO:     Started server process [13744]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
```
### Reload your browser. If you can select the "Algorithm" tab, it is working correctly!
Done!