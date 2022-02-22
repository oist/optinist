Table of Contents
=================

* [Installation](#installation)
* [0. GitHub SSH access settings](#0-github-ssh-access-settings)
* [1. Clone optinist repository](#1-clone-optinist-repository)
* [2. Make frontend environment](#2-make-frontend-environment)
      * [Install docker](#install-docker)
      * [Install docker-compose](#install-docker-compose)
      * [docker-compose buld frontend](#docker-compose-buld-frontend)
      * [docker-compose up frontend](#docker-compose-up-frontend)
      * [Launch browser.  <a href="http://localhost:3000" rel="nofollow">http://localhost:3000</a>](#launch-browser--httplocalhost3000)
* [3. Make backend environment](#3-make-backend-environment)
   * [For Linux Users](#for-linux-users)
      * [Install gcc, g++](#install-gcc-g)
      * [Install Anaconda](#install-anaconda)
      * [Create anaconda environment](#create-anaconda-environment)
      * [Install mamba](#install-mamba)
      * [Install library](#install-library)
      * [run backend](#run-backend)
      * [Reload your browser. If you can select the "Algorithm" tab, it is working correctly!](#reload-your-browser-if-you-can-select-the-algorithm-tab-it-is-working-correctly)

# Installation
We introduce how to install optinist.
We have developed optinist python(backend) and typescript(frontend), so you need to make both environment.
Please follow instructions below.

<br />

**CAUTION**: We confirmed them on Ubuntu 18.04 or 20.04.

# 0. GitHub SSH access settings
**You only need to do the following once.**

On your shell,
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
### Install docker
```
sudo apt update
sudo apt install apt-transport-https ca-certificates curl software-properties-common

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu focal stable"
sudo apt update

sudo apt install docker-ce
```
- To avoid using docker commands with sudo, you need to set $USER to docker user group.
```
sudo usermod -aG docker $USER
```
### Install docker-compose
```
wget https://github.com/docker/compose/releases/download/v2.2.3/docker-compose-linux-x86_64
mv docker-compose-linux-x86_64 docker-compose
sudo mv docker-compose /usr/local/bin/docker-compose
```
### docker-compose buld frontend
```
docker-compose build frontend
```
- `docker-compose build frontend` log is as blow:
```
$ docker-compose up frontend
Sending build context to Docker daemon    465kB
Step 1/8 : FROM node:14 as node
14: Pulling from library/node
a834d7c95167: Pull complete
57b3fa6f1b88: Pull complete
778df3ecaa0f: Pull complete
d353c340774e: Pull complete
6370e0bc373d: Pull complete
fb61153482cd: Pull complete
78fb5822e501: Pull complete
ba3577a691be: Pull complete
bd38fd0dd57b: Pull complete
Digest: sha256:b2c75df8c9706156c38b4f1f678d00e11cb2bfda09fc4ab6e36ec17ac9163865
Status: Downloaded newer image for node:14
 ---> 9cb3f042a684
Step 2/8 : COPY . ./
 ---> f688cc7411c1
Step 3/8 : ENV GENERATE_SOURCEMAP false
 ---> Running in bdcec15e8d6a
Removing intermediate container bdcec15e8d6a
 ---> 3667fcdd6fda
Step 4/8 : RUN yarn install
 ---> Running in 5f5656b1e45e
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
fatal: Not a git repository (or any of the parent directories): .git
Done in 27.26s.
Removing intermediate container 5f5656b1e45e
 ---> e2b5f0cd7788
Step 5/8 : RUN yarn build
 ---> Running in 4f2bd6935a86
yarn run v1.22.17
$ react-scripts build
Creating an optimized production build...
Compiled successfully.

File sizes after gzip:

  1.44 MB   build/static/js/2.7d54d872.chunk.js
  21.99 KB  build/static/js/main.d29046c0.chunk.js
  1.58 KB   build/static/js/3.f30dcf7f.chunk.js
  1.14 KB   build/static/css/2.17fb7f09.chunk.css
  1.12 KB   build/static/js/runtime-main.3ef6eb24.js
  565 B     build/static/css/main.021ef9c2.chunk.css

The bundle size is significantly larger than recommended.
Consider reducing it with code splitting: https://goo.gl/9VhYWB
You can also analyze the project dependencies: https://goo.gl/LeUzfb

The project was built assuming it is hosted at /.
You can control this with the homepage field in your package.json.

The build folder is ready to be deployed.
You may serve it with a static server:

  yarn global add serve
  serve -s build

Find out more about deployment here:

  https://cra.link/deployment

Done in 55.35s.
Removing intermediate container 4f2bd6935a86
 ---> 28912f086e58
Step 6/8 : FROM nginx:1.20.1
1.20.1: Pulling from library/nginx
b380bbd43752: Already exists
83acae5e2daa: Already exists
33715b419f9b: Already exists
eb08b4d557d8: Already exists
74d5bdecd955: Already exists
0820d7f25141: Already exists
Digest: sha256:a98c2360dcfe44e9987ed09d59421bb654cb6c4abe50a92ec9c912f252461483
Status: Downloaded newer image for nginx:1.20.1
 ---> c8d03f6b8b91
Step 7/8 : COPY --from=node ./build /usr/share/nginx/html
 ---> be1970e94e1d
Step 8/8 : CMD nginx -g "daemon off;"
 ---> Running in 53d799653e81
Removing intermediate container 53d799653e81
 ---> a85095edc139
Successfully built a85095edc139
Successfully tagged optinist_frontend:latest

Use 'docker scan' to run Snyk tests against images to find vulnerabilities and learn how to fix them
$
```
### docker-compose up frontend
```
docker-compose up frontend
```
- `docker-compose up frontend` log is as blow:
```
$ docker-compose up frontend
[+] Running 1/1
 ⠿ Container optinist-frontend-1  Recre...                                 0.0s
Attaching to optinist-frontend-1
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: using the "epoll" event method
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: nginx/1.20.1
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: built by gcc 8.3.0 (Debian 8.3.0-6)
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: OS: Linux 5.4.0-99-generic
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: getrlimit(RLIMIT_NOFILE): 1048576:1048576
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker processes
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 8
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 9
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 10
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 11
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 12
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 13
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 14
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 15
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 16
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 17
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 18
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 19
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 20
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 21
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 22
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 23
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 24
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 25
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 26
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 27
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 28
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 29
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 30
optinist-frontend-1  | 2022/02/14 07:17:09 [notice] 7#7: start worker process 31
```
### Launch browser.  http://localhost:3000
It opens correctly!

<br />

# 3. Make backend environment
## For Linux Users
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
