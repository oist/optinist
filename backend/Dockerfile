FROM python:3.9.7

RUN mkdir /app
WORKDIR /app

RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y git gcc g++ libgl1

RUN /usr/local/bin/python -m pip install --upgrade pip
COPY requirements.txt ./
RUN pip install numpy Cython
RUN pip install --no-cache-dir -r requirements.txt
# RUN export PIP_NO_CACHE_DIR=off
# RUN pip install pipenv
# COPY requirements.txt ./
# COPY Pipfile ./
# RUN PIP_NO_CACHE_DIR=off
# COPY Pipfile.lock ./
# RUN pipenv --python 3.9
# RUN /bin/bash -c "pipenv install -r requirements.txt"
# RUN set -ex && pipenv install tensorflow 
# RUN /bin/bash -c "pipenv install --skip-lock --deploy --ignore-pipfile"
# RUN /bin/bash -c "pipenv sync"

# ENV HOME /root

# # caiman install
# RUN git clone https://github.com/flatironinstitute/CaImAn /CaImAn
# RUN cd /CaImAn && pip install -e .

RUN apt update && \
    apt-get install -y libgl1-mesa-dev

COPY . .

ENTRYPOINT ["python", "main.py" ]