FROM --platform=linux/amd64 python:3.8.18-slim

WORKDIR /app

# install default packages
RUN apt-get --allow-releaseinfo-change update && \
    apt-get install --no-install-recommends -y sudo procps iproute2 iputils-ping less vim

# install conda & packages
RUN mkdir -p /opt/miniconda3 && \
    apt-get install --no-install-recommends -y gcc g++ libgl1 libgl1-mesa-dev libopencv-dev wget && \
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -O /opt/miniconda3/miniconda.sh && \
    apt-get purge wget -y && apt-get autoremove -y && apt-get clean && \
    bash /opt/miniconda3/miniconda.sh -b -u -p /opt/miniconda3 && \
    rm /opt/miniconda3/miniconda.sh && \
    export PATH="$PATH:/opt/miniconda3/bin" && \
    conda upgrade -y --all && \
    conda config --set channel_priority strict && \
    conda clean -y --tarballs
ENV PATH $PATH:/opt/miniconda3/bin

# setup user profile
RUN echo 'alias ll="ls -la --color=auto"' >> /root/.bashrc && \
    echo 'set nu ic hls nowrap ts=4 sw=4 | syntax on' >> /root/.vimrc

# setup optinist
COPY pyproject.toml poetry.lock ./
RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install poetry
RUN poetry install --no-root

# copy optinist files
COPY frontend/build /app/frontend/build
COPY studio /app/studio
COPY sample_data /app/sample_data
COPY main.py ./

EXPOSE 8000
