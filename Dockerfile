FROM --platform=linux/amd64 python:3.8.16-slim

COPY requirements.txt /app/requirements.txt

RUN apt-get --allow-releaseinfo-change update && \
    apt-get install --no-install-recommends -y gcc g++ libgl1 libgl1-mesa-dev && \
    pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir -r /app/requirements.txt && \
    apt-get clean && rm -rf /root/.cache/pip/*

RUN mkdir -p /root/miniconda3 && \
    apt-get install --no-install-recommends -y wget && \
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
    apt-get purge wget -y && apt-get autoremove -y && apt-get clean && \
    bash /root/miniconda3/miniconda.sh -b -u -p /root/miniconda3 && \
    rm /root/miniconda3/miniconda.sh && \
    export PATH="$PATH:/root/miniconda3/bin" && \
    conda upgrade -y --all && \
    conda config --set channel_priority strict && \
    conda clean -y --tarballs

COPY frontend/build /app/frontend/build
COPY optinist /app/optinist
COPY main.py /app/main.py

ENV PATH $PATH:/root/miniconda3/bin
WORKDIR /app
EXPOSE 8000

ENTRYPOINT ["python3", "main.py"]
CMD ["--host", "0.0.0.0"]