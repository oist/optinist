FROM python:3.9.7-slim

COPY requirements.txt /app/requirements.txt

RUN apt-get --allow-releaseinfo-change update && \
    apt-get install --no-install-recommends -y git gcc g++ libgl1 libgl1-mesa-dev && \
    pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir -r /app/requirements.txt && \
    pip3 install --no-cache-dir \
        cython==0.29.30 \
        holoviews==1.14.9 \
        ipyparallel==8.4.1 \
        matplotlib==3.5.2 \
        opencv-python==4.6.0.66 \
        scikit-image==0.18.0 \
        scikit-learn==1.1.1 \
        tensorflow==2.9.1 && \
        watershed==2.2.2 \
    pip3 install --no-cache-dir git+https://github.com/flatironinstitute/CaImAn.git@914324989443fac5d481ef32aad4f327701294a8#egg=caiman && \
    apt-get purge git -y && apt-get autoremove -y && apt-get clean && rm -rf /root/.cache/pip/*

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