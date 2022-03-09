FROM python:3.9.7

RUN mkdir /app
WORKDIR /app

RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y git gcc g++ libgl1

RUN /usr/local/bin/python -m pip install --upgrade pip
COPY backend/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# # # caiman install
# RUN git clone https://github.com/flatironinstitute/CaImAn /CaImAn
# RUN cd /CaImAn && pip install -e .

# install conda
# FROM --platform=linux/amd64
RUN mkdir -p ~/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
RUN bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
RUN rm -rf ~/miniconda3/miniconda.sh

ENV PATH $PATH:/root/miniconda3/bin
RUN conda upgrade -y --all
RUN conda install -n base -c conda-forge mamba

RUN apt update && apt-get install -y libgl1-mesa-dev

COPY . .

# ENTRYPOINT ["cd /app/backend & python main.py"]