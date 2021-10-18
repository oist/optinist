FROM continuumio/miniconda3

RUN mkdir /app
WORKDIR /app

RUN conda config --set always_yes yes
RUN conda update --yes conda
RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y git gcc g++ libgl1

RUN git clone https://github.com/flatironinstitute/CaImAn

RUN cd ./CaImAn && conda env create -n caiman -f environment.yml && conda install --override-channels -c conda-forge -n caiman pip

RUN /bin/bash -c "cd ./CaImAn && source activate caiman && pip install ."
RUN /bin/bash -c "cd ./CaImAn && source activate caiman && python setup.py build_ext -i && python caimanmanager.py install --inplace"

RUN apt update && \
    apt-get install -y libgl1-mesa-dev

COPY requirements.txt ./
RUN /bin/bash -c "source activate caiman && pip install --no-cache-dir -r requirements.txt"


SHELL ["conda", "run", "-n", "caiman", "/bin/bash", "-c"]

COPY . .

ENTRYPOINT [ "conda", "run", "--no-capture-output", "-n", "caiman", "python", "main.py" ]