FROM python:3.9.7

RUN mkdir /app
WORKDIR /app

RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y git gcc g++ libgl1

RUN /usr/local/bin/python -m pip install --upgrade pip
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# # caiman install
RUN git clone https://github.com/flatironinstitute/CaImAn /CaImAn
RUN cd /CaImAn && pip install -e .

RUN apt update && \
    apt-get install -y libgl1-mesa-dev

COPY . .

ENTRYPOINT ["python", "main.py"]