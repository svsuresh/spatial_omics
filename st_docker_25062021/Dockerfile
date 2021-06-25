FROM ubuntu:20.04

SHELL [ "/bin/bash", "--login", "-c" ]

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qqy \
    && apt-get install -qqy \
        apt-utils \
        wget \
        vim \
        sudo \
        git \
        build-essential

ARG MINICONDA_VERSION=latest
ARG USER=galaxy

RUN groupadd ${USER} \
    && useradd -m -g ${USER} ${USER} -d /home/${USER} -s /bin/bash -p "$(openssl passwd -1 galaxy)" \
    && usermod -aG sudo galaxy

WORKDIR /home/${USER}

USER ${USER}

RUN cd /tmp \
    && wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh \
    && bash Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -f -b -p /home/${USER}/.conda \
    && rm Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh 

ENV PATH /home/${USER}/.conda/bin:$PATH
COPY requirements.yml /tmp/requirements.yml
RUN conda env create -n "st_pipeline" python=3.6.3 -f /tmp/requirements.yml
RUN echo "source activate st_pipeline" > ~/.bashrc
ENV PATH /home/${USER}/.conda/envs/st_pipeline/bin:$PATH
RUN /bin/bash -c "source activate st_pipeline"
RUN conda install -c bioconda star samtools openssl -y
RUN pip install taggd
RUN git clone https://github.com/jfnavarro/st_pipeline.git \
	&& cd st_pipeline \
	&& python setup.py build \
	&& python setup.py install \
	&& python setup.py test

RUN echo $PATH

