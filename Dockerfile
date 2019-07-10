FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/crisprvar pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-crisprvar-1.0dev/bin:$PATH

# --- Begin copied from https://github.com/pinellolab/CRISPResso2/blob/master/Dockerfile --- #
#install ms fonts
RUN echo "deb http://httpredir.debian.org/debian jessie main contrib" > /etc/apt/sources.list \
  && echo "deb http://security.debian.org/ jessie/updates main contrib" >> /etc/apt/sources.list \
  && echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && apt-get update \
  && apt-get install -y ttf-mscorefonts-installer \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/*
# --- End copied from https://github.com/pinellolab/CRISPResso2/blob/master/Dockerfile --- #

RUN cd /home && \
    git clone https://github.com/pinellolab/CRISPResso2.git && \
    cd CRISPResso2 && \
    python setup.py install

RUN apt update && apt install -y dos2unix
