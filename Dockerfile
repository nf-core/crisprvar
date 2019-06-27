FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/crisprvar pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-crisprvar-1.0dev/bin:$PATH
RUN apt update && apt install -y dos2unix
