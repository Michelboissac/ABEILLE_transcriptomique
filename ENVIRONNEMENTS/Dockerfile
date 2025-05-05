FROM continuumio/miniconda3

COPY environment.yml .

RUN conda env update -n base -f environment.yml
ENV PATH=/opt/conda/bin:$PATH
