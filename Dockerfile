# syntax=docker/dockerfile:1

FROM continuumio/anaconda3

WORKDIR /src

COPY environment.yml .

RUN conda env create -f environment.yml && conda clean -afy

SHELL ["conda", "run", "-n", "bioenv", "/bin/bash", "-c"]

COPY . .

RUN conda list

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bioenv", "perl", "run_metaGPA.pl"]
CMD ["--help"]