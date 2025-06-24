# syntax=docker/dockerfile:1

FROM mambaforge/mambaforge:latest

WORKDIR /src

COPY environment.yml .

RUN mamba env create -f environment.yml && mamba clean --all --yes

SHELL ["mamba", "run", "-n", "bioenv", "/bin/bash", "-c"]

COPY . .

RUN mamba run -n bioenv conda list

ENTRYPOINT ["mamba", "run", "--no-capture-output", "-n", "bioenv", "perl", "run_metaGPA.pl"]
CMD ["--help"]