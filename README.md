# MetaGPA

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
  - [Using Conda/Mamba](#using-condamamba)
  - [Using Docker](#using-docker)
- [Usage](#usage)
  - [Quick Start](#quick-start)
  - [Command Line Options](#command-line-options)
  - [Example Workflow](#example-workflow)
- [Input/Output](#inputoutput)
- [Pipeline Structure](#pipeline-structure)
- [Configuration](#configuration)
- [Dependencies](#dependencies)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Overview
**MetaGPA (Metagenomics -GenomePhenome Association)** pipeline is used to link genetic information (e.g. genes encoding for enzymes) in metagenome with a dedicated functional phenotype (e.g. DNA cytosine modification).
The main MetaGPA pipeline identifies protein family domains (Pfam, from Pfam-A database, http://pfam.xfam.org/) which are significantly associated with the studied phenotype.
We also provide additional scripts for further analyses of any chosen interested Pfam identified from the main pipeline.
For more details, please refer to our paper: [A Genome-Phenome Association study in native microbiomes identifies a mechanism for cytosine modification in DNA and RNA.](https://elifesciences.org/articles/70021)

## Installation

Clone the package:

```bash
git clone https://github.com/ge-a/metaGPA.git
```

### Using Conda/Mamba
1. **Create the environment:**
    ```bash
    conda env create -f environment.yml
    # or, with mamba (faster):
    mamba env create -f environment.yml
    ```

2. **Activate the environment:**
    ```bash
    conda activate metaGPA
    ```

3. **Test the installation:**
    ```bash
    perl src/run_metaGPA.pl --help
    ```
### Using Docker
A Dockerfile is provided for containerized installation and reproducibility.
1. **Build the Docker image:**
    ```bash
    docker build -t metagpa .
    ```

2. **Run the pipeline in a container:**
    ```bash
    docker run --rm -it -v $(pwd):/data metagpa perl /data/src/run_metaGPA.pl --help
    ```

   - Replace `$(pwd)` with the path to your data if running from a different directory.
   - You can mount additional volumes or set environment variables as needed.

## Usage

### Quick Start
Here is a minimal example to get started with MetaGPA after installation:

1. **Prepare your input files:**  
   Ensure you have your control and selection FASTQ files and a configuration file (if needed).

2. **Run the pipeline:**  
   ```bash
   perl src/run_metaGPA.pl \
     --control_1 CONTROL_SAMPLE \
     --selection_1 SELECTION_SAMPLE \
     --output-dir results
   ```

   Replace `CONTROL_SAMPLE`, `SELECTION_SAMPLE`, and `EXPERIMENT_NAME` with your filenames or sample identifiers.

3. **Check the output:**  
   Results will be written to into the directory specified by `--output-dir` within the output dir (e.g., `output/results_phage/`).
---

### Command Line Options

You can see all available options by running:
```bash
perl src/run_metaGPA.pl --help
```

Common options include:
- `--control_1, -c1` : Name or path of the control 1 sample FASTQ file(s)
- `--selection_1, -s1` : Name or path of the selection 1 sample FASTQ file(s)
- `--output-dir, -o` : Directory to store results
- `--control_2, -c2` : (Optional) Name or path of the control 2 sample FASTQ file(s)
- `--selection_2, -s2` : (Optional) Name or path of the selection 2 sample FASTQ file(s)
- `--assembly, -A` : (Optional) Run assembly steps
- `--annotation, -AN` : (Optional) Run annotation steps
- `--mapping, -M` : (Optional) Run mapping steps
- `--enrichment, -E` : (Optional) Run enrichment analysis steps
- `--cutoff, -c` : (Optional) Enrichment cutoff value
- `--trim, -t` : (Optional) Do we want to trim files before running pipeline
- `--help` : Show help message

Note that if none of assembly, annotation, mapping, or enrichment are specified, then all steps will be ran. If a downstream step is ran without necessary preceding steps the program will assume the previous steps were ran in the previously output results and look for necessary data files accordingly.

Also note that while you may specify both control 1 and 2, as well as selection 1 and 2 when running the program, you may input only files for control 1 and selection 1 and the program will generate the file path for control 2 and selection 2 assuming they have the same base file name with a 2 instead of a 1 at the end.

---
### Example Workflow
Step-by-step example of a typical analysis.

## Input/Output
Description of required input files and generated output files.

## Pipeline Structure
Overview of the main scripts and their roles.

## Configuration
How to configure the pipeline for your data.

## Dependencies
List of required software and libraries.

## Troubleshooting
Common issues and solutions.

## Contact
How to get help or report issues.