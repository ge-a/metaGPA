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
The main MetaGPA pipeline identifies protein family domains (Pfam, from Pfam-A database, http://pfam.xfam.org/ and or TIGRFAM) which are significantly associated with the studied phenotype.
We also provide additional scripts for further analyses of any chosen interested Pfam identified from the main pipeline.
For more details, please refer to our paper: [A Genome-Phenome Association study in native microbiomes identifies a mechanism for cytosine modification in DNA and RNA.](https://elifesciences.org/articles/70021)

## Installation
Clone the repository using one of the following methods:

### Option 1: HTTPS 
```bash
git clone https://github.com/ge-a/metaGPA.git
```

### Option 1: SSH

```bash
git clone git@github.com:ge-a/metaGPA.git
```

### Using Conda
1. **Create the environment:**
    ```bash
    conda env create -f environment.yml
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
Please note that this workflow will create a directory called "output" which stores all the outputted data files from each step of the pipeline. The directory name you pass into the output-dir option will be stored inside of this generated "output" directory. When passing in your output directories please do not specify a path but simply the name of the directory you would like to create as the program begins looking from inside this generated "output" dir. For example, if you would like to create an output-dir "PF4", but pass in "output/PF4" to the output-dir flag, the folder that is generated will be stored at relative path output/output/PF4, as the program takes into account the output dir automatically.
### Quick Start
Here is a minimal example to get started with MetaGPA after installation:

1. **Prepare your input files:**  
   Ensure you have your control and selection FASTQ files.

2. **Run the pipeline:**  
   ```bash
   perl run_metaGPA.pl \
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
perl run_metaGPA.pl --help
```

Options include:
- `--control_1, -c1` :  \[Required\] Name or path of the control 1 sample FASTQ file
- `--selection_1, -s1` : \[Required\] Name or path of the selection 1 sample FASTQ file
- `--output-dir, -o` : \[Required\] Directory to store results
- `--control_2, -c2` : (Optional) Name or path of the control 2 sample FASTQ file
- `--selection_2, -s2` : (Optional) Name or path of the selection 2 sample FASTQ file
- `--assembly, -A` : (Optional) Run assembly steps
- `--annotation, -AN` : (Optional) Run annotation steps
- `--mapping, -M` : (Optional) Run mapping steps
- `--enrichment, -E` : (Optional) Run enrichment analysis steps
- `--cutoff, -c` : (Optional) Enrichment cutoff value
- `--trim, -t` : (Optional) Do we want to trim files before running pipeline
- `--threads` :  (Optional) Number of threads to run assembly with
- `--memory` :  (Optional) Memory allocation for assembly
- `--help` : Show help message

Note that if none of assembly, annotation, mapping, or enrichment are specified, then all steps will be ran. If a downstream step is ran without necessary preceding steps the program will assume the previous steps were ran in the previously output results and look for necessary data files accordingly.

Also note that while you may specify both control 1 and 2, as well as selection 1 and 2 when running the program, you may input only files for control 1 and selection 1 and the program will generate the file path for control 2 and selection 2 assuming they have the same base file name with a 2 instead of a 1 at the end.

### MultiSelection Usage

```bash
perl run_multitreatment.pl --help
```

Options include:
- `--control_1, -c1` : \[Required\] Name or path of the control 1 sample FASTQ file
- `--selection_1, -s1` : \[Required\] Name or path of the selection 1 sample FASTQ file(s)
- `--output-dir, -o` : \[Required\] Directory to store results
- `--control_2, -c2` : (Optional) Name or path of the control 2 sample FASTQ file(s)
- `--selection_2, -s2` : (Optional) Name or path of the selection 2 sample FASTQ file(s)
- `--assembly, -A` : (Optional) Run assembly steps
- `--annotation, -AN` : (Optional) Run annotation steps
- `--mapping, -M` : (Optional) Run mapping steps
- `--enrichment, -E` : (Optional) Run enrichment analysis steps
- `--cutoff, -c` : (Optional) Enrichment cutoff value
- `--trim, -t` : (Optional) Do we want to trim files before running pipeline
- `--threads` :  (Optional) Number of threads to run assembly with
- `--memory` :  (Optional) Memory allocation for assembly
- `--help` : Show help message

Note that to run multiple selections for the control you can simply write the s1 flag again and pass in the next selection sample you would like to run.

### Postpipeline Analysis Function Usage

```bash
perl run_analysis.pl ---command analysis_function \[options for analysis script\]
```
If you would like to see the help output for the the analysis function you can run this command
```bash
perl run_analysis.pl ---command analysis_function --subhelp
```
You can see all run options with command:
```bash
perl run_analysis.pl --help
```

Options include:
- `--control_1, -c1` : \[Required\]  Analysis type (e.g., create_tree_file, get_gff)
- `--subhelp` :Show help message for your specified command
- `--help` : Show help message

Note that to run multiple selections for the control you can simply write the s1 flag again and pass in the next selection sample you would like to run.

Available commands:
- `create_tree_file`
- `get_enriched_fasta_orf`
- `get_flanking_for_contig_vis`
- `get_gff`
- `parse_tree`
- `residue_analysis`


---
## Input/Output
### Input Files

FQ or GZ fasta files or compressed fasta files

### Output Files

#### Assembly
- Combined Control and Selected assembly fasta file
- Combined Control and Selected faa file
#### Annotation
- TIGRFAM Hit tab file
- PFAM Hit tab file
#### Mapping
- Control and Selected Reads bed file
#### Enrichmnet
- Enriched reads in Control and Selected groups txt file
- Enrichment Annotated PFAM file
- Enrichment Annotated TIGRFAM file

## Pipeline Structure

### Assembly
### Mapping
### Annotation
### Enrichment

## Dependencies
Please see environment.yml file for a list of required dependencies and versions. These dependencies can be installed directly with command: 
`conda env create -f environment.yml`

## Troubleshooting
Common issues and solutions.

## Contact
How to get help or report issues.
