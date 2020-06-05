# Requirements and installation

## Requirements

System requirements:
* [Python](https://www.python.org/) 3
* Minimap2 version 2.17-r941 (or higher)
* pbmm2 version 1.1.0 (or higher)
* LAST version last-1060
* Longshot version 0.3.5 (or higher)
* Flye version 2.7 (or higher)
* Miniasm version 0.3-179 (or higher)
* Wtdbg2 version 2.5 (or higher)
* Tandem-genotypes version 1.5.0 (or higher)
* Samtools version 1.9 (or higher)
* Bcftools version 1.9 (or higher)
* Bgzip version 1.9 (or higher)
* Gzip version 1.5 (or higher)
* Snakemake version 5.7.0 (or higher)

The dependencies are listed above and can be installed easily by using the `config/main_env.yaml` file.

## Install all the dependencies
In the following, I shall give you a walkthrough on how to set up a conda environment with all dependencies installed.
First, you've to create a new conda environment from the `config/main_env.yaml` file by using the following command:

    conda env create -f config/main_env.yaml

Subsequently, you can activate the environment with

    conda activate long-read_pipeline

To see if your environment is set up correctly, check whether or not the newly created environment `long-read pipeline` is present in the list of conda environments. You can view this list by using the command:

    conda env list

If the environment `long-read pipeline` is present, then you're set.

**Note**, that this conda environment has to be activated each time you want to run the pipelines.


## Getting started

This guide will help to get you started so that the pipelines can be used without any problems.

* First, install `snakemake`. For easy installation of `snakemake`, enter the following line in the command-line:

    conda install -c bioconda snakemake

* Second, clone the whole project to your work directory.
* Third, change the last line in both `config/main_env.yaml` and `config/py2_env.yaml` so that it specifies the path where you want your dependencies located.
* Fourth, install all the dependencies on your machine. An easy way to do this is by using the `config/py2_env.yaml` file to set up a new conda environment. Note, that this environment has to be activated before you run the pipeline.
* Fifth, perform the test run by entering the line below. The test run is completed successfully when it terminates without errors.

    snakemake -j 100 --latency-wait 180 --cluster-config cluster.json --use-conda --cluster "sbatch -p {cluster.partition} --qos {cluster.qos} --mem={cluster.mem} -t {cluster.time} --ntasks {cluster.ntasks} -c {cluster.cpus-per-task}"

* Sixth, if no errors occured during the testrun, change the global variable `WORK_PATH`, defined in `common.py`, to the path where all your samples are stored.

* Seventh, change `SAMPLE_FILES` in the Snakefile to a sample_id of interest, or uncomment that line for automatic sample selection.

* Lastly, uncomment the lines in `rule all` in the Snakefile as desired, as you can select which pipeline to activate in this way.

When everything is set and the above steps are successfully completed, then the pipelines should be working.

For further help, see the Snakefile which is documented thoroughly.
Alternatively, you can enter

    snakemake help

to show the help menu.

## Paper

For more details, please see: [Developing and validating bioinformatic pipelines that enable assessment of VNTR expansions across haploid human genomes using noisy long reads](docs/paper.pdf) by N Radunovic.
