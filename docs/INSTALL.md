# Requirements and installation

## Requirements

System requirements:
* [Python](https://www.python.org/) 3
* Minimap2 version 2.17-r941 (or higher)
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

The above dependencies are downloaded and installed automatically when the
pipeline is being used for the first time.

## Installation

* First, install `snakemake`. For easy installation of `snakemake`, enter the following line in the command-line:

    conda install -c bioconda snakemake

* Second, clone the whole project to your work directory.
* Third, change the last line in both `config/main_env.yaml` and `config/py2_env.yaml` so that it specifies the path where you want your dependencies located.
* Fourth, change global variable `WORK_PATH`, defined in `common.py`, to the path where all your samples are stored.
* Fifth, peform the testrun by entering the line below. The testrun is completed succesfully when it terminates without errors.

    snakemake -j 100 --latency-wait 180 --cluster-config cluster.json --use-conda --cluster "sbatch -p {cluster.partition} --qos {cluster.qos} --mem={cluster.mem} -t {cluster.time} --ntasks {cluster.ntasks} -c {cluster.cpus-per-task}"

* Sixth, when no errors occured when peforming the testrun, change `SAMPLE_FILES` in the Snakefile to a sample_id of interest, or uncomment that line for automatic sample selection.

* Lastly, uncomment the lines in `rule all` in the Snakefile as desired, as you can select which pipeline to activate in this way.

When everything is set and the above steps are succesfully completed, then the pipelines should be working.

For further help, see the Snakefile which is documented thoroughly.
Alternatively, you can enter

    snakemake help

to show the help menu.

## Paper

For more details, please see: [Developing and validating bioinformatic pipelines that enable assessment of VNTR expansions across haploid human genomes using noisy long reads](docs/paper.pdf) by N Radunovic.
