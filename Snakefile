# our long-read analysis pipeline

localrules: all, clean, help

rule all:
    """
    Run the whole pipeline.
    """
    input:
        "/tudelft.net/staff-bulk/ewi/insy/VUMC/pacbio/1_A01/m64050_190903_135620.subreads.fastq"

rule help:
    """
    ----------------------------------------------------------------------------
                        Bio-informatica TGS Pipeline

    Deze pipeline, bestaande uit Python scripts en Bash commands, verkrijgt
    de consensussequentie en de variantenlijst (VCF) van meegegeven fasta
    bestanden, waaronder de pair-end reads en het referentiegenoom.

    Hieronder is van iedere rule de beschrijving (=help) weergegeven, inclusief
    input- en outputbestanden van de desbetreffende rule.

    Om de gehele pipeline te runnen, voer onderstaande in in de command-line:

        $ snakemake --snakefile <snakefile name>

    Voor de volledige help, voer onderstaande regel in in de command-line:

        $ snakemake --snakefile <snakefile name> help

    ----------------------------------------------------------------------------
    """
    run:
        for rule in workflow.rules:
          print(rule.name)
          print(rule.docstring)

rule clean:
    shell:
        """
        rm /tudelft.net/staff-bulk/ewi/insy/VUMC/pacbio/1_A01/outputfiles-analysis/my_first_test.fasta.gz
        """

rule test_1:
    input:
        a="/tudelft.net/staff-bulk/ewi/insy/VUMC/pacbio/1_A01/m64050_190903_135620.subreads.bam"
    output:
        "/home/nfs/nradunovic/output-analysis/test.txt"
    script:
        "/home/nfs/nradunovic/output-analysis/start.py"

rule bam2fasta:
    input:
        "/tudelft.net/staff-bulk/ewi/insy/VUMC/pacbio/1_A01/m64050_190903_135620.subreads.bam"
    output:
        "/tudelft.net/staff-bulk/ewi/insy/VUMC/pacbio/1_A01/outputfiles-analysis/my_first_test.fasta.gz"
    shell:
        "bam2fasta -o /tudelft.net/staff-bulk/ewi/insy/VUMC/pacbio/1_A01/outputfiles-analysis/my_first_test {input}"





