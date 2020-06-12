from common import *
from multiprocessing.pool import ThreadPool as Pool
import concurrent.futures
import os

PY2_ENV = "config/py2_env.yaml" # Conda python 2 environment (don't touch)
SAMPLE_FILES = ["sample_for_testrun_CLR"] # COMMENT THIS LINE OUT WHEN NOT TESTING. PUT IN HERE YOUR SAMPLE_ID FOR TESTING.

# other samples:
# - m64050_190904_071400 -> Our centenarian samples
# - m64013_190124_221354 -> open-source PacBio dataset CLR
# - m64014_181209_091052 -> open-source PacBio dataset CCS

CHROMS = ["chr" + str(x) for x in range(23) if x != 0] # This is used for phasing with Longshot
CHROMS.append("chrX")
CHROMS.append("chrY")

wildcard_constraints:
    sample_id="[\w\d_]+",
    ref_genome="[\w\d_.]+"

localrules: all, clean, help


rule all:
    input:
        # THIS IS THE TESTRUN
        expand(os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/{samplefile}.aligned.sorted.bam'), samplefile=SAMPLE_FILES)

        # Uncomment for running the REFERENCE BASED approach of the HAPLOTYPE-RECONSTRUCTION pipeline
        #expand(os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{samplefile}/wtdbg2_asm_haplotypes_second_run/haplotype_2/{samplefile}.hap2.3nd_polishing_iteration.cns.fa'), samplefile=SAMPLE_FILES)

        # Uncomment for running the DE NOVO BASED approach of the HAPLOTYPE-RECONSTRUCTION pipeline
        #expand(os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{samplefile}/haplotype_2/assembly.fasta'), samplefile=SAMPLE_FILES)

        # Uncomment for running the STRUCTURAL VARIANT CALLING pipeline (not fully, properly working yet)
        #expand(os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/sv_calling_on_haplotypes/de_novo_asm_based_haplotypes/{samplefile}/haplotype_2/asm_haplotype_2.prg.refined.vcf'), samplefile=SAMPLE_FILES)


rule help:
    """
    To run the whole pipeline, enter the line below:

        $ snakemake -j 100 --latency-wait 180 --cluster-config cluster.json
        --use-conda --cluster "sbatch -p {cluster.partition} --qos {cluster.qos}
        --mem={cluster.mem} -t {cluster.time} --ntasks {cluster.ntasks} -c
        {cluster.cpus-per-task}"

    For help, enter the line below:

        $ snakemake [specific rule] help
    """
    run:
        for rule in workflow.rules:
          print(rule.name)
          print(rule.docstring)

rule clean:
    shell:
        """
        rm slurm*
        """


rule retrieve_files:
    """
    For manual bam-file retrieval. Moves every bam file contained in the
    sample_folder to /reads.
    """
    input:
        os.path.join(str(SCRIPTS_PATH), 'retrieve_files.sh')
    shell:
        """
        echo "Samples from which the subreads succesfully are retrieved\n" {SAMPLE_FILES}
        """


######## General rules #########################################################
################################################################################

rule bam2fasta:
    input:
        os.path.join(str(READS_DIR), '{sample_id}.subreads.bam')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/{sample_id}.subreads.fasta')
    shell:
        """
        samtools fasta {input} > {output}
        """

rule bam2fastq:
    input:
        os.path.join(str(READS_DIR), '{sample_id}.subreads.bam')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fastq_files/{sample_id}.subreads.fq')
    shell:
        """
        samtools fastq {input} > {output}
        """


rule index_reference:
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/{ref_genome}.fa')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/{ref_genome}.mmi'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/{ref_genome}.fai')
    shell:
        """
        pbmm2 index --preset 'SUBREAD' -k 19 -w 10 {input} {output[0]}
        samtools faidx {input} -o {output[1]}
        """


rule test_run:
    input:
        'test/{sample_id}.subreads.bam'
    output:
        os.path.join(str(READS_DIR), '{sample_id}.subreads.bam')
    shell:
        """
        cp {input} {output}
        """

rule retrieve_reference_genome:
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/hg38.fa')
    shell:
        """
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gunzip hg38.fa.gz
        mv hg38.fa {output}
        """


######## de-novo-based-approach, haplotype-reconstruction pipeline #############
################################################################################

rule bam2fasta_hap_de_novo_asm_based_pipeline:
    """
    Obtains fasta files from haplotype-specific bam files. Merges the unassigned
    reads with the reads of both haplotype 1 as haplotype 2.
    """
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.hap1.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.hap2.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.unassigned.bam')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_1.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_2.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.subreads_unassigned.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_1_incl_unassigned.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_2_incl_unassigned.fasta')
    shell:
        """
        samtools fasta {input[0]} > {output[0]}
        samtools fasta {input[1]} > {output[1]}
        samtools fasta {input[2]} > {output[2]}
        cat {input[0]} {output[2]} > {input[3]}
        cat {input[1]} {output[2]} > {input[4]}
        """


rule asm_to_indexed_asm:
    """
    Create an index file for the obtained diploid genome assembly.
    """
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/diploid_genome/assembly.fasta')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/diploid_genome/assembly.fasta.fai')
    shell:
        """
        samtools faidx {input} -o {output}
        """


rule map_to_asm: # Replace '--preset SUBREAD' with '--preset CCS' when analysing CCS/Hifi reads
    """
    Mapping the subreads against the obtained diploid genome assembly.
    """
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/{sample_id}.subreads.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/diploid_genome/assembly.fasta')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_asm/{sample_id}/{sample_id}.aligned.subreads_to_asm.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_asm/{sample_id}/{sample_id}.aligned.subreads_to_asm.sorted.bam')
    run:
        cmd = """
        pbmm2 align --preset SUBREADS --log-level DEBUG --unmapped {input[1]} {input[0]} {output[0]}
        samtools sort -o {output[1]}  {output[0]}
        samtools index -b {output[1]} {output[1]}.bai
        samtools stats {output[1]} > {output[1]}.stats
        """
        shell(cmd)


rule de_novo_asm_haplotypes_mapped_to_ref: # POSSIBLE EXTENTION: add option '-L FILE' to samtools view in order to only output alignments overlapping the input BED FILE [null].
    """
    Mapping the subreads against the obtained diploid genome assembly.
    """
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/hg38.fa'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/{sample_id}/haplotype_1/hap1_asm.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/{sample_id}/haplotype_2/hap2_asm.fasta')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/de_novo_asm_haplotypes_mapped_back_to_ref/{sample_id}/haplotype_1/{sample_id}.aligned.subreads_to_asm.sam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/de_novo_asm_haplotypes_mapped_back_to_ref/{sample_id}/haplotype_1/{sample_id}.aligned.subreads_to_asm.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/de_novo_asm_haplotypes_mapped_back_to_ref/{sample_id}/haplotype_1/{sample_id}.aligned.subreads_to_asm.sorted.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/de_novo_asm_haplotypes_mapped_back_to_ref/{sample_id}/haplotype_2/{sample_id}.aligned.subreads_to_asm.sam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/de_novo_asm_haplotypes_mapped_back_to_ref/{sample_id}/haplotype_2/{sample_id}.aligned.subreads_to_asm.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/de_novo_asm_haplotypes_mapped_back_to_ref/{sample_id}/haplotype_2/{sample_id}.aligned.subreads_to_asm.sorted.bam')
    run:
        threads_for_mapping = 19 # minimap2 uses INT + 1 threads
        cmd = """
        minimap2 --secondary=no -L -ax asm20 -t {threads_for_mapping} {input[0]} {input[1]} > {output[0]}
        samtools view -h -S -b -T {input[0]} {output[0]} > {output[1]}
        samtools sort -o {output[2]}  {output[1]}
        samtools index -b {output[2]} {output[2]}.bai

        minimap2 --secondary=no -L -ax asm20 -t {threads_for_mapping} {input[0]} {input[2]} > {output[3]}
        samtools view -h -S -b -T {input[0]} {output[3]} > {output[4]}
        samtools sort -o {output[5]}  {output[4]}
        samtools index -b {output[5]} {output[5]}.bai
        """
        shell(cmd)


def multi_threading_phasing_de_novo_based_pipeline(contig_and_id):
    contig = contig_and_id[0][1:]
    sample_id = contig_and_id[1]
    path_to_phased_bam = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/de_novo_asm_based_pipeline/' + sample_id)
    path_to_vcf_file = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/de_novo_asm_based_pipeline/' + sample_id)
    alignment = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_asm/' + sample_id + '/' + sample_id + '.aligned.subreads_to_asm.sorted.bam')
    asm = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/' + sample_id + '/assembly.fasta')
    cmd = """
    longshot -F -r {contig} -p {path_to_phased_bam}/{sample_id}.{contig} --bam {alignment} --ref {asm} --out {path_to_vcf_file}/{sample_id}.{contig}.vcf
    """
    shell(cmd)
    return "finished!"


rule phasing_de_novo_based_pipeline: # (bash-error: doesn't find headers.txt file)? TO-DO: vcf-files mergen en {output[0]} genereren.
    """
    Haplotype phasing the reads that mapped to the intermediate diploid genome
    assembly that was obtained in the previous step. Obtaining a bam file with
    reads corresponding to either haplotype 1, haplotype 2 or unassigned.
    """
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_asm/{sample_id}/{sample_id}.aligned.subreads_to_asm.sorted.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/diploid_genome/assembly.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/diploid_genome/assembly.fasta.fai')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.longshot.vcf'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.hap1.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.hap2.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.unassigned.bam')
    params:
        path_to_headers=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/headers.txt')
    run:
        path_to_phased_bam = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/de_novo_asm_based_pipeline/' + str(wildcards.sample_id))
        path_to_vcf_file = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/de_novo_asm_based_pipeline/' + str(wildcards.sample_id))
        cmd = """
        cat {input[1]} | grep '>' > {params.path_to_headers}
        """
        #shell(cmd)
        with open(params.path_to_headers) as f:
            content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        content = [x.rstrip() for x in content]
        executor = concurrent.futures.ProcessPoolExecutor(70)
        futures = [executor.submit(multi_threading_phasing_de_novo_based_pipeline, [contig, str(wildcards.sample_id)]) for contig in content]
        concurrent.futures.wait(futures)
        cmd = """

        touch {output[0]}

        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.hap1_1.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_1*.hap1.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_1*.hap1.bam
        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.hap2_1.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_1*.hap2.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_1*.hap2.bam
        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.unassigned_1.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_1*.unassigned.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_1*.unassigned.bam

        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.hap1_2.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_2*.hap1.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_2*.hap1.bam
        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.hap2_2.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_2*.hap2.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_2*.hap2.bam
        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.unassigned_2.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_2*.unassigned.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_2*.unassigned.bam

        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.hap1_3.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_3*.hap1.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_3*.hap1.bam
        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.hap2_3.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_3*.hap2.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_3*.hap2.bam
        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.unassigned_3.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_3*.unassigned.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_3*.unassigned.bam

        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.hap1_rest.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_[4-9]*.hap1.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_[4-9]*.hap1.bam
        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.hap2_rest.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_[4-9]*.hap2.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_[4-9]*.hap2.bam
        samtools merge -f {path_to_phased_bam}/{wildcards.sample_id}.unassigned_rest.bam --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.contig_[4-9]*.unassigned.bam {path_to_phased_bam}/{wildcards.sample_id}.scaffold_[4-9]*.unassigned.bam

        samtools merge -f {output[1]} --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.hap1_*
        samtools merge -f {output[2]} --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.hap2_*
        samtools merge -f {output[3]} --threads 70 {path_to_phased_bam}/{wildcards.sample_id}.unassigned_*

        rm {path_to_phased_bam}/{wildcards.sample_id}.hap1_* {path_to_phased_bam}/{wildcards.sample_id}.hap2_* {path_to_phased_bam}/{wildcards.sample_id}.unassigned_*
        """
        shell(cmd)
        cmd = """
        mkdir {path_to_phased_bam}/backup
        mv {path_to_phased_bam}/{wildcards.sample_id}.contig*.bam {path_to_phased_bam}/backup/
        mv {path_to_phased_bam}/{wildcards.sample_id}.scaffold*.bam {path_to_phased_bam}/backup/
        mv {path_to_vcf_file}/{wildcards.sample_id}.contig*.vcf* {path_to_vcf_file}/backup/
        mv {path_to_vcf_file}/{wildcards.sample_id}.scaffold*.vcf* {path_to_vcf_file}/backup/
        """
        shell(cmd)
        cmd = """
        rm {path_to_phased_bam}/{wildcards.sample_id}.ctg*.bam {path_to_vcf_file}/{wildcards.sample_id}.(contig|scaffold)*.vcf*
        """
        #shell(cmd) # Uncomment for automatic removal of intermediate bam- and vcf-files


rule flye_asm_de_novo_based: # replace '--pacbio-raw' for '--pacbio-hifi' when analysing hifi/CCS reads
    """
    De novo assembly on the subreads, obtaining a diploid genome assembly that
    represents the sample genome.
    """
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/{sample_id}.subreads.fasta')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/diploid_genome/assembly.fasta')
    params:
        path_to_out_dir1=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/diploid_genome')
    run:
        genome_size = 3000000000 # the human genome has 3 billion base pairs
        polishing_iterations = 3 # peforming three iterations is conventional
        number_of_parallel_threads = 36
        cmd = """
        flye --pacbio-raw {input} --out-dir {params.path_to_out_dir1} --genome-size {genome_size} --iterations {polishing_iterations} --threads {number_of_parallel_threads}
        """
        shell(cmd)


rule flye_asm_de_novo_based_haplotypes: # replace '--pacbio-raw' for '--pacbio-hifi' when analysing hifi/CCS reads
    """
    De novo assembly on the subreads corresponding to either of both haplotypes,
    obtaining a haplotype-specific assemblies.
    """
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_1_incl_unassigned.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/de_novo_asm_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_2_incl_unassigned.fasta')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/haplotype_1/assembly.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/haplotype_2/assembly.fasta')
    params:
        path_to_out_dir_hap1=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/haplotype_1'),
        path_to_out_dir_hap2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/de_novo_asm_based_pipeline/{sample_id}/haplotype_2')
    run:
        genome_size = 3000000000 # the human genome has 3 billion base pairs
        polishing_iterations = 3 # peforming three iterations is conventional
        number_of_parallel_threads = 36
        cmd = """
        flye --pacbio-raw {input[0]} --out-dir {params.path_to_out_dir_hap1} --genome-size {genome_size} --iterations {polishing_iterations} --threads {number_of_parallel_threads}
        flye --pacbio-raw {input[1]} --out-dir {params.path_to_out_dir_hap2} --genome-size {genome_size} --iterations {polishing_iterations} --threads {number_of_parallel_threads}
        """
        shell(cmd)


######## reference-based approach, haplotype-reconstruction pipeline ###########
################################################################################

rule map_to_reference:
    input:
        os.path.join(str(READS_DIR), '{sample_id}.subreads.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/hg38.mmi')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/{sample_id}.aligned.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/{sample_id}.aligned.sorted.bam')
    #    os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/{sample_id}.aligned.sorted.bam.bai'),
    #    os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/{sample_id}.aligned.sorted.bam.stats')
    shell: # pbmm2 uses the maximum number of threads as default
        """
        pbmm2 align --log-level DEBUG --unmapped {input[1]} {input[0]} {output[0]}
        samtools sort -o {output[1]}  {output[0]}
        samtools index -b {output[1]} {output[1]}.bai
        samtools stats {output[1]} > {output[1]}.stats
        """


rule get_unmapped:
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/{sample_id}.aligned.sorted.bam')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/{sample_id}.unmapped.bam')
    shell:
        """
        samtools view -b -f 4 {input} > {output}
        """


def multi_threading_phasing_reference_based_pipeline(contig_and_id):
    contig = contig_and_id[0]
    sample_id = contig_and_id[1]
    path_to_phased_bam = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/' + sample_id)
    path_to_vcf_file = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/reference_based_pipeline/' + sample_id)
    alignment = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/' + sample_id + '.aligned.sorted.bam')
    reference = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/hg38.fa')
    cmd = """
    longshot -F -r {contig} -p {path_to_phased_bam}/{sample_id}.{contig} --bam {alignment} --ref {reference} --out {path_to_vcf_file}/{sample_id}.{contig}.vcf
    bgzip {path_to_vcf_file}/{sample_id}.{contig}.vcf
    bcftools index {path_to_vcf_file}/{sample_id}.{contig}.vcf.gz
    """
    shell(cmd)
    return "finished!"


rule phasing_reference_based_pipeline:
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/{sample_id}.aligned.sorted.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/hg38.fa'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/hg38.fa.fai')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/reference_based_pipeline/{sample_id}/{sample_id}.longshot.vcf'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}/{sample_id}.hap1.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}/{sample_id}.hap2.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}/{sample_id}.unassigned.bam')
    run:
        path_to_phased_bam = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/' + str(wildcards.sample_id))
        path_to_vcf_file = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/reference_based_pipeline/' + str(wildcards.sample_id))
        executor = concurrent.futures.ProcessPoolExecutor(24)
        futures = [executor.submit(multi_threading_phasing_reference_based_pipeline, [contig, str(wildcards.sample_id)]) for contig in CHROMS]
        concurrent.futures.wait(futures)
        cmd = """
        samtools merge {output[1]} {path_to_phased_bam}/{wildcards.sample_id}.chr*.hap1.bam
        samtools merge {output[2]} {path_to_phased_bam}/{wildcards.sample_id}.chr*.hap2.bam
        samtools merge {output[3]} {path_to_phased_bam}/{wildcards.sample_id}.chr*.unassigned.bam
        bcftools merge --force-samples {path_to_vcf_file}/{wildcards.sample_id}.chr*.vcf.gz -o {output[0]}
        """
        shell(cmd)
        cmd = """
        rm {path_to_phased_bam}/{wildcards.sample_id}.chr*.bam {path_to_vcf_file}/{wildcards.sample_id}.chr*.vcf*
        """
        #shell(cmd) # Uncommand for automatic removal of intermediate bam- and vcf-files


rule bam2fasta_hap_reference_based_pipeline:
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}/{sample_id}.hap1.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}/{sample_id}.hap2.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}/{sample_id}.unassigned.bam')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_1.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_2.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_unassigned.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_1_incl_unassigned.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_2_incl_unassigned.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_1_incl_unassigned.fasta.gz'), # optional
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_2_incl_unassigned.fasta.gz') # optional
    shell:
        """
        samtools fasta {input[0]} > {output[0]}
        samtools fasta {input[1]} > {output[1]}
        samtools fasta {input[2]} > {output[2]}
        cat {output[0]} {output[2]} > {output[3]}
        cat {output[1]} {output[2]} > {output[4]}
        gzip -c {output[3]} > {output[5]}
        gzip -c {output[4]} > {output[6]}
        """


rule bam2fastq_hap_reference_based_pipeline:
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}.hap1.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}.hap2.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}.unassigned.bam')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fastq_files/reference_based_pipeline/{sample_id}.hap1.fq'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fastq_files/reference_based_pipeline/{sample_id}.hap2.fq'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fastq_files/reference_based_pipeline/{sample_id}.unassigned.fq'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fastq_files/reference_based_pipeline/{sample_id}.hap1_incl_unassigned.fq'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fastq_files/reference_based_pipeline/{sample_id}.hap2_incl_unassigned.fq')
    shell:
        """
        samtools fastq {input[0]} > {output[0]}
        samtools fastq {input[1]} > {output[1]}
        samtools fastq {input[2]} > {output[2]}
        cat {output[0]} {output[2]} > {output[3]}
        cat {output[1]} {output[2]} > {output[4]}
        """


rule bam2fasta_hap:
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}/{sample_id}.hap1.bam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'phased/reference_based_pipeline/{sample_id}/{sample_id}.hap2.bam')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.hap1.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.hap2.fasta')
    shell:
        """
        samtools fasta {input[0]} > {output[0]}
        samtools fasta {input[1]} > {output[1]}
        """


def multi_threading_map_to_hap(contig_and_id):
    contig = contig_and_id[0].split(' ', 1)[0][1:]
    sample_id = contig_and_id[1]
    alignment = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/' + sample_id + '/' + sample_id + '.aligned.subreads_to_hap')
    cmd = """
    samtools view -b {alignment}.sorted.bam "{contig}" > "{alignment}.{contig}.bam"
    samtools fasta -F 0x900 "{alignment}.{contig}.bam" > "{alignment}.{contig}.fasta"
    """
    shell(cmd)
    return "finished!"


rule map_to_hap: # In deze rule zat oorspronkelijk een bug. TO-DO: checken of de bug nu is verholpen!
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/{sample_id}.subreads.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/{sample_id}.wtdbg2.merged_haplotypes.cns.fa')
    output:
        #os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/{sample_id}/{sample_id}.aligned.subreads_to_hap.sam'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/{sample_id}/{sample_id}.extracted_subreads_after_map_to_hap.hap1.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/{sample_id}/{sample_id}.extracted_subreads_after_map_to_hap.hap2.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/{sample_id}/{sample_id}.extracted_subreads_after_map_to_hap.hap1.fasta.gz'), # possibly unnecessary
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/{sample_id}/{sample_id}.extracted_subreads_after_map_to_hap.hap2.fasta.gz') # possibly unnecessary
    params:
        path_to_headers=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/headers.txt'),
        split=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/{sample_id}/{sample_id}.aligned.subreads_to_hap')
    run:
        alignment = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes')
        threads_for_mapping = 49 # minimap2 uses INT + 1 threads
        cmd = """
        #minimap2 -L -ax map-pb --secondary=no -t {threads_for_mapping} {input[1]} {input[0]} > {params.split}.sam
        #samtools faidx {input[1]}
        #samtools view -Sbt {input[1]}.fai {params.split}.sam > {params.split}.bam
        #samtools sort -o {params.split}.sorted.bam {params.split}.bam
        #samtools index -b {params.split}.sorted.bam {params.split}.sorted.bam.bai
        cat {input[1]} | grep '>' > {params.path_to_headers}
        """
        shell(cmd)
        with open(params.path_to_headers) as f:
            content = f.readlines()
        content = [x.rstrip() for x in content]
        executor = concurrent.futures.ProcessPoolExecutor(50)
        futures = [executor.submit(multi_threading_map_to_hap, [contig, str(wildcards.sample_id)]) for contig in content]
        concurrent.futures.wait(futures)
        cmd = """
        cat {params.split}.hap1*.fasta > {output[1]}
        cat {params.split}.hap2*.fasta > {output[2]}
        gzip -c {output[1]} > {output[3]}
        gzip -c {output[2]} > {output[4]}
        rm {params.split}.hap1\|ctg* {params.split}.hap2\|ctg* {params.path_to_headers}
        """
        #shell(cmd) # TO-DO: original
        cmd = """
        cat {params.split}.hap1*.fasta > {output[0]}
        cat {params.split}.hap2*.fasta > {output[1]}
        gzip -c {output[0]} > {output[2]}
        gzip -c {output[1]} > {output[3]}
        #rm {params.split}.hap1\|ctg* {params.split}.hap2\|ctg* {params.path_to_headers}
        """
        shell(cmd)


rule wtdbg_asm_haplotypes_reference_based_first_run:
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_1_incl_unassigned.fasta.gz'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'fasta_files/reference_based_pipeline/{sample_id}/{sample_id}.subreads_haplotype_2_incl_unassigned.fasta.gz')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/haplotype_1/{sample_id}.hap1.renamed.3nd_polishing_iteration.cns.fa'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/haplotype_2/{sample_id}.hap2.renamed.3nd_polishing_iteration.cns.fa'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/{sample_id}.wtdbg2.merged_haplotypes.cns.fa')
    params:
        folder_hap1=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/haplotype_1'),
        folder_hap2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/haplotype_2'),
        prefix_haplotype_1=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/haplotype_1/{sample_id}.hap1'),
        prefix_haplotype_2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_first_run/haplotype_2/{sample_id}.hap2'),
        wtdbg2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'wtdbg2/wtdbg2'),
        wtpoa_cns=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'wtdbg2/wtpoa-cns')
    run:
        number_of_parallel_threads = 24
        number_of_parallel_threads_for_sorting = 8
        cmd = """
        cd {params.folder_hap1}
        ### Haplotype 1 ###
        # assemble long reads
        {params.wtdbg2} -x sequel -g 3g -t {number_of_parallel_threads} -i {input[0]} -fo {params.prefix_haplotype_1}

        # derive consensus
        {params.wtpoa_cns} -t {number_of_parallel_threads} -i {params.prefix_haplotype_1}.ctg.lay.gz -fo {params.prefix_haplotype_1}.raw.fa

        # 1st polishing iteration, not necessary if you want to polish the assemblies using other tools
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_1}.raw.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_1}.1st_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_1}.1st_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_1}.raw.fa -i - -fo {params.prefix_haplotype_1}.1st_polishing_iteration.cns.fa

        # 2nd polishing iteration
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_1}.1st_polishing_iteration.cns.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_1}.2nd_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_1}.2nd_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_1}.1st_polishing_iteration.cns.fa -i - -fo {params.prefix_haplotype_1}.2nd_polishing_iteration.cns.fa

        # 3nd polishing iteration
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_1}.2nd_polishing_iteration.cns.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_1}.3nd_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_1}.3nd_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_1}.2nd_polishing_iteration.cns.fa -i - -fo {params.prefix_haplotype_1}.3nd_polishing_iteration.cns.fa

        # rename headers (set identifier hap1)
        sed -E 's/(^>)(ctg*)/\1hap1|\2/' {params.prefix_haplotype_1}.3nd_polishing_iteration.cns.fa > {params.prefix_haplotype_1}.renamed.3nd_polishing_iteration.cns.fa

        cd {params.folder_hap2}
        ### Haplotype 2 ###
        # assemble long reads
        {params.wtdbg2} -x sequel -g 3g -t {number_of_parallel_threads} -i {input[0]} -fo {params.prefix_haplotype_2}

        # derive consensus
        {params.wtpoa_cns} -t {number_of_parallel_threads} -i {params.prefix_haplotype_2}.ctg.lay.gz -fo {params.prefix_haplotype_2}.raw.fa

        # 1st polishing iteration, not necessary if you want to polish the assemblies using other tools
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_2}.raw.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_2}.1st_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_2}.1st_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_2}.raw.fa -i - -fo {params.prefix_haplotype_2}.1st_polishing_iteration.cns.fa

        # 2nd polishing iteration
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_2}.1st_polishing_iteration.cns.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_2}.2nd_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_2}.2nd_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_2}.1st_polishing_iteration.cns.fa -i - -fo {params.prefix_haplotype_2}.2nd_polishing_iteration.cns.fa

        # 3nd polishing iteration
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_2}.2nd_polishing_iteration.cns.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_2}.3nd_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_2}.3nd_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_2}.2nd_polishing_iteration.cns.fa -i - -fo {params.prefix_haplotype_2}.3nd_polishing_iteration.cns.fa

        # rename headers (set identifier hap2)
        sed -E 's/(^>)(ctg*)/\1hap2|\2/' {params.prefix_haplotype_2}.3nd_polishing_iteration.cns.fa > {params.prefix_haplotype_2}.renamed.3nd_polishing_iteration.cns.fa

        # merge haplotypes
        cat {params.prefix_haplotype_1}.renamed.3nd_polishing_iteration.cns.fa {params.prefix_haplotype_2}.renamed.3nd_polishing_iteration.cns.fa > {output[2]}
        """
        shell(cmd)


rule wtdbg_asm_haplotypes_reference_based_second_run:
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/{sample_id}/{sample_id}.extracted_subreads_after_map_to_hap.hap1.fasta.gz'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'alignments/subreads_mapped_back_to_haplotypes/{sample_id}/{sample_id}.extracted_subreads_after_map_to_hap.hap2.fasta.gz')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_second_run/haplotype_1/{sample_id}.hap1.3nd_polishing_iteration.cns.fa'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_second_run/haplotype_2/{sample_id}.hap2.3nd_polishing_iteration.cns.fa')
    params:
        folder_hap1=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_second_run/haplotype_1'),
        folder_hap2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_second_run/haplotype_2'),
        prefix_haplotype_1=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_second_run/haplotype_1/{sample_id}.hap1'),
        prefix_haplotype_2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/reference_based_pipeline/{sample_id}/wtdbg2_asm_haplotypes_second_run/haplotype_2/{sample_id}.hap2'),
        wtdbg2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'wtdbg2/./wtdbg2'),
        wtpoa_cns=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'wtdbg2/./wtpoa-cns')
    run:
        number_of_parallel_threads = 24
        number_of_parallel_threads_for_sorting = 8
        cmd = """
        cd {params.folder_hap1}
        ### Haplotype 1 ###
        # assemble long reads
        {params.wtdbg2} -x sequel -g 3g -t {number_of_parallel_threads} -i {input[0]} -fo {params.prefix_haplotype_1}

        # derive consensus
        {params.wtpoa_cns} -t {number_of_parallel_threads} -i {params.prefix_haplotype_1}.ctg.lay.gz -fo {params.prefix_haplotype_1}.raw.fa

        # 1st polishing iteration, not necessary if you want to polish the assemblies using other tools
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_1}.raw.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_1}.1st_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_1}.1st_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_1}.raw.fa -i - -fo {params.prefix_haplotype_1}.1st_polishing_iteration.cns.fa

        # 2nd polishing iteration
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_1}.1st_polishing_iteration.cns.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_1}.2nd_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_1}.2nd_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_1}.1st_polishing_iteration.cns.fa -i - -fo {params.prefix_haplotype_1}.2nd_polishing_iteration.cns.fa

        # 3nd polishing iteration
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_1}.2nd_polishing_iteration.cns.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_1}.3nd_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_1}.3nd_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_1}.2nd_polishing_iteration.cns.fa -i - -fo {params.prefix_haplotype_1}.3nd_polishing_iteration.cns.fa

        cd {params.folder_hap2}
        ### Haplotype 2 ###
        # assemble long reads
        {params.wtdbg2} -x sequel -g 3g -t {number_of_parallel_threads} -i {input[0]} -fo {params.prefix_haplotype_2}

        # derive consensus
        {params.wtpoa_cns} -t {number_of_parallel_threads} -i {params.prefix_haplotype_2}.ctg.lay.gz -fo {params.prefix_haplotype_2}.raw.fa

        # 1st polishing iteration, not necessary if you want to polish the assemblies using other tools
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_2}.raw.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_2}.1st_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_2}.1st_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_2}.raw.fa -i - -fo {params.prefix_haplotype_2}.1st_polishing_iteration.cns.fa

        # 2nd polishing iteration
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_2}.1st_polishing_iteration.cns.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_2}.2nd_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_2}.2nd_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_2}.1st_polishing_iteration.cns.fa -i - -fo {params.prefix_haplotype_2}.2nd_polishing_iteration.cns.fa

        # 3nd polishing iteration
        minimap2 -t {number_of_parallel_threads} -ax map-pb -r2k {params.prefix_haplotype_2}.2nd_polishing_iteration.cns.fa {input[0]} | samtools sort -@ {number_of_parallel_threads_for_sorting} > {params.prefix_haplotype_2}.3nd_polishing_iteration.bam
        samtools view -F0x900 {params.prefix_haplotype_2}.3nd_polishing_iteration.bam | {params.wtpoa_cns} -t {number_of_parallel_threads} -d {params.prefix_haplotype_2}.2nd_polishing_iteration.cns.fa -i - -fo {params.prefix_haplotype_2}.3nd_polishing_iteration.cns.fa
        """
        shell(cmd)


######## Structural Variant Calling pipeline ###################################
################################################################################

rule reveal_sv_caller_de_novo_based_haplotypes: # Doesn't work properly yet. Contains some errors still.
    input:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'references/hg38.fa'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/{sample_id}/haplotype_1/hap1_asm.fasta'),
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'assemblies/{sample_id}/haplotype_2/hap2_asm.fasta')
    output:
        os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/sv_calling_on_haplotypes/de_novo_asm_based_haplotypes/{sample_id}/haplotype_2/asm_haplotype_2.prg.refined.vcf')
    params:
        prefix_hap1=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/sv_calling_on_haplotypes/de_novo_asm_based_haplotypes/{sample_id}/haplotype_1'),
        prefix_hap2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/sv_calling_on_haplotypes/de_novo_asm_based_haplotypes/{sample_id}/haplotype_2'),
        input1=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/sv_calling_on_haplotypes/de_novo_asm_based_haplotypes/{sample_id}/haplotype_1/asm_haplotype_1'),
        input2=os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'vcf/sv_calling_on_haplotypes/de_novo_asm_based_haplotypes/{sample_id}/haplotype_2/asm_haplotype_2'),
        number_of_threads=20
    conda:
        PY2_ENV # This specifies the conda env containing Python 2, because reveal uses explicitly Python 2.
    shell:
        """
        cd {params.prefix_hap1}
        #(1) Convert draft assemblies to graphs (address rearrangements)
        ##reveal transform --64 {input[0]} {input[1]} -o {params.input1}.gfa

        ##reveal split --64 {params.prefix_hap1}.gfa --nocycles

        #(2) Use REM to construct an anchor based alignment graph (brake down the problem)
        ##reveal rem --64 {params.input1}.gfa_chr1* -m20 --threads={params.number_of_threads} -o {params.input1}.merged1.gfa
        ##reveal rem --64 {params.input1}.gfa_chr2* -m20 --threads={params.number_of_threads} -o {params.input1}.merged2.gfa
        ##reveal rem --64 {params.input1}.gfa_chr3* {params.input1}.gfa_chr4* {params.input1}.gfa_chr5* {params.input1}.gfa_chr6* {params.input1}.gfa_chr7* {params.input1}.gfa_chr8* {params.input1}.gfa_chr9* {params.input1}.gfa_chrU* {params.input1}.gfa_chrX* {params.input1}.gfa_chrY* -m20 --threads={params.number_of_threads} -o {params.input1}.merged3.gfa
        ##reveal rem --64 {params.input1}.merged1.gfa {params.input1}.merged2.gfa {params.input1}.merged3.gfa -m20 --threads={params.number_of_threads} -o {params.input1}.merged.gfa
        ##reveal rem --64 {input[0]} {params.input1}.merged.gfa -m20 --threads={params.number_of_threads} -o {params.input1}.prg.gfa
        #(3) Unzip all bubbles in the graph
        ##reveal unzip --64 {params.input1}.prg.gfa -u10
        #(4) Refine all bubbles in the graph using MSA
        ##reveal refine --64 {params.input1}.prg.unzipped.gfa --nproc={params.number_of_threads} --all --maxsize=10000 --minsize=2 --mindiff=0 --minconf=90
        #(5) Output variants
        ##reveal variants --64 {params.input1}.prg.gfa --bed > {params.input1}.prg.anchored.bed
        ##reveal variants --64 {params.input1}.prg.unzipped.gfa --bed > {params.input1}.prg.unzipped.bed
        reveal variants --64 {params.input1}.prg.unzipped.realigned.gfa --vcf > {params.input1}.prg.refined.vcf

        cd {params.prefix_hap2}
        #(1) Convert draft assemblies to graphs (address rearrangements)
        ##reveal transform --64 {input[0]} {input[2]} -o {params.input2}.gfa

        ##reveal split --64 {params.input2}.gfa --nocycles

        #(2) Use REM to construct an anchor based alignment graph (brake down the problem)
        ##reveal rem --64 {params.input2}.gfa_chr1* -m20 --threads={params.number_of_threads} -o {params.input2}.merged1.gfa
        ##reveal rem --64 {params.input2}.gfa_chr2* -m20 --threads={params.number_of_threads} -o {params.input2}.merged2.gfa
        ##reveal rem --64 {params.input2}.gfa_chr3* {params.input2}.gfa_chr4* {params.input2}.gfa_chr5* {params.input2}.gfa_chr6* {params.input2}.gfa_chr7* {params.input2}.gfa_chr8* {params.input2}.gfa_chr9* {params.input2}.gfa_chrU* {params.input2}.gfa_chrX* {params.input2}.gfa_chrY* -m20 --threads={params.number_of_threads} -o {params.input2}.merged3.gfa
        ##reveal rem --64 {params.input2}.merged1.gfa {params.input2}.merged2.gfa {params.input2}.merged3.gfa -m20 --threads={params.number_of_threads} -o {params.input2}.merged.gfa
        ##reveal rem --64 {input[0]} {params.input2}.merged.gfa -m20 --threads={params.number_of_threads} -o {params.input2}.prg.gfa
        #(3) Unzip all bubbles in the graph
        ##reveal unzip --64 {params.input2}.prg.gfa -u10
        #(4) Refine all bubbles in the graph using MSA
        reveal refine --64 {params.input2}.prg.unzipped.gfa --nproc={params.number_of_threads} --all --maxsize=10000 --minsize=2 --mindiff=0 --minconf=90
        #(5) Output variants
        ##reveal variants --64 {params.input2}.prg.gfa --bed > {params.input2}.prg.anchored.bed
        ##reveal variants --64 {params.input2}.prg.unzipped.gfa --bed > {params.input2}.prg.unzipped.bed
        reveal variants --64 {params.input2}.prg.unzipped.realigned.gfa --vcf > {params.input2}.prg.refined.vcf
        """
