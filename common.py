"""
This script looks whether there are any samples that are yet to be processed.
If so, than the script gathers the subreads for each sample and puts them in
/reads. Subsequently, the sample_id's are collected, as they are used by the
pipelines to select the right subreads for each sample.
"""

import os

SCRIPTS_PATH="scripts"
WORK_PATH="/tudelft.net/staff-bulk/ewi/insy/VUMC/pacbio" # Change WORK_PATH to the path where all samples are stored
WORK_DIR_NAME="downstream_files"

READS_DIR = os.path.join(str(WORK_PATH), str(WORK_DIR_NAME), 'reads')

SAMPLE_FILES = []

def gather_files(path):
    cmd = "mkdir -p " + os.path.join(str(WORK_PATH), WORK_DIR_NAME)
    os.system(cmd)
    cmd = "mkdir -p " + READS_DIR
    os.system(cmd)
    for d in os.listdir(str(WORK_PATH)):
        d = os.path.join(str(WORK_PATH), d)
        for f in os.listdir(d):
            if f.endswith('subreads.bam') or f.endswith('subreads.bam.pbi'): #necessary raw long-read file
                os.replace(os.path.join(d, f), path + "/" + f)
    cmd = """find ${work_path} -mindepth 1 ! -regex "^${WORK_PATH}/${WORK_DIR_NAME}\(/.*\)?" -delete"""
    #os.system(cmd) # uncomment for automatic removal of files other than the actual subreads

#read in all long-read files in specific dir as samplefiles
gather_files(READS_DIR)
for f in os.listdir(READS_DIR):
    if f.endswith('subreads.bam'):
        SAMPLE_FILES.append(f.replace(".subreads.bam", ""))
SAMPLE_FILES.sort()
