#!/bin/bash

work_path=$1
work_dir_name=$2
samples=$3

mkdir -p ${work_path}/${work_dir_name}/reads

for dir in ${work_path}/*/ ; do
    cd ${dir};
    if [[ $( ls | egrep ${samples} | wc -l ) == 1 ]] ; then
        mv ${samples} ${work_path}/${work_dir_name}/reads;
    fi
done

######################################## statement below will delete every map except work_dir
#find ${work_path} -mindepth 1 ! -regex "^${work_path}/${work_dir_name}\(/.*\)?" -delete
