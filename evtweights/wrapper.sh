#!/bin/bash

year=$1
dataset=$2
home=$(pwd)

tar -xf utils.tar.gz

. utils/haa_setup.sh

cd /ceph/jhornung/MC_2018/$year/$dataset/
ls > scopes

while read s;
do 
    if [ "$s" != "scopes" ] && [ -d $s ];
    then 
        cd $s
        echo $s 
        rm $dataset.root
        hadd -f $dataset.root *.root
        python3 $home/utils/add_crosssection_weight.py $dataset.root $home/utils/crosssections.json
        cd ..
    fi
done < scopes

rm scopes