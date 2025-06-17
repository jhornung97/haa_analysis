#!/bin/bash

cd /ceph/jhornung/MC_2018/2016postVFP/
pwd

procs=("${@:1:$#-1}")

target=(${!#})

> tmp_processes
for proc in ${procs[@]}; do
    ls -d1 $proc* >> tmp_processes
done

if [ ! -d $target ]; 
then
    echo "Creating directory $target"
    mkdir $target
fi

for subdir in mm em ee; 
do
    if [ ! -d $target/$subdir ]; 
    then
        echo "Creating directory $target/$subdir"
        mkdir $target/$subdir
    fi
done

while read p;
do
    echo $p
    ls -1 $p > scopes
    while read s;
    do
        echo $s
        rm $target/$s/${p}.root
	    cp $p/$s/${p}.root $target/$s/
    done < scopes
    rm scopes
done < tmp_processes

rm -f tmp_processes

rm -f $target/mm/$target.root
rm -f $target/em/$target.root
rm -f $target/ee/$target.root

hadd $target/mm/$target.root $target/mm/*.root
hadd $target/em/$target.root $target/em/*.root
hadd $target/ee/$target.root $target/ee/*.root

cd $HAA