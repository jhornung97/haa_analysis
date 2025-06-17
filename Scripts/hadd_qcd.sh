#!/bin/bash

cd /ceph/jhornung/MC_2018/2018/
pwd

procs=("${@:1:$#-1}")

target=(${!#})

for proc in ${procs[@]}; do
    ls -d1 $proc* >> tmp_processes
done

while read p;
do
    ls -1 $p > scopes
    while read s;
    do
    rm $target/$s/${p}.root
	cp $p/$s/${p}.root $target/$s/
    done < scopes
    rm scopes
done < tmp_processes

rm tmp_processes

rm $target/mm/$target.root
rm $target/em/$target.root
rm $target/ee/$target.root

hadd $target/mm/$target.root $target/mm/*.root
hadd $target/em/$target.root $target/em/*.root
hadd $target/ee/$target.root $target/ee/*.root

cd $HAA