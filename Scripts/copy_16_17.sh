#!/bin/bash

MIT_dir=$(xrdfs root://xrootd.cmsaf.mit.edu:1094/ ls /store/user/paus/nanosu/A02/ | xargs basename --multiple)

dir_array=($MIT_dir)

touch copy.txt

for d in $(cat /work/jhornung/Haa/KingMaker/bkgmc.txt); do
    dir16=${d/18/16}
    dir17=${d/18/17}
    
    cut_dir16=$(echo $dir16 | awk '{split($0,a,"16"); print a[1] "16"}')
    cut_dir17=$(echo $dir17 | awk '{split($0,a,"17"); print a[1] "17"}')

    for dir in ${dir_array[@]}; do
        if [[ $dir == *$cut_dir16* ]]; then
            echo $dir >> copy.txt
        fi
        if [[ $dir == *$cut_dir17* ]]; then
            echo $dir >> copy.txt
        fi
    done

done

cat copy.txt | sort | uniq > copy_16_17.txt

cat copy_16_17.txt | xargs -n1 -P8 -i xrdcp -r root://xrootd.cmsaf.mit.edu:1094//store/user/paus/nanosu/A02/{} root://cmsxrootd-kit-disk.gridka.de:1094//store/user/jhornung/bkg_raw