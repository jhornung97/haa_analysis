#!/bin/bash

xrdfs root://cmsxrootd-kit-disk.gridka.de:1094/ ls /store/user/jhornung/bkg_raw/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/ | xargs basename --multiple > tmp_files

declare -a noMentions
declare -a multipleMentions

tag=$1

for f in $(cat tmp_files); do
    echo "Checking file: $f"
    shortName=${f:0:8}
    nMentions=$(grep -rl $shortName $HAA/KingMaker/data/logs/$tag/Output/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/ | wc -l)
    if [ $nMentions -eq 0 ]; then
        echo "File $f not found in logs"
    else
        echo "File $f found $nMentions times in logs"
        if [ $nMentions -ge 2 ]; then
            multipleMentions+=($f)
        fi      
    fi
done

echo "Files with no mentions:"
for f in ${noMentions[@]}; do
    echo $f
done

echo "Files with multiple mentions:"
for f in ${multipleMentions[@]}; do
    echo $f
done

rm tmp_files