#!/bin/bash

declare -A map

xrdfs root://cmsxrootd-kit-disk.gridka.de:1094/ ls /store/user/jhornung/bkg_raw/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/ | xargs basename --multiple > tmp_files

while read -r f; do
    shortName=${f:0:6}
    if [[ -n ${map[$shortName]} ]]; then
        echo "Overlap found: $shortName"
    else
        map[$shortName]=1
    fi
done < tmp_files