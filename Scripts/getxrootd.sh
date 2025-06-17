#!/bin/bash

DATASET="/DYJetsToMuMu_M-50_TuneCP5_ZptWeighted_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM"
dasgoclient -query="file dataset=${DATASET} | grep site.name, site.dataset_fraction" -json | python -c "
import sys, json
for item in json.load(sys.stdin):
    for site in item.get('site', []):
        if site.get('dataset_fraction', 0) > 0:
            print(site.get('name', ''))" | while read -r site; do
    echo "Site: ${site}"  # Debugging line
    dasgoclient -query="file site=${site} dataset=${DATASET}" -json | python -c "
import sys, json
for item in json.load(sys.stdin):
    for file in item.get('file', []):
        print(file.get('name', ''))" | while read -r file; do
        echo "File: ${file}"  # Debugging line
        echo "root://${site}/${file}"
    done
done
