#!/bin/bash

cd /ceph/jhornung/MC_2018/2018
ls -d -I "qcd" QCD_HT* > QCD_HT_dirs.txt

while read line; do
    if [ "$line" != "QCD_HT_dirs.txt" ]; then
        echo $line
        cd $line
        ls > scopes
        while read s; do
            if [ "$s" != "scopes" ]; then
                echo $s
                cd $s
                rm $line.root
                rm ${line}_final.root
                ls *.root > files
                while read f; do
                    if [ "$f" != "files" ]; then
                        echo $f
                        nevents=$(python3 $HAA/Scripts/check_file.py $f)
                        if [ $nevents -gt 0 ]; then
                            echo "file not empty"
                            hadd $line.root $f
                        fi                        
                    fi
                done < files
                rm files
                python3 $HAA/Scripts/add_crosssection_weight.py $line.root
                cp $line.root /ceph/jhornung/MC_2018/2018/qcd/$s
                cd ..
            fi
        
        done < scopes
        rm scopes
        cd ..
    fi

done < QCD_HT_dirs.txt
rm QCD_HT_dirs.txt

hadd -f /ceph/jhornung/MC_2018/2018/qcd/mm/qcd.root /ceph/jhornung/MC_2018/2018/qcd/mm/*.root
hadd -f /ceph/jhornung/MC_2018/2018/qcd/ee/qcd.root /ceph/jhornung/MC_2018/2018/qcd/ee/*.root
hadd -f /ceph/jhornung/MC_2018/2018/qcd/em/qcd.root /ceph/jhornung/MC_2018/2018/qcd/em/*.root

cd $HAA