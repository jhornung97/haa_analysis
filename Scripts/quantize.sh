#!/bin/bash

touch quantize_changes.txt

cd /ceph/jhornung

ls /ceph/jhornung/analysis_old/2018/ | grep -v "Single*" | grep -v "merged_mc" | grep -v "single_muon_data" > tmp

while read p;
do 
    echo "$p"
    echo "$p" >> $HAA/quantize_changes.txt
    ls analysis_old/2018/$p | xargs basename --multiple > scopes
    while read s;
    do 
	echo "$s"
	echo "$s" >> $HAA/quantize_changes.txt
	python3 $HAA/Scripts/read_cutflow.py analysis_old/2018/$p/$s/$p.root >> $HAA/quantize_changes.txt
	python3 $HAA/Scripts/read_cutflow.py analysis/2018/$p/$s/$p.root >> $HAA/quantize_changes.txt
    done < scopes
done < tmp

rm tmp

cd $HAA
