#!/bin/bash
#ls -I "HaaKKK*" -I "diboson" -I "single_top" -I "ttbar" -I "wh" /ceph/jhornung/MC_2018/ > bkg_processes

for f in $(cat /work/jhornung/Haa/data_17_18.txt); do
	echo "Checking directory: $f"
	xrdfs root://xrootd.cmsaf.mit.edu:1094/ ls /store/user/paus/nanosu/A02/$f | xargs basename --multiple > tmp_files
	
	# Check if the directory exists at GridKA, if not, create it
    if ! xrdfs "root://cmsdcache-kit-disk.gridka.de:1094/" stat "/store/user/jhornung/data_raw/$f" &> /dev/null; then
        echo "Creating directory $f at GridKA."
        xrdfs "root://cmsdcache-kit-disk.gridka.de:1094/" mkdir -p "/store/user/jhornung/data_raw/$f"
    fi
	
	rm -f $f.txt

	touch $f.txt

	dir_is_good=0

	while [ "$dir_is_good" -eq "0" ];
	do	
    	counter=0
    	while read p;
    	do
		if [[ "$p" != "tmp"* ]];
		then
	    	MIT_size=$(xrdfs "root://xrootd.cmsaf.mit.edu:1094/" ls -l "/store/user/paus/nanosu/A02/$f/$p" | awk '{print $4}')
	    	GridKA_size=$(xrdfs "root://cmsdcache-kit-disk.gridka.de:1094/" ls -l "/store/user/jhornung/data_raw/$f/$p" | awk '{print $4}')
	    	if [ "$MIT_size" != "$GridKA_size" ];
	    	then
			echo "File $p is corrupted at GridKA as sizes are not equal."
			echo "Size at MIT: $MIT_size"
			echo "Size at GridKA: $GridKA_size"
			echo "$p" >> $f.txt
			let counter++
		    fi
		fi
    	done < tmp_files
		
    	cat $f.txt | xargs -n1 -P8 -i xrdcp -rf root://xrootd.cmsaf.mit.edu:1094//store/user/paus/nanosu/A02/$f/{} root://cmsdcache-kit-disk.gridka.de:1094//store/user/jhornung/data_raw/$f/

    	rm $f.txt

    	echo "Last iteration found $counter corrupted files"
    	if [ "$counter" -eq "0" ];
    	then
		echo "Directory is good! :partyparrot:"
		dir_is_good=1
    	fi
	done
done
