#!/bin/bash

process_mc_dir() {
    if [ ! $1 == "tmp_processes" ];
    then
        cd $1
        echo $1
        ls > scopes
        while read s;
        do
            if [ "$s" != "scopes" ];
            then
                cd $s
                echo $s
                rm $1.root
                hadd -f $1.root *.root
                python3 $HAA/Scripts/add_crosssection_weight.py $1.root
                cd ..
            fi
        done < scopes
        rm scopes
        cd ..
    fi
    pwd
}

export -f process_mc_dir

if [ "$1" == "data" ];
then
    cd /ceph/jhornung/Data_2018/2016preVFP/
    ls -d1 SingleMuon* > tmp_processes
    
    while read p;
    do
	if [ "$p" != "tmp_processes" ];
	then
	    echo "Now processing $p"
	    cd $p
	    pwd
	    ls > scopes
	    while read s;
	    do
		if [ "$s" != "scopes" ];
		then
		    pwd
		    echo $s
		    cd $s
		    rm $p.root
		    hadd -f $p.root *.root
		    cp $p.root /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/$s
		    cd ..
		fi
	    done < scopes
	    rm scopes
            cd ..
	    pwd
	fi
    done < tmp_processes

    rm /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/ee/single_muon_data.root
    rm /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/mm/single_muon_data.root
    rm /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/em/single_muon_data.root

    hadd -f /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/ee/single_muon_data.root /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/ee/*.root
    hadd -f /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/mm/single_muon_data.root /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/mm/*.root
    hadd -f /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/em/single_muon_data.root /ceph/jhornung/Data_2018/2016preVFP/single_muon_data/em/*.root

    rm tmp_processes
fi

if [ "$1" == "mc" ];
then
    cd /ceph/jhornung/MC_2018/2017/
    #ls  > tmp_processes
#    ls -d > tmp_processes
    ls -I "diboson" -I "wh" -I "ttbar" -I "single_top" -I qcd -I "VH*" -I "QCD*" -I tmp_processes -I "full_bkg" > tmp_processes
  #ls -d ST* | grep -v "ST_tW_Dilept*" >  tmp_processes
#    while read p;
#    do
#        if [ "$p" != "tmp_processes" ];
#        then
#            #echo "Now processing $p"
#            cd $p
#            pwd
#            echo " "
#            ls > scopes
#            while read s;
#            do
#                if [ "$s" != "scopes" ];
#                then
#                    cd $s
#                    pwd
#                    echo $s
#		            rm $p.root
#		            rm ${p}_final.root
#                    hadd -f $p.root *.root
#		            python3 $HAA/Scripts/add_crosssection_weight.py $p.root 
#                    cd ..
#                fi
#            done < scopes
#            rm scopes
#            cd ..
#            pwd
#        fi
#    done < tmp_processes

    cat tmp_processes | xargs -I {} -P 1 bash -c 'process_mc_dir {}'

    rm tmp_processes
fi

cd $HAA
