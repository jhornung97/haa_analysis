#!/bin/bash

home=$(pwd)

tar -xzf $home/utils.tar.gz 

. $home/utils/haa_setup.sh

dirs=$(ls -d $home/utils/input/2* | xargs -n1 basename)

iso_mask=$1

mkdir $home/limits_$iso_mask

for dir in $dirs; do
    echo "Processing year: $dir"
    #scopes=$(ls -d $HAA/limits_$iso_mask/$dir/* | xargs -n1 basename)

    mkdir $home/limits_$iso_mask/$dir

    scopes=$(ls -d $home/utils/input/$dir/* | xargs -n1 basename | sort -u)
    scopes="$(echo "$scopes" | grep -v '^combined$'; echo "$scopes" | grep '^combined$')"

    for scope in $scopes; do
        mkdir $home/limits_$iso_mask/$dir/$scope
        if [ "$scope" != "combined" ]; then
            echo "Processing scope: $scope"
            python3 $home/utils/scripts/systematics_signal_dcb+gauss.py $scope $dir $iso_mask $home
            python3 $home/utils/scripts/bkg_fit.py $scope $dir $iso_mask $home

            cp $home/utils/input/$dir/$scope/$scope.txt $home/limits_$iso_mask/$dir/$scope/

            cd $home/limits_$iso_mask/$dir/$scope
            echo "Running combine for $scope in $dir"
            combine -M AsymptoticLimits $scope.txt -m 125 --saveWorkspace -n .limits_iso_cut_$iso_mask --run blind
        fi
        
        if [ "$scope" == "combined" ]; then
            cd $home/limits_$iso_mask/$dir/$scope
            cp $home/utils/input/$dir/$scope/combined.txt $home/limits_$iso_mask/$dir/$scope/
            echo "Running combine for combined scope in $dir"
            combine -M AsymptoticLimits $scope.txt -m 125 --saveWorkspace -n .limits_iso_cut_$iso_mask --run blind
            python3 $home/utils/scripts/grab_output.py $iso_mask $dir $home
        fi

        echo "Finished processing scope: $scope"
    done
    echo "Finished processing year: $dir"
done