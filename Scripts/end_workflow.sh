#!/bin/bash

. $HOME/haa_setup.sh

dirs=2018 #$(ls -d $HAA/limits/2* | xargs -n1 basename)

iso_mask=0.7

for dir in $dirs; do
    echo "Processing year: $dir"
    #scopes=$(ls -d $HAA/limits/$dir/* | xargs -n1 basename)

    scopes=$(ls -d $HAA/limits/$dir/* | xargs -n1 basename | sort -u)
    scopes="$(echo "$scopes" | grep -v '^combined$'; echo "$scopes" | grep '^combined$')"

    for scope in $scopes; do
        if [ "$scope" != "combined" ]; then
            echo "Processing scope: $scope"
            python3 $HAA/Scripts/systematics_signal_dcb+gauss.py $scope $dir $iso_mask
            python3 $HAA/Scripts/bkg_fit.py $scope $dir $iso_mask

            cd $HAA/limits/$dir/$scope
            echo "Running combine for $scope in $dir"
            combine -M AsymptoticLimits $scope.txt -m 125 --saveWorkspace -n .limits_iso_cut_$iso_mask --run blind
        fi
        
        if [ "$scope" == "combined" ]; then
            cd $HAA/limits/$dir/$scope
            echo "Running combine for combined scope in $dir"
            combine -M AsymptoticLimits $scope.txt -m 125 --saveWorkspace -n .limits_iso_cut_$iso_mask --run blind
        fi

        echo "Finished processing scope: $scope"
    done
    echo "Finished processing year: $dir"
done

#cd $HAA/limits/all_years_combined
#echo "Running combine for all years combined"
#combine -M AsymptoticLimits all_years_combined_with_2016.txt -m 125 --saveWorkspace -n .limits_iso_cut_$iso_mask --run blind
#combine -M AsymptoticLimits all_years_combined_only_2016.txt -m 125 --saveWorkspace -n .limits2016_iso_cut_$iso_mask --run blind
#
#mkdir -p /ceph/jhornung/limit_scan/all_years_combined
#mv $HAA/limits/all_years_combined/higgsCombine.limits_iso_cut_$iso_mask.AsymptoticLimits.mH125.root /ceph/jhornung/limit_scan/all_years_combined/
#mkdir -p /ceph/jhornung/limit_scan/all_years_combined_only_2016
#mv $HAA/limits/all_years_combined/higgsCombine.limits2016_iso_cut_$iso_mask.AsymptoticLimits.mH125.root /ceph/jhornung/limit_scan/all_years_combined_only_2016/