#/bin/bash

# streamline limit production

scope=$1

source /cvmfs/cms.cern.ch/cmsset_default.sh

# do fits and prep workspaces
source $HAA/alt_setup.sh
python3 /work/jhornung/Haa/Scripts/signal_fit.py $scope
python3 /work/jhornung/Haa/Scripts/bkg_fit.py $scope
deactivate
# make datacard

cd /work/jhornung/Haa/limits/CMSSW_11_3_4/src/
cmsenv
cd HiggsAnalysis/CombinedLimit

git fetch origin
git checkout v9.2.1
scramv1 b clean; scramv1 b # always make a clean build

cd /work/jhornung/Haa/limits
echo "---------------" > /work/jhornung/Haa/limits/$scope.txt
echo "imax 1" >> /work/jhornung/Haa/limits/$scope.txt
echo "jmax 1" >> /work/jhornung/Haa/limits/$scope.txt
echo "kmax *" >> /work/jhornung/Haa/limits/$scope.txt
echo "---------------" >> /work/jhornung/Haa/limits/$scope.txt
echo "shapes Haa H_mass workspace_sig_$scope.root workspace_sig_$scope:model_signal_H_mass" >> /work/jhornung/Haa/limits/$scope.txt
echo "shapes data_obs H_mass workspace_bkg_$scope.root workspace_bkg_$scope:model_data_H_mass" >> /work/jhornung/Haa/limits/$scope.txt
echo "shapes bkg_mass H_mass workspace_bkg_$scope.root workspace_bkg_$scope:model_bkg_H_mass" >> /work/jhornung/Haa/limits/$scope.txt
echo "---------------" >> /work/jhornung/Haa/limits/$scope.txt
echo "bin H_mass" >> /work/jhornung/Haa/limits/$scope.txt
echo "observation -1" >> /work/jhornung/Haa/limits/$scope.txt
echo "---------------" >> /work/jhornung/Haa/limits/$scope.txt
echo "bin H_mass H_mass" >> /work/jhornung/Haa/limits/$scope.txt
echo "process Haa bkg_mass" >> /work/jhornung/Haa/limits/$scope.txt
echo "process 0 1" >> /work/jhornung/Haa/limits/$scope.txt
echo "rate 59e3 1" >> /work/jhornung/Haa/limits/$scope.txt
echo "---------------" >> /work/jhornung/Haa/limits/$scope.txt

# run limits
cd /work/jhornung/Haa/limits/CMSSW_11_3_4/src/
cmsenv
cd HiggsAnalysis/CombinedLimit

git fetch origin
git checkout v9.2.1
scramv1 b clean; scramv1 b # always make a clean build

cd /work/jhornung/Haa/limits
combine -M AsymptoticLimits $scope.txt -m 125 --freezeParameters MH --saveWorkspace -n $scope.bestfit