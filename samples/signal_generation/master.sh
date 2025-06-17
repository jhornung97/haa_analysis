#!/bin/bash


rand=`shuf -i 0-5000000 -n1`

source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc700

### GEN
echo "###########"
echo "GEN step"
echo "###########"
scram p CMSSW CMSSW_10_6_30_patch1
cd CMSSW_10_6_30_patch1/src
eval `scram runtime -sh`
scram b -j1
cd ../../
cmsRun gen.py seed=$rand

### SIM
echo "###########"
echo "SIM step"
echo "###########"
scram p CMSSW CMSSW_10_6_17_patch1
cd CMSSW_10_6_17_patch1/src
eval `scram runtime -sh`
scram b -j1
cd ../../
cmsRun sim.py

### DIGI
echo "###########"
echo "DIGI step"
echo "###########"
cd CMSSW_10_6_17_patch1/src
eval `scram runtime -sh`
scram b -j1
cd ../../
cmsRun digi.py

### HLT
echo "###########"
echo "HLT step"
echo "###########"
scram p CMSSW CMSSW_10_2_16_UL
cd CMSSW_10_2_16_UL/src
eval `scram runtime -sh`
scram b -j1
cd ../../
cmsRun hlt.py

### RECO
echo "###########"
echo "RECO step"
echo "###########"
cd CMSSW_10_6_17_patch1/src
eval `scram runtime -sh`
scram b -j1
cd ../../
cmsRun reco.py

### MINIAOD
echo "###########"
echo "MINIAOD step"
echo "###########"
scram p CMSSW CMSSW_10_6_20
cd CMSSW_10_6_20/src
eval `scram runtime -sh`
scram b -j1
cd ../../
cmsRun miniaod.py

