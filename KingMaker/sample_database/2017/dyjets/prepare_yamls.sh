##xrdfs root://xrootd.cmsaf.mit.edu:1094/ ls /store/user/paus/nanosu/A02/ > tmp_processes
#
#while read p;
#do
#    file=$(basename $p)
#    process=${file%.*} 
#    echo dbs: ${process} > ${process}.yaml
#    echo era: 2018 >> ${process}.yaml
#    echo filelist: >> ${process}.yaml
#
#    xrdfs root://cmsdcache-kit-disk.gridka.de:1094/ ls /store/user/jhornung/bkg_raw/${process} | xargs -n1 -P1 -i echo - root://cmsdcache-kit-disk.gridka.de:1094/{} >> ${process}.yaml
#
#    echo generator_weight: 1.0 >> ${process}.yaml
#    echo nevents: 197600000 >> ${process}.yaml
#
#    n=$(gfal-ls root://cmsdcache-kit-disk.gridka.de:1094//store/user/jhornung/bkg_raw/${process} | wc -l)
#
#    echo nfiles: $n >> ${process}.yaml
#    echo nick: $process >> ${process}.yaml
#    echo sample_type: bkg >> ${process}.yaml
#    echo xsec: 1.0 >> ${process}.yaml
#done < $HAA/KingMaker/vh.txt
#

sample_name="/DYJetsToMuMu_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9_ext1-v1/NANOAODSIM"
nick="DYJetsToMuMu_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos"

echo "dbs: ${sample_name}" > ${nick}.yaml
echo "era: 2017" >> ${nick}.yaml
echo "filelist:" >> ${nick}.yaml

dasgoclient -query="file dataset=${sample_name}" | xargs  -n1 -P1 -i echo - root://cmsxrootd.fnal.gov/{} >> ${nick}.yaml

echo "generator_weight: 1.0" >> ${nick}.yaml
echo "nevents: 197600000" >> ${nick}.yaml

n=$(dasgoclient -query="file dataset=${sample_name}" | wc -l)

echo "nfiles: $n" >> ${nick}.yaml
echo "nick: $nick" >> ${nick}.yaml
echo "sample_type: dyjets" >> ${nick}.yaml
echo "xsec: 1.0" >> ${nick}.yaml