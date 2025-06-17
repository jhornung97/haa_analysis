ls /work/jhornung/Haa/KingMaker/sample_database/2018/data/EGamma+Run2018*.yaml > tmp_processes

while read p;
do
    file=$(basename $p)
    process=${file%.*}
    echo "$process"
    echo dbs: ${process} > ${process}.yaml
    echo era: 2018 >> ${process}.yaml
    echo filelist: >> ${process}.yaml

    xrdfs root://xrootd.cmsaf.mit.edu:1094/ ls /store/user/paus/nanosu/A02/${process} | xargs basename --multiple | xargs -n1 -P1 -i echo - root://cmsdcache-kit-disk.gridka.de:1094//store/user/jhornung/data_raw/${process}/{} >> ${process}.yaml

    echo generator_weight: 1.0 >> ${process}.yaml
    echo nevents: 197600000 >> ${process}.yaml

    n=$(xrdfs root://xrootd.cmsaf.mit.edu:1094/ ls -l /store/user/paus/nanosu/A02/${process} | wc -l)

    echo nfiles: $n >> ${process}.yaml
    echo nick: $process >> ${process}.yaml
    echo sample_type: data >> ${process}.yaml
    echo xsec: 1.0 >> ${process}.yaml
done < tmp_processes
