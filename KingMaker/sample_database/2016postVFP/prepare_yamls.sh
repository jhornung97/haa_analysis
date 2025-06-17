xrdfs root://cmsdcache-kit-disk.gridka.de:1094/ ls /store/user/jhornung/bkg_raw | grep 20UL16 | grep -v preVFP | xargs basename --multiple > tmp_processes

while read p;
do
    file=$(basename $p)
    process=${file%.*}
    echo "$process"
    echo dbs: ${process} > /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml
    echo era: 2016postVFP >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml
    echo filelist: >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml

    xrdfs root://cmsdcache-kit-disk.gridka.de:1094/ ls /store/user/jhornung/bkg_raw/${process} | xargs basename --multiple | xargs -n1 -P1 -i echo - root://cmsdcache-kit-disk.gridka.de:1094//store/user/jhornung/bkg_raw/${process}/{} >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml

    echo generator_weight: 1.0 >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml
    echo nevents: 197600000 >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml

    n=$(xrdfs root://cmsdcache-kit-disk.gridka.de:1094/ ls /store/user/jhornung/bkg_raw/${process} | wc -l)
    echo $n

    echo nfiles: $n >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml
    echo nick: $process >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml
    echo sample_type: bkg >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml
    echo xsec: 1.0 >> /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/bkg/${process}.yaml
done < tmp_processes
