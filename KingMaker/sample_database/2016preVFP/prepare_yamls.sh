xrdfs root://cmsdcache-kit-disk.gridka.de:1094/ ls /store/user/jhornung/data_raw | grep 16 | grep HIPM | xargs basename --multiple > tmp_processes

while read p;
do
    file=$(basename $p)
    process=${file%.*}
    echo "$process"
    echo dbs: ${process} > /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml
    echo era: 2016preVFP >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml
    echo filelist: >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml

    xrdfs root://cmsdcache-kit-disk.gridka.de:1094/ ls /store/user/jhornung/data_raw/${process} | xargs basename --multiple | xargs -n1 -P1 -i echo - root://cmsdcache-kit-disk.gridka.de:1094//store/user/jhornung/data_raw/${process}/{} >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml

    echo generator_weight: 1.0 >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml
    echo nevents: 197600000 >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml

    n=$(xrdfs root://cmsdcache-kit-disk.gridka.de:1094/ ls /store/user/jhornung/data_raw/${process} | wc -l)
    echo $n

    echo nfiles: $n >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml
    echo nick: $process >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml
    echo sample_type: data >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml
    echo xsec: 1.0 >> /work/jhornung/Haa/KingMaker/sample_database/2016preVFP/data/${process}.yaml
done < tmp_processes
