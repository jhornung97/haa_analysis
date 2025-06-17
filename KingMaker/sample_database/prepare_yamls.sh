ls /work/jhornung/Haa/KingMaker/sample_database/2016postVFP/data/*MINIAOD.yaml > tmp_processes

while read p;
do
    file=$(basename $p)
    process=${file%.*}
    echo $process
    echo "$process:" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
    echo "  dbs: ${process}" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
    echo "  era: 2016postVFP" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
    echo "  generator_weight: 1.0" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
    echo "  nevents: 197600000" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
    n=$(xrdfs root://cmsdcache-kit-disk.gridka.de:1094/ ls /store/user/jhornung/data_raw/${process} | wc -l)
    echo $n
    echo "  nfiles: $n" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
    echo "  nick: $process" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
    echo "  sample_type: data" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
    echo "  xsec: 1.0" >> /work/jhornung/Haa/KingMaker/sample_database/datasets.yaml
done < tmp_processes
