#!/bin/bash

# Überprüfen, ob zwei Verzeichnisse als Argumente übergeben wurden
#if [ $# -ne 2 ]; then
#    echo "Usage: $0 xrootd://host/dir1 xrootd://host/dir2"
#    exit 1
#fi

dir1="$1"
output_file="$1_diff.txt"

echo $(xrdfs "root://xrootd.cmsaf.mit.edu:1094/" ls -l "/store/user/paus/nanosu/A02/$1" | wc -l)
echo $(xrdfs "root://cmsxrootd-kit-disk.gridka.de:1094/" ls -l "/store/user/jhornung/data_raw/$1" | wc -l)

# Überprüfen, ob die angegebenen Verzeichnisse existieren
xrdfs "root://xrootd.cmsaf.mit.edu:1094/" stat "/store/user/paus/nanosu/A02/$1"
if [ $? -ne 0 ]; then
    echo "Directory $dir1 does not exist or cannot be accessed at MIT."
    exit 1
fi

xrdfs "root://cmsxrootd-kit-disk.gridka.de:1094/" stat "/store/user/jhornung/data_raw/$1"
if [ $? -ne 0 ]; then
    echo "Directory $dir1 does not exist or cannot be accessed at GridKA."
    exit 1
fi

# Dateien in Verzeichnis 1 finden
files_dir1=$(xrdfs "root://xrootd.cmsaf.mit.edu:1094/" ls -l "/store/user/paus/nanosu/A02/$1" | awk '{print $NF}' | awk -F '/' '{print $NF}')

# Dateien in Verzeichnis 2 finden
files_dir2=$(xrdfs "root://cmsxrootd-kit-disk.gridka.de:1094/" ls -l "/store/user/jhornung/data_raw/$1"| awk '{print $NF}' | awk -F '/' '{print $NF}')

# Unterschiede in den Dateien finden
diff_files=$(comm -13 <(echo "$files_dir2" | sort) <(echo "$files_dir1" | sort))

# Ausgabe in Textdatei schreiben
echo "$diff_files" > "$output_file"

echo "Die Differenz der Dateien zwischen $dir1 und $dir2 wurde in $output_file gespeichert."

