#!/bin/bash

home=$(pwd)

tar -xzf $home/utils.tar.gz 

. $home/utils/haa_setup.sh

python3 $home/utils/bkg_fit.py $1 $2