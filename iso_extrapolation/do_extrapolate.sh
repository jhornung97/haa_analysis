#!/bin/bash

home=$(pwd)

tar -xzf $home/utils.tar.gz

. $home/utils/haa_setup.sh

python3 $home/utils/extrapolation.py $1 $2