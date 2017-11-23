#!/usr/bin/env bash

################################################################################
## MIDAS tutorial, following
## https://github.com/snayfach/MIDAS/blob/master/docs/tutorial.md
##
## author: sankaran.kris@gmail.com
## date: 11/23/2017
################################################################################

## download and unzip reference and sample data
cd ../../applications/MIDAS/database/
wget http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz
tar -zxvf midas_db_v1.2.tar.gz

cd ../../programming/research/readings/midas/example/
wget http://lighthouse.ucsf.edu/MIDAS/example.tar.gz
tar -zxvf example.tar.gz

cd ..
export PYTHONPATH=$PYTHONPATH:MIDAS
export PATH=$PATH:MIDAS/scripts
export MIDAS_DB=midas_db_v1.2

## first step of MIDAS workflow (species and gene profiling)
mkdir output
run_midas.py species output/sample_1 -1 example/sample_1.fq.gz
run_midas.py species output/sample_2 -1 example/sample_2.fq.gz

run_midas.py genes output/sample_1 -1 example/sample_1.fq.gz
run_midas.py genes output/sample_2 -1 example/sample_2.fq.gz

## second step of MIDAS workflow (merging)
merge_midas.py species species -i output/sample_1,output/sample_2 -t list
merge_midas.py genes -i output/sample_1,output/sample_2 -t list
