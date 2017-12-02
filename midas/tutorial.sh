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

cd ../../programming/research/readings/midas/
wget http://lighthouse.ucsf.edu/MIDAS/example.tar.gz
tar -zxvf example.tar.gz

export PYTHONPATH=$PYTHONPATH:/scratch/users/kriss1/applications/MIDAS
export PATH=$PATH:/scratch/users/kriss1/applications/MIDAS/scripts
export MIDAS_DB=/scratch/users/kriss1/applications/MIDAS/database/midas_db_v1.2

## first step of MIDAS workflow (species and gene profiling)
mkdir output
run_midas.py species output/sample_1 -1 example/sample_1.fq.gz
run_midas.py species output/sample_2 -1 example/sample_2.fq.gz

module load biology
module load samtools/1.6
run_midas.py genes output/sample_1 -1 example/sample_1.fq.gz
run_midas.py genes output/sample_2 -1 example/sample_2.fq.gz

## second step of MIDAS workflow (merging)
merge_midas.py species output/species -i output/sample_1,output/sample_2 -t list
merge_midas.py genes output/genes -i output/sample_1,output/sample_2 -t list
