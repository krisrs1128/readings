#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=nmf-2
#################
#a file for job output, you can check job progress
#SBATCH --output=nmf-2.out
#################
# a file for errors from the job
#SBATCH --error=nmf-2.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=01:00:00
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=6000
#SBATCH --partition=normal,hns
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kriss1@stanford.edu
ml R/3.3.0
Rscript /scratch/users/kriss1/programming/readings/nmf/src/nmf_script.R /scratch/users/kriss1/programming/readings/nmf/batch/config.json 2
