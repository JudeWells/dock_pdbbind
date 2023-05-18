#!/bin/bash
#$ -l tmem=24G
#$ -l h_vmem=24G
#$ -l h_rt=16:0:0
#$ -S /bin/bash
#$ -N train_gpu
#$ -wd /SAN/orengolab/nsp13/dock_pdbbind/
#$ -t 1-12
#$ -o /SAN/orengolab/nsp13/dock_pdbbind/src/configuration/cluster_output/
# EXAMPLE SUBMISSION COMMAND:

#These are optional flags but you probably want them in all jobs
#$ -j y
hostname
date
nvidia-smi

cd /SAN/orengolab/nsp13/dock_pdbbind/
# UPDATE FOR YOUR ENVIRONMENT

export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH
source /share/apps/source_files/python/python-3.8.5.source

export PYTHONPATH=/share/apps/python-3.8.5-shared/bin
export PYTHONPATH=/SAN/orengolab/nsp13/dock_pdbbind:$PYTHONPATH
python3 src/train_nn.py ${SGE_TASK_ID}