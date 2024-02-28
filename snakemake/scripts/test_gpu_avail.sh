#!/bin/bash
#SBATCH -N 1 #find requested resources on a single node
#SBATCH -p gpu-el8 #queue where our gpu machines are 
#SBATCH -G 1 #ask for num gpus
#SBATCH --mem=10G #ask for memory
#SBATCH --exclude=gpu[38-39]
#SBATCH -n 2 #ask for num tasks per allocated gpu

# source ~/.bashrc
# conda init
# source activate scvi-env 

# script
python scripts/test_gpu_avail.py

