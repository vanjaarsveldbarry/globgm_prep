#!/bin/bash
#SBATCH --partition=fat_genoa
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --output=/home/bvjaarsv/projects/globgm_prep/calculate_gwRecharge/_getData/get30secRecharge.out

source ${HOME}/.bashrc
mamba activate globgm

# python /home/bvjaarsv/projects/globgm_prep/calculate_gwRecharge/_getData/get30secRecharge.py
python /projects/0/einf4705/gwRecharge_correction/test.py 