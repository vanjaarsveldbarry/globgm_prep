#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH -J gwRecharge_cf

source ${HOME}/.bashrc
mamba activate globgm

saveDirectory=$1  
sim=$2

predRechargeFile=/scratch/depfg/7006713/globgm_input/_tempGLOBGM/predicted_recharge.nc

scenarios=(
        "historical_natural" 
        # "historical" 
    )
for scen in "${scenarios[@]}"; do
    subSaveDir=$saveDirectory/$sim/$scen
    echo $subSaveDir $sim $predRechargeFile
    python -u /eejit/home/7006713/projects/globgm_prep/3_calculate_gwRecharge/calculate_cf.py $subSaveDir $sim $predRechargeFile
done