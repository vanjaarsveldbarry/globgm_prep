#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH -J gwRecharge_cf

source ${HOME}/.bashrc
mamba activate globgm

saveDirectory=$1  
sim=$2
iRodsPassword=$3

cd "$(dirname "$0")"

# python -u _download_recharge.py $saveDirectory $sim $iRodsPassword

scenarios=(
        "historical_natural" 
    )
for scen in "${scenarios[@]}"; do
    subSaveDir=$saveDirectory/$sim/$scen
    python -u 3_gwRech_correction.py $subSaveDir $sim
done