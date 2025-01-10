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
scen="historical"
cd /eejit/home/7006713/projects/globgm_prep/3_gwRech_correction

python -u ./_download_recharge.py $saveDirectory $sim $iRodsPassword
wait
cp $saveDirectory/$sim/historical/predicted_recharge.nc $saveDirectory/$sim/historical_natural
wait
python -u ./3_gwRech_correction.py $saveDirectory/$sim/"historical" $sim
wait
python -u ./3_gwRech_correction.py $saveDirectory/$sim/"historical_natural" $sim
wait

if [ "$sim" != "gswp3-w5e5" ]; then
    cp -r $saveDirectory/$sim/historical/gwRecharge_correction_factor.zarr $saveDirectory/$sim/ssp126
    cp -r $saveDirectory/$sim/historical/gwRecharge_correction_factor.zarr $saveDirectory/$sim/ssp370
    cp -r $saveDirectory/$sim/historical/gwRecharge_correction_factor.zarr $saveDirectory/$sim/ssp585
    wait
fi