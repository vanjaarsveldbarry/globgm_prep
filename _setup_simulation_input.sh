#!/bin/bash -l 

cd "$(dirname "$0")"
iRodsPassword="xy6mTpxsvDogAVHq0JVLzQz1_LMD5IEc"

#Run the pre pre stuff so that you dont wait so long
saveDirectory="/scratch/depfg/7006713/globgm_input/cmip6_input"
tempDir="/scratch/depfg/7006713/globgm_input/_tempGLOBGM"
simulations=(
    # "gswp3-w5e5" 
    # "gfdl-esm4"
    # "ipsl-cm6a-lr"
    # "mpi-esm1-2-hr"
    # "mri-esm2-0"
    # "ukesm1-0-ll"
)

# sbatch -o /eejit/home/7006713/projects/globgm_prep/_calculate_gwRecharge/slurmOut/_calculate_gwRecharge.out ./_calculate_gwRecharge/createCorrectedRecharge.slurm

# for sim in "${simulations[@]}"; do
#     if [ "$sim" == "gswp3-w5e5" ]; then
        # sbatch -o /eejit/home/7006713/projects/globgm_prep/_natural_runs/slurmOut/${sim}.out ./_natural_runs/config/run_natural.slurm $saveDirectory $sim $tempDir
    # fi
# done
# for sim in "${simulations[@]}"; do
# done

# for sim in "${simulations[@]}"; do
    # mkdir -p ./3_gwRech_correction/slurmOut ./4_saturatedAreaFraction/slurmOut 
    # bash ./1_download_cmip6_data/download_cmip6_data.sh $saveDirectory $sim #BUSY  screen -rd dw
    # for sim in "${simulations[@]}"; do
    #     if [ "$sim" == "gswp3-w5e5" ]; then
    #         python ./2_download_natural_runs/2_download_natural_runs.py $saveDirectory $sim $iRodsPassword
    #     fi
    # done
    # sbatch -o ./3_gwRech_correction/slurmOut/${sim}_cf.out ./3_gwRech_correction/3_gwRech_correction.slurm $saveDirectory $sim $iRodsPassword
    # sbatch -o ./4_saturatedAreaFraction/slurmOut/${sim}.out ./4_saturatedAreaFraction/calc_sat_area_fraction_globgm.slurm $saveDirectory $sim $tempDir
# done