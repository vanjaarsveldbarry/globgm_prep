#!/bin/bash -l 

cd "$(dirname "$0")"
iRodsPassword="tp7n_6mQru0CnuJZ9I55vgxF53rYeSdV"

#Run the pre pre stuff so that you dont wait so long
saveDirectory="/scratch/depfg/7006713/globgm_input/cmip6_input"
tempDir="/scratch/depfg/7006713/globgm_input/_tempGLOBGM"
simulations=(
    "gswp3-w5e5" #RUN DONE - UPLOAD DONE
    # "gfdl-esm4" RUN DONE - UPLOAD DONE
    # "ipsl-cm6a-lr" RUN DONE - UPLOAD DONE
    # "mpi-esm1-2-hr" RUN DONE - UPLOAD DONE
    # "mri-esm2-0"  RUN DONE  - UPLOAD DONE
    # "ukesm1-0-ll"  #RUN BUSY
)

# sbatch -o /eejit/home/7006713/projects/globgm_prep/_calculate_gwRecharge/slurmOut/_calculate_gwRecharge.out ./_calculate_gwRecharge/createCorrectedRecharge.slurm
# for sim in "${simulations[@]}"; do
#     sbatch -o /eejit/home/7006713/projects/globgm_prep/_natural_runs/slurmOut/${sim}.out ./_natural_runs/config/run_natural.slurm $saveDirectory $sim $tempDir
# done

for sim in "${simulations[@]}"; do
    mkdir -p ./4_saturatedAreaFraction/slurmOut ./3_gwRech_correction/slurmOut/
    # bash ./1_download_cmip6_data/download_cmip6_data.sh $saveDirectory $sim
    # python ./2_download_natural_runs/2_download_natural_runs.py $saveDirectory $sim $iRodsPassword
    # sbatch -o ./3_gwRech_correction/slurmOut/${sim}_cf.out ./3_gwRech_correction/3_gwRech_correction.slurm $saveDirectory $sim $iRodsPassword
    sbatch -o ./4_saturatedAreaFraction/slurmOut/${sim}.out ./4_saturatedAreaFraction/calc_sat_area_fraction_globgm.slurm $saveDirectory $sim $tempDir
done