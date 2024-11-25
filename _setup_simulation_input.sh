#!/bin/bash -l 

#Run the pre pre stuff so that you dont wait so long
saveDirectory="/scratch/depfg/7006713/globgm_input/cmip6_input"
tempDir="/scratch/depfg/7006713/globgm_input/_tempGLOBGM"
simulations=(
    # "gswp3-w5e5" RUN DONE - UPLOAD DONE
    # "gfdl-esm4" RUN DONE - UPLOAD DONE
    # "ipsl-cm6a-lr" RUN DONE - UPLOAD DONE
    # "mpi-esm1-2-hr" RUN BUSY
    # "mri-esm2-0"  RUN BUSY
    # "ukesm1-0-ll" TODO
)

cd "$(dirname "$0")"
for sim in "${simulations[@]}"; do
    sbatch -o /eejit/home/7006713/projects/globgm_prep/_natural_runs/slurmOut/${sim}.out ./_natural_runs/config/run_natural.slurm $saveDirectory $sim $tempDir
done

# saveDirectory="/scratch/depfg/7006713/globgm_input/cmip6_input"
# tempDir="/scratch/depfg/7006713/globgm_input/_tempGLOBGM"
# simulations=(
#     "gswp3-w5e5" 
#     # "gfdl-esm4" 
#     # "ipsl-cm6a-lr" 
#     # "mpi-esm1-2-hr"
#     # "mri-esm2-0" 
#     # "ukesm1-0-ll" 
# )

# cd "$(dirname "$0")"
# for sim in "${simulations[@]}"; do
#     # bash ./1_download_cmip6_data/download_cmip6_data.sh $saveDirectory $sim
#     # sbatch -o /eejit/home/7006713/projects/globgm_prep/2_natural_runs/slurmOut/${sim}.out ./2_natural_runs/config/run_natural.slurm $saveDirectory $sim $tempDir
#     bash ./2_natural_runs/config/run_natural.slurm $saveDirectory $sim $tempDir
#     # sbatch -o /eejit/home/7006713/projects/globgm_prep/3_calculate_gwRecharge/slurmOut/${sim}_cf.out ./3_calculate_gwRecharge/calculate_cf.slurm $saveDirectory $sim $tempDir
#     # sbatch -o /eejit/home/7006713/projects/globgm_prep/4_saturatedAreaFraction/slurmOut/${sim}.out ./4_saturatedAreaFraction/calc_sat_area_fraction_globgm.slurm $saveDirectory $sim $tempDir
# done