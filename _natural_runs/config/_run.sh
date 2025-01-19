#!/bin/bash -l 

cd "$(dirname "$0")"
iRodsPassword="xy6mTpxsvDogAVHq0JVLzQz1_LMD5IEc"

#Run the pre pre stuff so that you dont wait so long
saveDirectory="/scratch/depfg/7006713/globgm_input/cmip6_input"
tempDir="/scratch/depfg/7006713/globgm_input/_tempGLOBGM"
simulations=(
    "gswp3-w5e5" 
    "gfdl-esm4"
    "ipsl-cm6a-lr"
    "mpi-esm1-2-hr"
    "mri-esm2-0"
    "ukesm1-0-ll"
)

for sim in "${simulations[@]}"; do
    sbatch -o /eejit/home/7006713/projects/globgm_prep/_natural_runs/slurmOut/${sim}.out ./_natural_runs/config/run_natural.slurm $saveDirectory $sim $tempDir
done