#!/bin/bash
#SBATCH --partition=genoa
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --output=/home/bvjaarsv/projects/globgm_prep/calculate_gwRecharge/slurmOut/createCorrectedRecharge.out

source ${HOME}/.bashrc
mamba activate globgm

# python /home/bvjaarsv/projects/globgm_prep/calculate_gwRecharge/preprocessData.py
# wait 
# echo "Done with preprocessData.py"

# python /home/bvjaarsv/projects/globgm_prep/calculate_gwRecharge/createTrainTest.py
# wait 
# echo "Done with createTrainTest.py"

python -u /home/bvjaarsv/projects/globgm_prep/calculate_gwRecharge/finalModel.py
wait 
echo "Done with finalModel.py"


