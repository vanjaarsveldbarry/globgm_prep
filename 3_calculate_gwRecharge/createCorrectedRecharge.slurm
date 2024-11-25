#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --output=/eejit/home/7006713/projects/globgm_prep/3_calculate_gwRecharge/slurmOut/createCorrectedRecharge.out

source ${HOME}/.bashrc
mamba activate globgm

cd /eejit/home/7006713/projects/globgm_prep/3_calculate_gwRecharge

#Get data from various online/local sources
# python ./getData.py
# wait 
# echo "Done with getData.py"

# python ./preprocessData.py
# wait 
# echo "Done with preprocessData.py"

# python ./createTrainTest.py
# wait 
# echo "Done with createTrainTest.py"

python -u ./finalModel.py
wait 
echo "Done with finalModel.py"


