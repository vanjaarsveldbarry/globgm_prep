#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00

source ${HOME}/.bashrc
mamba activate globgm

cd ./_calculate_gwRecharge

Get data from various online/local sources
python -u ./getData.py
wait 
echo "Done with getData.py"

python ./preprocessData.py
wait 
echo "Done with preprocessData.py"

python ./createTrainTest.py
wait 
echo "Done with createTrainTest.py"

python -u ./finalModel.py
wait 
echo "Done with finalModel.py"

wait
rm -r /scratch/depfg/7006713/gwRecharge_corrected/_tempPredict
rm -r /scratch/depfg/7006713/gwRecharge_corrected/_temp


