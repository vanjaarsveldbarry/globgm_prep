#!/bin/bash

sim="ipsl-cm6a-lr"
dataFolder=/scratch/depfg/7006713/natural_runs
dagFolder=/nluu14p/home/deposit-pilot/globgm/preprocess/_natural_runs

echo $dataFolder/$sim $dagFolder
iput -rvP $dataFolder/$sim $dagFolder


# K20uIT50MLnzlCm60n3qPJthV9JjM4Wk