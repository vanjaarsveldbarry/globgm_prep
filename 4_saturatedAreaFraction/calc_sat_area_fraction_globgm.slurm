#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -J sat_area_frac
#SBATCH -t 24:00:00

source ${HOME}/.bashrc
mamba activate globgm

saveDirectory=$1
simulation=$2
tempDir=$3/sat_area_frac
mkdir -p $tempDir

if [ "$simulation" == "gswp3-w5e5" ]; then
    scenarios=(
        "historical" 
        "historical_natural" 
    )
else
    scenarios=(
        "historical_natural" 
        "historical" 
        "ssp126" 
        "ssp370" 
        "ssp585" 
    )
fi

for scen in "${scenarios[@]}"; do
    if [[ "$scen" == *"historical"* ]]; then
        if [ "$simulation" == "gswp3-w5e5" ]; then
            startDate=1960-01-01
            endDate=2019-12-31
        else
            startDate=1960-01-01
            endDate=2014-12-31
        fi
    else
        startDate=2015-01-01
        endDate=2100-12-31
    fi
    outDir=$tempDir/${simulation}/${scen}/

    dataInputDir=$saveDirectory/$simulation/$scen
    saveFolder=$saveDirectory/$simulation/$scen/

    monthly_s1_file=$(find $dataInputDir -type f -name '*storUppTotal_global_monthly-average*.nc' | head -n 1)
    monthly_s2_file=$(find $dataInputDir -type f -name '*storLowTotal_global_monthly-average*.nc' | head -n 1)
    iniFile="/eejit/home/7006713/projects/globgm_prep/4_saturatedAreaFraction/calc_sat_area_globgm.ini"
    cd /eejit/home/7006713/projects/globgm_prep/4_saturatedAreaFraction/pgb_sat_area_frac/scripts
    python deterministic_runner_for_calculating_saturated_area_fraction.py $iniFile $outDir $startDate $endDate $monthly_s1_file $monthly_s2_file
    cp $outDir/sat_area_fraction.nc  $saveFolder/sat_area_fraction_monthly.nc 
    cdo yearsum $saveFolder/sat_area_fraction_monthly.nc  $saveFolder/sat_area_fraction_annual_sum.nc 
    cdo timmean $saveFolder/sat_area_fraction_annual_sum.nc  $saveFolder/sat_area_fraction_average.nc.temp 
    ncks -O -C -x -v time,time_bnds $saveFolder/sat_area_fraction_average.nc.temp $saveFolder/sat_area_fraction_average.nc.temp1
    ncwa -O -a time $saveFolder/sat_area_fraction_average.nc.temp1 $saveFolder/sat_area_fraction_average.nc
    rm -r $outDir
    rm $saveFolder/sat_area_fraction_annual_sum.nc
    rm $saveFolder/sat_area_fraction_average.nc.temp 
    rm $saveFolder/sat_area_fraction_average.nc.temp1 
    wait
done
