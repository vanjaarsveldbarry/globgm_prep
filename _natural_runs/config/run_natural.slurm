#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J run_natural
#SBATCH -t 120:00:00

saveDir=$1
simulation=$2
tempDir=$3/natural_runs
tempSaveDir=/scratch/depfg/7006713/natural_runs/$simulation/historical_natural
finalSaveDir=$saveDir/$simulation/historical_natural
mkdir -p $finalSaveDir
mkdir -p $tempSaveDir/netcdf

# #configuration (.ini) file
INI_FILE="/eejit/home/7006713/projects/globgm_prep/_natural_runs/config/setup_05min_historical_version_2021-09-16_naturalized_no-water_no-reservoirs_two-natural-landcovers.ini"

# # # # starting and end dates
STARTING_DATE="1960-01-01"
END_DATE="2019-12-31"

# # # # location/folder, where you will store output files of your 
MAIN_OUTPUT_DIR="${tempDir}/${simulation}"

# # # meteorological forcing files
if [ "$simulation" == "gswp3-w5e5" ]; then
    RELATIVE_HUMIDITY_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3a_version_2021-09-XX/copied_on_2021-09-XX/GSWP3-W5E5/merged/gswp3-w5e5_obsclim_hurs_global_daily_1901_2019_version_2021-09-XX.nc"
    PRECIPITATION_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3a_version_2021-09-XX/copied_on_2021-09-XX/GSWP3-W5E5/merged/gswp3-w5e5_obsclim_pr_global_daily_1901_2019_version_2021-09-XX.nc"
    PRESSURE_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3a_version_2021-09-XX/copied_on_2021-09-XX/GSWP3-W5E5/merged/gswp3-w5e5_obsclim_ps_global_daily_1901_2019_version_2021-09-XX.nc"
    SHORTWAVE_RADIATION_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3a_version_2021-09-XX/copied_on_2021-09-XX/GSWP3-W5E5/merged/gswp3-w5e5_obsclim_rsds_global_daily_1901_2019_version_2021-09-XX.nc"
    WIND_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3a_version_2021-09-XX/copied_on_2021-09-XX/GSWP3-W5E5/merged/gswp3-w5e5_obsclim_sfcwind_global_daily_1901_2019_version_2021-09-XX.nc"
    TEMPERATURE_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3a_version_2021-09-XX/copied_on_2021-09-XX/GSWP3-W5E5/merged/gswp3-w5e5_obsclim_tas_global_daily_1901_2019_version_2021-09-XX.nc"

else
    RELATIVE_HUMIDITY_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3b_version_2021-05-XX/copied_on_2021-06-XX/merged/historical/${simulation}/${simulation}_w5e5_historical_hurs_global_daily_1850_2014.nc"
    PRECIPITATION_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3b_version_2021-05-XX/copied_on_2021-06-XX/merged/historical/${simulation}/${simulation}_w5e5_historical_pr_global_daily_1850_2014.nc"
    PRESSURE_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3b_version_2021-05-XX/copied_on_2021-06-XX/merged/historical/${simulation}/${simulation}_w5e5_historical_ps_global_daily_1850_2014.nc"
    SHORTWAVE_RADIATION_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3b_version_2021-05-XX/copied_on_2021-06-XX/merged/historical/${simulation}/${simulation}_w5e5_historical_rsds_global_daily_1850_2014.nc"
    WIND_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3b_version_2021-05-XX/copied_on_2021-06-XX/merged/historical/${simulation}/${simulation}_w5e5_historical_sfcwind_global_daily_1850_2014.nc"
    TEMPERATURE_FORCING_FILE="/scratch/depfg/sutan101/data/isimip_forcing/isimip3b_version_2021-05-XX/copied_on_2021-06-XX/merged/historical/${simulation}/${simulation}_w5e5_historical_tas_global_daily_1850_2014.nc"
fi

# # # # # initial conditions
MAIN_INITIAL_STATE_FOLDER="/scratch/depfg/sutan101/data/pcrglobwb_input_aqueduct/version_2021-09-16/initial_conditions"
DATE_FOR_INITIAL_STATES="2019-12-31"

# # # number of spinup years
NUMBER_OF_SPINUP_YEARS="25"

# # # location of your pcrglobwb model scripts
PCRGLOBWB_MODEL_SCRIPT_FOLDER=/eejit/home/7006713/projects/globgm_prep/_natural_runs/PCR-GLOBWB_model/model

# # # load the conda enviroment
source ${HOME}/.bashrc
conda activate glob

# # # unset pcraster working threads 
unset PCRASTER_NR_WORKER_THREADS

# # - you may have to activate the following
export OPENBLAS_NUM_THREADS=1

# # go to the folder that contain PCR-GLOBWB scripts
cd ${PCRGLOBWB_MODEL_SCRIPT_FOLDER}

# run the model for all clones, from 1 to 53
for i in {1..53};do
    CLONE_CODE=${i}
    python3 deterministic_runner_with_arguments.py ${INI_FILE} debug_parallel ${CLONE_CODE} -mod ${MAIN_OUTPUT_DIR} -sd ${STARTING_DATE} -ed ${END_DATE} -pff ${PRECIPITATION_FORCING_FILE} -tff ${TEMPERATURE_FORCING_FILE} -presff ${PRESSURE_FORCING_FILE} -windff ${WIND_FORCING_FILE} -swradff ${SHORTWAVE_RADIATION_FORCING_FILE} -relhumff ${RELATIVE_HUMIDITY_FORCING_FILE} -misd ${MAIN_INITIAL_STATE_FOLDER} -dfis ${DATE_FOR_INITIAL_STATES} -num_of_sp_years ${NUMBER_OF_SPINUP_YEARS} &
done


# # process for merging files at the global extent
python3 deterministic_runner_merging_with_arguments.py ${INI_FILE} parallel -mod ${MAIN_OUTPUT_DIR} -sd ${STARTING_DATE} -ed ${END_DATE} &

# # wait until process is finished
wait
echo "end of model runs (please check your results)"

inputDir=$MAIN_OUTPUT_DIR/global/netcdf

cd $inputDir

unique_files=()
for file in *; do 
    unique_file=$(echo "${file}" | rev | cut -d'_' -f4- | rev)
    if [[ ! " ${unique_files[@]} " =~ " ${unique_file} " ]]; then
        unique_files+=("$unique_file")
    fi
done

merge files
index=1
batch_size=15
counter=0
for unique_file in "${unique_files[@]}"; do
    unique_file=$(echo "${unique_file}" | rev | cut -d'_' -f2- | rev)
    varName=$(echo "${unique_file}" | rev | cut -d'_' -f2- | rev)
    accuTimeInfo=$(echo "${unique_file}" | rev | cut -d'_' -f1 | rev)
    accu=${accuTimeInfo: -3}
    timeStr="${accuTimeInfo::-3}"

    if [ "$accu" = "Tot" ]; then
        accu_string="total"
    fi

    if [ "$accu" = "Avg" ]; then
        accu_string="average"
    fi

    if [ "$accu" = "Max" ]; then
        accu_string="maximum"
    fi

    if [ "$timeStr" = "month" ]; then
        time_string="monthly"
    fi

    if [ "$timeStr" = "annua" ]; then
        time_string="yearly"
    fi

    if [ "$simulation" = "gswp3-w5e5" ]; then
        date1="1960"
        date2="2019"
    else
        date1="1960"
        date2="2014"
    fi
    name=pcrglobwb_cmip6-isimip3-${simulation}_image-aqueduct_historical-natural_${varName}_global_${time_string}-${accu_string}_${date1}_${date2}_basetier1
    cdo -L -f nc4 -z zip_5 -mergetime $unique_file*.nc $tempSaveDir/netcdf/$name.nc &
    echo "Processing file $index->79: $unique_file"
    index=$(($index+1))
    counter=$(($counter+1))

    if [ $counter -ge $batch_size ]; then
        wait
        counter=0
    fi

done
wait

# collect initial states
mv $MAIN_OUTPUT_DIR/global/states $tempSaveDir