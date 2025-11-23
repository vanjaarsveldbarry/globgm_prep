

#!/bin/bash
saveDirectory=$1
simulation=$2
sourceURL="https://geo.public.data.uu.nl/vault-pcrglobwb-cmip6/research-pcrglobwb-cmip6%5B1690540205%5D/original/hypflowsci6_v1.0/output/"

variables=(
    "gwRecharge_global_monthly-total"
    "storLowTotal_global_monthly-average"
    "storUppTotal_global_monthly-average"
    "totalGroundwaterAbstraction_global_monthly-total"
    "totalRunoff_global_monthly-total"
    "precipitation_global_monthly-total"
)

if [ "$simulation" == "gswp3-w5e5" ]; then
    scenarios=("historical")
else
    scenarios=(
        "historical" 
        "ssp126" 
        "ssp370" 
        "ssp585" 
    )
fi

for scen in "${scenarios[@]}"; do
    if [ "$scen" == "historical" ]; then
        startYear=1960
        endYear=2014
        scen_name=$scen

        if [ "$simulation" == "gswp3-w5e5" ]; then
            startYear=1960
            endYear=2019
            scen_name="historical-reference"
        fi
    else
        startYear=2015
        endYear=2100
        scen_name=$scen
    fi
    _saveDir=$saveDirectory/$simulation/$scen
    mkdir -p $_saveDir
    cd $_saveDir
    for var in "${variables[@]}"; do
        (
            dwldURL="${sourceURL}${simulation}/${scen_name}/pcrglobwb_cmip6-isimip3-${simulation}_image-aqueduct_${scen_name}_${var}_${startYear}_${endYear}_basetier1.nc"
            wget -P "$tempSaveDir" "$dwldURL"
        ) &
    done
    if [ "$simulation" == "gswp3-w5e5" ]; then
        break
    fi
    wait
done