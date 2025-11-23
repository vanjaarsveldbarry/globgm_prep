# GLOBGM Input Data Preparation

This repository contains scripts and workflows for preparing input data for GLOBGM (Global Groundwater Model) simulations using CMIP6 climate data and PCR-GLOBWB model outputs.

## Overview

The repository provides an automated workflow to:
1. Download CMIP6 climate data from online repositories
2. Run natural (pristine) groundwater recharge simulations
3. Apply regression-based corrections to groundwater recharge estimates
4. Calculate saturated area fractions for the model


## Repository Structure

### Main Workflow Steps

1. **`1_download_cmip6_data/`** - Downloads CMIP6 climate forcing data

2. **`2_download_natural_runs/`** - Downloads natural run outputs from iRODS

3. **`3_gwRech_correction/`** - Applies corrections to groundwater recharge data

4. **`4_saturatedAreaFraction/`** - Calculates saturated area fractions
