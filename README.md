# GLOBGM Input Data Preparation

This repository contains scripts and workflows for preparing input data for GLOBGM (Global Groundwater Model) simulations using CMIP6 climate data and PCR-GLOBWB model outputs.

## Overview

The repository provides an automated workflow to:
1. Download CMIP6 climate data from online repositories
2. Run natural (pristine) groundwater recharge simulations
3. Apply machine learning-based corrections to groundwater recharge estimates
4. Calculate saturated area fractions for the model

These preprocessing steps are essential for preparing high-quality input data for global groundwater modeling studies.

## Repository Structure

The workflow is organized into sequential numbered directories and supporting utilities:

### Main Workflow Steps

1. **`1_download_cmip6_data/`** - Downloads CMIP6 climate forcing data
   - Downloads various hydroclimatic variables (groundwater recharge, storage, abstraction, runoff, precipitation)
   - Supports multiple climate models (GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL)
   - Handles both historical and future scenarios (SSP126, SSP370, SSP585)

2. **`2_download_natural_runs/`** - Downloads natural run outputs from iRODS
   - Retrieves simulation outputs for pristine/natural conditions
   - Uses iRODS data management system for data access

3. **`3_gwRech_correction/`** - Applies corrections to groundwater recharge data
   - Uses machine learning predictions to correct recharge estimates
   - Interpolates and resamples data to required grid resolution
   - Calculates correction factors between predicted and original recharge

4. **`4_saturatedAreaFraction/`** - Calculates saturated area fractions
   - Computes saturation area fraction based on PCR-GLOBWB soil storage states
   - Reference: [van Beek & Bierkens (2009)](https://vanbeek.geo.uu.nl/suppinfo/vanbeekbierkens2009.pdf)
   - Processes both upper and lower soil storage layers
   - Outputs monthly, annual, and average saturation fractions

### Supporting Utilities

- **`_calculate_gwRecharge/`** - Machine learning workflow for recharge correction
  - `getData.py` - Retrieves training data from various sources
  - `preprocessData.py` - Preprocesses and cleans input data
  - `createTrainTest.py` - Creates training and testing datasets
  - `finalModel.py` - Trains the final machine learning model
  - `correct_final.py` - Applies corrections to the final product

- **`_natural_runs/`** - PCR-GLOBWB model configuration for natural runs
  - Contains PCR-GLOBWB model setup and configuration files
  - Scripts for running pristine groundwater simulations

- **`_setup_simulation_input.sh`** - Master orchestration script
  - Coordinates the entire workflow
  - Submits SLURM batch jobs for each processing step
  - Manages simulation loops for multiple climate models and scenarios

## Prerequisites

### Software Dependencies

- **Python 3.x** with packages:
  - `xarray` - Multi-dimensional data analysis
  - `pyinterp` - Spatial interpolation
  - `numpy` - Numerical computing
  - `h5netcdf` - NetCDF file handling
  - `tqdm` - Progress bars
  - `python-irodsclient` - iRODS data access
  
- **Climate Data Operators (CDO)** - For NetCDF file manipulation

- **SLURM** - Job scheduling system (for HPC environments)

- **Mamba/Conda** - Environment management (environment name: `globgm`)

### Data Access

- Access to CMIP6 data repository: `https://geo.public.data.uu.nl/vault-pcrglobwb-cmip6/`
- iRODS credentials for accessing natural run outputs
- Sufficient storage space (simulations can be several hundred GB)

## Usage

### Quick Start

The master script `_setup_simulation_input.sh` orchestrates the entire workflow:

```bash
# Edit the script to configure:
# - saveDirectory: Where to store output data
# - tempDir: Temporary processing directory
# - simulations: Which climate models to process
# - iRodsPassword: Your iRODS access credentials

bash _setup_simulation_input.sh
```

### Running Individual Steps

Each numbered directory can also be run independently:

#### 1. Download CMIP6 Data
```bash
cd 1_download_cmip6_data
bash download_cmip6_data.sh <saveDirectory> <simulation>
```

#### 2. Download Natural Runs
```bash
cd 2_download_natural_runs
python 2_download_natural_runs.py <saveDirectory> <simulation> <iRodsPassword>
```

#### 3. Apply Recharge Corrections
```bash
cd 3_gwRech_correction
sbatch 3_gwRech_correction.slurm <saveDirectory> <simulation> <iRodsPassword>
```

#### 4. Calculate Saturated Area Fraction
```bash
cd 4_saturatedAreaFraction
sbatch calc_sat_area_fraction_globgm.slurm <saveDirectory> <simulation> <tempDir>
```

## Supported Climate Models

- **GSWP3-W5E5** (historical only, 1960-2019)
- **GFDL-ESM4**
- **IPSL-CM6A-LR**
- **MPI-ESM1-2-HR**
- **MRI-ESM2-0**
- **UKESM1-0-LL**

## Output Data

The workflow produces several key datasets:

- **Corrected groundwater recharge** - ML-corrected recharge estimates
- **Correction factors** - Spatial correction fields
- **Saturated area fractions** - Monthly, annual sum, and time-averaged fractions
- **Climate forcing variables** - Precipitation, runoff, abstraction, soil storage

All outputs are in NetCDF format with standardized spatial grids.

## Notes

- The workflow is designed for HPC environments with SLURM job scheduling
- Processing times can be significant (up to 24 hours per step for some models)
- Paths in the scripts are configured for specific server environments and may need adjustment
- The workflow uses screen sessions for long-running download tasks

## References

- van Beek, L.P.H., & Bierkens, M.F.P. (2009). The Global Hydrological Model PCR-GLOBWB: Conceptualization, Parameterization and Verification. [Report](https://vanbeek.geo.uu.nl/suppinfo/vanbeekbierkens2009.pdf)
- CMIP6 climate data: [PCR-GLOBWB CMIP6 Archive](https://geo.public.data.uu.nl/vault-pcrglobwb-cmip6/)

## License

See LICENSE file for details.

## Contact

For questions or issues, please contact the repository maintainer or open an issue on GitHub.
