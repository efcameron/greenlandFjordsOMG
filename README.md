
# README: greenlandFjordsOMG
Author: Emma Cameron, University of St Andrews
Date: 08 Apr 2026

## Overview 
This repository contains scripts used to process and plot CTD data from the Oceans Melting Greenland (OMG) project (OMG, 2019,2020). The code reproduces the figures from the paper:

"Greenland fjord processes have a depth-dependent influence on predicted submarine melt rates"
Submitted to Geophysical Research Letters (April 2026).

The analysis focuses on comparing CTD profiles near glacier termini
to water properties at fjord mouths.

This allows investigation of:

1. Differences in water properties between near-glacier CTDs and fjord mouth CTDs
2. The resulting impact on predicted submarine melt rates.


## Project Folder structure

This repository is designed to sit within a larger working directory: 

Fjord_Shelf_Observations/
|-- Working_data
|   |-- OMG_data      # put downloaded CTD and AXCTD files here
|       |-- GSW_CTDs     # put TEOS-10 GSW converted netCDF files here
|       |-- greenland       # put csv files here
|           |-- klu_fjord_data.csv
|           |-- qim_fjord_data.csv 
|           |-- fjord_data.csv
|           |-- PDS_profiles.csv
|           |-- allFilenames.txt <-- list of all the AXCTD and CTD filenames
|-- scripts 
|   |--greenlandFjordsOMG    # <-- This repository
|      |-- README.md
|      |-- environment.yml 
|      |-- config.yaml
|      |-- config_loader.py
|      |-- gsw_transformation_of_nc_files.ipynb
|      |-- fig2_anomalyPlots.ipynb
|      |-- fig3_TFandPDS.ipynb
|      |-- fig_4_caseStudies.ipynb
|      |-- pds_histograms.ipynb
|-- Results
|   |-- thermalForcing
|   |-- deltaProfiles
|   |-- caseStudies_qimKlu

## Configuration

File paths are managed through config.yaml:

paths:
  omg_dir: "Working_data/OMG_data"
  csv_dir: "Working_data/OMG_data/greenland"
  nc_dir: "Working_data/OMG_data/GSW_CTDs"
  results_dir: "Results"

These paths are resolved relative to a base directory, defined in config_loader.py as:

base_dir = Path(__file__).resolve().parents[2]

This means:

The repository assumes it lives inside the folder structure shown above
Paths are constructed automatically from this base directory.

## IMPORTANT ## 
You will likely need to modify config.yaml if your folder structure differs from this

## Data requirements

This repository does not include raw data.

## netCDF data
CTD data is available as netCDF files at: https://doi.org/10.5067/OMGEV-AXCT1
AXCTD data is available as netCDF files at: https://doi.org/10.5067/OMGEV-CTDS1

You will need:

AXCTD and CTD NetCDF files from the OMG dataset should be downloaded and placed in: 
Working_data/OMG_data


The temperature and salinity values in these netCDF files should be converted to conservative temperature and absolute salinity using the Thermodynamic Equation of Seawater 2010 (TEOS-10), using the GSW-Python library (https://teos-10.github.io/GSW-Python/). The script to do this from a .txt file of the OMG filenames is "gsw_transformation_of_nc_files.ipynb".

This script will place updated copies of the netCDF files with absolute salinity and conservative temperature in the folder:
Working_data/OMG_data/GSW_CTDs_v3 

## CSV Files 
Preprocessed CSV files should be placed in:
Working_data/OMG_data/greenland/

Expected CSV files include:

klu_fjord_data.csv 
--> contains data on Kangerluluk Fjord

qim_fjord_data.csv 
--> contains data on Qarassap Imaa Fjord

fjord_data.csv 
--> contains data on all fjords in the study and matches CTD pairs 

PDS_profiles.csv 
--> long-form CSV that contains PDS (%) values for each fjord at each depth increment, for beta=1.2 and beta=1.6.

allFilenames.txt
--> Text file listing all the CTD and AXCTD filenames to allow them to be converted to absolute salinity and conservative temperature.


CSV and txt files are publicly accessible from: DOI for PURE repository XXX

## Usage

The analysis is done using Jupyter notebooks. 

1. gsw_transformation_of_nc_files.ipynb 
--> Takes the locally downloaded .nc files from the OMG AXCTDs and CTDs, converts in-situ temperature and practical salinity to conservatove temperature and absolute salinity and exports new .nc files with these variables to a separate directory, ready for analysis. This script is applied to all files, consecutive scripts isolate the data used for the study. 

2. fig2_anomalyPlots.ipynb 
--> calculates and plots the temperature-difference and salinity-difference plots from Figure 2 in the accompanying paper. 

3. fig3_TFandPDS.ipynb 
--> calculates and plots profiles of thermal forcing and the percentage difference in submarine melt rate (PDS) between the near-glacier CTDs and the fjord mouth CTDs. Also exports PDS profiles for all fjords in one csv file (PDS_profiles.csv)

4. fig_4_caseStudies.ipynb 
--> Plots the CTD data and temperature-salinity diagrams for two example fjords.

5. pds_histograms.ipynb 
--> Calculates the depth-averaged PDS for each profile and plots results on regional and greenland-wide histograms to understand the overall impact of fjord processes on melt rates.

## Example workflow:

1. Ensure data is in the correct directories
2. Update config.yaml if needed
3. Launch Jupyter
4. Set up a new environment by running --> conda env create -f environment.yml
5. Open and run notebooks in order as required

## Dependencies: (see environment.yml)

In the terminal run the line:

conda env create -f environment.yml

## References for the OMG dataset:

OMG. 2020. OMG CTD Conductivity Temperature Depth. Ver. 1. PO.DAAC, CA, USA. Dataset accessed [2024-09-26] at https://doi.org/10.5067/OMGEV-CTDS1

OMG. 2019. OMG AXCTD Profiles. Ver. 1. PO.DAAC, CA, USA. Dataset accessed [2024-09-26] at https://doi.org/10.5067/OMGEV-AXCT1

## Reference for the sill and grounding line depth data
Mas e Braga, M. (2025) “Supplementary datasets for Mas e Braga et al., "Controls on fjord temperature throughout Greenland in a reduced-physics model"”. Zenodo. doi:10.5281/zenodo.15880691.

with any additional measurements sourced from:

Morlighem, M. et al. (2022). IceBridge BedMachine Greenland. (IDBMG4, Version 5). [Data Set]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/GMEVBWFLWA7X. Date Accessed 05-06-2024.

## Contact:

Emma Cameron, University of St Andrews
efc4@st-andrews.ac.uk

