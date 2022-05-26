# LepPheno_BestPractices

Scripts and data associated with manuscript titled "Phenological research based on natural history collections: practical guidelines and a Lepidopteran case study". 

If you download or clone this repository, you can fully replicate LMM and PLMM results and figures. The entire workflow can also be reproduced if you download and unzip the raw occurrence records from our OSF project (https://osf.io/wdzay/). These zip files should be unzipped into the /data subdirectory for scripts to run properly. 

## Overview of repository

## data

### climate
Contains csv of climate data used at our cell basis for modeling.

### FFP
Contains .tif file that has the mean beginning date of the frost free period across north america. Used to caclulate our seasonality trait.

### LMM_Data
Contains "model dataframes" that were used in statistical analyses. 
- mdf.csv is a csv with estimated phenology metrics but outliers remain.
- mdf_removeOutliersResiduals.csv has the outliers removed from mdf.csv.
- mdf_removeOutliersResiduals_wSeasonality.csv has the outliers removed and the seasonality trait calculated. This is the dataset used in our statistical analyses.

### phylogeny
Contains .csv and .tre files used to generate phylogy used in analses. insect_tree_wBranches.tre is the tree used in analyses.

### traits
Contains .csvs used to generate species list and associated traits. 

## scripts

### buildingPhylogenies
Script used to build phylogeny. Note script will not run on Windows. Successfully run on a linux machine.

### dataCleaning
Scripts used to clean data prior to analyses, listed in sequential order. Between Step 03 and 04, phenoestimates were made using the scripts/phenoEstimates/generate_phenometrics_try2.R script.

### figures
Scripts used to make manuscript figures. Interaction figures were made by first making individual interaction plots using the /scripts/figures/interactionScripts_revision scripts. 

### phenoEstimates
Script to generate phenological estimates. Note this would take a long time to run locally. I ran this script on a 64GB & 40 core linux cluster.

### statisticalAnalyses
Scripts used to generate LMMs and PLMMs used in manuscript. 

## outputs
Csv that saves the output of the script used to generate phenometrics.

## singleInteractionFigures
.Rdata of single interaction figures. Running the scripts to make single interaction figures takes some time, so figures with interactions were made by calling in these .Rdata files. 
