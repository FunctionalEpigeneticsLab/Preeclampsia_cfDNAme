# Cell-free DNA Methylome Analysis for Early Preeclampsia Prediction

The scripts provide details for the analyses performed in the manuscript entitled 'Cell-free DNA Methylome Analysis for Early Preeclampsia Prediction'.

> This is not the release of a software package. All scripts and binaries are provided as is, without any warranty. We are only providing this information and code in addition to the description of methods for making it easier to reproduce the analyses.

## Versions
```
Python 3.7.6
R version 3.6.2
```

## Bioinformatics processing

TrimGalore (version 0.6.4_dev) was used to trim adapters and low quality bases. Sequencing reads alignment, deduplication, and methylation calling were performed using Bismark suite (v0.22.3). Reads were aligned to human reference GRCh37. 

Sequencing depth was calculated by GATK3.8 `DepthOfCoverage`. Sequencing reads on targets were examined by GATK3.8 `DiagnoseTargets` to identify regions with bad coverage and mapping. 

## Data availability

Bisulfite sequencing data can be downloaded from European Genome-Phenome Archive under EGAS00001007071. Captured regions can be found in the folder `data`.

## Analysis

For differential methylation related analyses, please see the folder `diffmethyl`.

For preeclampsia prediction related analyses, please see the folder `classification`. 
