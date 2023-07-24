# SPEEDI

## Table of Contents

- [Overview of SPEEDI](#overview-of-speedi)
- [Using the SPEEDI Website](#using-the-speedi-website)
- [Running SPEEDI Locally](#running-speedi-locally)
- [Citing SPEEDI](#citing-speedi)

## Overview of SPEEDI

*Important Note*: SPEEDI is currently a work in progress. If you encounter any issues, feel free to contact William (wat2@princeton.edu).

Single-cell Pipeline for End to End Data Integration (SPEEDI) is a fully automated, end-to-end pipeline that facilitates single cell data analysis and improves robustness and reproducibility. SPEEDI computationally infers batch labels and automates the application of state of the art processing and analysis tools. Additionally, SPEEDI implements a reference-based cell type annotation method coupled with a majority-vote system. SPEEDI takes raw count feature-by-barcode single cell data matrices as input and outputs an integrated and annotated single-cell object, a log file with auto-selected analysis parameters and a set of preliminary analyses.

## Using the SPEEDI Website

We are currently developing a website where users can upload their single cell datasets for processing through SPEEDI and then view and download results. Coming soon! 

## Running SPEEDI Locally

To install SPEEDI locally, you can use `devtools` and `BiocManager`:

```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github('FunctionLab/SPEEDI', repos = BiocManager::repositories())
```

All R-related dependencies should be installed automatically. After installing, the easiest way to use SPEEDI is with the `run_SPEEDI()` wrapper function:

```
library(SPEEDI)
# Learn more about the SPEEDI pipeline wrapper function
?run_SPEEDI
# Example parameters for run_SPEEDI - note that some optional parameters 
# (metadata_df, reference_file_name, analysis_name, and sample_id_list) were not used
reference_tissue <- "PBMC"
data_type <- "RNA"
species <- "human"
data_path <- "~/test_input"
reference_dir <- "~/references"
output_dir <- "~/test_output"
record_doublets <- FALSE
# Run SPEEDI pipeline
run_SPEEDI(reference_tissue = reference_tissue, data_type = data_type, species = species, data_path = data_path, reference_dir = reference_dir, 
           output_dir = output_dir, record_doublets = record_doublets)
```

You can also run the individual steps of the pipeline separately. Read through the documentation for `run_SPEEDI()` to learn more!

## Citing SPEEDI

The SPEEDI manuscript is currently in preparation.
