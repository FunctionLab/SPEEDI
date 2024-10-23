# SPEEDI

## Table of Contents

- [Overview of SPEEDI](#overview-of-speedi)
- [Using the SPEEDI Website](#using-the-speedi-website)
- [Running SPEEDI Locally](#running-speedi-locally)
- [Citing SPEEDI](#citing-speedi)
- [Need Help?](#need-help)

## Overview of SPEEDI

![Overview of SPEEDI](https://github.com/FunctionLab/SPEEDI/blob/main/SPEEDI_overview.png?raw=true)

Single-cell Pipeline for End to End Data Integration (SPEEDI) is a fully automated, end-to-end pipeline that facilitates single cell data analysis and improves robustness and reproducibility. SPEEDI computationally infers batch labels and automates the application of state of the art processing and analysis tools. Additionally, SPEEDI implements a reference-based cell type annotation method coupled with a majority-vote system. SPEEDI takes raw count feature-by-barcode single cell data matrices as input and outputs an integrated and annotated single-cell object, a log file with auto-selected analysis parameters, and a set of preliminary analyses.

## Using the SPEEDI Website

The [SPEEDI Website](https://speedi.princeton.edu/) allows users to upload their single cell datasets to our server for processing. Users can then view and download results once processing completes. Please visit the website to learn more!

## Running SPEEDI Locally

To install the SPEEDI R package locally, you can use `devtools` and `BiocManager`:

```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github('FunctionLab/SPEEDI', repos = BiocManager::repositories())
```

All R-related dependencies should be installed automatically. Note that [RTools](https://cran.r-project.org/bin/windows/Rtools/) is required to install the SPEEDI R package in Windows. To learn how to use the SPEEDI R package, please view the [SPEEDI vignette](https://speedi.princeton.edu/vignette).

## Citing SPEEDI

The SPEEDI paper was published in Cell Systems on October 16, 2024. 

**Citation**: Wang Y, Thistlethwaite W, Tadych A, Ruf-Zamojski F, Bernard DJ, Cappuccio A, Zaslavsky E, Chen X, Sealfon SC, Troyanskaya OG. Automated single-cell omics end-to-end framework with data-driven batch inference. Cell Systems. 2024 Oct 16;15(10):982-990.e5.

## Need Help?

If you encounter any issues using SPEEDI, feel free to contact a SPEEDI administrator (speedi@genomics.princeton.edu).
