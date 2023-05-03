# SPEEDI

## Table of Contents

- [Overview of SPEEDI](#overview-of-speedi)
- [Using the SPEEDI Website](#using-the-speedi-website)
- [Running SPEEDI Locally](#running-speedi-locally)
- [Citing SPEEDI](#citing-speedi)

## Overview of SPEEDI

Single-cell Pipeline for End to End Data Integration (SPEEDI) is a fully automated, end-to-end pipeline that facilitates single cell data analysis and improves robustness and reproducibility. SPEEDI computationally infers batch labels and automates the application of state of the art processing and analysis tools. Additionally, SPEEDI implements a reference-based cell type annotation method coupled with a majority-vote system. SPEEDI takes raw count feature-by-barcode single cell data matrices as input and outputs an integrated and annotated single-cell object, a log file with auto-selected analysis parameters and a set of preliminary analyses.

## Using the SPEEDI Website

We are currently developing a website where users can upload their single cell datasets for processing through SPEEDI and then view and download results. Coming soon! 

## Running SPEEDI Locally

To install SPEEDI locally, you can use `devtools`:

```
library(devtools)
devtools::install_github('FunctionLab/SPEEDI')
```

All dependencies will be installed automatically.

## Citing SPEEDI

The SPEEDI manuscript is currently in preparation.
