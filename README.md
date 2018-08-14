# Overview

Welcome to the GitHub repository for following publication:
["Systematic Analysis of the Determinants of Gene Expression Noise in Embryonic Stem Cells" (Faure AJ, Schmiedel JM, Lehner B, Cell Systems, doi:10.1016/j.cels.2017.10.003)](https://www.cell.com/cell-systems/abstract/S2405-4712(17)30440-4)

Here you'll find all the scripts needed to reproduce the main text figures and results from the computational analyses described in the paper.

The pipeline relies heavily on the [SCDE](https://github.com/hms-dbmi/scde) package which implements a set of statistical methods for analyzing single-cell RNA-seq data.

# Pipeline

To run this pipeline, you will first need to download matrices of read/molecule counts in single cells as well as custom gene annotation files: 

* **Read counts** for each experiment ([chromatin_noise_paper_raw.zip](https://www.dropbox.com/s/13a3ddosp4ojopc/chromatin_noise_paper_raw.zip?dl=0))
* **Supplementary data** including custom gene annotations ([chromatin_noise_paper_misc.zip](https://www.dropbox.com/s/5q4ch3y7a5szgjb/chromatin_noise_paper_misc.zip?dl=0))

The first three scripts perform variance normalisation, cell cycle correction and gene set over-dispersion analysis using the PAGODA framework (contained in the [SCDE](https://github.com/hms-dbmi/scde) package):

## 1. Variance normalisation

[01\_pagoda\_varnorm.R](./01_pagoda_varnorm.R) performs read count filtering, construction of cell-specific error models and normalisation of gene expression variances relative to transcriptome-wide expectations. Commandline options specific to each dataset are provided in [01\_pagoda\_varnorm.sh](./01_pagoda_varnorm.sh).







