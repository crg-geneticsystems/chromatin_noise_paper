# Overview

Welcome to the GitHub repository for following publication:
["Systematic Analysis of the Determinants of Gene Expression Noise in Embryonic Stem Cells" (Faure AJ, Schmiedel JM, Lehner B, Cell Systems, doi:10.1016/j.cels.2017.10.003)](https://www.cell.com/cell-systems/abstract/S2405-4712(17)30440-4)

Here you'll find all the scripts needed to reproduce the main text figures and results from the computational analyses described in the paper.

The pipeline relies heavily on the [SCDE](https://github.com/hms-dbmi/scde) package which implements a set of statistical methods for analyzing single-cell RNA-seq data.

# Pipeline

The first three stages (01 - 03) perform variance normalisation, cell cycle correction and gene set over-dispersion analysis using the PAGODA framework (contained in the [SCDE](https://github.com/hms-dbmi/scde) package).

To run this pipeline, you will first need to download matrices of read/molecule counts in single cells as well as custom gene annotation files: 

* **Read counts** for each experiment ([chromatin_noise_paper_raw.zip](https://www.dropbox.com/s/13a3ddosp4ojopc/chromatin_noise_paper_raw.zip?dl=0))
* **Supplementary data** including custom gene annotations ([chromatin_noise_paper_misc.zip](https://www.dropbox.com/s/5q4ch3y7a5szgjb/chromatin_noise_paper_misc.zip?dl=0))

Once unzipped, their contents should be transferred to repository base path ("chromatin_noise_paper").

The second three stages (04 - 06) produce plots appearing in the above publication. To run these scripts, you will once again need the Supplementary data files ([chromatin_noise_paper_misc.zip](https://www.dropbox.com/s/5q4ch3y7a5szgjb/chromatin_noise_paper_misc.zip?dl=0)) as well as the output of stages 01 - 03.

Alternatively, to skip the first three stages and reproduce the plots only, you can simply download the required output of stages 01 - 03:

* **Processed data** for each experiment ([chromatin_noise_paper_processed.zip](https://www.dropbox.com/s/wurwzs5ihwngwhk/chromatin_noise_paper_processed.zip?dl=0))

Finally, you can also simply download the resulting plots produced by stages 04 - 06:

* **Final plots** summarising results from all analyses ([chromatin_noise_paper_plots.zip](https://www.dropbox.com/s/3fyxd38owur5d3u/chromatin_noise_paper_plots.zip?dl=0))

## 1. Variance normalisation

[01\_pagoda\_varnorm.R](./01_pagoda_varnorm.R) performs read count filtering, construction of cell-specific error models and normalisation of gene expression variances relative to transcriptome-wide expectations. Command-line options specific to each dataset are provided in [01\_pagoda\_varnorm.sh](./01_pagoda_varnorm.sh).

## 2. Cell cycle control

[02\_pagoda\_cellcyclecontrol.R](./02_pagoda_cellcyclecontrol.R) controls for cell cycle-related variation using the set of genes annotated to the GO biological process term "cell cycle". This procedure is repeated iteratively until cell cycle gene covariability is indistinguishable from that expected by chance. Command-line options specific to each dataset are provided in [02\_pagoda\_cellcyclecontrol.sh](./02_pagoda_cellcyclecontrol.sh).

## 3. Gene Set Over-Dispersion and Total Noise Level Bias Analysis

[03\_pagoda\_adjvariance\_bias.R](./03_pagoda_adjvariance_bias.R) performs chromatin/promoter type gene set and pathway covariability as well as total noise level bias analysis. Command-line options specific to each dataset are provided in [03\_pagoda\_adjvariance\_bias.sh](./03_pagoda_adjvariance_bias.sh).

## 4. Total Noise Level Bias Plots

[04\_adjvariance\_bias\_plots.R](./04_adjvariance_bias_plots.R) produces chromatin/promoter type gene set total noise level bias and binomial smooth plots.

## 5. Integrated Noise Model Analysis and Plots

[05\_integratedmodels\_plots.R](./05_integratedmodels_plots.R) performs all integrated model analyses and produces summary plots of model coefficients.

## 6. mRNA Stability Meta-analysis and Plots

[06\_mRNAstability\_plots.R](./06_mRNAstability_plots.R) performs mRNA stability meta-analysis and produces associated plots.



