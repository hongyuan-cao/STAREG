# STAREG

Spatial transcriptomic analysis of replicable expressed genes

## Overview

STAREG is an empirical Bayesian method for identifying replicable spatially variable genes in data generated from various spatially resolved transcriptomic techniques. It provides effective control of false discovery rate and higher power for the replicability analysis by borrowing information across genes and across different studies, and is scalable to datasets with tens of thousands of genes measured on tens of thousands of spatial locations.

## Installation

```R
## Install dependency packages if necessary
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")
install.packages(c("Iso"))

## Install STAREG
install.packages("devtools")
devtools::install_github("hongyuan-cao/STAREG")
# An alternative installation method
install.packages("STAREG")

## Load STAREG
library(STAREG)
```

## An numeric example

We illustrate the usage of STAREG for replicability analysis of two large-scale multiple testing problem using simulated data.

```R
## Pre-specify the number of hypotheses, m, the prior probabilities of the joint hidden states, xi's, and the alternative settings
m = 10000; xi00 = 0.9; xi01 = 0.025; xi10 = 0.025; xi11 = 0.05
mu1 = 2; mu2 = 2.5; sigma1 = 1; sigma2 = 1

## Generate the hidden states and corresponding p-values in two studies 
h = sample(0:3, m, replace = TRUE, prob = c(xi00, xi01, xi10, xi11))
states1 = rep(0, m)
states1[c(which(h==2), which(h==3))] = 1
states2 = rep(0, m)
states2[c(which(h==1), which(h==3))] = 1

stat1 = rnorm(m, states1*mu1, sigma1)
stat2 = rnorm(m, states2*mu2, sigma2)

p1 = 1 - pnorm(stat1, mean = 0, sd = sigma1)
p2 = 1 - pnorm(stat2, mean = 0, sd = sigma2)

## Replicability analysis
library(STAREG)
alpha <- 0.05
rep.obj <- stareg(p1, p2)
rep.genes <- which(rep.obj$fdr <= alpha)
```

## Data and reproducibility

All the R functions and scripts to reproduce the numeric simulation results in the manuscript are contained in “simulation” folder.

The SRT data used in the real data analysis can be downloaded from the links provided in the *Data availability* section in the manuscript, and an example of the replicability analysis of mouse olfactory bulb data is provided in the “real_data” folder for illustration of the usage. 

More details on real data analysis can be found [here](https://github.com/YanLi15/STAREG-Analysis).
