# BVCMM4GT: Bayesian Varying Coefficients Mixed Models for Group Testing data
This repository contains R/Rcpp programs for the article: **“Bayesian Varying Coefficients Mixed Models for Group Testing Data”** by Yizeng Li and Dewei Wang, which has been submitted for publication. 

## Getting Started

### Dependencies
A list of required packages is available:
- BayesLogit"
- geoR
- ltsa
- mvtnorm
- Matrix
- hdf5r
- Rcpp
- glue
- Hmisc

Please refer to the [install and upgrade documentation on R software official website](https://www.r-project.org/) for all available installation methods.

### Installing

[UNIX/MAC] Open terminal and install it via git
```sh
~$ git clone yizenglistat/BVCMM4GT
~$ cd BVCMM4GT
~$ pwd # expected to see USER_SAVED_PATH/BVCMM4GT
```

[WIN] Download the zip file from this repository

### Reproduce simulation

- Step 0: Make sure your current path is `USER_SAVED_PATH/BVCMM4GT`
```sh
~$ pwd
```

- Step 1: Run simulation.R locally (not recommendded)
```r
# this will run 500 repetitions for all possible scenarios in setting1
# and it is very computationally heavy (>24 hours)
> source("simulation.R") # all saved .RData files will be saved under output folders
```

- Step 1 (alternative): Run setting.R on cluster (recommend)
Suppose we had 125 nodes available, identified id as `1,2,...,125`. Within each node, we only need to run 4 repetitions for all possible scenariors in setting1. In total, we still have completed 500 runs. 
```sh
~$ Rscript simulation.R NODE_ID
``` 

- Step 2: Generate figures
After gathering all saved `.RData` files from Step 1, the figures are generated in R by Python3
```sh
# first argument is string, "known" or "unknown", for Se, Sp
# second argument is sample size, 3000 or 5000, for number of individuals
# third argument is name of model, "m1" or "m2" for different sets of functions in simulations

# example 1: generate graph for known Se Sp, N=5000, M1 (with DT, AT, or IT and cj=5 or 10) 
~$ python3 figures.py known 5000 m1

# example 2: generate graph for unknown Se Sp, N=3000, M2 (with DT, AT, or IT and cj=5 or 10) 
~$ python3 figures.py unknown 3000 m2
```

- Step 3: Generate tables
After gathering all saved `.RData` files from Step 1, the tables are rendered to `.tex` files: `tab_unknown.tex` and `tab_known.tex`

## Authors

[Yizeng Li](https://yizengli.com)
[Dewei Wang](https://sites.google.com/view/deweiwang)
