# BVCMM4GT: Bayesian Varying Coefficients Mixed Models for Group Testing Data

This repository contains R/Rcpp programs for the article: **“Bayesian Varying Coefficients Mixed Models for Group Testing Data”** by Yizeng Li and Dewei Wang, which has been submitted for publication. 

## Getting Started

### Dependencies
A list of required packages is available:
- BayesLogit
- geoR
- ltsa
- mvtnorm
- Matrix
- hdf5r
- Rcpp
- glue
- Hmisc

Please refer to the [install and upgrade documentation on the R software official website](https://www.r-project.org/) for all available installation methods.

### Installing

#### UNIX/MAC
Open the terminal and install it via git
```sh
~$ git clone yizenglistat/BVCMM4GT
~$ cd BVCMM4GT
~$ pwd # Expected to see USER_SAVED_PATH/BVCMM4GT
```

#### Windows
Download the zip file from this repository.

### Reproduce Simulation

- **Step 0:** Ensure your current path is `USER_SAVED_PATH/BVCMM4GT`
```sh
~$ pwd
```

- **Step 1:** Run `main.R` locally (not recommended)
```r
# This will run 500 repetitions for all possible scenarios in Setting1
# It is computationally heavy (>24 hours)
> source("main.R") # All saved .RData files will be saved under output folders
```

- **Step 1 (alternative):** Run `setting.R` on the cluster (recommended)
Assuming 125 nodes are available, identified as `1,2,...,125`. Within each node, only 4 repetitions are needed for all possible scenarios in Setting1. In total, 500 runs are completed.
```sh
~$ Rscript main.R NODE_ID
``` 

- **Step 2:** Generate figures
After gathering all saved `.RData` files from Step 1, the figures are generated in R by Python3
```sh
# The first argument is a string, "known" or "unknown", for Se, Sp
# The second argument is the sample size, 3000 or 5000, for the number of individuals
# The third argument is the name of the model, "m1" or "m2" for different sets of functions in simulations

# Example 1: Generate a graph for known Se Sp, N=5000, M1 (with DT, AT, or IT and cj=5 or 10) 
~$ python3 figures.py known 5000 m1

# Example 2: Generate a graph for unknown Se Sp, N=3000, M2 (with DT, AT, or IT and cj=5 or 10) 
~$ python3 figures.py unknown 3000 m2
```

- **Step 3:** Generate tables
After gathering all saved `.RData` files from Step 1, the tables are rendered to `.tex` files: `tab_unknown.tex` and `tab_known.tex`

## Authors

- [Yizeng Li](https://yizengli.com)
- [Dewei Wang](https://sites.google.com/view/deweiwang)
