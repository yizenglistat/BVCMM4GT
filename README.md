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

### Key arguments
- `task_id`: integer, node id, default `1`.
- `nreps`: integer, number of repetitions, default `500`.
- `knowns`: boolean or a vector of boolean, TRUE or FALSE or both, default `c(TRUE, FALSE)`. 
- `Ns`: integer or a vector of integers, 3000 or 5000 or both, default `c(3000,5000)`.
- `pool_sizes`: integer or a vector of integers, 5 or 10 or both, default `c(5, 10)`.
- `model_names`: string or a vector of strings, m1 or m2 or both, default `c("m1", "m2")`.
- `testings`: string or a vector of strings, DT or AT or IT or any of them, default `c("DT", "AT", "IT")`.
  
### Reproduce simulations

- **Step 0:** Ensure your current path is `USER_SAVED_PATH/BVCMM4GT`
```sh
~$ pwd
```

- **Step 1:** Run `main.R` locally (not recommended)
```r
# This will run 500 repetitions for all possible scenarios
# It is computationally heavy (>24 hours)
> source("main.R") # All saved .RData files will be saved under output folders
```

- **Step 1 (alternative):** Run `main.R` on the cluster (recommended)
  
Assuming 125 nodes are available, identified as `1,2,...,125`. Within each node, only 4 repetitions are needed for all possible scenarios. In total, 500 runs are completed.
```sh
~$ Rscript main.R NODE_ID
``` 

- **Step 2:** Generate figures

Given all saved `.RData` files from Step 1, the figures are generated by Python3 after running the following function to generate required dataframes (automatically saved under `output/` folder)
```r
# load required libraries in main.R before generate figures
> figures(folder="output", 
	pool_sizes=c(5,10), 
	model_names=c("m1", "m2"),
	testings=c("IT", "DT","AT"),
	knowns=c(TRUE, FALSE), 
	Ns=c(3000,5000),
	sigma=0.5,
	nreps=500)
```
Visualizing the posterior estimates with 95% pointwise credible intervals by
```sh
# The first argument is a string, "known" or "unknown", for Se, Sp
# The second argument is the sample size, 3000 or 5000, for the number of individuals
# The third argument is the name of the model, "m1" or "m2" for different sets of functions in simulations

# Example 1: Generate a graph for known Se Sp, N=5000, M1 (with DT, AT, or IT and cj=5 or 10) 
~$ python3 src/figures.py known 5000 m1

# Example 2: Generate a graph for unknown Se Sp, N=3000, M2 (with DT, AT, or IT and cj=5 or 10) 
~$ python3 src/figures.py unknown 3000 m2
```

#### Figure demo

**<p align="center"> N=5000, M1, and known Se, Sp </p>**
![known_m1_wo_MPT](https://github.com/yizenglistat/BVCMM4GT/assets/43308957/273d095a-85e9-465e-a9c0-5cc14885860e)

**<p align="center"> N=5000, M2, and unknown Se, Sp </p>**
![unknown_m2_wo_MPT](https://github.com/yizenglistat/BVCMM4GT/assets/43308957/5a4d7b81-92a6-4890-a768-27b6d42b2391)

- **Step 3:** Generate tables
  
After gathering all saved `.RData` files from Step 1, the tables are rendered to `.tex` files: `tab_unknown.tex` and `tab_known.tex` by
```r
# load required libraries in main.R before generate tables
# see details in the main.R 

# generate the table when Se and Sp are known
> tables(folder="output", 
	pool_sizes=c(5,10), 
	model_names=c("m1", "m2"),
	testings=c("IT", "DT","AT"),
	isknown=TRUE,
	Ns=c(3000,5000),
	sigma=0.5,
	nreps=500)

# generate the table when Se and Sp are unknown
> tables(folder="output",
	pool_sizes=c(5,10),
	model_names=c("m1", "m2"),
	testings=c("IT", "DT","AT"),
	isknown=FALSE,
	Ns=c(3000,5000),
	sigma=0.5,
	nreps=500)
```
##### Table demo (results from N=5000 and unknown Se, Sp case)

<img width="701" alt="demo_tab" src="https://github.com/yizenglistat/BVCMM4GT/assets/43308957/b0065f47-edde-4aaa-abe4-b6ee7655cab6">

### Data Application
We simulated fake data, similar to realistic data, to conduct real data analysis through running the following
```r
> source("application_fake.R")
```

##### Figure demo (results from real data analysis)
![app_vcm](https://github.com/yizenglistat/BVCMM4GT/assets/43308957/73cb97f4-d1e2-4958-8d37-518fbd7eb117)

## Authors

[Yizeng Li](https://yizengli.com) and [Dewei Wang](https://sites.google.com/view/deweiwang)
