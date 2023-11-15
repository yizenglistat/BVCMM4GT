rm(list=ls(all=TRUE))
graphics.off()
cluster <- FALSE
#cluster <- TRUE

if(!cluster){ 
	#setwd("~/Projects/gpvcm")
	task_id <- 1
	nreps <- 2
	knowns <- c(TRUE, FALSE)
	Ns <- c(3000, 5000)
	pool_sizes <- c(5, 10)
	model_names <- c("m1", "m2")
	testings <- c("DT", "AT", "IT")
	known <- FALSE
	N_test <- 600
	sigma <- 0.5
	folder <- 'output'
}else{
	task_id <- as.integer(commandArgs(trailingOnly = TRUE))
	nreps <- 4
	pool_sizes <- c(5, 10)
	Ns <- c(3000,5000)
	model_names <- c("m1", "m2")
	testings <- c("DT", "AT", "IT")
	knowns <- c(TRUE, FALSE)
	N_test <- 600
	folder <- 'output'
}

packages <- c("BayesLogit", "geoR", "ltsa", "mvtnorm", "Matrix","hdf5r","Rcpp","glue", "Hmisc")
temp <- suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
temp <- sapply(list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE), function(x) source(x))
temp <- sapply(list.files("src", pattern="*.cpp$", full.names=TRUE, ignore.case=TRUE), function(x) sourceCpp(x))
rm(temp)

for(rep in 1:nreps){
	for(known in knowns){
		for(N in Ns){
			for(model_name in model_names){
				for(pool_size in pool_sizes){

					#set.seed(4455)

					sim_options <- syn_options(
						N=N, 
						pool_size=pool_size, 
						u_lower=-3, 
						u_upper=3, 
						N_test=N_test, 
						nsites=64, 
						sigma=sigma, 
						se=c(0.95,0.98), 
						sp=c(0.98,0.99)
					)

					full_data <- synthetic(
						model_name=model_name, 
						options=sim_options
					)

					
					for(testing in testings){
						data <- full_data[[testing]]
						print(data$prevalence)
						options <- mcmc_options(
							task_id=nreps*(task_id-1)+rep,
							model_name=model_name, 
							outdir=glue('output/{ifelse(known, "known", "unknown")}/{N}/{model_name}/cj{pool_size}/{testing}/'), 
							nchain=1, 
							nburn=2, 
							nkeep=2, 
							nmem=2, 
							nknots=100,
							nthin=1, 
							ndisp=100, 
							a_se=0.5, b_se=0.5, 
							a_sp=0.5, b_sp=0.5,
							a_sigma2=2, b_sigma2=1, 
							known=known, 
							phi_sd=0.1, kappa=2,
							delete=TRUE, 
							seed=FALSE
						)

						fit <- try(gpp_estimate(data, options), silent=TRUE)
						
						if(all(class(fit)=="try-error")){
							rdata_file <- glue(options$outdir, "fitted_task", options$task_id, ".ERROR")
							if(!file.exists(rdata_file)) file.create(rdata_file)
						}else{
							rdata_file <- glue(options$outdir, "fitted_task", options$task_id, ".RData")
							saveRDS(fit, file=rdata_file)
						}

					}
				}
			}
		}
	}
}
