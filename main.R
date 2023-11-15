rm(list=ls(all=TRUE))
graphics.off()

cluster <- FALSE
#cluster <- TRUE

if(!cluster){ 
	# setwd("YOURLOCALPATH")
	task_id <- 1
	nreps <- 1
	pool_sizes <- c(5)
	model_names <- c("m2")
	testings <- c("DT")
	known <- TRUE
	N_test <- 600
}else{
	task_id <- as.integer(commandArgs(trailingOnly = TRUE))
	nreps <- 4
	pool_sizes <- c(5)
	model_names <- c("m1")
	testings <- c("DT", "AT", "IT", "MPT")
	known <- TRUE
}

packages <- c("BayesLogit", "geoR", "ltsa", "mvtnorm", "Matrix", "hdf5r", "Rcpp", "glue")
temp <- suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
temp <- sapply(list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE), function(x) source(x))
temp <- sapply(list.files("src", pattern="*.cpp$", full.names=TRUE, ignore.case=TRUE), function(x) sourceCpp(x))
rm(temp)

for(rep in 1:nreps){
	for(model_name in model_names){
		for(pool_size in pool_sizes){

			#set.seed(4455)

			sim_options <- syn_options(
				N=5000, 
				pool_size=pool_size, 
				u_lower=-3, 
				u_upper=3, 
				N_test=N_test, 
				nsites=64, 
				sigma=0.5, 
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
					outdir=glue('output/{model_name}/cj{pool_size}/{testing}/'), 
					nchain=1, 
					nburn=2000, 
					nkeep=5000, 
					nmem=1000, 
					nknots=100,
					nthin=1, 
					ndisp=1000, 
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
