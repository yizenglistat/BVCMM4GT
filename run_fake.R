rm(list=ls(all=TRUE))
graphics.off()
cluster <- FALSE
#cluster <- TRUE

if(!cluster){
	task_id <- 1
	nreps <- 1
	N_test <- 600
	folder <- 'output'
	model_name <- 'fake'
}else{
	task_id <- as.integer(commandArgs(trailingOnly = TRUE))
	nreps <- 1
	N_test <- 600
	folder <- 'output'
}

packages <- c("BayesLogit", "geoR", "ltsa", "mvtnorm", "Matrix","hdf5r","Rcpp","glue", "Hmisc")
temp <- suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
temp <- sapply(list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE), function(x) source(x))
temp <- sapply(list.files("src", pattern="*.cpp$", full.names=TRUE, ignore.case=TRUE), function(x) sourceCpp(x))
rm(temp)


data <- preprocess_fake(folder='./data', vars=1:6, linear=NULL, twoway=FALSE, N_test=N_test, reformat=FALSE)

# mcmc options
options <- mcmc_options(
	task_id=task_id,
	model_name=model_name, 
	outdir=glue('output/{model_name}/'), 
	nchain=2, 
	nburn=5000, 
	nkeep=2500, 
	nmem=2500, 
	nknots=100,
	nthin=10, 
	ndisp=500, 
	a_se=0.5, b_se=0.5, 
	a_sp=0.5, b_sp=0.5,
	a_sigma2=2, b_sigma2=1, 
	known=FALSE,
	phi_sd=0.1, kappa=1.5,
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

rdata_file <- glue("output/{model_name}/fitted_task1.RData")
fitted <- readRDS(rdata_file)

knots <- rep(seq(-0.8, 1.25, length.out=fit$N_test), times=length(fitted$beta_test_hat_med)/fit$N_test)
knots <- do_unnormalization(knots, log=TRUE)
med <- fitted$beta_test_hat_med
lower <- fitted$beta_test_hat_lower
upper <- fitted$beta_test_hat_upper
label <- rep(c("Intercept", "Race", "New", "Multiple", "Contact", "Symptoms"), each=fitted$N_test)

df <- data.frame(knots=knots, med=med, lower=lower, upper=upper, label=label)
saveRDS(df, glue("output/{model_name}/df_vcm.RData"))
