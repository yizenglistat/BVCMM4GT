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

# posterior processing

# ------------------------------------- figures
# pool_sizes <- c(5, 10)
# folder <- "output"
# model_names <- c('m1', 'm2')
# sigma <- 0.5
# nreps <- 5
# testings <- c("IT", "DT", "AT")
# knowns <- c(TRUE, FALSE)
# Ns <- c(3000, 5000)

for(known in knowns){
	for(N in Ns){
		for(model_name in model_names){
			df <- c()
			for(pool_size in pool_sizes){
				for(testing in testings){
					print(glue("processing posterior {testing} data ..."))
					res <- post_summary(folder=folder,
						known=ifelse(known, "known", "unknown"),
						sample_size=N, 
						model_name=model_name, 
						pool_size=pool_size, testing=testing, 
						sigma=sigma, nreps=nreps, 
						cache=TRUE, replace=FALSE)
					is_pts <- is.na(res$levs)
					knots <- res$levs[!is_pts]
					med <- res$med[!is_pts]
					true <- res$true[!is_pts]
					lower <- res$lower[!is_pts]
					upper <- res$upper[!is_pts]
					df_testing <- data.frame(
						knots=knots,
						med=med,
						true=true,
						lower=lower,
						upper=upper,
						label=rep(c("f0","f1","f2"),each=length(knots)/3),
						testing=testing,
						pool=pool_size
					)
					df <- rbind(df, df_testing)
					print(glue("{testing} data done!"))
				}
			}
			rdata_file <- glue("{folder}/{ifelse(known, 'known', 'unknown')}/{N}/{model_name}/df_pool.RData")
			saveRDS(df, rdata_file)
		}
	}
}



# ------------------------------------- tables (known)
# pool_sizes <- c(5, 10)
# folder <- "output"
# model_names <- c('m1', 'm2')
# sigma <- 0.5
# nreps <- 5
# testings <- c("IT", "DT", "AT")
# knowns <- c(TRUE, FALSE)
# Ns <- c(3000, 5000)

df <- c()
dataset <- matrix(NA, ncol=5, nrow=47) # data part
colnames(dataset) <- c("IT", "DT", "AT", "DT", "AT")
rm_names <- function(x) {x[c(FALSE, x[-1]==x[-length(x)])] <- ""; x} 
front_cols <- apply(expand.grid(
	Summary=c("Bias (CP95)", "SSD (ESE)", "", "Average Number of Tests (savings in %)", "Average Time (in minutes)", ""), 
	Parameter=c("$\\sigma$", ""),
	Model=c("M1", "M2"),
	Samplesize=c("$3000$", "$5000$"))[, 4:1], 2, rm_names)
front_cols <- head(front_cols,-1)
colnames(front_cols) <- c("Sample Size", "Model", "Parameter", "Summary")
tabset <- cbind(front_cols, dataset)
tabset <- tabset[-c(4:9, 18:23, 30:35, 42:47), ]

for(N in Ns){
	for(model_name in model_names){
		for(pool_size in pool_sizes){
			for(testing in testings){
				print(glue("processing posterior {testing} data ..."))
				#if(testing=='IT') nreps <- 35 else nreps <- 136
				res <- post_summary(folder=folder,
						known="known",
						sample_size=N, 
						model_name=model_name, 
						pool_size=pool_size, testing=testing, 
						sigma=sigma, nreps=nreps, 
						cache=TRUE, replace=FALSE)
				
				is_pts 	<- is.na(res$levs)
				
				#table
				bias 	<- res$bias[is_pts][1]
				ese	 	<- res$ese[is_pts][1]
				ssd 	<- res$ssd[is_pts][1]
				eci 	<- apply(res$param_eci, 1, mean)[1]

				ext 	<- apply(res$param_ext, 1, mean)
				avg_time <- round(ext[4]/60,2)
				if(testing == "IT"){
					save_pct <- substr(glue_tab_item(round(100*(1-ext[5]/N),2),0),1,4)
				}else{
					save_pct <- substr(glue_tab_item(round(100*(1-ext[5]/N),2),0),1,5)
				}
				avg_tests <- glue("{round(ext[4],2)}({save_pct}%)")


				if(N==5000){
					if(model_name=='m1'){
						if(testing=='IT'){
							tabset[1, 5] <- glue_tab_item(bias, eci)
							tabset[2, 5] <- glue_tab_item(ssd, ese)
							tabset[4, 5] <- avg_tests
							tabset[5, 5] <- avg_time
						}
						if(testing=="DT"&pool_size==5){
							tabset[1, 6] <- glue_tab_item(bias, eci)
							tabset[2, 6] <- glue_tab_item(ssd, ese)
							tabset[4, 6] <- avg_tests
							tabset[5, 6] <- avg_time
						}
						if(testing=="AT"&pool_size==10){
							tabset[1, 7] <- glue_tab_item(bias, eci)
							tabset[2, 7] <- glue_tab_item(ssd, ese)
							tabset[4, 7] <- avg_tests
							tabset[5, 7] <- avg_time
						}
						if(testing=="DT"&pool_size==5){
							tabset[1, 8] <- glue_tab_item(bias, eci)
							tabset[2, 8] <- glue_tab_item(ssd, ese)
							tabset[4, 8] <- avg_tests
							tabset[5, 8] <- avg_time
						}
						if(testing=="AT"&pool_size==10){
							tabset[1, 9] <- glue_tab_item(bias, eci)
							tabset[2, 9] <- glue_tab_item(ssd, ese)
							tabset[4, 9] <- avg_tests
							tabset[5, 9] <- avg_time
						}
					}
					if(model_name=='m2'){
						if(testing=='IT'){
							tabset[7, 5] <- glue_tab_item(bias, eci)
							tabset[8, 5] <- glue_tab_item(ssd, ese)
							tabset[10, 5] <- avg_tests
							tabset[11, 5] <- avg_time
						}
						if(testing=="DT"&pool_size==5){
							tabset[7, 6] <- glue_tab_item(bias, eci)
							tabset[8, 6] <- glue_tab_item(ssd, ese)
							tabset[10, 6] <- avg_tests
							tabset[11, 6] <- avg_time
						}
						if(testing=="AT"&pool_size==10){
							tabset[7, 7] <- glue_tab_item(bias, eci)
							tabset[8, 7] <- glue_tab_item(ssd, ese)
							tabset[10, 7] <- avg_tests
							tabset[11, 7] <- avg_time
						}
						if(testing=="DT"&pool_size==5){
							tabset[7, 8] <- glue_tab_item(bias, eci)
							tabset[8, 8] <- glue_tab_item(ssd, ese)
							tabset[10, 8] <- avg_tests
							tabset[11, 8] <- avg_time
						}
						if(testing=="AT"&pool_size==10){
							tabset[7, 9] <- glue_tab_item(bias, eci)
							tabset[8, 9] <- glue_tab_item(ssd, ese)
							tabset[10, 9] <- avg_tests
							tabset[11, 9] <- avg_time
						}
					}
				}

				if(N==5000){
					if(model_name=='m1'){
						if(testing=='IT'){
							tabset[13, 5] <- glue_tab_item(bias, eci)
							tabset[14, 5] <- glue_tab_item(ssd, ese)
							tabset[16, 5] <- avg_tests
							tabset[17, 5] <- avg_time
						}
						if(testing=="DT"&pool_size==5){
							tabset[13, 6] <- glue_tab_item(bias, eci)
							tabset[14, 6] <- glue_tab_item(ssd, ese)
							tabset[16, 6] <- avg_tests
							tabset[17, 6] <- avg_time
						}
						if(testing=="AT"&pool_size==10){
							tabset[13, 7] <- glue_tab_item(bias, eci)
							tabset[14, 7] <- glue_tab_item(ssd, ese)
							tabset[16, 7] <- avg_tests
							tabset[17, 7] <- avg_time
						}
						if(testing=="DT"&pool_size==5){
							tabset[13, 8] <- glue_tab_item(bias, eci)
							tabset[14, 8] <- glue_tab_item(ssd, ese)
							tabset[16, 8] <- avg_tests
							tabset[17, 8] <- avg_time
						}
						if(testing=="AT"&pool_size==10){
							tabset[13, 9] <- glue_tab_item(bias, eci)
							tabset[14, 9] <- glue_tab_item(ssd, ese)
							tabset[16, 9] <- avg_tests
							tabset[17, 9] <- avg_time
						}
					}
					if(model_name=='m2'){
						if(testing=='IT'){
							tabset[19, 5] <- glue_tab_item(bias, eci)
							tabset[20, 5] <- glue_tab_item(ssd, ese)
							tabset[22, 5] <- avg_tests
							tabset[23, 5] <- avg_time
						}
						if(testing=="DT"&pool_size==5){
							tabset[19, 6] <- glue_tab_item(bias, eci)
							tabset[20, 6] <- glue_tab_item(ssd, ese)
							tabset[22, 6] <- avg_tests
							tabset[23, 6] <- avg_time
						}
						if(testing=="AT"&pool_size==10){
							tabset[19, 7] <- glue_tab_item(bias, eci)
							tabset[20, 7] <- glue_tab_item(ssd, ese)
							tabset[22, 7] <- avg_tests
							tabset[23, 7] <- avg_time
						}
						if(testing=="DT"&pool_size==5){
							tabset[19, 8] <- glue_tab_item(bias, eci)
							tabset[20, 8] <- glue_tab_item(ssd, ese)
							tabset[22, 8] <- avg_tests
							tabset[23, 8] <- avg_time
						}
						if(testing=="AT"&pool_size==10){
							tabset[19, 9] <- glue_tab_item(bias, eci)
							tabset[20, 9] <- glue_tab_item(ssd, ese)
							tabset[22, 9] <- avg_tests
							tabset[23, 9] <- avg_time
						}
					}
				}
				print(glue("{testing} data done!"))
			}
		}
	}
}

latex(tabset,
      file="tab_known.tex",
      cgroup=c("","c=1", "c=5", "c=10"),
      n.cgroup=c(4, 1, 2, 2),
      na.blank=TRUE,
      #rowlabel=c("Sample Size", "Model", "Parameter", "Summary"),
      booktabs=TRUE,
      rownames=NULL,
      collabel.just=rep("c", 9),
      col.just=c(rep("c",4), rep("r", 8)),
      where="!htbp",
      caption="Simulation results for models M1 and M2 with sample sizes $N\\in\\{3000,5000\\}$, $c\\in\\{5,10\\}$ under both DT and AT protocols, when assay sensitivities and specificities are known, respectively. 
      Note the summary statistics for $\\sigma$ include average bias (Bias), sample standard deviation (SSD) of the 500 posterior median estimates, the average of the 500 estimates of the posterior standard deviation (ESE), 
      and empirical coverage probability of 95\\% credible intervals (CP95).",
      bold.header=FALSE,
)

# ------------------------------------- tables (unknown)
# pool_sizes <- c(5, 10)
# folder <- "output"
# model_names <- c('m1', 'm2')
# sigma <- 0.5
# nreps <- 5
# testings <- c("IT", "DT", "AT")
# knowns <- c(TRUE, FALSE)
# Ns <- c(3000, 5000)

df <- c()
dataset <- matrix(NA, ncol=5, nrow=119) # data part
colnames(dataset) <- c("IT", "DT", "AT", "DT", "AT")
rm_names <- function(x) {x[c(FALSE, x[-1]==x[-length(x)])] <- ""; x} 
front_cols <- apply(expand.grid(
	Summary=c("Bias (CP95)", "SSD (ESE)", "", "Average Number of Tests (savings in %)", "Average Time (in minutes)", ""), 
	Parameter=c("$\\sigma$", "$S_{e(1)}$", "$S_{e(2)}$", "$S_{p(1)}$", "$S_{p(2)}$"),
	Model=c("M1", "M2"),
	Samplesize=c("$3000$", "$5000$"))[, 4:1], 2, rm_names)
front_cols <- head(front_cols,-1)
colnames(front_cols) <- c("Sample Size", "Model", "Parameter", "Summary")
tabset <- cbind(front_cols, dataset)
tabset <- tabset[-c(4:6, 10:12, 16:18, 22:24,
	34:36, 40:42, 46:48, 52:54,
	64:66, 70:72, 76:78, 82:84,
	94:96, 100:102, 106:108, 112:114), ]

for(N in Ns){
	for(model_name in model_names){
		for(pool_size in pool_sizes){
			for(testing in testings){
				print(glue("processing posterior {testing} data ..."))
				#if(testing=='IT') nreps <- 35 else nreps <- 136
				res <- post_summary(folder=folder,
						known="unknown",
						sample_size=N, 
						model_name=model_name, 
						pool_size=pool_size, testing=testing, 
						sigma=sigma, nreps=nreps, 
						cache=TRUE, replace=FALSE)
				
				is_pts 	<- is.na(res$levs)
				
				#table

				if(testing %in% c("DT", "AT") ){
					# sigma, Se, Sp
					bias 	<- res$bias[is_pts]
					ese	 	<- res$ese[is_pts]
					ssd 	<- res$ssd[is_pts]
					eci 	<- apply(res$param_eci, 1, mean)
				}else{
					bias 	<- ""
					ese 	<- ""
					ssd 	<- ""
					eci 	<- ""
				}

				ext 	<- apply(res$param_ext, 1, mean)
				avg_time <- round(ext[4]/60,2)
				if(testing == "IT"){
					save_pct <- substr(glue_tab_item(round(100*(1-ext[5]/N),2),0),1,4)
				}else{
					save_pct <- substr(glue_tab_item(round(100*(1-ext[5]/N),2),0),1,5)
				}
				avg_tests <- glue("{round(ext[4],2)}({save_pct}%)")


				if(N==3000){
					if(model_name=='m1'){
						if(testing=='IT'){
							tabset[1:17, 5] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="DT"&pool_size==5){
							tabset[1:17, 6] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="AT"&pool_size==10){
							tabset[1:17, 7] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="DT"&pool_size==5){
							tabset[1:17, 8] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="AT"&pool_size==10){
							tabset[1:17, 9] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
					}
					if(model_name=='m2'){
						if(testing=='IT'){
							tabset[19:35, 5] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="DT"&pool_size==5){
							tabset[19:35, 6] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="AT"&pool_size==10){
							tabset[19:35, 7] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="DT"&pool_size==5){
							tabset[19:35, 8] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="AT"&pool_size==10){
							tabset[19:35, 9] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
					}
				}

				if(N==5000){
					if(model_name=='m1'){
						if(testing=='IT'){
							tabset[37:53, 5] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="DT"&pool_size==5){
							tabset[37:53, 6] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="AT"&pool_size==10){
							tabset[37:53, 7] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="DT"&pool_size==5){
							tabset[37:53, 8] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="AT"&pool_size==10){
							tabset[37:53, 9] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
					}
					if(model_name=='m2'){
						if(testing=='IT'){
							tabset[55:71, 5] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="DT"&pool_size==5){
							tabset[55:71, 6] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="AT"&pool_size==10){
							tabset[55:71, 7] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="DT"&pool_size==5){
							tabset[55:71, 8] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
						if(testing=="AT"&pool_size==10){
							tabset[55:71, 9] <- join_tab_item(c(glue_tab_item(bias, eci), avg_tests), c(glue_tab_item(ssd, ese), avg_time))
						}
					}
				}
				print(glue("{testing} data done!"))
			}
		}
	}
}

latex(tabset,
      file="tab_unknown.tex",
      cgroup=c("","c=1", "c=5", "c=10"),
      n.cgroup=c(4, 1, 2, 2),
      na.blank=TRUE,
      #rowlabel=c("Sample Size", "Model", "Parameter", "Summary"),
      booktabs=TRUE,
      rownames=NULL,
      collabel.just=rep("c", 9),
      col.just=c(rep("c",4), rep("r", 8)),
      where="!htbp",
      caption="Simulation results for models M1 and M2 with sample sizes $N\\in\\{3000,5000\\}$, $c\\in\\{5,10\\}$ under both DT and AT protocols, when assay sensitivities and specificities are unknown, respectively. 
      Note the summary statistics for $\\sigma$ include average bias (Bias), sample standard deviation (SSD) of the 500 posterior median estimates, the average of the 500 estimates of the posterior standard deviation (ESE), 
      and empirical coverage probability of 95\\% credible intervals (CP95).",
      bold.header=FALSE,
)
