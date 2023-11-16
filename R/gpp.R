# nmem must be divisible by nkeep
mcmc_options <- function(
	task_id=NULL, model_name='m1', outdir='output/', 
	nchain=1, nburn=500, nkeep=500, nmem=500, nthin=2, ndisp=100, nknots=100,
	a_se=0.5, b_se=0.5, a_sp=0.5, b_sp=0.5,
	a_sigma2=2, b_sigma2=1, known=FALSE, 
	phi_sd=0.1, kappa=2,
	delete=TRUE, seed=FALSE)
{

	if(nmem>nkeep|!(nkeep%%nmem)) nmem <- nkeep
	
	return(list(task_id=task_id, model_name=model_name, outdir=outdir, 
		nchain=nchain, nburn=nburn, nkeep=nkeep, nmem=nmem, nthin=nthin, ndisp=ndisp, nknots=nknots, 
		known=known, a_se=a_se, b_se=b_se, a_sp=a_sp, b_sp=b_sp, a_sigma2=a_sigma2, b_sigma2=b_sigma2,
		phi_sd=phi_sd, kappa=kappa,
		delete=delete, seed=seed))
}

gpp_mcmc <- function(chain_id=1, data, options=mcmc_options()){

	# data
	name 	<- data$name
	testing <- data$testing
	u 		<- data$u
	X 		<- data$X
	X_lin 	<- if(is.null(data$X_lin)) NULL else data$X_lin
	u_seq 	<- data$u_seq
	S 		<- data$S
	Z 		<- data$Z
	Y		<- data$Y
	pools 	<- data$pools
	npools 	<- data$npools
	nsites 	<- data$nsites
	nassay 	<- length(unique(Z[,3]))
	nbeta 	<- ncol(X)
	N 		<- data$N
	N_test 	<- data$N_test
	pool_size <- data$pool_size
	
	SE 		<- data$se
	SP 		<- data$sp

	# mcmc options
	nburn 	<- options$nburn
	nkeep 	<- options$nkeep
	nthin 	<- options$nthin
	nknots 	<- options$nknots
	ndisp 	<- options$ndisp
	nmem 	<- options$nmem
	known 	<- if(testing %in% c("DT", "AT")) options$known else TRUE
	phi_sd 	<- options$phi_sd
	kappa 	<- options$kappa
	a_se 	<- options$a_se
	b_se 	<- options$b_se
	a_sp 	<- options$a_sp
	b_sp 	<- options$b_sp
	a_sigma2 <- options$a_sigma2
	b_sigma2 <- options$b_sigma2 

	task_id <- options$task_id
	outdir 	<- options$outdir
	model_name <- options$model_name

	# config gpp
	config 			<- gpp_config(t=u, 					# t sequence
								  t_new=u_seq,
							 	  nknots=nknots, 		# number of selected knots
							 	  nbeta=nbeta,			# number of beta including intercept
							 	  phi_sd=phi_sd,		# hyperparameter phi_sd, default is 0.1
							 	  kappa=kappa)			# hyperparameter kappa, default is 2

	nnew 			<- config$nnew 					# *note* number of t_new = number of t_unique 

	t_unique 		<- config$t_unique 				# unique t values in t sequence
	nunique 		<- config$nunique 				# number of unique t values
	t_knots  		<- config$t_knots 				# initial selected knots
	nknots 			<- config$nknots				# number of selected knots
	
	# in order to calculate correlation matrix of GPP easily, compute distance matrix in advance.
	unique_dist 	<- config$unique_dist			# nunique x nunique, L1 dist matrix
	knots_dist  	<- config$knots_dist			# nknot x nknot, L1 dist matrix
	cross_dist  	<- config$cross_dist 			# nunique x nknot, L1 dist matrix
	cross_dist_new 	<- config$cross_dist_new 		# nnew x nknot, L1 dist matrix

	# ---------------------- beta_config matrix description ----------------------- #
			# each columns means a beta function											#
			# row1: gamma distribution shape, prior distribution for tau 					#
			# row2: gamma distribution scale, prior distribution for tau 					#
			# row3: sampled tau from gamma prior distribution based on row1 and row2 		#
			# row4: uniform distribution lower bound for phi 								#
			# row5: uniform distribution upper bound for phi 								#
			# row6: phi averaged value based on uniform distribution based on row4 and row5	#
			# row7: proposed phi for M-H sampler, which made a transformation on row6 		#
			# row8: hyperparameter phi_sd, default 0.1										#
			# row9: hyperparameter kappa, default 2											#
			# ----------------------------------------------------------------------------- #

	beta_config 	<- config$beta_config 			# matrix of beta configuration for gpp 

	R_knots_list 	<- config$R_knots_list 			# list of correlation matrix, nknots x nknots
	R_knots_inv_list<- config$R_knots_inv_list		# list of inverse correlation matrix, nknots x nknots
	R_cross_list 	<- config$R_cross_list 			# list of cross correlation matrix, nunique x nknots
	Q_list 			<- config$Q_list 				# convert beta(t_knots) vector back to beta(t_unique) vector

	beta			<- matrix(0, N, nbeta)				# initial beta(t) matrix
	beta_unique 	<- matrix(0, nunique, nbeta)		# initial beta(t_unique) matrix
	beta_knots 	 	<- matrix(0, nknots, nbeta)			# initial beta(t_knots) matrix

	# initialization
	sigma <- 0.5
	gamma <- rnorm(nsites,0,sigma)

	if(testing=='IT'){
		Y[,1] <- Z[,1]
		
	}else if(testing=='MPT'){
		for(idx in 1:npools){
	  		Y[seq((idx-1)*pool_size+1,idx*pool_size),1] <- Z[idx,1]
		}

	}else if(testing=='DT'){
		if(model_name!="app"){
			Y[(Y[,2]==1),1] <- 0
			Y[(Y[,2]==2),1] <- Z[Y[(Y[,2]==2),4]]
		}
	}else if(testing=='AT'){
		Y[Y[,2]==2,1]	<- 0
		Y[(Y[,2]==3),1]	<- Z[Y[(Y[,2]==3),5]]
	}
	
	omega <- rpg(N, 1, 0)
	y_latent <- Y[,1]
	h <- (y_latent - 0.5)/omega
	rands <- 0
	nalpha <- if(is.null(X_lin)) 1 else ncol(X_lin)
	alpha <- rep(0, nalpha)
	lins <- if(is.null(X_lin)) rep(0, N) else X_lin %*% alpha

	# storage
	if(!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings = FALSE)

	# create h5 chain to save matrices	
	file_name <- paste0(outdir, name, 'chain', chain_id, '_task', task_id, '.h5')
	if(file.exists(file_name)) file.remove(file_name)
	chain <- H5File$new(file_name, mode = "a")	

	# create y_store
	y_store <- chain$create_dataset(name="y_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep,N), maxdims=c(nkeep,N)),chunk_dims=c(nkeep,1))
	y_store_tmp <- matrix(NA, nmem, N)

	# create beta_store
	beta_test_store <- chain$create_dataset(name="beta_test_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep, N_test*nbeta), maxdims=c(nkeep,N_test*nbeta)),chunk_dims=c(nkeep,1))
	beta_test_store_tmp <- matrix(NA, nmem, N_test*nbeta)

	# create gamma_store
	gamma_store <- chain$create_dataset(name="gamma_store", dtype=h5types$H5T_NATIVE_FLOAT,
		space=H5S$new(dims=c(nkeep,nsites), maxdims=c(nkeep,nsites)),chunk_dims=c(nkeep,1))
	gamma_store_tmp <- matrix(NA, nmem, nsites)

	# create alpha_store
	alpha_store <- chain$create_dataset(name="alpha_store", dtype=h5types$H5T_NATIVE_FLOAT,
		space=H5S$new(dims=c(nkeep,nalpha), maxdims=c(nkeep,nalpha)),chunk_dims=c(nkeep,1))
	alpha_store_tmp <- matrix(NA, nmem, nalpha)

	# create se and sp store
	se_store <- chain$create_dataset(name="se_store", dtype=h5types$H5T_NATIVE_FLOAT,
		space=H5S$new(dims=c(nkeep,nassay), maxdims=c(nkeep,nassay)),chunk_dims=c(nkeep,1))
	se_store_tmp <- matrix(NA, nmem, nassay)

	sp_store <- chain$create_dataset(name="sp_store", dtype=h5types$H5T_NATIVE_FLOAT,
		space=H5S$new(dims=c(nkeep,nassay), maxdims=c(nkeep,nassay)),chunk_dims=c(nkeep,1))
	sp_store_tmp <- matrix(NA, nmem, nassay)

	# create sigma_store
	sigma_store <- chain$create_dataset(name="sigma_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep,1), maxdims=c(nkeep,1)),chunk_dims=c(nkeep,1))
	sigma_store_tmp <- matrix(NA, nmem, 1)

	# ----- summary statistics store required ------

	# create prevalence store
	prevalence_store <- chain$create_dataset(name="prevalence_store", dtype=h5types$H5T_NATIVE_FLOAT, 
	space=H5S$new(dims=c(nkeep,1), maxdims=c(nkeep,1)),chunk_dims=c(nkeep,1))
	prevalence_store_tmp <- matrix(NA, nmem, 1)

	# create number of tests stores
	ntests_store <- chain$create_dataset(name="ntests_store", dtype=h5types$H5T_NATIVE_FLOAT, 
	space=H5S$new(dims=c(nkeep,1), maxdims=c(nkeep,1)),chunk_dims=c(nkeep,1))
	ntests_store_tmp <- matrix(NA, nmem, 1)
	
	niter 	<- nburn + nkeep * nthin

	# mcmc
	ikeep <- 0
	slice_end <- 0
	sample_state <- 'burn in'

	for(iter in 1:niter){

		# update beta
		for (d in 1:nbeta){
			Xd				<- X[,d]														
			fixs_d 			<- do_rowsums(matrix(X[,-d],N,nbeta-1), matrix(beta[,-d],N,nbeta-1))

			h_beta			<- h - fixs_d - rands - lins									# construct h_beta = h - h without beta_d
			E_mat 			<- config$E_mat 												# beta_d indicator mat, from t_unique to t seq
			Q_mat 			<- Q_list[[d]] 													# beta_d transform mat, from t_knots to t_unique
			R_knots_mat 	<- R_knots_list[[d]] 											# beta_d knots correlation matrix 
			R_knots_inv_mat <- R_knots_inv_list[[d]]										# beta_d inverse knots correlation matrix
			R_cross_mat 	<- R_cross_list[[d]]											# beta_d nunique x nknots cross corr matrix
			beta_d_config   <- beta_config[,d]												# beta_d config matrix
			
			# update beta_d fun and pour ouput into a list
			beta_d_output 	<- update_beta(center=FALSE, 
											Xd=Xd, omega=omega, h=h, h_beta=h_beta, 
											E_mat=E_mat, Q_mat=Q_mat, 
											knots_dist=knots_dist, cross_dist=cross_dist, 
											R_knots_mat=R_knots_mat, 
											R_knots_inv_mat=R_knots_inv_mat, 
											R_cross_mat=R_cross_mat,
											beta_d_config=beta_d_config)

			beta[,d] 				<- beta_d_output$beta_d 								# update beta_d(t)
			beta_unique[,d] 		<- beta_d_output$beta_d_unique							# update beta_d(t_unique)
			beta_knots[,d] 			<- beta_d_output$beta_d_knots 							# update beta_d(t_knots)
			beta_config[3,d] 		<- beta_d_output$tau 									# update tau in beta_config (d_th)
			beta_config[6,d] 		<- beta_d_output$phi 									# update phi in beta_config (d_th)
			beta_config[7,d] 		<- beta_d_output$trphi									# update trphi in beta_config (d_th)
			R_knots_list[[d]] 		<- beta_d_output$R_knots_mat 							# update knots corr mat (d_th)
			R_knots_inv_list[[d]] 	<- beta_d_output$R_knots_inv_mat						# update inverse knots corr mat (d_th)
			R_cross_list[[d]] 		<- beta_d_output$R_cross_mat							# update nunique x nknots cross corr mat (d_th)
			Q_list[[d]] 			<- beta_d_output$Q_mat 									# update transform mat (d_th)
 		
 		}

 		fixs <- do_rowsums(X, beta)

		# update Se, Sp
		if(!known){
			sesp_config	<- matrix(0,nrow=nassay,ncol=4)										# col1, col2 updated a_Se, b_Se; col3, col4 updated a_Sp, b_Sp					
			sesp_output <- update_sesp(N, npools, Y, Z, sesp_config, nassay) 				# update error prior
			updated_se 	<- matrix(sesp_output[,1:2],nassay, 2)								# updated part of Se prior
			updated_sp 	<- matrix(sesp_output[,3:4],nassay, 2)								# updated part of Sp prior
			se 			<- rbeta(nassay, a_se+updated_se[,1], b_se+updated_se[,2])
			sp 			<- rbeta(nassay, a_sp+updated_sp[,1], b_sp+updated_sp[,2])
		}else{
			se 			<- if(testing=="MPT") SE[1] else if(testing=="IT") SE[2] else SE
			sp 			<- if(testing=="MPT") SP[1] else if(testing=="IT") SP[2] else SP
		}

		# print(glue(
		# 	"iter:{iter} | sigma:{round(sigma,2)} | alpha:{round(alpha,2)}
		# 	se1:{round(se[1],4)} | se2:{round(se[2],4)} | se3:{round(se[3],4)}
		# 	sp1:{round(sp[1],4)} | sp2:{round(sp[2],4)} | sp3:{round(sp[3],4)}"))

		# update auxiliary h (or omega) 
		fixs 			<- do_rowsums(X, beta)
		rands 			<- sapply(S, function(idx) gamma[idx])
		omega 			<- rpg(N, 1, fixs + rands + lins)
		h 				<- (y_latent - 0.5)/omega 	

		# update gamma
		for(site in 1:nsites){
			gamma_sd 	<- sqrt((sigma^2)/(1+(sigma^2)*sum(omega[(S==site)]))) 
			gamma_mu 	<- (gamma_sd^2)*sum(omega[(S==site)]*(h[S==site]-fixs[S==site]-lins[S==site]))
			gamma[site] <- rnorm(1,gamma_mu,gamma_sd)
		}
		rands 			<- sapply(S, function(idx) gamma[idx])

		# update sigma 
		sigma 			<- sqrt(1/rgamma(1,a_sigma2+nsites/2, b_sigma2+sum(gamma^2)/2))

		# update y_latent
		prob  		<- do_logit_inv(fixs+rands+lins) 					
		uniform_var <- runif(N)
		newY 		<- rep(0, N)
		
		Y[,1] 		<- update_y(N, prob, Y, Z, newY, uniform_var, se, sp)
		y_latent 	<- Y[,1]

		# update alpha
		if(nalpha>1){
	 		Sigma_alpha <- solve(t(X_lin)%*%Diagonal(x=omega)%*%X_lin+solve(diag(50/4,nalpha,nalpha)))
	 		mu_alpha 	<- Sigma_alpha%*%t(X_lin)%*%Diagonal(x=omega)%*%(h-fixs-rands)
	 		alpha		<- as.vector(mvtnorm::rmvnorm(1,mu_alpha, as.matrix(Sigma_alpha), method = "svd"))

 			lins 		<- X_lin%*%alpha 
 		}
 		
		# display
		if(iter %% ndisp==0){
			if(iter > nburn){
				sample_state = 'sampling'
			}
			
			verbose = paste0('iteration ',iter, ' (', sample_state,')')
			print(verbose)
		}

		# save
		if(iter > nburn){
			if(iter %% nthin==0){
				ikeep <- ikeep + 1

				beta_test <- c()
				for (d in 1:nbeta){											
					R_knots_mat_new 				<- matern(cross_dist_new, phi=beta_config[6,d], kappa=beta_config[9,d])
					Q_mat_new 						<- R_knots_mat_new%*%R_knots_inv_list[[d]]
					Sigma_mat 						<- rep(1,nnew)-diag(tcrossprod(Q_mat_new,R_knots_mat_new))
					Sigma_mat[Sigma_mat<=0] 		<- 0
					beta_test 						<- c(beta_test, as.vector(rnorm(nnew, Q_mat_new%*%beta_knots[,d], Sigma_mat/beta_config[3,d])))
				}


				y_store_tmp[ikeep, ] 			<- y_latent
				beta_test_store_tmp[ikeep, ]	<- beta_test
				gamma_store_tmp[ikeep, ]		<- gamma
				alpha_store_tmp[ikeep, ]		<- alpha
				sigma_store_tmp[ikeep, ] 		<- sigma
				prevalence_store_tmp[ikeep, ]	<- mean(y_latent)
				ntests_store_tmp[ikeep, ]		<- nrow(Z)

				se_store_tmp[ikeep, ] 			<- se
				sp_store_tmp[ikeep, ] 			<- sp


			}
			if(ikeep==nmem){
				slice_start <- slice_end + 1
				slice_end <- slice_start + nmem - 1
				slice <- slice_start:slice_end
				print(paste0('Storing task ', options$task_id, ' chain ', chain_id, ': ', round(100*slice_end/nkeep,2),'%'))

				# hdf5r
				y_store[slice, ] 			<- y_store_tmp
				gamma_store[slice, ] 		<- gamma_store_tmp
				alpha_store[slice, ]		<- alpha_store_tmp
				sigma_store[slice, ] 		<- sigma_store_tmp
				prevalence_store[slice, ]	<- prevalence_store_tmp
				ntests_store[slice, ]		<- ntests_store_tmp 
				beta_test_store[slice, ]	<- beta_test_store_tmp
				se_store[slice, ] 			<- se_store_tmp
				sp_store[slice, ] 			<- sp_store_tmp

				ikeep <- 0
			}
		}
	}
	chain$close_all()
}

post_mean <- function(data_name=NULL, options=mcmc_options(), param_name, nparam, binary=0){
	
	draws <- array(NA,c(options$nkeep, nparam, options$nchain))

	for(chain_id in 1:options$nchain){
		file_name <- paste0(options$outdir, data_name, 'chain', chain_id, '_task', options$task_id, '.h5')
		chain <- H5File$new(file_name, mode="r+")
		#draws[,,chain_id] <- h5read(file_name, paste0(param_name,'_store'))
		draws[,,chain_id] <- chain[[paste0(param_name,'_store')]][,]
		if(binary>0){
			for(binary_idx in 1:binary){
				zero <- 2*binary_idx-1
				ones <- 2*binary_idx
				difference <- draws[,ones,chain_id] - draws[,zero,chain_id]
				draws[,zero,chain_id] <- difference
				draws[,ones,chain_id] <- difference
			}
		}
		chain$close_all()
	}

	post_mean <- apply(draws, 2, mean)
	post_med <- apply(draws, 2, median)
	post_std <- apply(draws, 2, sd)
	post_lower <- apply(draws, 2, quantile, 0.025, type=1)
	post_upper <- apply(draws, 2, quantile, 0.975, type=1)
	
	return(list(mean=post_mean,median=post_med,std=post_std,lower=post_lower,upper=post_upper))
}

gpp_estimate <- function(data, options=mcmc_options()){
	
	if(options$seed){
		set.seed(options$seed)
	}

	tic <- proc.time()

	# paralell computing
	for(chain_id in 1:options$nchain){
		gpp_mcmc(chain_id, data, options)
	}

	toc <- as.vector(proc.time() - tic)[3]

	print(glue('TASK [{options$task_id}]; estimation time [min]: {round(toc/60,2)}'))

	post_y_hat <- post_mean(data$name, options, 'y', data$N)
	y_hat <- post_y_hat$mean
	y_hat_med <- post_y_hat$median
	y_hat_std <- post_y_hat$std
	y_hat_lower <- post_y_hat$lower
	y_hat_upper <- post_y_hat$upper

	post_beta_test_hat <- post_mean(data$name, options, 'beta_test', data$N_test*ncol(data$X))
	beta_test_hat <- post_beta_test_hat$mean
	beta_test_hat_med <- post_beta_test_hat$median
	beta_test_hat_std <- post_beta_test_hat$std
	beta_test_hat_lower <- post_beta_test_hat$lower
	beta_test_hat_upper <- post_beta_test_hat$upper

	post_se_hat <- post_mean(data$name, options, 'se', length(unique(data$Z[,3])))
	se_hat <- post_se_hat$mean
	se_hat_med <- post_se_hat$median
	se_hat_std <- post_se_hat$std
	se_hat_lower <- post_se_hat$lower
	se_hat_upper <- post_se_hat$upper

	post_sp_hat <- post_mean(data$name, options, 'sp', length(unique(data$Z[,3])))
	sp_hat <- post_sp_hat$mean
	sp_hat_med <- post_sp_hat$median
	sp_hat_std <- post_sp_hat$std
	sp_hat_lower <- post_sp_hat$lower
	sp_hat_upper <- post_sp_hat$upper

	post_sigma_hat <- post_mean(data$name, options, 'sigma', 1)
	sigma_hat <- post_sigma_hat$mean
	sigma_hat_med <- post_sigma_hat$median
	sigma_hat_std <- post_sigma_hat$std
	sigma_hat_lower <- post_sigma_hat$lower
	sigma_hat_upper <- post_sigma_hat$upper

	post_alpha_hat <- post_mean(data$name, options, 'alpha', if(is.null(data$X_lin)) 1 else ncol(data$X_lin))
	alpha_hat <- post_alpha_hat$mean
	alpha_hat_med <- post_alpha_hat$median
	alpha_hat_std <- post_alpha_hat$std
	alpha_hat_lower <- post_alpha_hat$lower
	alpha_hat_upper <- post_alpha_hat$upper

	post_gamma_hat <- post_mean(data$name, options, 'gamma', data$nsites)
	gamma_hat <- post_gamma_hat$mean
	gamma_hat_med <- post_gamma_hat$median
	gamma_hat_std <- post_gamma_hat$std
	gamma_hat_lower <- post_gamma_hat$lower
	gamma_hat_upper <- post_gamma_hat$upper

	post_prevalence_hat <- post_mean(data$name, options, 'prevalence', 1)
	prevalence_hat <- post_prevalence_hat$mean
	prevalence_hat_med <- post_prevalence_hat$median
	prevalence_hat_std <- post_prevalence_hat$std
	prevalence_hat_lower <- post_prevalence_hat$lower
	prevalence_hat_upper <- post_prevalence_hat$upper

	post_ntests_hat <- post_mean(data$name, options, 'ntests', 1)
	ntests_hat <- post_ntests_hat$mean
	ntests_hat_med <- post_ntests_hat$median
	ntests_hat_std <- post_ntests_hat$std
	ntests_hat_lower <- post_ntests_hat$lower
	ntests_hat_upper <- post_ntests_hat$upper

	if(options$delete){
		for(chain_id in 1:options$nchain){
			file_name <- paste0(options$outdir, data$name, 'chain', chain_id, '_task', options$task_id, '.h5')
			file.remove(file_name)
		}
	}

	return(list(time=toc/options$nchain,N_test=data$N_test,
		y_hat=y_hat, y_hat_med=y_hat_med, y_hat_std=y_hat_std, y_hat_lower=y_hat_lower, y_hat_upper=y_hat_upper,
		sigma_hat=sigma_hat, sigma_hat_med=sigma_hat_med, sigma_hat_std=sigma_hat_std, sigma_hat_lower=sigma_hat_lower, sigma_hat_upper=sigma_hat_upper,
		gamma_hat=gamma_hat, gamma_hat_med=gamma_hat_med, gamma_hat_std=gamma_hat_std, gamma_hat_lower=gamma_hat_lower, gamma_hat_upper=gamma_hat_upper,
		alpha_hat=alpha_hat, alpha_hat_med=alpha_hat_med, alpha_hat_std=alpha_hat_std, alpha_hat_lower=alpha_hat_lower, alpha_hat_upper=alpha_hat_upper,
		beta_test_hat=beta_test_hat, beta_test_hat_med=beta_test_hat_med, beta_test_hat_std=beta_test_hat_std, beta_test_hat_lower=beta_test_hat_lower, beta_test_hat_upper=beta_test_hat_upper,
		se_hat=se_hat, se_hat_med=se_hat_med, se_hat_std=se_hat_std, se_hat_lower=se_hat_lower, se_hat_upper=se_hat_upper,
		sp_hat=sp_hat, sp_hat_med=sp_hat_med, sp_hat_std=sp_hat_std, sp_hat_lower=sp_hat_lower, sp_hat_upper=sp_hat_upper,
		prevalence_hat=prevalence_hat, prevalence_hat_med=prevalence_hat_med, prevalence_hat_std=prevalence_hat_std, prevalence_hat_lower=prevalence_hat_lower, prevalence_hat_upper=prevalence_hat_upper,
		ntests_hat=ntests_hat, ntests_hat_med=ntests_hat_med, ntests_hat_std=ntests_hat_std, ntests_hat_lower=ntests_hat_lower, ntests_hat_upper=ntests_hat_upper
		))
}















